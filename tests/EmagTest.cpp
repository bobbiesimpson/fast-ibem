#include <iostream>
#include "HAssembly.h"
#include "Functor.h"
#include "Kernel.h"
#include "Point3D.h"
#include "Geometry.h"
#include "HConformingForest.h"
#include "MaterialParam.h"
#include "OutputVTK.h"
#include "hlib.hh"

using namespace nurbs;

int main(int argc, char* argv[])
{
    // initiate HLib library
    HLIB::INIT();
    HLIB::CFG::set_verbosity( 3 );
    
    if(argc < 3)
    {
        std::cerr << "Please run as ./emagtest <input_file> <wavenumber> <optional: hrefinement>\n";
        return EXIT_FAILURE;
    }
    
    // TODO: bundle all the input tasks into a single Input class possibly
    // returning a std::unique_ptr<MultiForest> or std::unique_ptr<Forest>
    
    // First read in the geometry file
    std::cout << "Trying to open hbs input file....\n";
    std::ifstream ifs(argv[1]);
    if(!ifs)
        error("Cannot open file for reading\n");
    Geometry g;
    if(!g.loadHBSFile(ifs))
        error("Failed to load geometry from hbs data");
    
    // Construct the necessary forests
    HDivForest multiforest(g);
    
    // Apply hrefinement
    uint refine = 0;
    const uint max_refine = 10;
    if(argc > 3) {
        auto input = std::atoi(argv[3]);
        if(input < 0)
            std::cout << "Cannot supplied negative refinement integer. Carrying on with no h-refinement";
        if(input > max_refine) {
            std::cout << "Truncating refinement to " << max_refine << " levels.";
            refine = max_refine;
        }
        else {
            std::cout << "Applying " << input << " levels of h-refinement\n";
            refine = input;
        }
    }
    multiforest.hrefine(refine);
    
    // Some hardcoded input parameters
    const double k = std::atof(argv[2]);
    
    // Polarisation vector and wavevector
    const Point3D pvec(1.0, 0.0, 0.0);
    const Point3D kvec(0.0, 0.0, k);
    
    const double omega = k;
    const double mu = 1.0;
    
    // Create class for assembling Hmatrix and force vector
    fastibem::HAssemblyEmag hassembly(multiforest,
                                      kvec,
                                      pvec,
                                      mu,
                                      omega);
    
    auto A = hassembly.assembleHmatrix();
    
    
    std::cout << "kvalues ....\n\n";
    for(size_t i = 0; i < 10; ++i)
    {
        for(size_t j = 0; j < 10; ++j)
            std::cout << A->centry(i, j) << "\t";
        std::cout << "\n";
    }
    
    // Force vector
    auto f = A->row_vector();
    hassembly.assembleForceVector(f.get());
    
    // Now solve!
    HLIB::TTimer                    timer( HLIB::WALL_TIME );
    HLIB::TConsoleProgressBar       progress;
    
    // Preconditioning
    auto  B = A->copy();
    timer.start();
    auto  A_inv = HLIB::factorise_inv( B.get(), hassembly.trunAccInstance(), & progress );
    timer.pause();
    
    std::cout << " done in " << timer << std::endl;
    std::cout << "    size of LU factor = " << HLIB::Mem::to_string( B->byte_size() ) << std::endl;
    std::cout << "    inversion error   = " << std::scientific << std::setprecision( 4 )
    << HLIB::inv_approx_2( A.get(), A_inv.get() ) << std::endl;
    
    // GMRES solve
    std::cout << std::endl << "━━ solving system" << std::endl;
    HLIB::TAutoSolver     solver( 1000 );
    HLIB::TSolver::TInfo  solve_info( false, HLIB::verbose( 2 ) );
    
    auto x = A->col_vector();
    
    timer.start();
    solver.solve( A.get(), x.get(), f.get(), A_inv.get(), &solve_info );
    
    if ( solve_info.converged() ) {
        std::cout << "  converged in " << timer << " and "
        << solve_info.n_iter() << " steps with rate " << solve_info.conv_rate()
        << ", |r| = " << solve_info.res_norm() << std::endl;
        
        hassembly.clusterTree()->perm_i2e()->permute( x.get() );
        
        // create output vector and set to VTK file
        const auto n = multiforest.globalDofN();
        std::vector<std::complex<double>> solnvec(n);
        
        std::cout << "solution\n\n";
        for(size_t i = 0; i < n; ++i)
        {
            const auto entry = x->centry(i);
            solnvec[i] = std::complex<double>(entry.re(), entry.im());
            std::cout << entry.re() << "," << entry.im() << "\n";
        }
        
        // Output solution
        nurbs::OutputVTK output("sphere_test");
        //output.outputComplexVectorField(multiforest, "surface current", solnvec);
        //output.outputComplexAnalysisField(forest, "acoustic potential", solnvec);
    }
    else
        std::cout << "  not converged in " << timer << " and "
        << solve_info.n_iter() << " steps " << std::endl;

    
    return EXIT_SUCCESS;
}