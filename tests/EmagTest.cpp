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
#include "NURBSCache.h"

using namespace nurbs;

int main(int argc, char* argv[])
{
    // initiate HLib library
    HLIB::INIT();
    HLIB::CFG::set_verbosity( 3 );
    
//    HLIB::CFG::set_nthreads(1);
    
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
    
    //g.normalise();
//        g.rescale(2.0/34.4);
    
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

//    multiforest.prefine(1);
    multiforest.hrefine(refine);
    
    std::cout << "Performing emag scattering analysis on multiforest with "
              << multiforest.elemN() << " elements, "
              << multiforest.globalDofN() << " dof.....\n";
    
    // Some hardcoded input parameters
    const double k = std::atof(argv[2]);
    
//    // Polarisation vector and wavevector for sphere problem
//    const Point3D pvec(0.0, 0.0, 1.0);
//    const Point3D kvec(k, 0.0, 0.0);
//    
//    const double omega = k;
//    const double mu = 1.0;
    
    const double mu = 1.257e-6;
    const double epsilon = 8.854e-12;
    const double omega = k / std::sqrt(mu * epsilon);
    
    const Point3D kvec(k, 0.0, 0.0);
    const Point3D pvec(0.0, 1.0, 0.0);

    // Create class for assembling Hmatrix and force vector
    fastibem::HAssemblyEmag hassembly(multiforest,
                                      kvec,
                                      pvec,
                                      mu,
                                      omega);
    
    auto A = hassembly.assembleHmatrix();

    // Force vector
    auto f = A->row_vector();
    hassembly.assembleForceVector(f.get());

	// clear the cache
	nurbshelper::NURBSCache::Instance().clear();
	
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
        std::vector<std::complex<double>> solnvec;
        
        std::cout << "solution\n\n";
        for(size_t i = 0; i < x->size(); ++i)
        {
            const auto entry = x->centry(i);
            solnvec.push_back(std::complex<double>(entry.re(), entry.im()));
            //solnvec.push_back(std::complex<double>(1.0, 0.0));
            std::cout << entry.re() << "," << entry.im() << "\n";
        }
        
        // Output solution
        nurbs::OutputVTK output("emagtest");
        output.outputComplexVectorField(multiforest, "surface_current", solnvec);
    }
    else
        std::cout << "  not converged in " << timer << " and "
        << solve_info.n_iter() << " steps " << std::endl;

    
    return EXIT_SUCCESS;
}
