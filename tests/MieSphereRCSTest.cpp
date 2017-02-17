#include <iostream>
#include <string>

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
#include "Norm.h"

#include <boost/filesystem.hpp>

using namespace boost::filesystem;

using namespace nurbs;

int main(int argc, char* argv[])
{
    // initiate HLib library
    HLIB::INIT();
    HLIB::CFG::set_verbosity( 3 );
    
    if(argc < 3)
    {
        std::cerr << "Please run as ./miespherercstest <input_file> <wavenumber min> <wavenumber max> <nsample> <optional: hrefinement>\n";
        return EXIT_FAILURE;
    }
    
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
    
    // Get wavenumber inputs
    const double kmin = std::atof(argv[2]);
    const double kmax = std::atof(argv[3]);
    if(kmax < kmin)
        throw std::runtime_error("Kmax cannot be less that kmin");
    
    const int nsample = std::atoi(argv[4]);
    
    // Apply hrefinement
    uint refine = 0;
    const uint max_refine = 10;
    if(argc > 5)
    {
        auto input = std::atoi(argv[5]);
        if(input < 0)
            std::cout << "Cannot supplied negative refinement integer. Carrying on with no h-refinement";
        if(input > max_refine)
        {
            std::cout << "Truncating refinement to " << max_refine << " levels.";
            refine = max_refine;
        }
        else
        {
            std::cout << "Applying " << input << " levels of h-refinement\n";
            refine = input;
        }
    }
    
    multiforest.hrefine(refine);
    
    std::cout << "Performing emag scattering analysis on multiforest with "
    << multiforest.elemN() << " elements, "
    << multiforest.globalDofN() << " dof.....\n";
    
    // Some hardcoded input parameters
    const Point3D observe(0.0,0.0,-1.0e3);
    const Point3D rhat(0.0,0.0,-1.0);
    const std::vector<std::complex<double>> pvec
    {
        std::complex<double>(-1.0, 0.0),
        std::complex<double>(0.0, 0.0),
        std::complex<double>(0.0, 0.0)
    };

    const double mu = 1.0;
    
    std::string filename("mie-sphere-rcs.dat");
    std::ofstream ofs( filename.c_str() );
    if( !ofs )
        error( "Error opening file" );
    ofs.precision( 18 );
    std::cout.precision( 18 );
    
    double k = kmin;
    const double kinc = (kmax - kmin) / nsample;
    
    
    // loop over wavenumber samples
    while(k < kmax)
    {
        const double omega = k;
        Point3D kvec(0.0, 0.0, k);
        
        // Create class for assembling Hmatrix and force vector
        fastibem::HAssemblyEmag hassembly(multiforest,
                                          kvec,
                                          pvec,
                                          mu,
                                          k);

        std::unique_ptr<HLIB::TMatrix> A;
        
        // attempt to read matrix
        std::string filename(argv[1]);
        size_t lastindex = filename.find_last_of(".");
        std::string rawname = filename.substr(0, lastindex);
        rawname = rawname + "_h" + std::to_string(refine) + "_k" + std::to_string(k) + "_dof" + std::to_string(multiforest.globalDofN());
        {
            HLIB::THLibMatrixIO io;
            
            if(exists(rawname))
            {
                std::cout << "Trying to read existing Hmatrix file: " << rawname << "....\n";
                A = io.read(rawname);
            }
            else
            {
                A = hassembly.assembleHmatrix();
                io.write(A.get(), rawname);
            }
        }
        
        //auto A = hassembly.assembleHmatrix();
        
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
            }
            
//            // Output solution
            const std::string fname = "miesphere-" + std::to_string(k);
            nurbs::OutputVTK output(fname, 2);
            output.outputComplexVectorField(multiforest, "surface_current", solnvec);
            
            // Uncomment to output analytical mie solution
            //output.outputAnalyticalMieComplexVectorField(multiforest, "exact mie", k);
            
            // Get Radar cross section data
            double rcs = output.computeRCS(multiforest, observe, k, rhat, mu , omega , solnvec);
            std::cout<< "k = "<< k << "\t" << "rcs = "<< rcs/PI <<"\n";
            ofs << k << "\t" << rcs/PI << "\n";
            
            // Calculate surface norm
            std::cout << "L2 graph norm\t" << nurbs::L2graphNormMieSphere(multiforest, k, solnvec) << "\n";
            
        }
        else
            std::cout << "  not converged in " << timer << " and " << solve_info.n_iter() << " steps " << std::endl;

        nurbshelper::NURBSCache::Instance().clear();
        
        k += kinc;
        
    }
    return EXIT_SUCCESS;
}