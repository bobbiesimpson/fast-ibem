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
        std::cerr << "Please run as ./emagtest <input_file> <wavenumber> <optional: hrefinement>\n";
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
    
    // hardcoded rotation and scaling for stealth model
    //g.rotate(nurbs::CartesianComponent::X, nurbs::PI * 0.5);
    
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
    //multiforest.graded_hrefine(refine, 0.25);
    
    std::cout << "Performing emag scattering analysis on multiforest with "
    << multiforest.elemN() << " elements, "
    << multiforest.globalDofN() << " dof.....\n";
    
    // Some hardcoded input parameters
    const double k = std::atof(argv[2]);
    
    const double mu = 1.25663706e-6;
    const double epsilon = 8.85418782e-12;

//    const double mu = 1.0;
//    const double epsilon = 1.0;
    
    const double omega = k / std::sqrt(mu * epsilon);
    
    const Point3D kvec(k, 0.0, 0.0);
    const std::vector<std::complex<double>> pvec{   std::complex<double>(0.0, 0.0),
                                                    std::complex<double>(0.0, 0.0),
                                                    std::complex<double>(0.0, 0.0)};
    
    // Create class for assembling Hmatrix and force vector
    fastibem::HAssemblyEmag hassembly(multiforest,
                                      kvec,
                                      pvec,
                                      mu,
                                      omega);
    
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
    
    std::cout << "now computing RCS data...\n";
    
    // Open file for writing RCS data
    std::string rcsfilename("rcs.dat");
    std::ofstream ofs( rcsfilename.c_str() );
    if( !ofs )
        error( "Error opening file" );
    ofs.precision( 18 );
    std::cout.precision( 18 );
    
    const int nseg = 200;                            // number of points for sampling RCS
    const int noutput = 5;
    
    const double delta = nurbs::PI / nseg;          // theta increment
    const double rfar = 1.0e4;                      // distance of far-field points from origin
    
//    const double start = nurbs::PI * 2.0/3.0;
//    const double end = nurbs::PI;
//    const double inc = (end - start)/nseg;
//    std::vector<double> tvals;
//    for(size_t i = 0; i < nseg; ++i)
//        tvals.push_back(start + i * inc);
    
    // loop over RCS points
    for (uint isample = 0; isample < nseg+1; ++isample)
    {
        std::cout << "Computing RCS sample point " << isample + 1 << " / " << nseg + 1 << "....\n";
        
        const double theta = isample * delta;       // current theta

//        const double theta = tvals[isample];
        
        Point3D sample_pt(-rfar * cos(theta),
                          -rfar * sin(theta),
                          0.0);                     // coords of sample point
        Point3D newkvec(k * cos(theta),
                        k * sin(theta),
                        0.0);                       // wave vector
        Point3D rhat(-cos(theta),
                     -sin(theta),
                     0.0);                          // unit vector that defines direction of the plane wave
        
        hassembly.setWaveVector(newkvec);
        
        const std::vector<std::complex<double>> polarvec
        {
            -sin(theta)* std::complex<double>(0.0, 0.0),
            cos(theta) * std::complex<double>(0.0, 0.0),
                         std::complex<double>(1.0, 0.0)
        };

        // is this necessary
        hassembly.setPolarVector(polarvec);
            
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
        auto  A_inv = HLIB::factorise_inv( B.get(), hassembly.trunAccInstance(), & progress );
        
        // GMRES solve
        HLIB::TAutoSolver     solver( 1000 );
        HLIB::TSolver::TInfo  solve_info( false, HLIB::verbose( 2 ) );
        
        auto x = A->col_vector();
        
        timer.start();
        solver.solve( A.get(), x.get(), f.get(), A_inv.get(), &solve_info );
        
        if ( solve_info.converged() ) {
            
            hassembly.clusterTree()->perm_i2e()->permute( x.get() );
            
            // create output vector and set to VTK file
            std::vector<std::complex<double>> solnvec;
            for(size_t i = 0; i < x->size(); ++i)
            {
                const auto entry = x->centry(i);
                solnvec.push_back(std::complex<double>(entry.re(), entry.im()));
            }
            
            // Output solution
            nurbs::OutputVTK output("emagtest_VV_" + std::to_string(isample));
            
            if(isample % (nseg / noutput) == 0)
                output.outputComplexVectorField(multiforest, "surface_current", solnvec);
            
            // Write rcs data
            const double raw_rcs = output.computeRCS(multiforest, sample_pt, k, rhat, mu, omega, solnvec);
            const double normalised_rcs = 10.0 * log10(raw_rcs);
            
            std::cout<< "theta = " << theta << " " << "rcs = "<< normalised_rcs << "\n";
            ofs << 180.0 - (180.0 * theta/PI) << "\t" << normalised_rcs << "\n";
        }
        else
            std::cout << "  not converged in " << timer << " and "
            << solve_info.n_iter() << " steps " << std::endl;
    }
    
    return EXIT_SUCCESS;
}