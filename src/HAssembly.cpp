#include "HAssembly.h"

namespace fastibem {
    
    std::unique_ptr<HLIB::TMatrix> HAssemblyEmag::generateHmatrix()
    {
        std::cout << "━━ building H-matrix ( eps = " << precision() << " )" << std::endl;
        
        HLIB::TTimer timer( HLIB::WALL_TIME );
        HLIB::TConsoleProgressBar progress;
        
        EmagCoeffFn coefffn(this, clusterTree()->perm_i2e(), clusterTree()->perm_i2e());
        
        HLIB::TACAPlus<HLIB::complex> aca(&coefffn);
        HLIB::TDenseMBuilder<HLIB::complex> h_builder(&coefffn, &aca);
        HLIB::TPSMatrixVis mvis;
        h_builder.set_coarsening(false);
        
        timer.start();
        auto A = h_builder.build(blockClusterTree(), trunAccInstance(), &progress);
        timer.pause();
        
        return A;
    }
    
    std::unique_ptr<HLIB::TVector> HAssemblyEmag::generateForceVector()
    {
        std::unique_ptr<HLIB::TVector> fvec;
        
        //clusterTree()->perm_e2i()->permute( fvec.get() );
        return fvec;
    }
    
    HAssemblyEmag::MatrixType HAssemblyEmag::evalSubmatrix(const std::vector<uint>& rows,
                                                           const std::vector<uint>& cols)
    {
        // Evaluate emag BE terms through Galerkin quadrature routines
        git 
        
        
        
        
        
        
        
        return MatrixType();
    }
    
}