#include "HAssembly.h"

namespace fastibem {
    
    std::unique_ptr<HLIB::TMatrix> HAssemblyEmag::assembleHmatrix() const
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
    
    void HAssemblyEmag::assembleForceVector(HLIB::TVector* f) const
    {
        assert(f->size() == forest().globalDofN());
        clusterTree()->perm_e2i()->permute(f);
    }
    
    HAssemblyEmag::MatrixType HAssemblyEmag::evalSubmatrix(const std::vector<uint>& rows,
                                                           const std::vector<uint>& cols) const
    {
        // Evaluate emag BE terms through Galerkin quadrature routines
        
        MatrixType matrix(rows.size());
        for(size_t i = 0; i < rows.size(); ++i)
            matrix[i].resize(cols.size());
        return matrix;
    }
    
}