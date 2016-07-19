#ifndef FASTIBEM_HASSEMBLY_H
#define FASTIBEM_HASSEMBLY_H

#include <stdexcept>
#include "hlib.hh"

#include "MultiForest.h"
#include "Forest.h"
#include "Kernel.h"

namespace fastibem {
    
    ///
    /// This class encapsulates the assembly of the global matrix and force
    /// vector using H-matrices.
    /// All details of the Hmatrix construction are hidden within this class.
    ///
    class HAssembly {
        
    public:
        
        /// Function that generates the LHS matrix
        virtual std::unique_ptr<HLIB::TMatrix> generateHmatrix() = 0;
        
        /// Function that generates RHS matrix
        virtual std::unique_ptr<HLIB::TMatrix> generateGmatrix()
        {
            throw std::runtime_error("generateGMatrix() function has not been implemented for this class");
        }
        
        /// Generate the RHS force vector
        virtual std::unique_ptr<HLIB::TVector> generateForceVector() = 0;
        
        /// Get the block cluster tree
        const HLIB::TBlockClusterTree* blockClusterTree() const
        { return mBlockClusterTree.get(); }
        
        /// Get the cluster tree
        const HLIB::TClusterTree* clusterTree() const
        { return mClusterTree.get(); }
        
        const HLIB::TCoordinate* coords() const
        { return mCoords.get(); }
        
        const HLIB::TTruncAcc& trunAccInstance() const
        {
            return mTruncAccuracy;
        }
        
        /// Are we caching terms?
        bool usingCache() const { return mUseCache; }
        
        const double precision() const { return mPrecision; }
        
        const unsigned minBlockSize() const { return mMinBlockSize; }
        
        const double admissibilityCondition() const { return mAdmissibleCondition; }
        
        
    protected:
        
        /// Construct with multiforest
        HAssembly(const nurbs::MultiForest& f,
                  bool cache = false,
                  const double eps = 1.0e-4,
                  const unsigned nmin = 50,
                  const double admissible = 2.0)
        :
        mUseCache(cache),
        mPrecision(eps),
        mMinBlockSize(nmin),
        mAdmissibleCondition(admissible),
        mTruncAccuracy(eps, 0.0)
        { initClusterData(f.nodalForest()); }
        
        /// construct with forest
        HAssembly(const nurbs::Forest& f,
                  bool cache = false,
                  const double eps = 1.0e-4,
                  const unsigned nmin = 50,
                  const double admissible = 2.0)
        :
        mUseCache(cache),
        mPrecision(eps),
        mMinBlockSize(nmin),
        mAdmissibleCondition(admissible),
        mTruncAccuracy(eps, 0.0)
        { initClusterData(f); }
        
    private:
        
        /// Set up the block cluster tree using the mesh
        void initClusterData(const nurbs::Forest& f);
        
        /// Cache flag
        bool mUseCache;
        
        /// Precision of Hmatrix
        const double mPrecision;
        
        /// Min number of blocks required during computation
        const unsigned mMinBlockSize;
        
        /// Admissibility condition which defines far-field terms
        const double mAdmissibleCondition;
        
        /// Coordinate data
        std::unique_ptr<HLIB::TCoordinate> mCoords;
        
        /// The block cluster tree
        std::unique_ptr<HLIB::TBlockClusterTree> mBlockClusterTree;
        
        /// The cluster tree
        std::unique_ptr<HLIB::TClusterTree> mClusterTree;
        
        /// Hlib accuracy instance
        HLIB::TTruncAcc mTruncAccuracy;
        
    };
    
    /// Create a specialisation for electromagnetic analysis
    class HAssemblyEmag : public HAssembly {
        
    public:
        
        HAssemblyEmag(const nurbs::MultiForest& f,
                      const nurbs::Point3D& wvec,
                      const nurbs::Point3D& pvec,
                      const double mu,
                      const double omega,
                      bool cache = false)
        :
        HAssembly(f, cache),
        mWaveVec(wvec),
        mPolarVec(pvec),
        mMu(mu),
        mOmega(omega),
        mKernel(wvec.length())
        {}
        
        /// Function that generates the LHS matrix
        virtual std::unique_ptr<HLIB::TMatrix> generateHmatrix() override
        {
            std::unique_ptr<HLIB::TMatrix> hmat;
            return hmat;
        }

        
        /// Generate the RHS force vector
        virtual std::unique_ptr<HLIB::TVector> generateForceVector() override
        {
            std::unique_ptr<HLIB::TVector> fvec;
            
            //clusterTree()->perm_e2i()->permute( fvec.get() );
            return fvec;
        }
        
        
    private:
        
        /// The wavevector accessor
        const nurbs::Point3D& wavevector() const { return mWaveVec; }
        
        /// The polarisation vector accessor
        const nurbs::Point3D& polarvector() const { return mPolarVec; }
        
        /// Mu material parameter accessor
        const double mu() const { return mMu; }
        
        /// Omega (angular frequency) accessor
        const double omega() const { return mOmega; }
        
        /// The wavevector
        const nurbs::Point3D mWaveVec;
        
        /// The polarisation vector
        const nurbs::Point3D mPolarVec;
        
        /// Mu material parameter
        const double mMu;
        
        /// Omega (angular frequency)
        const double mOmega;
        
        /// The kernel instance
        fastibem::EmagKernel mKernel;
    };
    
}

#endif