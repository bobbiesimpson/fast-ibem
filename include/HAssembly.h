#ifndef FASTIBEM_HASSEMBLY_H
#define FASTIBEM_HASSEMBLY_H

#include <stdexcept>
#include "hlib.hh"

#include "MultiForest.h"
#include "BoundingBoxIterator.h"
#include "Forest.h"
#include "Kernel.h"
#include "OutputVTK.h"

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
        template<typename F>
        HAssembly(const F& f,
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
        template<typename F>
        void initClusterData(const F& forest)
        {
            //
            // Bounding box setup
            //
            const uint n = forest.collocPtN();
            std::vector<nurbs::Point3D> vertices;
            std::vector<nurbs::Point3D> bbmin;
            std::vector<nurbs::Point3D> bbmax;
            
            std::vector<double*> p_vertices(n);
            std::vector<double*> p_bbmin(n);
            std::vector<double*> p_bbmax(n);
            
            // Iterate over all bounding boxes
            for(nurbs::BoundingBoxIterator it(forest); !it.isDone(); ++it)
            {
                
                const uint icurrent = it.currentIndex();
                
                // insert point data
                vertices.push_back(it.currentPt());
                const nurbs::Point3D p = it.currentPt();
                p_vertices[icurrent] = new double[3];
                for(uint i = 0; i < 3; ++i)
                    p_vertices[icurrent][i] = p.getCoord(i);
                
                bbmin.push_back(it.currentLowerBound());
                const nurbs::Point3D bbminpt = it.currentLowerBound();
                p_bbmin[icurrent] = new double[3];
                for(uint i = 0; i < 3; ++i)
                    p_bbmin[icurrent][i] = bbminpt.getCoord(i);
                
                bbmax.push_back(it.currentUpperBound());
                const nurbs::Point3D bbmaxpt = it.currentUpperBound();
                p_bbmax[icurrent] = new double[3];
                for(uint i = 0; i < 3; ++i)
                    p_bbmax[icurrent][i] = bbmaxpt.getCoord(i);
            }
            
            // now push data into a vector and create TCoordinate vector
            std::vector<std::pair<nurbs::Point3D, nurbs::Point3D>> bbdata;
            for(uint i = 0; i < n; ++i)
                bbdata.push_back(std::make_pair(bbmin[i], bbmax[i]));
            
            // output bounding box data
            nurbs::OutputVTK output("boundingbox");
            output.outputBoundingBoxSet(bbdata);
            
            mCoords = nurbs::make_unique<HLIB::TCoordinate>(p_vertices, 3, p_bbmin, p_bbmax);
            
            HLIB::TAutoBSPPartStrat  part_strat;
            HLIB::TBSPCTBuilder      ct_builder(&part_strat, minBlockSize());
            mClusterTree           = ct_builder.build(mCoords.get());
            HLIB::TStdGeomAdmCond    adm_cond(admissibilityCondition());
            HLIB::TBCBuilder         bct_builder;
            mBlockClusterTree = bct_builder.build(mClusterTree.get(), mClusterTree.get(), &adm_cond);
            
            // print out cluster visualisations
            if(HLIB::verbose(2))
            {
                HLIB::TPSClusterVis        c_vis;
                HLIB::TPSBlockClusterVis   bc_vis;
                
                c_vis.print( clusterTree()->root(), "multiforest_ct" );
                bc_vis.print( blockClusterTree()->root(), "multiforest_bct" );
            }
        }
        
        
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
        
        typedef std::complex<double> ReturnType;
        
        typedef std::vector<std::vector<ReturnType>> MatrixType;
        
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
        virtual std::unique_ptr<HLIB::TMatrix> generateHmatrix() override;

        
        /// Generate the RHS force vector
        virtual std::unique_ptr<HLIB::TVector> generateForceVector() override;
        
        
    private:
        
        /// Evaluate the submatrix for the given row and column indices
        MatrixType evalSubmatrix(const std::vector<uint>& rows,
                                 const std::vector<uint>& cols);
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
        
        /// Nested class that represents emag Hmatrix coefficients
        class EmagCoeffFn : public HLIB::TPermCoeffFn<HLIB::complex> {
            
        public:
            
            /// Typedef for assembly type
            typedef HAssemblyEmag AssemblyType;
            
            // Constructor
            EmagCoeffFn(AssemblyType* a,
                        const HLIB::TPermutation* row_perm,
                        const HLIB::TPermutation* col_perm)
            :
            TPermCoeffFn<HLIB::complex>(row_perm, col_perm),
            mpAssembly(a) {}
            
            /// Override the evaluation function
            void eval(const std::vector<HLIB::idx_t>& rowidxs,
                      const std::vector<HLIB::idx_t>& colidxs,
                      HLIB::complex* matrix) const override
            {
                const size_t  n = rowidxs.size();
                const size_t  m = colidxs.size();
                
                std::vector<uint> gcolloc_vec;
                for(auto irow = 0; irow < n; ++irow)
                    gcolloc_vec.push_back(rowidxs[irow]);
                    
                std::vector<uint> gbasis_vec;
                for(auto icol = 0; icol < m; ++icol)
                    gbasis_vec.push_back(colidxs[icol]);
                    
                auto rmat = assemblyInstance()->evalSubmatrix(gcolloc_vec, gbasis_vec);
                    
                    for(size_t j = 0; j < m; ++j)
                        for (size_t i = 0; i < n; ++i) {
                            const auto cval = rmat[i][j];
                            matrix[j*n +i] = HLIB::complex(cval.real(),
                                                       cval.imag());
                        }
                
            }
            
            using HLIB::TPermCoeffFn<HLIB::complex>::eval;
            
            /// Collocation BE matrix is non-symmetric
            virtual HLIB::matform_t matrix_format() const override { return HLIB::MATFORM_NONSYM; }
            
        private:
            
            /// assembly instance getter
            AssemblyType* assemblyInstance() const { return mpAssembly; }
            
            /// Pointer to non-owning assembly instance
            AssemblyType* mpAssembly;
        };
        
    };
    
}

#endif