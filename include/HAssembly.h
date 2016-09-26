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
        virtual std::unique_ptr<HLIB::TMatrix> assembleHmatrix() const = 0;
        
        /// Function that generates RHS matrix
        virtual std::unique_ptr<HLIB::TMatrix> assembleGmatrix() const
        {
            throw std::runtime_error("generateGMatrix() function has not been implemented for this class");
        }
        
        /// Generate the RHS force vector
        virtual void assembleForceVector(HLIB::TVector* fvec) const = 0;
        
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
        
        /// Accessor for quadrature orders
        const std::vector<unsigned>& defaultQuadratureOrder() const { return mQuadratureOrder; }
        
        /// Get the precision specified
        const double precision() const { return mPrecision; }
        
        /// Get the minimum H-matrix block size
        const unsigned minBlockSize() const { return mMinBlockSize; }
        
        /// Get the admissilibity condition
        const double admissibilityCondition() const { return mAdmissibleCondition; }

        
    protected:
        
        /// Construct with multiforest
        template<typename F>
        HAssembly(const F& f,
                  const std::vector<unsigned>& qorder = {2,2},
                  bool cache = false,
                  const double eps = 1.0e-4,
                  const unsigned nmin = 50,
                  const double admissible = 2.0)
        :
        mQuadratureOrder(qorder),
        mUseCache(cache),
        mPrecision(eps),
        mMinBlockSize(nmin),
        mAdmissibleCondition(admissible),
        mTruncAccuracy(eps, 0.0)
        { init(f); }

        /// GEt the set of elements connected to the given global basis index
        const std::vector<unsigned>& connectedEls(const unsigned ibasis) const
        {
            return mConnectedEls[ibasis];
        }
        
    private:
        
        /// Set up the block cluster tree using the mesh
        template<typename F>
        void init(const F& forest)
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
            
            // now set up the set of elements connected to each basis function
            mConnectedEls.clear();
            mConnectedEls.resize(forest.globalDofN());
            
            for(unsigned ielem = 0; ielem < forest.elemN(); ++ielem)
            {
                const auto el = forest.bezierElement(ielem);
                const auto gbasisvec = el->signedGlobalBasisFuncI();
                for(const auto& gindex : gbasisvec)
                {
                    if(-1 == gindex)
                        continue;
                    mConnectedEls[gindex].push_back(ielem);
                }
            }
        }
        
        /// Order of quadrature in each parametric dirction
        std::vector<unsigned> mQuadratureOrder;
        
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
        
        /// The set of element indices connected to each basis function
        std::vector<std::vector<unsigned>> mConnectedEls;
        
    };
    
    /// Create a specialisation for electromagnetic analysis
    class HAssemblyEmag : public HAssembly {
        
    public:
        
        typedef std::complex<double> ReturnType;
        
        typedef std::vector<std::vector<ReturnType>> MatrixType;
        
        /// Construcotr
        HAssemblyEmag(const nurbs::MultiForest& f,
                      const nurbs::Point3D& wvec,
                      const std::vector<std::complex<double>>& pvec,
                      const double mu,
                      const double omega,
                      const std::vector<unsigned>& qorder = {4,4},
                      const std::vector<unsigned>& nsubcells = {1,1},
                      bool cache = false)
        :
        HAssembly(f, qorder, cache),
        mWaveVec(wvec),
        mPolarVec(pvec),
        mMu(mu),
        mOmega(omega),
        mWavenumber(wvec.length()),
        mDefaultSubcellN(nsubcells),
        mMultiForest(f)
        {}
        
        /// Function that generates the LHS matrix
        virtual std::unique_ptr<HLIB::TMatrix> assembleHmatrix() const override;

        
        /// Generate the RHS force vector
        virtual void assembleForceVector(HLIB::TVector* f) const override;
        
        

        
        /// Evaluate the submatrix for the given row and column indices
        MatrixType evalSubmatrix(const std::vector<uint>& rows,
                                 const std::vector<uint>& cols) const;
        
            private:
        
        /// For the given source and field element with an edge
        /// singularity, evaluate the emag kernel and assemble
        /// terms into the given matrix
        void evalEdgeSingularity(const unsigned isrcel,
                                 const unsigned ifieldel,
                                 const nurbs::Edge e1,
                                 const nurbs::Edge e2,
                                 const std::map<int, int>& g2locals,
                                 const std::map<int, int>& g2localf,
                                 MatrixType& mat) const;
        
        /// For the given source and field element with an vertex
        /// singularity, evaluate the emag kernel and assemble
        /// terms into the given matrix
        void evalVertexSingularity(const unsigned isrcel,
                                   const unsigned ifieldel,
                                   const nurbs::Vertex v2,
                                   const std::map<int, int>& g2locals,
                                   const std::map<int, int>& g2localf,
                                   MatrixType& mat) const;
        
        /// For the given source and field element which are
        /// coincident, evaluate the emag kernel and assemble
        /// terms into the given matrix
        void evalCoincidentSingularity(const unsigned iel,
                                       const std::map<int, int>& g2locals,
                                       const std::map<int, int>& g2localf,
                                       MatrixType& mat) const;
        
        /// The wavevector accessor
        const nurbs::Point3D& wavevector() const { return mWaveVec; }
        
        /// Get the wavenumber
        const double wavenumber() const { return mWavenumber; }
        
        /// The polarisation vector accessor
        const std::vector<std::complex<double>>& polarvector() const { return mPolarVec; }
        
        /// Mu material parameter accessor
        const double mu() const { return mMu; }
        
        /// Omega (angular frequency) accessor
        const double omega() const { return mOmega; }
        
        /// Multiforest accessor
        const nurbs::MultiForest& forest() const { return mMultiForest; }
        
        /// Subcell number accessor
        const nurbs::UIntVec& defaultSubcellN() const { return mDefaultSubcellN; }
        
        /// The wavevector
        const nurbs::Point3D mWaveVec;
        
        /// The polarisation vector
        const std::vector<std::complex<double>> mPolarVec;
        
        /// Mu material parameter
        const double mMu;
        
        /// Omega (angular frequency)
        const double mOmega;
        
        /// The wavenumber (i.e. k)
        const double mWavenumber;
        
        /// Default number of subcells used for polar integration
        const nurbs::UIntVec mDefaultSubcellN;

        /// Reference to multiforest
        const nurbs::MultiForest& mMultiForest;
        
        /// Nested class that represents emag Hmatrix coefficients
        class EmagCoeffFn : public HLIB::TPermCoeffFn<HLIB::complex> {
            
        public:
            
            /// Typedef for assembly type
            typedef HAssemblyEmag AssemblyType;
            
            // Constructor
            EmagCoeffFn(const AssemblyType* a,
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
                    for (size_t i = 0; i < n; ++i)
                    {
                        const auto cval = rmat[i][j];
                        matrix[j*n +i] = HLIB::complex(cval.real(),
                                                       cval.imag());
                    }
                
            }
            
            using HLIB::TPermCoeffFn<HLIB::complex>::eval;
            
            /// Collocation BE matrix is non-symmetric
            virtual HLIB::matform_t matrix_format() const override { return HLIB::MATFORM_NONSYM; }
            
            virtual bool is_complex() const override { return true; }
            
        private:
            

            /// assembly instance getter
            const AssemblyType* assemblyInstance() const { return mpAssembly; }
            
            /// Pointer to non-owning assembly instance
            const AssemblyType* mpAssembly;
        };
        
    };

}

#endif