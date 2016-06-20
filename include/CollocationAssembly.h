#ifndef FASTIBEM_HELMHOLTZ_ASSEMBLY_H
#define FASTIBEM_HELMHOLTZ_ASSEMBLY_H


#include <thread>
#include <vector>
#include <tuple>

#include "Forest.h"
#include "Kernel.h"
#include "IPolarIntegrate.h"

namespace fastibem {
    
    /// A handy helper function for calculating thread bounds
    std::vector<long int> calculateThreadBounds(long int parts, long int mem);
    
    /// A class repsonsible for computing the discrete components (row,column) for
    /// collocation BE analysis. The primary use of this class is for computing
    /// H-matrix terms for fast isogeometric BE analysis.
    
    /// Parameterised by the kernel type K
    
    template<typename K>
    class CollocationAssembly {
        
    private:
        
        /// Typedef of the kernel type
        typedef K KernelType;
        
        /// Typedef of the kernel return type
        typedef typename K::ReturnType DataType;

    public:
        
        /// Default constructor
        CollocationAssembly()
        :
        CollocationAssembly(nullptr) {}
        
        /// Construct with given mesh
        CollocationAssembly(const nurbs::Forest* f,
                            const K& k = K(),
                            bool precompute = false)
        :
        mMesh(f),
        mKernel(k) {  if(f != nullptr) init(precompute); }
        
        /// Const forest getter
        const nurbs::Forest* forest() const { return mMesh; }
        
        /// The important evaluation function
        DataType eval(const uint icolloc,
                      const uint ibasis) const
        {
        
            return DataType();
        
        }
        
        /// kernel getter
        const K& kernel() const { return mKernel; }
        
    private:
        
        void clear()
        {
            mCache.clear();
            mConnectedEls.clear();
        }
        
        /// Init function for precomputation.
        void init(bool precompute)
        {
            const nurbs::Forest* f = forest();
            
            clear();
            mConnectedEls.resize(f->globalDofN());
            
            // Set up map that specifies what elements lie in the span of each basis function
            for(uint ielem = 0; ielem < f->elemN(); ++ielem) {
                const auto el = f->element(ielem);
                const auto gbasisvec = el->globalBasisFuncI();
                for(const auto& gindex : gbasisvec)
                    mConnectedEls[gindex].push_back(ielem);
            }
            if(precompute) {
                 // calculate all singular terms
                const auto nthreads = std::thread::hardware_concurrency();
                const auto bounds = calculateThreadBounds(nthreads, forest()->elemN());
                for(const auto& val : bounds)
                    std::cout << val << "\n";
                std::vector<std::thread> threads;
                for(uint ithread = 0; ithread < nthreads; ++ithread)
                    threads.push_back(std::thread(&CollocationAssembly<K>::multithreadSingularWorker,
                                                  this,
                                                  bounds[ithread],
                                                  bounds[ithread+1]));
                for(auto& t : threads)
                    t.join();
            }
            // Now the cache contains all singular terms. Any requests for non-singular
            // terms will compute the relevant non-singular integrals on demand.
        }
        
        /// Multithread function for computing singular element contributions.
        /// We specifiy the lower and upper element indices of the range that we
        /// want to compute
        void multithreadSingularWorker(const uint ilower,
                                       const uint iupper)
        {
            std::map<std::pair<uint, uint>, DataType> result; // Store results for present range
            
            std::vector<std::tuple<uint, uint, DataType>> data;
            
            // loop over elements in range
            for(uint ielem = ilower; ielem < iupper; ++ielem) {

                const auto el = forest()->element(ielem);
                const auto basisvec = el->globalBasisFuncI(); // global basis func. vector
                std::vector<DataType> qvec(basisvec.size(), 0.0); // local element vector
                
                // loop over collocation points that lie within this element
                for(uint icolloc = 0; icolloc < el->collocPtN(); ++icolloc) {
                    const uint gcolloc_i = el->globalCollocI(icolloc);
                    const auto cpt_parent = el->collocParentCoord(icolloc);
                    const nurbs::Point3D xs = el->eval(cpt_parent);
                    for(nurbs::IPolarIntegrate igpt(cpt_parent, el->integrationOrder(4)); !igpt.isDone(); ++igpt) {
                        const nurbs::GPt2D gpt = igpt.get();
                        const nurbs::Point3D xf = el->eval(gpt);
                        const auto basis = el->basis(gpt.s, gpt.t);
                        for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis) {
                            qvec[ibasis] += kernel().eval(xs, xf) * basis[ibasis] * el->jacDet(gpt) * igpt.getWeight();
                        }
                    }
                    for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis) {
                        result.insert(std::make_pair(std::make_pair(gcolloc_i, basisvec[ibasis]), qvec[ibasis]));

                        data.emplace_back(std::make_tuple(gcolloc_i, basisvec[ibasis], qvec[ibasis]));
                    }
                }
            }
            
            // now add the result map for this set of elements to the global cache maps
            std::lock_guard<std::mutex> lock(mMutex);
//            for(const auto& data : result) {
//                std::cout << data.first.first << "," << data.first.second << "\n";
//                auto search = mCache.find(data.first);
//                if(search != mCache.end())
//                    mCache[search->first] += data.second;
//                else
//                    mCache.insert(std::make_pair(data.first, data.second));
//            }
            for(const auto& d : data) {
                const auto p = std::make_pair(std::get<0>(d), std::get<1>(d));
                std::cout << p.first << "," << p.second << "\n";
                auto search = mCache.find(p);
                const auto val = std::get<2>(d);
                if(search != mCache.end())
                    mCache[search->first] += val;
                else
                    mCache.insert(std::make_pair(p,val));
            }
        }
        
        /// Non -const forest getter
        const nurbs::Forest* forest() { return mMesh; }
        
        /// Pointer to NURBS mesh
        const nurbs::Forest* mMesh;
        
        /// Map from (cpt,global basis) index pairing to cached value
        std::map<std::pair<uint, uint>, DataType> mCache;
        
        /// Map from a global basis index to the connected elements
        std::vector<std::vector<uint>> mConnectedEls;
        
        /// The kernel instance
        const KernelType mKernel;
        
        /// The mutex for this class
        std::mutex mMutex;
        
        
    };
   

}
#endif