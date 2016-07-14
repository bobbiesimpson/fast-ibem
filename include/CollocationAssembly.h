#ifndef FASTIBEM_HELMHOLTZ_ASSEMBLY_H
#define FASTIBEM_HELMHOLTZ_ASSEMBLY_H


#include <thread>
#include <vector>
#include <tuple>
#include <set>
#include <cmath>
#include "Geometry.h"
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
                            const unsigned& nt = 1/*std::thread::hardware_concurrency()*/)
        :
        mMesh(f),
        mKernel(k),
        mMaxThreads(nt)
        {  if(f != nullptr) init(); }
        
        /// Const forest getter
        const nurbs::Forest* forest() const { return mMesh; }
        
        /// Perform integratino without any caching
        std::vector<std::vector<DataType>> evalWithoutCache(const std::vector<uint>& rows,
                                                            const std::vector<uint>& cols)
        {
            // Default quadrature order
            const std::vector<unsigned> qorder{2,2};
            
            // Resize return matrix with zeros
            std::vector<std::vector<DataType>> rmat(rows.size(), std::vector<DataType>(cols.size()));
            
            // Vectors of global collocation point and element indices we must compute
            const auto& igcolloc_vec = rows;
            std::vector<uint> igel_vec;
            
            // We must also define a map from a global basis index to the local row index
            // of the returned matrix. Any basis function not present is denoted by -1.
            std::map<uint, int> g2localmap;
            
            // A map which indicates for a given collocation point
            // the set of (element index, local cpt index) pairs which require singular integration
            std::map<unsigned, std::vector<std::pair<unsigned, unsigned>>> singular_cpt_el_map;
            
            // Now fill up the vector of required elements
            for(uint ibasis = 0; ibasis < cols.size(); ++ibasis)
            {
                const uint gb_i = cols[ibasis]; // global basis index
                g2localmap[gb_i] = ibasis;
                
                const auto elvec = connectedEls(gb_i);
                
                // Add all elements in the span of this basis function to the requested list
                for(const auto& e : elvec)
                {
                    igel_vec.push_back(e);
                    const auto el = forest()->bezierElement(e);
                    const auto& gbasisvec = el->globalBasisFuncI();
                    
                    // Set all basis function indices not in the required list equal to -1
                    for(const auto& igbasis : gbasisvec)
                    {
                        if(g2localmap.find(igbasis) == g2localmap.end())
                            g2localmap[igbasis] = -1;
                    }
                }
            }
            
            // now remove duplicate entries
            std::sort(igel_vec.begin(), igel_vec.end());
            auto last = std::unique(igel_vec.begin(), igel_vec.end());
            igel_vec.erase(last, igel_vec.end());
            
            // determine what elements correspond to singular integration for each
            // collocation point
            std::vector<std::vector<unsigned>> cpt_singular_els(igcolloc_vec.size());
            std::vector<std::vector<unsigned>> cpt_local_indices(igcolloc_vec.size());
            
            for(const auto& iel : igel_vec)
            {
                const auto el = forest()->bezierElement(iel);
                const auto& colloc_conn = el->globalCollocConn();
                
                for(size_t i = 0; i < igcolloc_vec.size(); ++i)
                {
                    const auto igcpt = igcolloc_vec[i];
                    auto find = std::find(colloc_conn.begin(), colloc_conn.end(), igcpt);
                    if(find != colloc_conn.end())
                    {
                        const auto lcindex = find - colloc_conn.begin();
                        cpt_singular_els[i].push_back(iel);
                        cpt_local_indices[i].push_back(lcindex);
                    }
                }
            }
            
            // generate a list of the global coordinate points
            std::vector<nurbs::Point3D> cpt_coords(igcolloc_vec.size());
            
            for(unsigned i = 0; i < igcolloc_vec.size(); ++i)
            {
                auto colloc_conn = localCollocEntry(igcolloc_vec[i]);
                if(!colloc_conn.second)
                    throw std::runtime_error("Bad collocation mapping in CollocationAssembly");
                const uint icel = colloc_conn.first.first;
                const uint iclocal = colloc_conn.first.second;
                const auto cel = forest()->bezierElement(icel);
                const nurbs::GPt2D s_parent = cel->collocParentCoord(iclocal);
                cpt_coords[i] = cel->eval(s_parent);
            }
            
            // Let's compute nonsingular integrals!
            for(nurbs::IElemIntegrate igpt(qorder); !igpt.isDone(); ++igpt)
            {
                const nurbs::GPt2D& gpt = igpt.get();
                const auto& w = igpt.getWeight();
                
                for(const auto& iel : igel_vec)
                {
                    const auto el = forest()->bezierElement(iel);
                    const auto gcolloc_conn = el->globalCollocConn();
                    
                    // Get relevant element parameters
                    const nurbs::Point3D& xf = el->eval(gpt);
                    const auto& t1 = el->tangent(gpt, nurbs::ParamDir::S);
                    const auto& t2 = el->tangent(gpt, nurbs::ParamDir::T);
                    const auto m = cross(t1,t2);
                    auto normal = m.asNormal();
                    if(forest()->geometry()->normalsFlipped())
                        normal *= -1.0;
                    
                    const auto jdet = m.length() * el->jacDetParam(gpt.s, gpt.t);
                    
                    const auto& basis = el->basis(gpt.s, gpt.t);
                    const auto& gbasisvec = el->globalBasisFuncI();
                    
                    // loop over requested collocation points
                    for(size_t lcindex = 0; lcindex < igcolloc_vec.size(); ++lcindex)
                    {
                        // if this is a singular integral, continue
                        const auto& singular_els = cpt_singular_els[lcindex];
                        auto find = std::find(singular_els.begin(), singular_els.end(), iel);
                        if(find != singular_els.end())
                            continue;
                        
                        const auto& xs = cpt_coords[lcindex];
                        
                        for(size_t ibasis = 0; ibasis < gbasisvec.size(); ++ibasis)
                        {
                            const auto& ilbasis = g2localmap[gbasisvec[ibasis]];
                            if(-1 != ilbasis)
                                rmat[lcindex][ilbasis] += kernel().evalDLP(xs, xf, normal) * basis[ibasis] * jdet * w;
                        }
                    }
                }
            }
            
            // now perform singular integration
            for(size_t icolloc = 0; icolloc < cpt_singular_els.size(); ++icolloc)
            {
                const auto& elvec = cpt_singular_els[icolloc];
                const auto& ilocal = cpt_local_indices[icolloc];
                
                // loop over singular elements
                for(size_t ielem = 0; ielem < elvec.size(); ++ielem)
                {
                    const auto igel = cpt_singular_els[icolloc][ielem];
                    
                    const auto el = forest()->bezierElement(igel);
                    const auto& gbasisvec = el->globalBasisFuncI();
                    
                    const auto cpt_parent = el->collocParentCoord(ilocal[ielem]);
                    const nurbs::Point3D& xs = cpt_coords[icolloc];
                
                    for(nurbs::IPolarIntegrate igpt(cpt_parent, qorder); !igpt.isDone(); ++igpt)
                    {
                        const nurbs::GPt2D& gpt = igpt.get();
                        const double w = igpt.getWeight();
                        
                        const nurbs::Point3D xf = el->eval(gpt);
                        const auto basis = el->basis(gpt.s, gpt.t);
                        
                        const auto t1 = el->tangent(gpt, nurbs::ParamDir::S);
                        const auto t2 = el->tangent(gpt, nurbs::ParamDir::T);
                        const auto m = cross(t1,t2);
                        
                        auto normal = m.asNormal();
                        if(forest()->geometry()->normalsFlipped())
                            normal *= -1.0;
                        const auto jdet = m.length() * el->jacDetParam(gpt.s, gpt.t);
                        
                        for(size_t ibasis = 0; ibasis < gbasisvec.size(); ++ibasis)
                        {
                            const auto& ilbasis = g2localmap[gbasisvec[ibasis]];
                            if(-1 != ilbasis)
                                rmat[icolloc][ilbasis] += kernel().evalDLP(xs, xf, normal) * basis[ibasis] * jdet * w;
                        }
                    }
                }
                
            }
            
            // and add jump terms using the first element in the singular element lists
            for(size_t icolloc = 0; icolloc < cpt_singular_els.size(); ++icolloc)
            {
                if(cpt_singular_els[icolloc].size() == 0)
                    continue;
                
                const auto iel = cpt_singular_els[icolloc][0];
                const auto ilocal = cpt_local_indices[icolloc][0];
                const double jval = (forest()->geometry()->normalsFlipped()) ? -0.5 : 0.5;
                
                const auto el = forest()->bezierElement(iel);
                const nurbs::GPt2D parent_cpt = el->collocParentCoord(ilocal);
                const auto basis = el->basis(parent_cpt.s, parent_cpt.t);
                const auto gbasis_ivec = el->globalBasisFuncI();
                
                for(uint ibasis = 0; ibasis < basis.size(); ++ibasis)
                {
                    const auto jterm = -jval * basis[ibasis];
                    const auto& ilbasis = g2localmap[gbasis_ivec[ibasis]];
                    if(ilbasis != -1)
                        rmat[icolloc][ilbasis] += jterm;
                }
            }
            return rmat;
        }

        /// Given a set of collocation (row) indices and global basis function (column)
        /// indices, return a matrix of the computed values
        std::vector<std::vector<DataType>> eval(const std::vector<uint>& cindices,
                                                const std::vector<uint>& bindices)
        {
            // First resize return matrix with zeros
            std::vector<std::vector<DataType>> rmat(cindices.size(), std::vector<DataType>(bindices.size()));
            
            // Now fill in cached entries and determine which colloc pt, element index
            // pairings need to be computed.
            std::vector<std::pair<uint, uint>> required_epairs; // (cpt, element) index pairs required
            //std::vector<std::pair<uint, uint>> required_cbpairs;// (cpt, basis) index pairs required
            
            for(uint icolloc = 0; icolloc < cindices.size(); ++icolloc) {
                const uint gc_i = cindices[icolloc]; // global colloc. index
                for(uint ibasis = 0; ibasis < bindices.size(); ++ibasis) {
                    const uint gb_i = bindices[ibasis]; // global basis index
                    const auto p = isCached(gc_i, gb_i);
                    if(!p.second) {
                        //required_cbpairs.push_back(std::make_pair(gc_i, gb_i));
                        const auto elvec = connectedEls(gb_i);
                        for(const auto& e : elvec)
                            required_epairs.push_back(std::make_pair(gc_i, e));
                    }
                }
            }
            // now remove duplicate entries
            std::sort(required_epairs.begin(), required_epairs.end());
            auto last = std::unique(required_epairs.begin(), required_epairs.end());
            required_epairs.erase(last, required_epairs.end());
            
            /// Only compute when we have non-zero entries
            if(required_epairs.size() != 0) {
                
                multithreadNonSingularWorker(std::vector<std::pair<uint,uint>>(required_epairs.begin(), required_epairs.begin() + required_epairs.size()));
                                             
                // Otherwise, now do the actual computation for the required pairings
//                const auto availablethread_n = maxAvailableThreads();
//                const auto bounds = calculateThreadBounds(availablethread_n, required_epairs.size());
//                const auto nthreads = bounds.size() - 1; // the actual number of threads we will use
//                
//                std::vector<std::thread> threads;
//                
//                // create the threads and start the work by calling join() on each
//                auto it = required_epairs.begin();
//                for(uint ithread = 0; ithread < nthreads; ++ithread) {
//                    const uint ilower = bounds[ithread];
//                    const uint iupper = bounds[ithread + 1];
//                    const uint nterms = iupper - ilower;
//                    
//                    threads.push_back(std::thread(&CollocationAssembly::multithreadNonSingularWorker,
//                                                  this,
//                                                  std::vector<std::pair<uint,uint>>(it, it + nterms)));
//                    it += nterms;
//                }
//                
//                for(auto& t : threads)
//                    t.join();
            }
            
            std::lock_guard<std::mutex> lock(mMutex);
            // And finally put the newly created entries into the return matrix
            for(uint irow = 0; irow < cindices.size(); ++irow) {
                const uint igcolloc = cindices[irow];
                for(uint icol = 0; icol < bindices.size(); ++icol) {
                    const uint igbasis = bindices[icol];
                    const auto cachedval = isCached(igcolloc, igbasis);
                    if(!cachedval.second)
                        throw std::runtime_error("Error: Requesting result for a (collocation,basis)"
                                                 " index pairing that was expected to be cached");
                    rmat[irow][icol] = cachedval.first;
                    //eraseCache(igcolloc, igbasis);
                }
            }
            return rmat;
        }
        
        
        /// kernel getter
        const K& kernel() const { return mKernel; }
        
        void clear()
        {
            mCache.clear();
            mConnectedEls.clear();
            mGlobalToLocalCollocMap.clear();
            mTempCache.clear();
            mTempComputedEls.clear();
            mJumpCache.clear();
        }
        
    private:
        

        
        /// Init function for precomputation. Precompute all singular terms.
        void init()
        {
            const nurbs::Forest* f = forest();
            clear();
            
            // Set up map that specifies what elements lie in the span of each basis function
            mConnectedEls.resize(f->globalDofN());
            
            for(uint ielem = 0; ielem < f->elemN(); ++ielem) {
                const auto el = f->bezierElement(ielem);
                const auto gbasisvec = el->globalBasisFuncI();
                for(const auto& gindex : gbasisvec) {
                    mConnectedEls[gindex].push_back(ielem);
                }
                
                // and insert appropriate entries into the global to local colloc map
                for(uint icolloc = 0; icolloc < el->collocPtN(); ++icolloc)
                    insertGlobalToLocalCollocEntry(el->globalCollocI(icolloc), std::make_pair(ielem, icolloc));
            }
            
            
            // compute jump terms and add to cache
            
            const double jval = (f->geometry()->normalsFlipped()) ? -0.5 : 0.5;
            for(uint ielem = 0; ielem < f->elemN(); ++ielem) {
                const auto el = f->bezierElement(ielem);
                for(uint icolloc = 0; icolloc < el->collocPtN(); ++icolloc) {
                    const uint gcolloc_i = el->globalCollocI(icolloc);
                    const nurbs::GPt2D parent_cpt = el->collocParentCoord(icolloc);
                    const auto basis = el->basis(parent_cpt.s, parent_cpt.t);
                    const auto gbasis_ivec = el->globalBasisFuncI();
                    for(uint ibasis = 0; ibasis < basis.size(); ++ibasis) {
                        const uint gbasis_i = gbasis_ivec[ibasis];
                        DataType jterm = -jval * basis[ibasis];
                        auto p = std::make_pair(gcolloc_i, gbasis_i);
                        auto find = mJumpCache.find(p);
                        if(find == mJumpCache.end())
                            mJumpCache[p] = jterm;
                    }
                }
            }
            
            // calculate all singular terms
            const auto availablethread_n = maxAvailableThreads();
            const auto bounds = calculateThreadBounds(availablethread_n, f->elemN());
            const auto nthreads = bounds.size() - 1; // the actual number of threads we will use
            
            std::vector<std::thread> threads;
            for(uint ithread = 0; ithread < nthreads; ++ithread)
                threads.push_back(std::thread(&CollocationAssembly<K>::multithreadSingularWorker,
                                              this,
                                              bounds[ithread],
                                              bounds[ithread+1]));
            
            //std::cout << "Starting singular quadrature computations with " << nthreads << " threads...\n";
            for(auto& t : threads)
                t.join();
            
        }
        
        /// Multithread function for computing singular element contributions.
        /// We specifiy the lower and upper element indices of the range that we
        /// want to compute
        void multithreadSingularWorker(const uint ilower,
                                       const uint iupper)
        {
            std::vector<std::tuple<uint, uint, uint, DataType>> data; // we store the final quadrature result for each cpt and global basis index here
            
            // loop over elements in range
            for(uint ielem = ilower; ielem < iupper; ++ielem) {
                
                const auto el = forest()->bezierElement(ielem); // pointer to element
//                const auto degree_vec = el->degree();
                const std::vector<unsigned> degree_vec{2,2};
                const auto basisvec = el->globalBasisFuncI(); // global basis func. vector
                
                // loop over collocation points that lie within this element
                for(uint icolloc = 0; icolloc < el->collocPtN(); ++icolloc) {
                    
                    std::vector<DataType> qvec(basisvec.size(), 0.0); // local element vector
                    const uint gcolloc_i = el->globalCollocI(icolloc);
                    const auto cpt_parent = el->collocParentCoord(icolloc);
                    const nurbs::Point3D xs = el->eval(cpt_parent);
                    
                    for(nurbs::IPolarIntegrate igpt(cpt_parent, degree_vec); !igpt.isDone(); ++igpt) {
                        const nurbs::GPt2D gpt = igpt.get();
                        const nurbs::Point3D xf = el->eval(gpt);
                        const auto basis = el->basis(gpt.s, gpt.t);
                        
                        const auto t1 = el->tangent(gpt, nurbs::ParamDir::S);
                        const auto t2 = el->tangent(gpt, nurbs::ParamDir::T);
                        const auto m = cross(t1,t2);
                        
                        auto normal = m.asNormal();
                        if(forest()->geometry()->normalsFlipped())
                            normal *= -1.0;
                        const auto jdet = m.length() * el->jacDetParam(gpt.s, gpt.t);
                    
                        for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis) {
                            qvec[ibasis] += kernel().evalDLP(xs, xf, normal) * basis[ibasis] * jdet * igpt.getWeight();
                            
//                            const auto temp = kernel().evalDLP(xs, xf, normal) * basis[ibasis] * el->jacDet(gpt) * igpt.getWeight();
//                            if(std::isnan(temp.real()) || std::isnan(temp.imag()))
//                                std::cout << "bad singular basis\n";
                        }
                    }
                    // now assemble the local vector into the data map
                    const auto cpt_basis = el->basis(cpt_parent.s, cpt_parent.t);
                    for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis) {
                        const auto gbasis_i = basisvec[ibasis];
                        data.emplace_back(std::make_tuple(gcolloc_i, gbasis_i, ielem, qvec[ibasis]));
                    }
                }
            }
            
            // now add the result map for this set of elements to the global cache maps
            std::lock_guard<std::mutex> lock(mMutex);
            for(const auto& d : data) {
                const uint igcolloc = std::get<0>(d);
                const uint igbasis = std::get<1>(d);
                const uint iel = std::get<2>(d);
                const auto val = std::get<3>(d);
                insertCacheVal(igcolloc, igbasis, iel, val);
            }
        }
        
        /// Given a set of (cpt, element) index pairs, compute the BE matrix entries for
        /// each and store in the cache.
        void multithreadNonSingularWorker(const std::vector<std::pair<uint, uint>>& pvec)
        {
            std::vector<std::tuple<uint, uint, uint, DataType>> data; // we put the computed entries here before caching
            
            for(const auto& p : pvec) {
            
                const uint igcolloc = p.first;
                const uint iel = p.second;
                const auto el = forest()->bezierElement(iel);
                const std::vector<unsigned> degree_vec{2,2};//el->degree();
                const auto& gbasis_vec = el->globalBasisFuncI();
                
                // Now, through the collocation map, generate the collocation point physical coordinate
                auto colloc_conn = localCollocEntry(igcolloc);
                if(!colloc_conn.second)
                    throw std::runtime_error("Bad collocation mapping in CollocationAssembly");
                const uint icel = colloc_conn.first.first;
                const uint iclocal = colloc_conn.first.second;
                const auto cel = forest()->bezierElement(icel);
                const nurbs::GPt2D s_parent = cel->collocParentCoord(iclocal);
                const nurbs::Point3D xs = cel->eval(s_parent);
                
                std::vector<DataType> qvec(gbasis_vec.size(), 0.0);
                
                for(nurbs::IElemIntegrate igpt(degree_vec); !igpt.isDone(); ++igpt) {
                    
                    const nurbs::GPt2D& gpt = igpt.get();
                    const nurbs::Point3D& xf = el->eval(gpt);
                    const auto& basis = el->basis(gpt.s, gpt.t);
                    
                    const auto& t1 = el->tangent(gpt, nurbs::ParamDir::S);
                    const auto& t2 = el->tangent(gpt, nurbs::ParamDir::T);
                    const auto m = cross(t1,t2);
                    auto normal = m.asNormal();
                    if(forest()->geometry()->normalsFlipped())
                        normal *= -1.0;
                    const auto jdet = m.length() * el->jacDetParam(gpt.s, gpt.t);
                    
                    for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis)
                        qvec[ibasis] += kernel().evalDLP(xs, xf, normal) * basis[ibasis] * jdet * igpt.getWeight();
                }
                
                // insert the element matrix into the data map
                for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis)
                    data.emplace_back(std::make_tuple(igcolloc, gbasis_vec[ibasis], iel, qvec[ibasis]));
            }
            
            // finally, cache the values using a mutex to prevent thread locking
            std::lock_guard<std::mutex> mutex(mMutex);
            for(const auto& t : data)
            {
                insertCacheVal(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t));
            }

        }
        
        /// Determine whether the value for a collocation index, basis index pairing
        /// has been cached.
        std::pair<DataType, bool> isCached(const uint igcolloc,
                                           const uint igbasis) const
        {
            auto search = mCache.find(std::make_pair(igcolloc, igbasis));
            if(search != mCache.end())
                return std::make_pair(search->second, true);
            else
                return std::make_pair(DataType(), false);
        }
        
        /// Erase cache entry
        void eraseCache(const uint igcolloc,
                        const uint igbasis)
        {
            mCache.erase(std::make_pair(igcolloc, igbasis));
        }
        
        /// Number of connected elements
        uint connnectedElN(const uint igcolloc) const
        {
            return connectedEls(igcolloc).size();
        }
        
        /// Get the elements that cover the span of this basis function.
        std::vector<uint> connectedEls(const uint igbasis) const
        {
            assert(igbasis < mConnectedEls.size());
            return mConnectedEls.at(igbasis);
        }
        
        /// insert an entry into the cache (wrapper function)
        /// returns false if already cached/insertion failed.
        void insertCachedEntry(const std::pair<uint, uint> indices,
                               const DataType val)
        {
            mCache[indices] += val;
        }
        
        /// Insert an entry into the global to local collocation map.
        /// Will only insert when an entry does not already exist.
        void insertGlobalToLocalCollocEntry(const uint igcolloc,
                                            const std::pair<uint, uint> entry)
        {
            auto find = mGlobalToLocalCollocMap.find(igcolloc);
            if(find == mGlobalToLocalCollocMap.end())
                mGlobalToLocalCollocMap.insert(std::make_pair(igcolloc, entry));
//            else
//                std::cout << "Entry already exists for collocation index: "
//                            << igcolloc << ". Carrying on....\n";
        }
        
        std::pair<std::pair<uint, uint>, bool> localCollocEntry(const uint igcolloc) const
        {
            auto search = mGlobalToLocalCollocMap.find(igcolloc);
            if(search != mGlobalToLocalCollocMap.end())
                return std::make_pair(search->second, true);
            else
                return std::make_pair(std::make_pair(nurbs::INVALID_UINT, nurbs::INVALID_UINT), false);
        }
        
        /// for a given global collocation index, global basis function index and global
        /// element index, insert a cached value.
        void insertCacheVal(const uint igcolloc,
                            const uint igbasis,
                            const uint iel,
                            const DataType val)
        {
            /// Do not insert if this entry has been computed already
            auto cached = isCached(igcolloc, igbasis);
            if(cached.second) {
//                std::cout << "The entry for this collocation point index and "
//                          << "basis function index has already been calculated\n";
                return;
            }
            
            // Output error message to check if we are computing terms more than once.
            const auto p = std::make_pair(igcolloc, igbasis);
            auto elvec = mTempComputedEls[p];
            auto find = std::find(elvec.begin(), elvec.end(), iel);
            
            if(find != elvec.end()) {
                //std::cout << "Value already computed for this collocation point, basis index and element index.\n";
                return;
            }
            else {
                
                mTempCache[p] += val;
                mTempComputedEls[p].push_back(iel);
            }
            
            // finally, if all element contributions are in place. Set the value into the final cache
            auto elconn = connectedEls(igbasis);
            
            // first construct a vectof of all the computed element indices
            std::vector<uint> tels;
            for(const auto& iel : mTempComputedEls[p])
                tels.push_back(iel);
            
            // sort this vector and the element connectivity vector
            std::sort(tels.begin(), tels.end());
            std::sort(elconn.begin(), elconn.end());
            
            /// insert final result into cache if all elements computed, erase temp cache
            // and insert jump term if this exists in cache.
            if(tels == elconn) {
                insertCachedEntry(p, mTempCache[p]);
                auto find = mTempComputedEls.find(p);
                mTempComputedEls.erase(find);
                
                auto find2 = mTempCache.find(p);
                mTempCache.erase(find2);
                
                // jump term
                auto jsearch = mJumpCache.find(p);
                if(jsearch != mJumpCache.end()) {
                    insertCachedEntry(p, jsearch->second);
                    mJumpCache.erase(jsearch);
                }

            }
        }

        
        /// Non -const forest getter
        const nurbs::Forest* forest() { return mMesh; }
        
        //// max number of threads available
        const unsigned& maxAvailableThreads() const { return mMaxThreads; }
        
        /// Pointer to NURBS mesh
        const nurbs::Forest* mMesh;
        
        /// Map from (cpt,global basis) index pairing to cached value
        std::map<std::pair<uint, uint>, DataType> mCache;
        
        /// Map from a global basis index to the connected elements
        std::vector<std::vector<uint>> mConnectedEls;
        
        /// Map from a global cpt, global basis index pairing to
        /// the global element indices computed and cached in temp store
        std::map<std::pair<uint,uint>, std::vector<uint>> mTempComputedEls;
        
        /// We store partially computed values for a collocation point,
        /// basis index pairing here
        std::map<std::pair<uint, uint>, DataType> mTempCache;
        
        /// Map from a global collocation index to a connected element and local
        /// collocation index
        std::map<uint, std::pair<uint, uint>> mGlobalToLocalCollocMap;
        
        /// Cache for storing jump terms to be assembled in future.
        std::map<std::pair<uint, uint>, DataType> mJumpCache;
        
        /// The kernel instance
        const KernelType mKernel;
        
        /// Maximum number of threads
        const unsigned mMaxThreads;
        
        /// The mutex for this class
        std::mutex mMutex;
        
        
    };
   

}
#endif