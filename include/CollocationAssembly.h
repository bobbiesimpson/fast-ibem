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
                            const K& k = K())
        :
        mMesh(f),
        mKernel(k)
        {  if(f != nullptr) init(); }
        
        /// Const forest getter
        const nurbs::Forest* forest() const { return mMesh; }
        
        /// Perform integratino without any caching
        std::vector<std::vector<DataType>> evalHSubmatrix(const std::vector<uint>& rows,
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
        
        /// kernel getter
        const K& kernel() const { return mKernel; }
        
        void clear()
        {
            mConnectedEls.clear();
            mGlobalToLocalCollocMap.clear();
        }
        
    private:
        
        /// Init function for precomputation. Precompute all singular terms.
        void init()
        {
            const nurbs::Forest* f = forest();
            clear();
            
            // Set up map that specifies what elements lie in the span of each basis function
            mConnectedEls.resize(f->globalDofN());
            
            for(uint ielem = 0; ielem < f->elemN(); ++ielem)
            {
                const auto el = f->bezierElement(ielem);
                const auto gbasisvec = el->globalBasisFuncI();
                for(const auto& gindex : gbasisvec)
                    mConnectedEls[gindex].push_back(ielem);
                
                // and insert appropriate entries into the global to local colloc map
                for(uint icolloc = 0; icolloc < el->collocPtN(); ++icolloc)
                    insertGlobalToLocalCollocEntry(el->globalCollocI(icolloc), std::make_pair(ielem, icolloc));
            }
        }
        
        /// Insert an entry into the global to local collocation map.
        /// Will only insert when an entry does not already exist.
        void insertGlobalToLocalCollocEntry(const uint igcolloc,
                                            const std::pair<uint, uint> entry)
        {
            auto find = mGlobalToLocalCollocMap.find(igcolloc);
            if(find == mGlobalToLocalCollocMap.end())
                mGlobalToLocalCollocMap.insert(std::make_pair(igcolloc, entry));
        }
        
        
        /// Get the local collocation index
        std::pair<std::pair<uint, uint>, bool> localCollocEntry(const uint igcolloc) const
        {
            auto search = mGlobalToLocalCollocMap.find(igcolloc);
            if(search != mGlobalToLocalCollocMap.end())
                return std::make_pair(search->second, true);
            else
                return std::make_pair(std::make_pair(nurbs::INVALID_UINT, nurbs::INVALID_UINT), false);
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
        
        /// Non -const forest getter
        const nurbs::Forest* forest() { return mMesh; }
        
        /// Pointer to NURBS mesh
        const nurbs::Forest* mMesh;
        
        /// Map from a global basis index to the connected elements
        std::vector<std::vector<uint>> mConnectedEls;
        
        /// Map from a global collocation index to a connected element and local
        /// collocation index
        std::map<uint, std::pair<uint, uint>> mGlobalToLocalCollocMap;
        
        /// The kernel instance
        const KernelType mKernel;
        
        
    };
   

}
#endif