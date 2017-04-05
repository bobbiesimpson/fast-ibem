#include <complex>

#include "HAssembly.h"
#include "AnalysisElement.h"
#include "IGalerkinIntegrate.h"
#include "IEdgeQuadrature.h"
#include "IVertexQuadrature.h"
#include "IEqualQuadrature.h"
#include "IRegularQuadrature.h"
#include "IPolarIntegrate.h"
#include "Functor.h"
#include "base.h"
#include "ISubElemIntegrate.h"
#include "NURBSCache.h"
#include "IEqualQuadratureTri.h"
#include "IPolarDegenerate.h"

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
        
        // temp force vector
        std::vector<std::complex<double>> ftemp(f->size());
        
        EmagPlaneWave pw_functor(wavevector(), polarvector());
        
        uint offset = 2;
        
        for(unsigned ielem = 0; ielem < forest().elemN(); ++ielem)
        {
            
            const auto p_el = forest().bezierElement(ielem);
            
            const auto econn = p_el->signedGlobalBasisFuncI();
            const auto& qorder = p_el->equalIntegrationOrder(offset);
            
            for(nurbs::IElemIntegrate igpt(qorder); !igpt.isDone(); ++igpt)
            {
                const auto gpt = igpt.get();
                const auto w = igpt.getWeight();
                
                // physical coordinate
                const auto x = p_el->eval(gpt);
                
                // emag plane wave at x
                const auto pw = pw_functor(x);
                
                // tangent vectors
                const auto& t1 = p_el->tangent(gpt.s, gpt.t, nurbs::ParamDir::S);
                const auto& t2 = p_el->tangent(gpt.s, gpt.t, nurbs::ParamDir::T);
                
                // basis funcs
                const auto basis = p_el->basis(gpt.s, gpt.t, t1, t2);
                const auto jdet = p_el->jacDet(gpt);
                
                for(size_t ibasis = 0; ibasis < econn.size(); ++ibasis)
                {
                    // global basis index
                    const auto igbasis = econn[ibasis];
                    if(-1 == igbasis)
                        continue; // degenerate dof
                    
                    for(unsigned j = 0; j < 3; ++j)
                        ftemp[igbasis] += 1.0 / (std::complex<double>(0.0, mu() * omega())) * basis[ibasis][j] * pw[j] * w * jdet;
                    //                        ftemp[igbasis] -= 1.0 / (std::complex<double>(0.0, mu() * omega()))  * basis[ibasis][j] * pw[j] * w * jdet;
                }
            }
        }
        
        for(size_t i = 0; i < ftemp.size(); ++i)
            f->set_centry(i, HLIB::complex(ftemp[i].real(), ftemp[i].imag()));
        
        clusterTree()->perm_e2i()->permute(f);
    }
    
    HAssemblyEmag::MatrixType HAssemblyEmag::evalSubmatrix(const std::vector<uint>& rows,
                                                           const std::vector<uint>& cols) const
    {
        
        // First resize matrix
        MatrixType matrix(rows.size());
        for(size_t i = 0; i < rows.size(); ++i)
            matrix[i].resize(cols.size());
        
        const auto k = wavenumber();
        const double minDistRatio = 2.0;                // the minimum d / h ratio for adaptive quadrature
        
        const nurbs::UIntVec regularorder{3,3};         // default quadrature order for regular integrals
        const nurbs::UIntVec degen_order{4,4};
        
        const nurbs::UIntVec nearsingular_sorder{4,4};
        const nurbs::UIntVec nearsingular_forder{4,4};
        
        // Vectors of source (test) and field (trial) elements we must compute
        std::vector<unsigned> isrc_els;
        std::vector<unsigned> ifield_els;
        
        // maps from global basis indices to local row or column indices.
        // a local index of -1 denotes that the term is ignored.
        std::map<int, int> g2local_field;
        std::map<int, int> g2local_src;
        
        // first fill up the source element vector and map
        for(size_t ibasis = 0; ibasis < rows.size(); ++ibasis)
        {
            const uint igbasis = rows[ibasis]; // global basis index
            g2local_src[igbasis] = ibasis;
            
            // Add all elements in the span of this basis function to the requested list
            for(const auto& e : connectedEls(igbasis))
            {
                isrc_els.push_back(e);
                const auto el = forest().bezierElement(e);
                
                // Set all basis function indices not in the required list equal to -1
                for(const auto& i : el->signedGlobalBasisFuncI())
                {
                    if(-1 == i) // degenerate point
                        continue;
                    
                    if(g2local_src.find(i) == g2local_src.end())
                        g2local_src[i] = -1;
                }
            }
        }
        
        // and the same for field elements
        for(size_t ibasis = 0; ibasis < cols.size(); ++ibasis)
        {
            const uint igbasis = cols[ibasis]; // global basis index
            g2local_field[igbasis] = ibasis;
            
            // Add all elements in the span of this basis function to the requested list
            for(const auto& e : connectedEls(igbasis))
            {
                ifield_els.push_back(e);
                const auto el = forest().bezierElement(e);
                
                // Set all basis function indices not in the required list equal to -1
                for(const auto& i : el->signedGlobalBasisFuncI())
                {
                    if(-1 == i)
                        continue;
                    
                    if(g2local_field.find(i) == g2local_field.end())
                        g2local_field[i] = -1;
                }
            }
        }
        
        // remove any duplicate entries
        std::sort(isrc_els.begin(), isrc_els.end());
        auto last_src = std::unique(isrc_els.begin(), isrc_els.end());
        isrc_els.erase(last_src, isrc_els.end());
        
        std::sort(ifield_els.begin(), ifield_els.end());
        auto last_field = std::unique(ifield_els.begin(), ifield_els.end());
        ifield_els.erase(last_field, ifield_els.end());
        
        // now compute the sets of coincident, edge adjacent
        // and vertex adjacent integrals along with the required edges/vertices.
        std::vector<std::tuple<unsigned, unsigned, nurbs::Vertex, nurbs::Vertex>> vertex_integrals;
        std::vector<std::tuple<unsigned, unsigned, nurbs::Edge, nurbs::Edge>> edge_integrals;
        std::vector<std::tuple<unsigned, unsigned, nurbs::Edge, nurbs::Edge>> degen_edge_integrals;
        std::vector<unsigned> coincident_integrals;
        std::vector<std::tuple<unsigned, unsigned>> degen_vertex_integrals;
        
        // The sets of regular source and field elements
        std::map<unsigned, std::vector<unsigned>> regular_elmap;
        std::map<unsigned, std::vector<unsigned>> regular_degen_elmap;
        std::map<unsigned, std::vector<unsigned>> nearsingular_elmap;
        
        nurbs::Edge e1, e2;
        nurbs::Vertex v1, v2;
        
        // Now loop over source and field elements and determine singularity cases
        // and thus fill up the vectors vertex_integrals, edge_integrals etc.
        for(size_t isrcel = 0; isrcel < isrc_els.size(); ++isrcel)
        {
            const auto igsrcel = isrc_els[isrcel];
            const auto p_srcel = forest().bezierElement(igsrcel);
            
            for(size_t ifieldel = 0; ifieldel < ifield_els.size(); ++ifieldel)
            {
                const auto igfieldel = ifield_els[ifieldel];
                const auto p_fieldel = forest().bezierElement(igfieldel);
                
                if(igsrcel == igfieldel)
                    coincident_integrals.push_back(igsrcel);
                
                // edge singularity
                else if(nurbs::edgeConnected(*p_srcel, *p_fieldel, e1, e2))
                {
                    // seperate degenerate and non-degenerate cases
                    /*if(p_srcel->degenerate() && p_fieldel->degenerate())
                        degen_edge_integrals.push_back(std::make_tuple(igsrcel, igfieldel, e1, e2));
                    else*/
                        edge_integrals.push_back(std::make_tuple(igsrcel, igfieldel, e1, e2));
                }
//                else if(nurbs::connectedAtDegeneratePt(*p_srcel, *p_fieldel))
//                    degen_vertex_integrals.push_back(std::make_tuple(igsrcel, igfieldel));

                // vertex singularity
//                else if(nurbs::vertexConnected(*p_srcel, *p_fieldel, v1, v2))
//                    vertex_integrals.push_back(std::make_tuple(igsrcel, igfieldel, v1, v2));
                
                // increase quadrature for all degenerate, non-singular elements
//                else if(p_fieldel->degenerate())
//                    degen_vertex_integrals.push_back(std::make_tuple(igsrcel, igfieldel));
                
                // nearly singular integral
//                else if(nurbs::dist(p_srcel->eval(0.0, 0.0), p_fieldel->eval(0.0, 0.0)) < p_fieldel->size() * minDistRatio)
//                {
//                    /*if(p_srcel->degenerate() || p_fieldel->degenerate())
//                        degen_vertex_integrals.push_back(std::make_tuple(igsrcel, igfieldel));
//                    else*/
//                        nearsingular_elmap[igsrcel].push_back(igfieldel);
//                }
                
                // regular integral
                else
                {
//                    if(p_srcel->degenerate() || p_fieldel->degenerate())
//                        regular_degen_elmap[igsrcel].push_back(igfieldel);
//                    else
                        regular_elmap[igsrcel].push_back(igfieldel);
                }
            }
        }
        
        // evaluate edge singularity elements
        for(const auto& tuple : edge_integrals)
            evalEdgeSingularity(std::get<0>(tuple),
                                std::get<1>(tuple),
                                std::get<2>(tuple),
                                std::get<3>(tuple),
                                g2local_src,
                                g2local_field,
                                matrix);
        
//        // evaluate vertex elements
//        for(const auto& tuple : vertex_integrals)
//            evalVertexSingularity(std::get<0>(tuple),
//                                  std::get<1>(tuple),
//                                  std::get<2>(tuple),
//                                  std::get<3>(tuple),
//                                  g2local_src,
//                                  g2local_field,
//                                  matrix);
        
//        // eval vertex connected degenerate elements
//        for(const auto& tuple: degen_vertex_integrals)
//            evalVertexDegenerateSingularity(std::get<0>(tuple),
//                                            std::get<1>(tuple),
//                                            g2local_src,
//                                            g2local_field,
//                                            matrix);
        
        // evaluate coincidnet elements
        for(const auto& iel : coincident_integrals)
            evalCoincidentSingularity(iel,
                                      g2local_src,
                                      g2local_field,
                                      matrix);
        
        
        // cache source element terms
        std::map<std::pair<unsigned, unsigned>, nurbs::DoubleVecVec> sbasis_map;
        std::map<std::pair<unsigned, unsigned>, double> sjdet_map;
        std::map<std::pair<unsigned, unsigned>, nurbs::DoubleVec> sdiv_map;
        std::map<std::pair<unsigned, unsigned>, nurbs::Point3D> spoint_map;
        std::map<unsigned, nurbs::IntVec> sconn_map;
        
        for(const auto& isel : isrc_els)
        {
            const auto p_sel = forest().bezierElement(isel);
            
            const auto& sorder = regularorder;
            
            // compute set of basis func indices required
            nurbs::IntVec sset;
            const auto& conn = p_sel->signedGlobalBasisFuncI();
            for(uint ibasis = 0; ibasis < p_sel->basisFuncN(); ++ibasis)
            {
                if(g2local_src[conn[ibasis]] != -1)
                    sset.push_back(ibasis);
            }
            sconn_map[isel] = sset;
            
            // loop over source element quadrature points
            for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
            {
                const auto sparent = igpt_s.get();
                const auto index = igpt_s.currentIndex();
                const auto pair = std::make_pair(isel, index);
                
                // source element parameters
                spoint_map[pair] = p_sel->eval(sparent);
                
                // tangent vectors
                const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
                
                const double jpiola_s = nurbs::cross(t1, t2).length();
                
                // physical coordinate
                sbasis_map[pair] = p_sel->basis(sparent.s, sparent.t, t1, t2);
                const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
                const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
                
                // divergence operator
                nurbs::DoubleVec div_vec(p_sel->basisFuncN(), 0.0);
                
                for(unsigned ibasis = 0; ibasis < p_sel->basisFuncN(); ++ibasis)
                    div_vec[ibasis] = 1./jpiola_s * (ds_s[ibasis][0] + dt_s[ibasis][1]);
                sdiv_map[pair] = div_vec;
                
                // jacobian determinant
                sjdet_map[pair] = p_sel->jacDet(sparent.s, sparent.t, t1, t2);
            }
        }
        
        // Now perform computation of regular elements
        for(const auto& pair: regular_elmap)
        {
            const auto& igsrcel = pair.first;
            const auto p_sel = forest().bezierElement(igsrcel);
            const auto& sconn = p_sel->signedGlobalBasisFuncI();
            
            const auto sorder = regularorder;
            
            // loop over source element quadrature points
            for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
            {
                const auto sw = igpt_s.getWeight();
                
                const auto index = igpt_s.currentIndex();
                const auto cache_pair = std::make_pair(igsrcel, index);
                
                // fetch cached terms
                const auto& basis_s = sbasis_map[cache_pair];
                const auto& jdet_s = sjdet_map[cache_pair];
                const auto& divvec_s = sdiv_map[cache_pair];
                const auto& x = spoint_map[cache_pair];

                // loop over field elements
                for(const auto& igfieldel : pair.second)
                {
                    const auto p_fel = forest().bezierElement(igfieldel);
                    const auto& fconn = p_fel->signedGlobalBasisFuncI();
                    
                    const auto& forder = regularorder;
                    
                    // integrate over field elements
                    for(nurbs::IElemIntegrate igpt_f(forder); !igpt_f.isDone(); ++igpt_f)
                    {
                        const auto fparent = igpt_f.get();
                        const auto fw = igpt_f.getWeight();
                        
                        const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                        const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                        const auto y = p_fel->eval(fparent);
                        const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
                        const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                        const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                        const double jdet_f = p_fel->jacDet(fparent, t1, t2);
                        const double jpiola_f = nurbs::cross(t1, t2).length();
                        
                        // kernel
                        const double r = dist(x, y);
                        const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                        
                        // now loop over test and trial functions
                        for(size_t itest = 0; itest < sconn.size(); ++itest)
                        {
                            const auto igbasis_s = sconn[itest];
                            if(igbasis_s == -1)
                                continue;
                            
                            const auto ilocal_s = g2local_src[igbasis_s];
                            if(ilocal_s == -1)
                                continue;
                            
                            // divergence (source)
//                            const double div_s = 1./jpiola_s *(ds_s[itest][0] + dt_s[itest][1]);
                            const double div_s = divvec_s[itest];
                            
                            for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                            {
                                const auto igbasis_f = fconn[itrial];
                                if(igbasis_f == -1)
                                    continue;
                                
                                const auto ilocal_f = g2local_field[igbasis_f];
                                if(ilocal_f == -1)
                                    continue;
                                
                                // divergence (field)
                                const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                                
                                for(unsigned i = 0; i < 3; ++i)
                                    matrix[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                                matrix[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * sw * fw;
                                
                            }
                        }
                    }
                }
            }
        }
        
//        // clear cache
//        sbasis_map.clear();
//        sjdet_map.clear();
//        sdiv_map.clear();
//        spoint_map.clear();
//        
//        // Compute source terms for degnerate elements
//        for(const auto& isel : isrc_els)
//        {
//            const auto p_sel = forest().bezierElement(isel);
//            
//            const auto& sorder = degen_order;
//            
//            // loop over source element quadrature points
//            for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
//            {
//                const auto sparent = igpt_s.get();
//                const auto index = igpt_s.currentIndex();
//                const auto pair = std::make_pair(isel, index);
//                
//                // source element parameters
//                spoint_map[pair] = p_sel->eval(sparent);
//                
//                // tangent vectors
//                const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
//                const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
//                
//                const double jpiola_s = nurbs::cross(t1, t2).length();
//                
//                // physical coordinate
//                sbasis_map[pair] = p_sel->basis(sparent.s, sparent.t, t1, t2);
//                const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
//                const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
//                
//                // divergence operator
//                nurbs::DoubleVec div_vec(p_sel->basisFuncN(), 0.0);
//                
//                for(unsigned ibasis = 0; ibasis < p_sel->basisFuncN(); ++ibasis)
//                    div_vec[ibasis] = 1./jpiola_s * (ds_s[ibasis][0] + dt_s[ibasis][1]);
//                sdiv_map[pair] = div_vec;
//                
//                // jacobian determinant
//                sjdet_map[pair] = p_sel->jacDet(sparent.s, sparent.t, t1, t2);
//            }
//        }
//        
//        // Perform computation of degenerate, regular elements
//        for(const auto& pair: regular_degen_elmap)
//        {
//            const auto& igsrcel = pair.first;
//            const auto p_sel = forest().bezierElement(igsrcel);
//            const auto& sconn = p_sel->signedGlobalBasisFuncI();
//            
//            const auto sorder = degen_order;
//            
//            // loop over source element quadrature points
//            for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
//            {
//                const auto sw = igpt_s.getWeight();
//                
//                const auto index = igpt_s.currentIndex();
//                const auto cache_pair = std::make_pair(igsrcel, index);
//                
//                // fetch cached terms
//                const auto& basis_s = sbasis_map[cache_pair];
//                const auto& jdet_s = sjdet_map[cache_pair];
//                const auto& divvec_s = sdiv_map[cache_pair];
//                const auto& x = spoint_map[cache_pair];
//                
//                // loop over field elements
//                for(const auto& igfieldel : pair.second)
//                {
//                    const auto p_fel = forest().bezierElement(igfieldel);
//                    const auto& fconn = p_fel->signedGlobalBasisFuncI();
//                    
//                    const auto& forder = degen_order;
//                    
//                    // integrate over field elements
//                    for(nurbs::IElemIntegrate igpt_f(forder); !igpt_f.isDone(); ++igpt_f)
//                    {
//                        const auto fparent = igpt_f.get();
//                        const auto fw = igpt_f.getWeight();
//                        
//                        const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
//                        const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
//                        const auto y = p_fel->eval(fparent);
//                        const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
//                        const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
//                        const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
//                        const double jdet_f = p_fel->jacDet(fparent, t1, t2);
//                        const double jpiola_f = nurbs::cross(t1, t2).length();
//                        
//                        // kernel
//                        const double r = dist(x, y);
//                        const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
//                        
//                        // now loop over test and trial functions
//                        for(size_t itest = 0; itest < sconn.size(); ++itest)
//                        {
//                            const auto igbasis_s = sconn[itest];
//                            if(igbasis_s == -1)
//                                continue;
//                            
//                            const auto ilocal_s = g2local_src[igbasis_s];
//                            if(ilocal_s == -1)
//                                continue;
//                            
//                            // divergence (source)
//                            const double div_s = divvec_s[itest];
//                            
//                            for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
//                            {
//                                const auto igbasis_f = fconn[itrial];
//                                if(igbasis_f == -1)
//                                    continue;
//                                
//                                const auto ilocal_f = g2local_field[igbasis_f];
//                                if(ilocal_f == -1)
//                                    continue;
//                                
//                                // divergence (field)
//                                const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
//                                
//                                for(unsigned i = 0; i < 3; ++i)
//                                    matrix[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
//                                matrix[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * sw * fw;
//                                
//                            }
//                        }
//                    }
//                }
//            }
//        }
        
        return matrix;
    }
    
    void HAssemblyEmag::evalVertexDegenerateSingularity(const unsigned isrcel,
                                                        const unsigned ifieldel,
                                                        const std::map<int, int>& g2locals,
                                                        const std::map<int, int>& g2localf,
                                                        MatrixType& mat) const
    {
//        std::cout << "Evaluating vertex degenerate integral for source element: " << isrcel << " field element: " << ifieldel << "\n";
        
        const nurbs::UIntVec sorder{8,8};
        const nurbs::UIntVec forder{4,4};
        
        const auto p_sel = forest().bezierElement(isrcel);
        const auto p_fel = forest().bezierElement(ifieldel);
        
        const auto& sconn = p_sel->signedGlobalBasisFuncI();
        const auto& fconn = p_fel->signedGlobalBasisFuncI();
        
        const auto& k = wavenumber();
        
        // loop over source element quadrature points
        for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        {
            const auto sparent = igpt_s.get();
            const auto sw = igpt_s.getWeight();
            
            // source element parameters
            const auto x = p_sel->eval(sparent);
            
            const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1, t2);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const double jdet_s = p_sel->jacDet(sparent.s, sparent.t, t1, t2);
            const double jpiola_s = nurbs::cross(t1, t2).length();
            
            
            // integrate over field elements
            for(nurbs::IElemIntegrate igpt_f(forder); !igpt_f.isDone(); ++igpt_f)
            {
                const auto fparent = igpt_f.get();
                const auto fw = igpt_f.getWeight();
                
                const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                const auto y = p_fel->eval(fparent);
                const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
                const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                const double jdet_f = p_fel->jacDet(fparent, t1, t2);
                const double jpiola_f = nurbs::cross(t1, t2).length();
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < sconn.size(); ++itest)
                {
                    const auto igbasis_s = sconn[itest];
                    if(igbasis_s == -1)
                        continue;
                    
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = 1./jpiola_s *(ds_s[itest][0] + dt_s[itest][1]);
                    //                            const double div_s = divvec_s[itest];
                    
                    for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                    {
                        const auto igbasis_f = fconn[itrial];
                        if(igbasis_f == -1)
                            continue;
                        
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * sw * fw;
                        
                    }
                }
            }
            
        }
    }
    
    void HAssemblyEmag::evalEdgeSingularity(const unsigned isrcel,
                                            const unsigned ifieldel,
                                            const nurbs::Edge e1,
                                            const nurbs::Edge e2,
                                            const std::map<int, int>& g2locals,
                                            const std::map<int, int>& g2localf,
                                            MatrixType& mat) const
    {
        // Sauter and Schwab edge singularity quadrature
        
        // Quadrature orders determined from testing

        // Pointers to elements
        const auto p_sel = forest().bezierElement(isrcel);
        const auto p_fel = forest().bezierElement(ifieldel);
        
        nurbs::UIntVec sorder{5,5};
        nurbs::UIntVec forder{4,4};
        
        // Connectivities
        const auto& sconn = p_sel->signedGlobalBasisFuncI();
        const auto& fconn = p_fel->signedGlobalBasisFuncI();
        
        const auto& k = wavenumber();
        
        for(nurbs::IEdgeQuadrature igpt(sorder, forder, e1, e2); !igpt.isDone(); ++igpt)
        {
            const auto gpt4d = igpt.get();
            const auto sparent = gpt4d.srcPt();
            const auto fparent = gpt4d.fieldPt();
            
            const auto w = igpt.getWeight();
            
            // source element parameters
            const auto& t1_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            const auto& x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1_s, t2_s);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const auto& jdet_s = p_sel->jacDet(sparent, t1_s, t2_s);
            const auto& jpiola_s = nurbs::cross(t1_s, t2_s).length();
            
            // field point terms
            const auto& t1_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
            const auto& t2_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
            const auto y = p_fel->eval(fparent);
            const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1_f, t2_f);
            const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
            const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
            const double jdet_f = p_fel->jacDet(fparent, t1_f, t2_f);
            const double jpiola_f = nurbs::cross(t1_f, t2_f).length();
            
            // kernel
            const double r = dist(x, y);
            const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
            
            // now loop over test and trial functions
            for(size_t itest = 0; itest < sconn.size(); ++itest)
            {
                const auto igbasis_s = sconn[itest];
                if(igbasis_s == -1)
                    continue;
                
                const auto ilocal_s = g2locals.at(igbasis_s);
                if(ilocal_s == -1)
                    continue;
                
                // divergence (source)
                const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                
                for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                {
                    const auto igbasis_f = fconn[itrial];
                    if(igbasis_f == -1)
                        continue;
                    
                    const auto ilocal_f = g2localf.at(igbasis_f);
                    if(ilocal_f == -1)
                        continue;
                    
                    // divergence (field)
                    const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                    
                    for(unsigned i = 0; i < 3; ++i)
                        mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * w;
                    mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * w;
                    
                }
            }
        }


    }
    
    void HAssemblyEmag::evalEdgeSingularityPolar(const unsigned isrcel,
                                                 const unsigned ifieldel,
                                                 const nurbs::Edge e1,
                                                 const nurbs::Edge e2,
                                                 const std::map<int, int>& g2locals,
                                                 const std::map<int, int>& g2localf,
                                                 MatrixType& mat) const
    {
        const double k = wavenumber();
        
        // source element and connectivity
        const auto p_sel = forest().bezierElement(isrcel);
        const auto& sconn = p_sel->signedGlobalBasisFuncI();
        
        // quadrature orders
        const nurbs::UIntVec sorder{5,5};
        const nurbs::UIntVec forder{3,3};
        
        // field element and connectivity
        const auto p_fel = forest().bezierElement(ifieldel);
        const auto& fconn = p_fel->signedGlobalBasisFuncI();
        
        const auto degenerate_edgepair = p_fel->degenerateEdge();

        for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        {
            const auto sparent = igpt_s.get();
            const auto sw = igpt_s.getWeight();
            
            // cached tangent entries
            const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            
            // source element parameters
            const auto x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1, t2);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const double jdet_s = p_sel->jacDet(sparent, t1, t2);
            const double jpiola_s = nurbs::cross(t1, t2).length();
            
            // integrate over field elements
            for(nurbs::IPolarDegenerate igpt_f(nurbs::projectPt(sparent, e1, e2), degenerate_edgepair.second, forder); !igpt_f.isDone(); ++igpt_f)
            {
                const auto fparent = igpt_f.get();
                const auto fw = igpt_f.getWeight();
                
                // cached tangent entries
                const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                
                // field element parameters
                const auto y = p_fel->eval(fparent);
                const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
                const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                const double jdet_f = p_fel->jacDet(fparent, t1, t2);
                const double jpiola_f = nurbs::cross(t1, t2).length();
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < sconn.size(); ++itest)
                {
                    const auto igbasis_s = sconn[itest];
                    if(igbasis_s == -1)
                        continue;
                    
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                    
                    for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                    {
                        const auto igbasis_f = fconn[itrial];
                        if(igbasis_f == -1)
                            continue;
                        
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * sw * fw;
                        
                    }
                }
            }
        }
    }
    
    void HAssemblyEmag::evalVertexSingularity(const unsigned isrcel,
                                              const unsigned ifieldel,
                                              const nurbs::Vertex v1,
                                              const nurbs::Vertex v2,
                                              const std::map<int, int>& g2locals,
                                              const std::map<int, int>& g2localf,
                                              MatrixType& mat) const
    {
        // Sauter and Schwab vertex singularity quadrature
        
        // Quadrature orders determined from testing
        const nurbs::UIntVec sorder{6,6};
        const nurbs::UIntVec forder{3,3};
        
        // Pointers to elements
        const auto p_sel = forest().bezierElement(isrcel);
        const auto p_fel = forest().bezierElement(ifieldel);
        
        // Connectivities
        const auto& sconn = p_sel->signedGlobalBasisFuncI();
        const auto& fconn = p_fel->signedGlobalBasisFuncI();
        
        const auto& k = wavenumber();
        
        for(nurbs::IVertexQuadrature igpt(sorder, forder, v1, v2); !igpt.isDone(); ++igpt)
        {
            const auto gpt4d = igpt.get();
            const auto sparent = gpt4d.srcPt();
            const auto fparent = gpt4d.fieldPt();
            
            const auto w = igpt.getWeight();
            
            // source element parameters
            const auto& t1_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            const auto& x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1_s, t2_s);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const auto& jdet_s = p_sel->jacDet(sparent, t1_s, t2_s);
            const auto& jpiola_s = nurbs::cross(t1_s, t2_s).length();
            
            // field point terms
            const auto& t1_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
            const auto& t2_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
            const auto y = p_fel->eval(fparent);
            const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1_f, t2_f);
            const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
            const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
            const double jdet_f = p_fel->jacDet(fparent, t1_f, t2_f);
            const double jpiola_f = nurbs::cross(t1_f, t2_f).length();
            
            // kernel
            const double r = dist(x, y);
            const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
            
            // now loop over test and trial functions
            for(size_t itest = 0; itest < sconn.size(); ++itest)
            {
                const auto igbasis_s = sconn[itest];
                if(igbasis_s == -1)
                    continue;
                
                const auto ilocal_s = g2locals.at(igbasis_s);
                if(ilocal_s == -1)
                    continue;
                
                // divergence (source)
                const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                
                for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                {
                    const auto igbasis_f = fconn[itrial];
                    if(igbasis_f == -1)
                        continue;
                    
                    const auto ilocal_f = g2localf.at(igbasis_f);
                    if(ilocal_f == -1)
                        continue;
                    
                    // divergence (field)
                    const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                    
                    for(unsigned i = 0; i < 3; ++i)
                        mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * w;
                    mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * w;
                    
                }
            }
        }
        
    }
    
    void HAssemblyEmag::evalVertexSingularityPolar(const unsigned isrcel,
                                                   const unsigned ifieldel,
                                                   const nurbs::Vertex v1,
                                                   const nurbs::Vertex v2,
                                                   const std::map<int, int>& g2locals,
                                                   const std::map<int, int>& g2localf,
                                                   MatrixType& mat) const
    {
        const auto& nsubcells = defaultSubcellN();
        
        const double k = wavenumber();
        
        // source element and connectivity
        const auto p_sel = forest().bezierElement(isrcel);
        const auto& sconn = p_sel->signedGlobalBasisFuncI();
//        const auto& sorder = p_sel->equalIntegrationOrder();
        
        // quadrature orders
        const nurbs::UIntVec sorder{4,4};
        const nurbs::UIntVec forder{4,4};
        
        // field element and connectivity
        const auto p_fel = forest().bezierElement(ifieldel);
        const auto& fconn = p_fel->signedGlobalBasisFuncI();
        
//        uint offset = 0;
//        if(p_fel->degenerate())
//            offset += 2;
//        
//        const auto& forder = p_fel->equalIntegrationOrder(offset);
        
        // and finally loop over all regular integrals
        for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        {
            const auto sparent = igpt_s.get();
            const auto sw = igpt_s.getWeight();
            
            // cached tangent entries
            const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            
            // source element parameters
            const auto x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1, t2);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const double jdet_s = p_sel->jacDet(sparent, t1, t2);
            const double jpiola_s = nurbs::cross(t1, t2).length();
            
            // integrate over field elements
            for(nurbs::IPolarIntegrate igpt_f(nurbs::paramPt(v2), forder, nsubcells); !igpt_f.isDone(); ++igpt_f)
            {
                const auto fparent = igpt_f.get();
                const auto fw = igpt_f.getWeight();
                
                const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                
                // field element parameters
                const auto y = p_fel->eval(fparent);
                const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
                const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                const double jdet_f = p_fel->jacDet(fparent, t1, t2);
                const double jpiola_f = nurbs::cross(t1, t2).length();
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < sconn.size(); ++itest)
                {
                    const auto igbasis_s = sconn[itest];
                    if(igbasis_s == -1)
                        continue;
                    
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                    
                    for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                    {
                        const auto igbasis_f = fconn[itrial];
                        if(igbasis_f == -1)
                            continue;
                        
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * sw * fw;
                        
                    }
                }
            }
        }
    }

    void HAssemblyEmag::evalCoincidentSingularity(const unsigned iel,
                                                  const std::map<int, int>& g2locals,
                                                  const std::map<int, int>& g2localf,
                                                  MatrixType& mat) const
    {
        

        
        // Pointers to elements
        const auto p_sel = forest().bezierElement(iel);
        const auto p_fel = forest().bezierElement(iel);
        
        nurbs::UIntVec sorder{5,5};
        nurbs::UIntVec forder{4,4};
        
//        if(p_fel->degenerate())
//        {
//            std::for_each(sorder.begin(), sorder.end(), [](uint& v){ v *= 2; });
//            std::for_each(forder.begin(), forder.end(), [](uint& v){ v *= 2; });
//        }
        
        // Connectivities
        const auto& sconn = p_sel->signedGlobalBasisFuncI();
        const auto& fconn = p_fel->signedGlobalBasisFuncI();
        
        const auto& k = wavenumber();
        
        for(nurbs::IEqualQuadrature igpt(sorder, forder); !igpt.isDone(); ++igpt)
        {
            const auto gpt4d = igpt.get();
            const auto sparent = gpt4d.srcPt();
            const auto fparent = gpt4d.fieldPt();
            
            const auto w = igpt.getWeight();
            
            // source element parameters
            const auto& t1_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            const auto& x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1_s, t2_s);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const auto& jdet_s = p_sel->jacDet(sparent, t1_s, t2_s);
            const auto& jpiola_s = nurbs::cross(t1_s, t2_s).length();
            
            // field point terms
            const auto& t1_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
            const auto& t2_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
            const auto y = p_fel->eval(fparent);
            const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1_f, t2_f);
            const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
            const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
            const double jdet_f = p_fel->jacDet(fparent, t1_f, t2_f);
            const double jpiola_f = nurbs::cross(t1_f, t2_f).length();
            
            // kernel
            const double r = dist(x, y);
            const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
            
            // now loop over test and trial functions
            for(size_t itest = 0; itest < sconn.size(); ++itest)
            {
                const auto igbasis_s = sconn[itest];
                if(igbasis_s == -1)
                    continue;
                
                const auto ilocal_s = g2locals.at(igbasis_s);
                if(ilocal_s == -1)
                    continue;
                
                // divergence (source)
                const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                
                for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                {
                    const auto igbasis_f = fconn[itrial];
                    if(igbasis_f == -1)
                        continue;
                    
                    const auto ilocal_f = g2localf.at(igbasis_f);
                    if(ilocal_f == -1)
                        continue;
                    
                    // divergence (field)
                    const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                    
                    for(unsigned i = 0; i < 3; ++i)
                        mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * w;
                    mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * jdet_s * jdet_f * w;
                    
                }
            }
        }
    }
    
    void HAssemblyEmag::evalCoincidentSingularityPolar(const unsigned iel,
                                                       const std::map<int, int>& g2locals,
                                                       const std::map<int, int>& g2localf,
                                                       MatrixType& mat) const
    {
        
        
        const double k = wavenumber();
        const std::complex<double> iconst(0.0, 1.0);
        
        const auto p_el = forest().bezierElement(iel);
        const auto& conn = p_el->signedGlobalBasisFuncI();
        
        uint offset = 0;
        if(p_el->degenerate())
            offset += 2;
        const auto& forder = p_el->equalIntegrationOrder(offset);
        const auto& sorder = p_el->equalIntegrationOrder(offset);
        const auto& nsubcells = defaultSubcellN();

        for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        {
            const auto sparent = igpt_s.get();
            const auto sw = igpt_s.getWeight();
            
            // cached tangent entries
            const auto& t1 = p_el->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2 = p_el->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            
            // source element parameters
            const auto x = p_el->eval(sparent);
            const auto& basis_s = p_el->basis(sparent.s, sparent.t, t1, t2);
            const auto& ds_s = p_el->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_el->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const double jdet_s = p_el->jacDet(sparent, t1, t2);
            const double jpiola_s = nurbs::cross(t1, t2).length();
            
            // integrate over field elements
            for(nurbs::IPolarIntegrate igpt_f(sparent, forder, nsubcells); !igpt_f.isDone(); ++igpt_f)
            {
                const auto fparent = igpt_f.get();
                const auto fw = igpt_f.getWeight();
                
                // cached tangent entries
                const auto& t1 = p_el->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_el->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                
                // field element parameters
                const auto y = p_el->eval(fparent);
                const auto& basis_f = p_el->basis(fparent.s, fparent.t, t1, t2);
                const auto& ds_f = p_el->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                const auto& dt_f = p_el->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                const double jdet_f = p_el->jacDet(fparent, t1, t2);
                const double jpiola_f = nurbs::cross(t1, t2).length();
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < conn.size(); ++itest)
                {
                    const auto igbasis_s = conn[itest];
                    if(-1 == igbasis_s)
                        continue;
                    
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                    
                    for(size_t itrial = 0; itrial < conn.size(); ++itrial)
                    {
                        const auto igbasis_f = conn[itrial];
                        if(-1 == igbasis_f)
                            continue;
                        
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw * jdet_s * jdet_f;
                        
                    }
                }
            }
        }
    }
    
    void HAssemblyEmag::computeElCache(const nurbs::UIntVec& qorder,
                                       const std::vector<unsigned>& els,
                                       const std::map<int, int>& g2local,
                                       std::map<std::pair<unsigned, unsigned>, nurbs::DoubleVecVec>& basis,
                                       std::map<std::pair<unsigned, unsigned>, double>& jdet,
                                       std::map<std::pair<unsigned, unsigned>, nurbs::DoubleVec>& div,
                                       std::map<std::pair<unsigned, unsigned>, nurbs::Point3D>& pts,
                                       std::map<unsigned, nurbs::IntVec>& connmap) const
    {
        // first clear maps
        basis.clear();
        jdet.clear();
        div.clear();
        pts.clear();
        connmap.clear();
        
        // Precomputation of source element parameters.
        for(const auto& isel : els)
        {
            const auto p_sel = forest().bezierElement(isel);
            
            // compute set of basis func indices required
            nurbs::IntVec sset;
            const auto& conn = p_sel->signedGlobalBasisFuncI();
            for(uint ibasis = 0; ibasis < p_sel->basisFuncN(); ++ibasis)
            {
                if(g2local.at(conn[ibasis]) != -1)
                    sset.push_back(ibasis);
            }
            connmap[isel] = sset;
            
            // loop over source element quadrature points
            for(nurbs::IElemIntegrate igpt_s(qorder); !igpt_s.isDone(); ++igpt_s)
            {
                const auto sparent = igpt_s.get();
                const auto index = igpt_s.currentIndex();
                const auto pair = std::make_pair(isel, index);
                
                // source element parameters
                pts[pair] = p_sel->eval(sparent);
                
                // tangent vectors
                const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
                
                const double jpiola_s = nurbs::cross(t1, t2).length();
                
                // physical coordinate
                basis[pair] = p_sel->basis(sparent.s, sparent.t, t1, t2);
                const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
                const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
                
                // divergence operator
                nurbs::DoubleVec div_vec(p_sel->basisFuncN(), 0.0);
                
                for(unsigned ibasis = 0; ibasis < p_sel->basisFuncN(); ++ibasis)
                    div_vec[ibasis] = 1./jpiola_s * (ds_s[ibasis][0] + dt_s[ibasis][1]);
                div[pair] = div_vec;
                
                // jacobian determinant
                jdet[pair] = p_sel->jacDet(sparent.s, sparent.t, t1, t2);
            }
        }
        
    }
}
