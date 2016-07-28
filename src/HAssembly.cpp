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
        
        for(unsigned ielem = 0; ielem < forest().elemN(); ++ielem)
        {
            const auto p_el = forest().bezierElement(ielem);
            const auto econn = p_el->globalBasisFuncI();
            const auto& qorder = p_el->equalIntegrationOrder();
            
            for(nurbs::IElemIntegrate igpt(qorder); !igpt.isDone(); ++igpt)
            {
                const auto gpt = igpt.get();
                const auto w = igpt.getWeight();
                
                // physical coordinate
                const auto x = p_el->eval(gpt);
                
                // emag plane wave at x
                const auto pw = pw_functor(x);
                
                // basis funcs
                const auto basis = p_el->localBasis(gpt.s, gpt.t);
                
                // jacobian
                const auto jacob = p_el->jacob(gpt.s, gpt.t);
                
                for(size_t ibasis = 0; ibasis < econn.size(); ++ibasis)
                {
                    // global basis index
                    const auto igbasis = econn[ibasis];
                    
                    for(unsigned i = 0; i < 2; ++i)
                        for(unsigned j = 0; j < 3; ++j)
                            ftemp[igbasis] += -1.0 / std::complex<double>(0.0, omega() * mu()) * basis[ibasis][i] * jacob[i][j]
                            * pw[j] * w;
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
        
        // Vectors of source (test) and field (trial) elements we must compute
        std::vector<unsigned> isrc_els;
        std::vector<unsigned> ifield_els;
        
        // maps from global basis indices to local row or column indices.
        // a local index of -1 denotes that the term is ignored.
        std::map<unsigned, int> g2local_field;
        std::map<unsigned, int> g2local_src;
        
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
                for(const auto& i : el->globalBasisFuncI())
                {
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
                for(const auto& i : el->globalBasisFuncI())
                {
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
        std::vector<std::tuple<unsigned, unsigned, nurbs::Vertex>> vertex_integrals;
        std::vector<std::tuple<unsigned, unsigned, nurbs::Edge, nurbs::Edge>> edge_integrals;
        std::vector<unsigned> coincident_integrals;
        
        // The sets of regular source and field elements
        std::map<unsigned, std::vector<unsigned>> regular_elmap;
        
        nurbs::Edge e1, e2;
        nurbs::Vertex v1, v2;
        
        for(size_t isrcel = 0; isrcel < isrc_els.size(); ++isrcel)
        {
            const auto igsrcel = isrc_els[isrcel];
            const auto p_srcel = forest().bezierElement(igsrcel);
            
            for(size_t ifieldel = 0; ifieldel < ifield_els.size(); ++ifieldel)
            {
                const auto igfieldel = ifield_els[ifieldel];
                const auto p_fieldel = forest().bezierElement(igfieldel);
                
                // edge singularity
                /*if(nurbs::edgeConnected(*p_srcel, *p_fieldel, e1, e2))
                    edge_integrals.push_back(std::make_tuple(igsrcel, igfieldel, e1, e2));

                // vertex singularity
                else if(nurbs::vertexConnected(*p_srcel, *p_fieldel, v1, v2))
                    vertex_integrals.push_back(std::make_tuple(igsrcel, igfieldel, v2));
                
                // coincident singularity
                else */if(igsrcel == igfieldel)
                    coincident_integrals.push_back(igsrcel);
                
                // regular integral
                else
                    regular_elmap[igsrcel].push_back(igfieldel);
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
        
        
        // evaluate vertex elements
        for(const auto& tuple : vertex_integrals)
            evalVertexSingularity(std::get<0>(tuple),
                                  std::get<1>(tuple),
                                  std::get<2>(tuple),
                                  g2local_src,
                                  g2local_field,
                                  matrix);
                                  
        // evaluate coincidnet elements
        for(const auto& iel : coincident_integrals)
            evalCoincidentSingularity(iel,
                                      g2local_src,
                                      g2local_field,
                                      matrix);
        
        
        // and finally loop over all regular integrals
        for(const auto& pair: regular_elmap)
        {
            const auto& igsrcel = pair.first;
            const auto p_sel = forest().bezierElement(igsrcel);
            const auto& sconn = p_sel->globalBasisFuncI();
            
            const auto& sorder = p_sel->equalIntegrationOrder();
            
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
                
                // loop over field elements
                for(const auto& igfieldel : pair.second)
                {
                    const auto p_fel = forest().bezierElement(igfieldel);
                    const auto& fconn = p_fel->globalBasisFuncI();
                    
                    nurbs::UIntVec forder;
                    if(nurbs::dist(x, p_fel->eval(0.0, 0.0)) > p_fel->size() * 2.0)
                        forder = nurbs::UIntVec{2,2};
                    else
                        forder = p_fel->equalIntegrationOrder();
                    
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
                        
                        // kernel
                        const double r = dist(x, y);
                        const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                        
                        // now loop over test and trial functions
                        for(size_t itest = 0; itest < sconn.size(); ++itest)
                        {
                            const auto igbasis_s = sconn[itest];
                            const auto ilocal_s = g2local_src[igbasis_s];
                            if(ilocal_s == -1)
                                continue;
                            
                            // divergence (source)
                            const double div_s = (ds_s[itest][0] + dt_s[itest][1]);
                            
                            for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                            {
                                const auto igbasis_f = fconn[itrial];
                                const auto ilocal_f = g2local_field[igbasis_f];
                                if(ilocal_f == -1)
                                    continue;
                                
                                // divergence (field)
                                const double div_f = (ds_f[itrial][0] + dt_f[itrial][1]);
                                
                                for(unsigned i = 0; i < 3; ++i)
                                    matrix[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                                matrix[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw;
                                
                            }
                        }
                    }
                }
            }
        }
        return matrix;
    }
    
    void HAssemblyEmag::evalEdgeSingularity(const unsigned isrcel,
                                            const unsigned ifieldel,
                                            const nurbs::Edge e1,
                                            const nurbs::Edge e2,
                                            const std::map<unsigned, int>& g2locals,
                                            const std::map<unsigned, int>& g2localf,
                                            MatrixType& mat) const
    {
        const auto& nsubcells = defaultSubcellN();
        
        const double k = wavenumber();
        
        // source element and connectivity
        const auto p_sel = forest().bezierElement(isrcel);
        const auto& sconn = p_sel->globalBasisFuncI();
        const auto& sorder = p_sel->equalIntegrationOrder();
        
        // field element and connectivity
        const auto p_fel = forest().bezierElement(ifieldel);
        const auto& fconn = p_fel->globalBasisFuncI();
        const auto& forder = p_fel->equalIntegrationOrder();
        
        // and finally loop over all regular integrals
        for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        {
            const auto sparent = igpt_s.get();
            const auto sw = igpt_s.getWeight();
            
            // source element parameters
            const auto x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const double jdet_s = p_sel->jacDet(sparent);
            
            // integrate over field elements
            for(nurbs::IPolarIntegrate igpt_f(nurbs::projectPt(sparent, e1, e2), forder, nsubcells); !igpt_f.isDone(); ++igpt_f)
            {
                const auto fparent = igpt_f.get();
                const auto fw = igpt_f.getWeight();

                // field element parameters
                const auto y = p_fel->eval(fparent);
                const auto& basis_f = p_fel->basis(fparent.s, fparent.t);
                const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                const double jdet_f = p_fel->jacDet(fparent);
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < sconn.size(); ++itest)
                {
                    const auto igbasis_s = sconn[itest];
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = (ds_s[itest][0] + dt_s[itest][1]);
                    
                    for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                    {
                        const auto igbasis_f = fconn[itrial];
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw;
                        
                    }
                }
            }
        }
    }
    
    void HAssemblyEmag::evalVertexSingularity(const unsigned isrcel,
                                              const unsigned ifieldel,
                                              const nurbs::Vertex v2,
                                              const std::map<unsigned, int>& g2locals,
                                              const std::map<unsigned, int>& g2localf,
                                              MatrixType& mat) const
    {
        const auto& nsubcells = defaultSubcellN();
        
        const double k = wavenumber();
        
        // source element and connectivity
        const auto p_sel = forest().bezierElement(isrcel);
        const auto& sconn = p_sel->globalBasisFuncI();
        const auto& sorder = p_sel->equalIntegrationOrder();
        
        // field element and connectivity
        const auto p_fel = forest().bezierElement(ifieldel);
        const auto& fconn = p_fel->globalBasisFuncI();
        const auto& forder = p_fel->equalIntegrationOrder();
        
        // and finally loop over all regular integrals
        for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        {
            const auto sparent = igpt_s.get();
            const auto sw = igpt_s.getWeight();
            
            // source element parameters
            const auto x = p_sel->eval(sparent);
            const auto& basis_s = p_sel->basis(sparent.s, sparent.t);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const double jdet_s = p_sel->jacDet(sparent);
            
            // integrate over field elements
            for(nurbs::IPolarIntegrate igpt_f(nurbs::paramPt(v2), forder, nsubcells); !igpt_f.isDone(); ++igpt_f)
            {
                const auto fparent = igpt_f.get();
                const auto fw = igpt_f.getWeight();
                
                // field element parameters
                const auto y = p_fel->eval(fparent);
                const auto& basis_f = p_fel->basis(fparent.s, fparent.t);
                const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                const double jdet_f = p_fel->jacDet(fparent);
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < sconn.size(); ++itest)
                {
                    const auto igbasis_s = sconn[itest];
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = (ds_s[itest][0] + dt_s[itest][1]);
                    
                    for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                    {
                        const auto igbasis_f = fconn[itrial];
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw;
                        
                    }
                }
            }
        }
    }
    
    void HAssemblyEmag::evalCoincidentSingularity(const unsigned iel,
                                                  const std::map<unsigned, int>& g2locals,
                                                  const std::map<unsigned, int>& g2localf,
                                                  MatrixType& mat) const
    {

        
        const double k = wavenumber();
        
        const auto p_el = forest().bezierElement(iel);
        const auto& conn = p_el->globalBasisFuncI();
        
        const auto& forder = p_el->equalIntegrationOrder();
        const auto& sorder = p_el->equalIntegrationOrder();
        const auto& nsubcells = defaultSubcellN();
        
        // and finally loop over all regular integrals
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
                
                // kernel
                const double r = dist(x, y);
                const auto ekernel = std::exp(std::complex<double>(0.0, -k * r)) / (4.0 * nurbs::PI * r);
                
                // now loop over test and trial functions
                for(size_t itest = 0; itest < conn.size(); ++itest)
                {
                    const auto igbasis_s = conn[itest];
                    const auto ilocal_s = g2locals.at(igbasis_s);
                    if(ilocal_s == -1)
                        continue;
                    
                    // divergence (source)
                    const double div_s = (ds_s[itest][0] + dt_s[itest][1]);
                    
                    for(size_t itrial = 0; itrial < conn.size(); ++itrial)
                    {
                        const auto igbasis_f = conn[itrial];
                        const auto ilocal_f = g2localf.at(igbasis_f);
                        if(ilocal_f == -1)
                            continue;
                        
                        // divergence (field)
                        const double div_f = (ds_f[itrial][0] + dt_f[itrial][1]);
                        
                        for(unsigned i = 0; i < 3; ++i)
                            mat[ilocal_s][ilocal_f] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                        mat[ilocal_s][ilocal_f] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw;
                        
                    }
                }
            }
        }
    }

}