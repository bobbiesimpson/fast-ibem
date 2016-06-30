#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include "Forest.h"
#include "Geometry.h"
#include "base.h"
#include "IPolarIntegrate.h"
#include "OutputVTK.h"

int main(int argc, char* argv[])
{
    const double k = 1.0;
    const nurbs::Point3D k_dir(1.0, 0.0, 0.0);
    //
    // Geometry and forest setup
    //
    
    std::ifstream ifs(argv[1]);
    if(!ifs)
        nurbs::error("Cannot open file for reading\n");
    
    uint refine = 0;
    if(argc > 2)
        refine = std::atoi(argv[2]);
    
    nurbs::Geometry g;
    if(!g.loadHBSFile(ifs))
        nurbs::error("Failed to load geometry from hbs data");
    g.flipNormals(true);
    nurbs::Forest forest(g);
    forest.hrefine(refine);
    
    //
    // now assemble entries into A matrix
    //
    Eigen::MatrixXcf A(forest.globalDofN(), forest.globalDofN());

    std::map<uint, std::vector<uint>> cpt_elcache;
    std::map<uint, nurbs::Point3D> cpt_cache;
    std::map<uint, std::vector<double>> cpt_basis;
    std::map<uint, std::vector<uint>> cpt_ibasis;
    
    // first evaluate all singular integrals
    for(uint icelem = 0; icelem < forest.elemN(); ++icelem) {
        
        const auto cel = forest.element(icelem);
        const auto ibasisvec = cel->globalBasisFuncI();
        
        for(uint icolloc = 0; icolloc < cel->collocPtN(); ++icolloc) {
            
            const uint igcolloc = cel->globalCollocI(icolloc);
            cpt_elcache[igcolloc].push_back(icelem);
            const auto s_parent = cel->collocParentCoord(icolloc);
            const auto xs = cel->eval(s_parent.s, s_parent.t);
            cpt_cache[igcolloc] = xs;
            
            auto find = cpt_basis.find(igcolloc);
            if(find == cpt_basis.end()) {
                cpt_basis[igcolloc] = cel->basis(s_parent.s, s_parent.t);
                cpt_ibasis[igcolloc] = cel->globalBasisFuncI();
            }
            
            for(nurbs::IPolarIntegrate igpt(s_parent, cel->integrationOrder()); !igpt.isDone(); ++igpt) {
                
                const nurbs::GPt2D gpt = igpt.get();
                const nurbs::Point3D xf = cel->eval(gpt);
                const auto basis = cel->basis(gpt.s, gpt.t);
                const auto n = cel->normal(gpt);
                
                const double r = dist(xs, xf);
                const std::complex<double> const1(0.0, k * r);
                const double drdn = nurbs::dot((xf-xs) / r, n);
                const auto kernel = std::exp(const1) / (4.0 * nurbs::PI * r * r ) * (1.0 - const1) * drdn;
                
                for(uint ibasis = 0; ibasis < cel->basisFuncN(); ++ibasis)
                    A(igcolloc, ibasisvec[ibasis]) += kernel * basis[ibasis] * cel->jacDet(gpt) * igpt.getWeight();
                
            }
        }
    }
    
    //
    // and now regular integration
    //
    for(uint igcolloc = 0; igcolloc < forest.collocPtN(); ++igcolloc) {
        
        const auto xs = cpt_cache[igcolloc];
        
        for(uint ielem = 0; ielem < forest.elemN(); ++ielem) {
            
            auto elvec = cpt_elcache[igcolloc];
            
            auto find = std::find(elvec.begin(), elvec.end(), ielem);
            if(find != elvec.end()) // if it exists, already computed
                continue;

            const auto el = forest.element(ielem);
            const auto ibasisvec = el->globalBasisFuncI();
            
            for(nurbs::IElemIntegrate igpt(el->integrationOrder()); !igpt.isDone(); ++igpt) {
                
                const nurbs::GPt2D gpt = igpt.get();
                const nurbs::Point3D xf = el->eval(gpt);
                const auto basis = el->basis(gpt.s, gpt.t);
                const auto n = el->normal(gpt);
                
                const double r = dist(xs, xf);
                const std::complex<double> const1(0.0, k * r);
                const double drdn = nurbs::dot((xf-xs) / r, n);
                const auto kernel = std::exp(const1) / (4.0 * nurbs::PI * r * r ) * (1.0 - const1) * drdn;
                
                for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis)
                    A(igcolloc, ibasisvec[ibasis]) += kernel * basis[ibasis] * el->jacDet(gpt) * igpt.getWeight();
            }
        }
    }
    
    //
    // jump terms
    //
    for(uint igcolloc = 0; igcolloc < forest.globalDofN(); ++igcolloc) {
        const auto basis = cpt_basis[igcolloc];
        const auto basis_ivec = cpt_ibasis[igcolloc];
        assert(basis.size() == basis_ivec.size());
        for(size_t ibasis = 0; ibasis < basis.size(); ++ibasis)
            A(igcolloc, basis_ivec[ibasis]) += 0.5 * basis[ibasis];
    }
    
//    for(auto irow = 0; irow < A.rows(); ++irow) {
//        auto sum = A(0,0);
//        sum *= 0.0; // set to zero
//        for(auto icol = 0; icol < A.cols(); ++icol) {
//            sum += A(irow, icol);
//        }
//        std::cout << sum << "\n";
//    }
    
    //
    // fill RHS
    //
    Eigen::VectorXcf b(forest.globalDofN());
    b.setZero();
    for(uint igcolloc = 0; igcolloc < forest.globalDofN(); ++igcolloc) {
        const auto xs = cpt_cache[igcolloc];
        b(igcolloc) += std::exp(std::complex<double>(0.0, k * dot(xs, k_dir)));
    }

    std::cout << A << "\n";
    std::cout << b << "\n";
    //auto x = A.inverse() * b;
    Eigen::VectorXcf x = A.colPivHouseholderQr().solve(b);
    std::cout << x << "\n";
    std::cout << x.coeffRef(0, 0) << "\n";
    
    
    nurbs::OutputVTK output("direct_helmholtz");
    std::vector<std::complex<double>> soln(forest.globalDofN());
    for(uint i = 0; i < forest.globalDofN(); ++i)
        soln[i] = x(i);
    
    output.outputComplexAnalysisField(forest, "acoustic potential", soln);
    
    
    return EXIT_SUCCESS;
}