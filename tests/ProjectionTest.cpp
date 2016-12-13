#include <iostream>
#include <fstream>
#include <chrono>
#include "Forest.h"
#include "Geometry.h"
#include "BSplineSpace.h"
#include "BezierNodalElement.h"
#include "NodalElement.h"
#include "IElemIntegrate.h"
#include "OutputVTK.h"
#include "HConformingForest.h"
#include "Functor.h"
#include "Norm.h"
#include "NURBSCache.h"

#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>

using namespace nurbs;
using namespace fastibem;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Triplet;

int main(int argc, char* argv[]) {
    
    try {
        
        
        std::cout << "Running projection test.....\n";
        
        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            error("Cannot open file for reading\n");
        Geometry g;
        if(!g.loadHBSFile(ifs))
            error("Failed to load geometry from hbs data");
        
//        g.rescale(2.0/64.0);
        EmagPlaneWave pw_functor(Point3D(146.70915, 0.0, 0.0),
                                 {0.0, 0.0, 1.0});
        
        //SinusoidalFunctor s_functor;
        
        Forest forest(g);
        HDivForest divforest(g);
        
        // refinement
        uint refine = 0;
        const uint max_refine = 10;
        if(argc > 2) {
            auto input = std::atoi(argv[2]);
            if(input < 0)
                std::cout << "Cannot supplied negative refinement integer. Carrying on with no h-refinement";
            if(input > max_refine) {
                std::cout << "Truncating refinement to " << max_refine << " levels.";
                refine = max_refine;
            }
            else {
                std::cout << "Applying " << input << " levels of h-refinement\n";
                refine = input;
            }
        }

        divforest.hrefine(refine);
        
        // assembly
        const unsigned ndof = divforest.globalDofN();
        std::cout << "performing projection test with " << ndof << " degrees of freedom and " << divforest.elemN() << " elements\n";
        
        Eigen::VectorXd freal(ndof);
        Eigen::VectorXd fimag(ndof);
        
        for(size_t i = 0; i < ndof; ++i)
        {
            freal(i) = 0.0;
            fimag(i) = 0.0;
        }
        
        std::vector<Triplet> coefficients;
        
        SpMat M(ndof, ndof);
        
        for(uint ielem = 0; ielem < divforest.elemN(); ++ielem)
        {
//            std::cout << "Element: " << ielem << "\n";
            const auto el = divforest.bezierElement(ielem);
            
            const auto& conn = el->signedGlobalBasisFuncI();
//            std::cout << conn << "\n";
            
            std::vector<std::vector<double>> submatrix;
            for(size_t i = 0; i < conn.size(); ++i)
                submatrix.push_back(std::vector<double>(conn.size(), 0.0));
            
            for(IElemIntegrate igpt(el->integrationOrder(2)); !igpt.isDone(); ++igpt)
            {
                const auto gpt = igpt.get();
                const auto weight = igpt.getWeight();
                
                const auto& basis = el->basis(gpt.s, gpt.t);
                const auto& x = el->eval(gpt.s, gpt.t);
                const auto& jdet = el->jacDet(gpt);
                
                const auto pw = pw_functor(x);
//                const auto func = s_functor(x);
                
                for(size_t itest = 0; itest < conn.size(); ++itest)
                {
                    const auto& gtest_i = conn[itest];
                    if(-1 == gtest_i)
                        continue;
                    for(size_t itrial = 0; itrial < conn.size(); ++itrial)
                    {
                        const auto& gtrial_i = conn[itrial];
                        if(-1 == gtrial_i)
                            continue;
                        
                        double val = 0.0;
                        for(uint i = 0; i < 3; ++i)
                            val += basis[itest][i] * basis[itrial][i];
                        
                        val *= weight * jdet;
                        submatrix[itest][itrial] += val;
                    }
                    // force vector
                    
                    for(unsigned j = 0; j < 3; ++j)
                    {
//                        freal(gtest_i) += basis[itest][j] * func[j] * weight * jdet;
                        freal(gtest_i) += basis[itest][j] * pw[j].real() * weight * jdet;
                        fimag(gtest_i) += basis[itest][j] * pw[j].imag() * weight * jdet;
                    }
//                    for(unsigned i = 0; i < 2; ++i)
//                        for(unsigned j = 0; j < 3; ++j)
//                        {
//                            freal(gtest_i) += localbasis[itest][i] * jacob[i][j] * pw[j].real() * weight;
//                            fimag(gtest_i) += localbasis[itest][i] * jacob[i][j] * pw[j].imag() * weight;
//                        }

                }
            }
            
            for(size_t itest = 0; itest < conn.size(); ++itest)
            {
                const auto& gtest_i = conn[itest];
                if(-1 == gtest_i)
                    continue;
                for(size_t itrial = 0; itrial < conn.size(); ++itrial)
                {
                    const auto& gtrial_i = conn[itrial];
                    if(-1 == gtrial_i)
                        continue;
                    coefficients.push_back(Triplet(gtest_i, gtrial_i, submatrix[itest][itrial]));
                }
            }
        }
        
        M.setFromTriplets(coefficients.begin(), coefficients.end());
        
//        const auto& sspace = divforest.space(0, ParamDir::S);
//        const auto& tspace = divforest.space(0, ParamDir::T);
//        
//        for(uint iedge = 0; iedge < NEDGES; ++iedge)
//        {
//            for(uint iparam = 0; iparam < 2; ++iparam)
//            {
//                const auto indices = localBasisIVec(edgeType(iedge), ParamDirType(iparam), std::make_pair(sspace, tspace));
//                
//                for(const auto& index : indices)
//                {
//                    M.coeffRef(index, index) = 1.0;
//                    freal(index) = 0.0;
//                    fimag(index) = 0.0;
//                }
//            }
//        }
        

        Eigen::SimplicialCholesky<SpMat> chol(M);  // performs a Cholesky factorization of A
        Eigen::VectorXd xreal = chol.solve(freal);
        Eigen::VectorXd ximag = chol.solve(fimag);
        
        std::vector<std::complex<double>> solnvec;
        std::vector<double> real_solnvec;
        
        for(size_t i = 0; i < ndof; ++i)
        {
            std::complex<double> centry(xreal(i), ximag(i));
//            std::complex<double> centry = (87==i) ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
//            std::cout << centry << "\n";
            solnvec.push_back(centry);
            real_solnvec.push_back(xreal(i));
        }
        
        nurbshelper::NURBSCache::Instance().clear();
        
        // Output solution
        std::string filename("projectiontest");
        filename.append(std::to_string(refine));
        nurbs::OutputVTK output(filename, 30);
        output.outputComplexVectorField(divforest, "no_name", solnvec);
        
        nurbshelper::NURBSCache::Instance().clear();
        
        // Norm
        std::cout << "L2 graph norm: " << nurbs::L2graphNorm(divforest, solnvec) << "\n";
        
        
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e) {
        
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
        
    }
    catch(...) {
        
        std::cerr << "Unknown exception occured\n";
        return EXIT_FAILURE;
        
    }
}
