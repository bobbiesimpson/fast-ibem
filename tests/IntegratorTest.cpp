#include <iostream>
#include <numeric>
#include <hlib.hh>
#include "Forest.h"
#include "Geometry.h"
#include "base.h"
#include "CollocationAssembly.h"
#include "BoundingBoxIterator.h"
#include "OutputVTK.h"

using namespace HLIB;

using complex_t = HLIB::complex;
using real_t = HLIB::real;

///
/// A class to represent generation of the spline-based
/// coefficients for IGA BE analysis
///
class HelmholtzCoeffFn : public TPermCoeffFn<complex_t> {
    
public:
    
    /// Typedef for assembly type
    typedef fastibem::CollocationAssembly<fastibem::HelmholtzKernel> AssemblyType;
    
    // Constructor
    HelmholtzCoeffFn(AssemblyType* a,
                     const TPermutation * row_perm,
                     const TPermutation * col_perm)
    :TPermCoeffFn<complex_t>(row_perm, col_perm),
     mpAssembly(a) {}
    
    /// Override the evaluation function
    void eval(const std::vector<idx_t>& rowidxs,
              const std::vector<idx_t>& colidxs,
              complex_t* matrix ) const override
    {
        const size_t  n = rowidxs.size();
        const size_t  m = colidxs.size();
        
        //std::cout << "computing block of size: " << n << "x" << m << "\n";
        
        std::vector<uint> gcolloc_vec;
        for(auto irow = 0; irow < n; ++irow)
            gcolloc_vec.push_back(rowidxs[irow]);
        
        std::vector<uint> gbasis_vec;
        for(auto icol = 0; icol < m; ++icol)
            gbasis_vec.push_back(colidxs[icol]);
        
        auto rmat = assemblyInstance()->eval(gcolloc_vec, gbasis_vec);
        
        for(size_t j = 0; j < m; ++j)
            for (size_t i = 0; i < n; ++i) {
                const auto cval = rmat[i][j];
                matrix[j*n +i] = complex_t(cval.real(),
                                           cval.imag());
            }
        
    }
    
    using TPermCoeffFn<complex_t>::eval;
    
    /// Collocation BE matrix is non-symmetric
    virtual matform_t matrix_format() const override { return MATFORM_NONSYM; }
    
private:
    
    /// assembly instance getter
    AssemblyType* assemblyInstance() const { return mpAssembly; }
    
    /// Pointer to non-owning assembly instance
    AssemblyType* mpAssembly;
};

int main(const int argc, const char* argv[])
{
    try {
        
        // solver parameters
        const size_t  nmin = 20;                // min block size
        const real_t  eps  = real_t(1e-4);      // H-matrix precision
        
        // Helmholtz parameters
        const double k = 1.0;                  // wavenumber
        const nurbs::Point3D d(1.0, 0.0, 0.0);  // direction of plane wave
        
        INIT();
        CFG::set_verbosity( 3 );
        
        //
        // Geometry and forest setup
        //
        std::ifstream ifs(argv[1]);
        if(!ifs)
            nurbs::error("Cannot open file for reading\n");
        nurbs::Geometry g;
        if(!g.loadHBSFile(ifs))
            nurbs::error("Failed to load geometry from hbs data");
        g.flipNormals(true);
        nurbs::Forest forest(g);
        
        //
        // h-refinement
        //
        uint refine = 0;
        const uint max_refine = 5;
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
        //forest.degreeReduce(3); // reduce to linears
        forest.hrefine(refine);
        std::cout << "Performing BE analysis with " << forest.globalDofN() << " dof and " << forest.elemN() << " elements\n";
    
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
        for(nurbs::BoundingBoxIterator it(forest); !it.isDone(); ++it) {
            
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
        std::unique_ptr<TCoordinate> coord(nurbs::make_unique<TCoordinate>(p_vertices, 3, p_bbmin, p_bbmax));
        
        // Output bounding box data
        nurbs::OutputVTK output("sphere_test");
        output.outputGeometry(forest);
        output.outputBoundingBoxSet(bbdata);

        //
        // Kernel and assembly
        //
        fastibem::HelmholtzKernel hkernel(std::make_pair("k", k));
        fastibem::CollocationAssembly<fastibem::HelmholtzKernel> assembly(&forest,
                                                                          hkernel,
                                                                          true);

        //
        // Hlibpro cluster tree setup
        //
        TAutoBSPPartStrat  part_strat;
        TBSPCTBuilder      ct_builder(&part_strat, nmin);
        auto               ct = ct_builder.build(coord.get());
        TStdGeomAdmCond    adm_cond(2.0);
        TBCBuilder         bct_builder;
        auto               bct = bct_builder.build(ct.get(), ct.get(), &adm_cond);
        
        if(verbose(2)) {
            TPSClusterVis        c_vis;
            TPSBlockClusterVis   bc_vis;
            
            c_vis.print( ct->root(), "bem1d_ct" );
            bc_vis.print( bct->root(), "bem1d_bct" );
        }
        
        //
        // Hlibpro Hmatrix setup
        //
        std::cout << "━━ building H-matrix ( eps = " << eps << " )" << std::endl;
        
        TTimer                    timer( WALL_TIME );
        TConsoleProgressBar       progress;
        TTruncAcc                 acc(eps, 0.0 );
        HelmholtzCoeffFn          coefffn(&assembly, ct->perm_i2e(), ct->perm_i2e());
        TACAPlus<complex_t>       aca(&coefffn);
        TDenseMBuilder<complex_t> h_builder(&coefffn, &aca);
        TPSMatrixVis              mvis;
        h_builder.set_coarsening(false);
        
        timer.start();
        auto  A = h_builder.build(bct.get(), acc, &progress);
        timer.pause();
        
        std::cout << "    done in " << timer << std::endl;
        std::cout << "    size of H-matrix = " << Mem::to_string( A->byte_size() ) << std::endl;
        
        if(verbose(2)) {
            mvis.svd( true );
            mvis.print( A.get(), "bem1d_A" );
        }

        for(size_t j = 0; j < n; ++j)
            std::cout << A->centry(0, j) << " ";
        std::cout << "\n";
        
//        std::cout << "manual result\n\n";
//        std::vector<complex_t> result(n);
//        const auto xs = forest.collocPt(0, 0);
//        for(uint ielem = 0; ielem < forest.elemN(); ++ielem) {
//            const auto el = forest.element(ielem);
//            const auto gbasis_vec = el->globalBasisFuncI();
//            
//            for(nurbs::IElemIntegrate igpt(el->integrationOrder()); !igpt.isDone(); ++igpt) {
//                const nurbs::GPt2D gpt = igpt.get();
//                const nurbs::Point3D xf = el->eval(gpt);
//                const auto basis = el->basis(gpt.s, gpt.t);
//                const auto normal = el->normal(gpt);
//                for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis) {
//                    const auto kernel = hkernel.evalDLP(xs, xf, normal);
//                    result[gbasis_vec[ibasis]] += /*complex_t(kernel.real(),kernel.imag()) * */basis[ibasis] * el->jacDet(gpt) * igpt.getWeight();
//                }
//            }
//        }
//        
//        complex_t summer;
//        for(const auto& v : result) {
//            summer += v;
//            std::cout << "(" << v << ")";
//        }
        //
        // RHS vector setup (planewave)
        //
        std::map<uint, complex_t> rhsmap;
        for(uint ispace = 0; ispace < forest.spaceN(); ++ispace) {
            const auto space = forest.space(ispace);
            for(uint icolloc = 0; icolloc < space.grevilleAbscissaPtN(); ++icolloc) {
                const uint gindex = forest.globalCollocI(ispace, icolloc);
                auto search = rhsmap.find(gindex);
                if(search == rhsmap.end()) {
                    const nurbs::Point3D xc = forest.collocPt(ispace, icolloc);
                    std::cout << xc.length() << "\n";
                    const auto phi_i = std::exp(std::complex<double>(0.0, k * nurbs::dot(d, xc)));
                    const auto cval = complex_t(phi_i.real(), phi_i.imag());
                    rhsmap[gindex] = cval;
                }
            }
        }
    
        auto b = A->row_vector();
        for(size_t i = 0; i < n; i++)
            b->set_centry(i, rhsmap[i]);
        ct->perm_e2i()->permute( b.get() ); 
        
//        for(size_t i = 0; i < n; i++) {
//            const auto entry = b->centry(i);
//            std::cout << entry << "\t" << entry.abs() << "\n";
//        }

        //
        // Preconditioning
        //
        auto  B = A->copy();
        timer.start();
        auto  A_inv = factorise_inv( B.get(), acc, & progress );
        timer.pause();
        std::cout << " done in " << timer << std::endl;
        std::cout << "    size of LU factor = " << Mem::to_string( B->byte_size() ) << std::endl;
        std::cout << "    inversion error   = " << std::scientific << std::setprecision( 4 )
        << inv_approx_2( A.get(), A_inv.get() ) << std::endl;
        
        if( verbose( 2 ) )
            mvis.print( B.get(), "bem1d_LU" );

        //
        // GMRES solver
        //
        std::cout << std::endl << "━━ solving system" << std::endl;
        
//        TGMRES solver(100);
//        TSolver::TInfo  solve_info(false, verbose( 2 ));
        
        TAutoSolver     solver( 1000 );
        TSolver::TInfo  solve_info( false, verbose( 2 ) );
        
        auto x = A->col_vector();
        
        timer.start();
        solver.solve( A.get(), x.get(), b.get(), A_inv.get(), & solve_info );
        
        if ( solve_info.converged() ) {
            std::cout << "  converged in " << timer << " and "
            << solve_info.n_iter() << " steps with rate " << solve_info.conv_rate()
            << ", |r| = " << solve_info.res_norm() << std::endl;
            
            ct->perm_i2e()->permute( x.get() );
            
            // create output vector and set to VTK file
            std::vector<std::complex<double>> solnvec(n);
            for(size_t i = 0; i < n; ++i) {
                const auto entry = x->centry(i);
                solnvec[i] = std::complex<double>(entry.re(), entry.im());
                std::cout << entry.re() << "," << entry.im() << "\n";
            }
            output.outputComplexAnalysisField(forest, solnvec);
        }
        else
            std::cout << "  not converged in " << timer << " and "
            << solve_info.n_iter() << " steps " << std::endl;
        
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