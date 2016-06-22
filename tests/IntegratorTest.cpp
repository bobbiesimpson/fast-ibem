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


/// A class to represent generation of the spline-based
/// coefficients for IGA BE analysis

class HelmholtzCoeffFn : public TPermCoeffFn<complex_t> {
    
public:
    
    typedef fastibem::CollocationAssembly<fastibem::HelmholtzKernel> AssemblyType;
    
    // constructor
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
    
    virtual matform_t matrix_format() const { return MATFORM_NONSYM; }
    
private:
    
    /// assembly instance getter
    AssemblyType* assemblyInstance() const { return mpAssembly; }
    
    /// Pointer to non-owning assembly instance
    AssemblyType* mpAssembly;
};

int main(const int argc, const char* argv[])
{
    try {
        
        const size_t  nmin = 20;
        const real_t  eps  = real_t(1e-4);
        
        INIT();
        
        CFG::set_verbosity( 3 );
        
        // Read in the input and create the geometry object
        std::ifstream ifs(argv[1]);
        if(!ifs)
            nurbs::error("Cannot open file for reading\n");
        nurbs::Geometry g;
        if(!g.loadHBSFile(ifs))
            nurbs::error("Failed to load geometry from hbs data");
        nurbs::Forest forest(g);
        forest.hrefine(1);
        std::cout << "Performing BE analysis with " << forest.globalDofN() << " dof\n";
    
        // create vectors that manage the memory of bounding box data
        const uint n = forest.collocPtN();
        std::vector<nurbs::Point3D> vertices;
        std::vector<nurbs::Point3D> bbmin;
        std::vector<nurbs::Point3D> bbmax;
        
        std::vector<double*> p_vertices(n);
        std::vector<double*> p_bbmin(n);
        std::vector<double*> p_bbmax(n);

        for(nurbs::BoundingBoxIterator it(forest); !it.isDone(); ++it) {
            
            const uint icurrent = it.currentIndex();
            
            // insert point data
            vertices.push_back(it.currentPt());
            //p_vertices.push_back(vertices.back().data());
            const nurbs::Point3D p = it.currentPt();
            p_vertices[icurrent] = new double[3];
            for(uint i = 0; i < 3; ++i)
                p_vertices[icurrent][i] = p.getCoord(i);
            
            bbmin.push_back(it.currentLowerBound());
//            p_bbmin.push_back(bbmin.back().data());
            const nurbs::Point3D bbminpt = it.currentLowerBound();
            p_bbmin[icurrent] = new double[3];
            for(uint i = 0; i < 3; ++i)
                p_bbmin[icurrent][i] = bbminpt.getCoord(i);
            
            bbmax.push_back(it.currentUpperBound());
//            p_bbmax.push_back(bbmax.back().data());
            const nurbs::Point3D bbmaxpt = it.currentUpperBound();
            p_bbmax[icurrent] = new double[3];
            for(uint i = 0; i < 3; ++i)
                p_bbmax[icurrent][i] = bbmaxpt.getCoord(i);
        }
        
        // output bounding box data
        std::vector<std::pair<nurbs::Point3D, nurbs::Point3D>> bbdata;
        for(uint i = 0; i < n; ++i) {
//            std::cout << "vertex: " << vertices[i] << "\n";
//            std::cout << "minbb: " << bbmin[i] << "\n";
//            std::cout << "maxbb: " << bbmax[i] << "\n";
            bbdata.push_back(std::make_pair(bbmin[i], bbmax[i]));
        }
        nurbs::OutputVTK output("sphere_test");
        output.outputGeometry(forest);
        output.outputBoundingBoxSet(bbdata);
        
        
        std::unique_ptr<TCoordinate> coord(nurbs::make_unique<TCoordinate>(p_vertices, 3, p_bbmin, p_bbmax));

        // create kernel and assembly instance
        fastibem::HelmholtzKernel hkernel(std::make_pair("k", 0.0));
        fastibem::CollocationAssembly<fastibem::HelmholtzKernel> assembly(&forest,
                                                                          hkernel,
                                                                          true);
        
//        std::vector<uint> cvec{1};
//        std::vector<uint> bvec{3};
//        auto result = assembly.eval(cvec, bvec);
//        std::cout << "result = " << result[0][0] << "\n";

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
        }// if
        
        // now build the matrix
        std::cout << "━━ building H-matrix ( eps = " << eps << " )" << std::endl;
        
        TTimer                    timer( WALL_TIME );
        TConsoleProgressBar       progress;
        TTruncAcc                 acc(eps, 0.0 );
        HelmholtzCoeffFn          coefffn(&assembly, ct->perm_i2e(), ct->perm_i2e());
        TACAPlus<complex_t>       aca(&coefffn);
        TDenseMBuilder<complex_t> h_builder(&coefffn, &aca);
        TPSMatrixVis              mvis;
        
        // enable coarsening during construction
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
        if(A->is_complex())
            std::cout << "A is complex\n";
        
//        for(uint i = 0; i < 1; ++i) {
//            complex_t result;
//            for(uint j = 0; j < 20; ++j)
//                std::cout << A->centry(i, j) << "\t";
//            std::cout << "\n";
//        }

//        std::vector<uint> cindices{0};
//        std::vector<uint> bindices(20);
//        std::iota(bindices.begin(), bindices.end(), 0);
        
        std::vector<uint> cindices(forest.globalDofN());
        std::iota(cindices.begin(), cindices.end(), 0);
        std::vector<uint> bindices(forest.globalDofN());
        std::iota(bindices.begin(), bindices.end(), 0);
        
        auto hmat = assembly.eval(cindices, bindices);
        for(uint irow = 0; irow < cindices.size(); ++irow) {
            std::complex<double> sum(0.0, 0.0);
            for(uint icol = 0; icol < bindices.size(); ++icol) {
                const auto correct = hmat[irow][icol];
                const std::complex<double> approx(A->centry(irow, icol).re(), A->centry(irow, icol).im());
                sum += correct;
                if(irow == 1)
                    std::cout << std::setprecision(10) << "(" << correct << ")\t";
                //std::cout << std::abs(hmat[irow][icol] - approx) << "\t";
            }
            std::cout << "sum = " << sum << "\n";
        }
        
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