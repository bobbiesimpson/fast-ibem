#include "HAssembly.h"
#include "BoundingBoxIterator.h"
#include "OutputVTK.h"

namespace fastibem {
    
    /// Set up the block cluster tree using the mesh
    void HAssembly::initClusterData(const nurbs::Forest& forest)
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
    }
    
    
    
    
    
    
}