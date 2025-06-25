//geometry-central includes 
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/meshio.h"

//polyscope includes 
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include <vector> 
#include <optional>
#include <iostream>
#include <fstream>


using namespace geometrycentral;
using namespace geometrycentral::surface;


struct knitGraphVertex{
    int id = -1;//id of the vertex
    Vector3 position;//position embedded in R^3 (if position information is available)
    int row_in = -1;
    int row_out = -1;
    int col_in[2] = {-1, -1};
    int col_out[2] = {-1, -1};
};

struct KnitGraphHalfedge{
    int id = -1;//id of the halfedge
    bool isRowOut = false;
    bool isRowIn = false;
    bool isWaleOut = false;
    bool isWaleIn = false;
    bool isBoundary = false;
    int tip = -1;//id of vertex at the tip
    int tail = -1;//id of vertex at the tail
    KnitGraphHalfedge* twin;
    KnitGraphHalfedge() : twin(nullptr){};
    //if the halfedge has been visited in the face extraction procedure 
    bool isVisited = false;//has the halfedge been visited
    
    
};

struct KnitGraphEdge{

    int v1 = -1;
    int v2 = -1;
    bool isCourse = false;
    bool isWale = false;
    bool isBoundary = true;

};

struct KnitGraph{

    std::vector<knitGraphVertex> vertices;
    std::vector<KnitGraphHalfedge*> halfedges;

    std::vector<Vector3> embeddedVertices;
    std::vector<KnitGraphEdge> embeddedEdges;

    void renderGraph();

    void traceFaces();

    void traceShortRows();

    //for saving your graph as a line element obj
    void writeLineElementObj();

};
