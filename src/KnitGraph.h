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

    //row id and col id for DAG construction
    int rowId = -1;
    int colId = -1;

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

    bool isDisk = true;
    
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


//-----------------------------------------//
// Build DAG for integer coordinates for rows/cols for models of disk topology
//-----------------------------------------//

struct PairHash {
  size_t operator()(const std::pair<int,int>& p) const noexcept {
    // simple hash combine
    return (static_cast<size_t>(static_cast<uint32_t>(p.first)) << 32) ^
           static_cast<uint32_t>(p.second);
  }
};


struct Chain {
  int start;                 // start vertex id
  std::vector<int> verts;    // vertex ids in order
};

struct DAG {
  int n = 0;
  std::vector<std::vector<int>> adj;
  std::vector<int> indeg;
};

//collect rows of knit graph
std::vector<Chain> collectRows(const std::vector<knitGraphVertex>& vertices);

//collect cols of knit graph
std::vector<Chain> collectCols(const std::vector<knitGraphVertex>& vertices);


// Build vertex->chainId map from chains (rows or cols)
std::vector<int> buildVertexToChainMap(int nVerts, const std::vector<Chain>& chains);

    
//Build row DAG
DAG buildRowDAG(const std::vector<knitGraphVertex>& vertices,
                const std::vector<Chain>& rows,
                const std::vector<int>& v2row);


//Build col DAG 
// Column DAG edges: use course links v -> row_out
DAG buildColDAG(const std::vector<knitGraphVertex>& vertices,
                const std::vector<Chain>& cols,
                const std::vector<int>& v2col);


//perform a topological sort over a DAG
std::vector<int> topoSort(const DAG& g);

// Given: topo order vector `order` of size N, build rank map such that rank[node] = topoIndex
static std::vector<int> buildTopoRank(const std::vector<int>& order, int N) {
  if ((int)order.size() != N) throw std::runtime_error("buildTopoRank: size mismatch");
  std::vector<int> rank(N, -1);
  for (int i = 0; i < N; ++i) rank[order[i]] = i;
  return rank;
}
