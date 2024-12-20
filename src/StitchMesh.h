#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>



using namespace std;

struct FancyHalfedge; // Forward declaration

// FancyVertex structure
struct FancyVertex {
    int id;                     // Vertex ID
    FancyHalfedge* outgoing;    // Pointer to an outgoing halfedge

    FancyVertex(int id) : id(id), outgoing(nullptr) {}
};

// FancyHalfedge structure
struct FancyHalfedge {
    FancyVertex* origin;        // Origin vertex
    FancyHalfedge* twin;        // Twin halfedge
    FancyHalfedge* next;        // Next halfedge in the face
    FancyHalfedge* prev;        // Previous halfedge in the face
    int face_id;                // ID of the face to which this halfedge belongs
    bool isCourse = false;      // if the halfedge is a course halfedge
    bool isWale = false;        // if the halfedge is a wale halfedge

    FancyHalfedge() : origin(nullptr), twin(nullptr), next(nullptr), prev(nullptr), face_id(-1) {}
};

// FancyFace structure
struct FancyFace {
    int id;
    FancyHalfedge* edge;        // Pointer to one of its boundary halfedges

    FancyFace(int id) : id(id), edge(nullptr) {}
};

// Custom hash function for std::pair<int, int>
struct pair_hash {
    template <typename T1, typename T2>
    size_t operator ()(const pair<T1, T2>& p) const {
        auto h1 = hash<T1>{}(p.first);
        auto h2 = hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);  // Combine the two hash values
    }
};

void makeMesh(vector<Vector3> graphVertices, vector<KnitGraphEdge> edges) {
    
    // Step 1: Create vertices
    unordered_map<int, FancyVertex*> vertices;
    for (const auto& e : edges) {
        if (vertices.find(e.v1) == vertices.end())
            vertices[e.v1] = new FancyVertex(e.v1);
        if (vertices.find(e.v2) == vertices.end())
            vertices[e.v2] = new FancyVertex(e.v2);
    }

    // Step 2: Create halfedges
    vector<FancyHalfedge*> halfedges;
    unordered_map<pair<int, int>, FancyHalfedge*, pair_hash> edge_map;

    for (const auto& e : edges) {
        int v1 = e.v1, v2 = e.v2;

        if (!e.isBoundary){//not a boundary edge

            // Create two halfedges
            FancyHalfedge* he1 = new FancyHalfedge();
            FancyHalfedge* he2 = new FancyHalfedge();

            //set course and wale flags 
            if (e.isCourse){
                he1->isCourse = true;
                he2->isCourse = true;
            }
            if (e.isWale){
                he1->isWale = true;
                he2->isWale = true;
            }

            he1->origin = vertices[v1];
            he2->origin = vertices[v2];

            he1->twin = he2;
            he2->twin = he1;

            vertices[v1]->outgoing = he1; // Assign an outgoing halfedge
            vertices[v2]->outgoing = he2;

            edge_map[{v1, v2}] = he1;
            edge_map[{v2, v1}] = he2;

            halfedges.push_back(he1);
            halfedges.push_back(he2);
        }
        else{//hit a boundary edge, don't make the twin halfedge
            FancyHalfedge* he1 = new FancyHalfedge();
            if (e.isCourse)
                he1->isCourse = true;
            he1->origin = vertices[v1];
            vertices[v1]->outgoing = he1; 
            edge_map[{v1, v2}] = he1;
            halfedges.push_back(he1);

        }
    }

    for (const auto& he : halfedges){
        cout << "tip vertex = " << he->origin->id << std::endl;
        cout << "tail vertex = " << he->twin->origin->id << std::endl;
        if (he->isCourse)
            cout << "course edge " << std::endl;
        if (he->isWale)
            cout << "wale edge " << std::endl;
        std::cout << "-------------------" << std::endl;
    }

}
