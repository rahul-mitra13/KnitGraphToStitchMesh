#include "KnitGraph.h"
#include <fstream>

//render the graph
void KnitGraph::renderGraph(){

    std::vector<std::array<int, 2>> edges;
    //visualize the knit graph vertices with the virtual connections
    for (auto &v : vertices){
        embeddedVertices.push_back(v.position);
        if (v.row_out != -1){
            edges.push_back({v.id, v.row_out});
            KnitGraphEdge e;
            e.v1 = v.id;
            e.v2 = v.row_out;
            e.isCourse = true;
            //original halfedge
            KnitGraphHalfedge* he1 = new KnitGraphHalfedge();
            he1->tail = v.id;
            he1->tip = v.row_out;
            he1->isRowOut = true;
            //twin halfedge
            KnitGraphHalfedge* he2 = new KnitGraphHalfedge();
            he2->tail = v.row_out;
            he2->tip = v.id;
            he2->isRowIn = true;
            he1->twin = he2;
            he2->twin = he1;
            halfedges.push_back(he1);
            halfedges.push_back(he2);
            e.isBoundary = true;
            embeddedEdges.push_back(e);
        }
        if (v.col_out[0] != -1){
            edges.push_back({v.id, v.col_out[0]});
            KnitGraphEdge e;
            e.v1 = v.id;
            e.v2 = v.col_out[0];
            e.isWale = true;
            embeddedEdges.push_back(e);

            //original halfedge
            KnitGraphHalfedge* he1 = new KnitGraphHalfedge();
            he1->tail = v.id;
            he1->tip = v.col_out[0];
            he1->isWaleOut = true;
            //twin halfedge
            KnitGraphHalfedge* he2 = new KnitGraphHalfedge();
            he2->tail = v.col_out[0];
            he2->tip = v.id;
            he2->isWaleIn = true;
            he1->twin = he2;
            he2->twin = he1;
            halfedges.push_back(he1);
            halfedges.push_back(he2);

        }
        if (v.col_out[1] != -1){
            edges.push_back({v.id, v.col_out[1]});
            KnitGraphEdge e;
            e.v1 = v.id;
            e.v2 = v.col_out[1];
            e.isWale = true;
            embeddedEdges.push_back(e);

            //original halfedge
            KnitGraphHalfedge* he1 = new KnitGraphHalfedge();
            he1->tail = v.id;
            he1->tip = v.col_out[1];
            he1->isWaleOut = true;
            //twin halfedge
            KnitGraphHalfedge* he2 = new KnitGraphHalfedge();
            he2->tail = v.col_out[1];
            he2->tip = v.id;
            he2->isWaleIn = true;
            he1->twin = he2;
            he2->twin = he1;
            halfedges.push_back(he1);
            halfedges.push_back(he2);
        }
    }

    auto graphReal = polyscope::registerCurveNetwork("knit graph", embeddedVertices, edges);
    graphReal -> setRadius(0.001);
    graphReal -> setEnabled(true);

    for (int i = 0; i < edges.size(); i++){
        for (int j = 0; j < edges.size(); j++){
            if (i == j) continue;
            if (edges[i][0] == edges[j][0] && edges[i][1] == edges[j][1]){
                std::cout << "duplicate edges " << std::endl;
                exit(1);
            }
        }
    }


}

void KnitGraph::traceFaces(){

    
    //there's a 1-to-1 mapping between the face list and the edge list 
    std::vector<std::vector<size_t>> faces;
    
    std::map<std::pair<int, int>, KnitGraphHalfedge*> vPair2He;
    for (const auto& he : halfedges){
        vPair2He[std::make_pair(he->tail, he->tip)] = he;
    }

    //now set all the boundary halfedges to visited 
    for (const auto& he : halfedges){
        if (he->isBoundary) he->isVisited = true;
    }

    //now trace the faces 
    for (auto &he : halfedges){
        if (he->isVisited) continue;//skip halfedges we've already visited
        if (he -> isBoundary) continue;//skip boundary halfedges
    
        //need to come up with a way to trace bottom coure, wale, top course, wale
        std::vector<size_t> currFace;
        KnitGraphHalfedge * currHe = he;
        currHe -> isVisited = true;
        int startVertex = currHe->tail;
        do{
            std::cout << "looping through edges in a face " << std::endl;
            currFace.push_back(currHe->tail);
            if (currHe -> isRowOut){
                if (vertices[currHe->tip].col_in[0] != -1){//we're not at the bottom most row
                    //go to the next wale in
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_in[0])];
                    currHe -> isVisited = true;
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
            }
            if (currHe -> isRowIn){
                if (vertices[currHe->tip].col_out[0] != -1){//we're not at the top most row
                    if (vertices[currHe->tip].col_out[1] != -1){//we want to go to the second wale out (some extra work for row ins)
                        currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_out[1])];
                        currHe -> isVisited = true;
                        if (currHe -> tail == startVertex) break;
                        else{
                            continue;
                        }
                    }
                    else{
                        //go to the first wale out
                        currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_out[0])];
                        currHe -> isVisited = true;
                        if (currHe -> tail == startVertex) break;
                        else{
                            continue;
                        }
                    }
                    
                }
            }
            if (currHe -> isWaleIn){
                if ((vertices[currHe->tip].col_out[1] != -1) && (currHe -> tail == vertices[currHe->tip].col_out[1])){//maybe a bit redundant but okay
                    //go to the first wale out 
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_out[0])];
                    currHe -> isVisited = true;
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
                else if (vertices[currHe->tip].row_in != -1){//we're not at a short row 
                    //go to the next row in 
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].row_in)];
                    currHe -> isVisited = true;
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
                else{
                    //hit a short row, take the next wale in 
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_in[0])];
                    currHe -> isVisited = true;
                    //std::cout << currHe->tail << " ";
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
            }
            if (currHe -> isWaleOut){
                if ((vertices[currHe->tip].col_in[1] != -1) && (currHe -> tail == vertices[currHe->tip].col_in[0])){//maybe a bit redundant but okay
                    //go to second wale in 
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_in[1])];
                    currHe -> isVisited = true;
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
                else if (vertices[currHe->tip].row_out != -1){//we're not at a short row
                    //go to the next row out
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].row_out)];
                    currHe -> isVisited = true;
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
                else{
                    //hit a short row, take the next wale out 
                    currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_out[0])];
                    currHe -> isVisited = true;
                    if (currHe -> tail == startVertex) break;
                    else{
                        continue;
                    }
                }
            }
        }while(currHe -> tail != startVertex);
        if (currFace.size() > 1)
            faces.push_back(currFace);
    }
    
    Eigen::MatrixXd primalVertexPositions(vertices.size(), 3);
    for (const auto& v : vertices){
        primalVertexPositions(v.id, 0) = v.position[0];
        primalVertexPositions(v.id, 1) = v.position[1];
        primalVertexPositions(v.id, 2) = v.position[2];
    }

    ManifoldSurfaceMesh * primalMesh = new ManifoldSurfaceMesh(faces);
    VertexPositionGeometry * primalGeometry = new VertexPositionGeometry(*primalMesh, primalVertexPositions);
    // std::cout << "Number of faces in the primal mesh = " << primalMesh->nFaces() << std::endl;
    // std::cout << "Number of vertices in the primal mesh = " << primalMesh->nVertices() << std::endl;
    // std::cout << "Number of edges in the primal mesh = " << primalMesh->nEdges() << std::endl;
    // std::cout << "Vertices in the knit graph = " << vertices.size() << std::endl;
    // std::cout << "Is primalMesh manifold? = " << primalMesh -> isManifold() << std::endl;
    // std::cout << "Is primalMesh edge manifold? = " << primalMesh -> isEdgeManifold() << std::endl;
    // //can't be oriented if it's not manifold
    // std::cout << "Is primalMesh oriented? = " << primalMesh -> isOriented() << std::endl;
    auto psMesh = polyscope::registerSurfaceMesh("primal mesh", primalVertexPositions, faces);
    
    
    //std::vector<Vector3> dualVertexPositions;
    Eigen::MatrixXd dualVertexPositions(primalMesh->nFaces(), 3);
    std::unordered_map<Face, size_t> dualVertexIndices; // Map primal face -> dual vertex index
    for (Face f : primalMesh -> faces()){
        Vector3 avgPosition{0., 0., 0.};
        for (Vertex v : f.adjacentVertices()){
            avgPosition += primalGeometry->vertexPositions[v];
        }
        avgPosition /= f.degree();
        dualVertexIndices[f] = f.getIndex();
        dualVertexPositions(f.getIndex(), 0) = avgPosition[0];
        dualVertexPositions(f.getIndex(), 1) = avgPosition[1];
        dualVertexPositions(f.getIndex(), 2) = avgPosition[2];
    }
    //polyscope::registerPointCloud("dual vertex positions", dualVertexPositions);

    
    std::vector<std::vector<size_t>> dualFaces;
    std::vector<std::vector<int>> edgeLabels;
    for (Vertex v : primalMesh->vertices()){
        if (v.isBoundary()) continue;
        std::vector<size_t> dualFace;
        std::vector<int> currLabel;
        //Traverse the faces incident on this vertex
        Halfedge he = v.halfedge();
        do {
            Face f = he.face();
            if (!f.isBoundaryLoop()) {
                //dualFace.push_back(dualVertexIndices[f]);
                dualFace.push_back(f.getIndex());
                if (vPair2He[std::make_pair(he.tailVertex().getIndex(), he.tipVertex().getIndex())] -> isRowIn
                || vPair2He[std::make_pair(he.tailVertex().getIndex(), he.tipVertex().getIndex())] -> isRowOut){
                    currLabel.push_back(1);
                }
                else if (vPair2He[std::make_pair(he.tailVertex().getIndex(), he.tipVertex().getIndex())] -> isWaleIn){
                    currLabel.push_back(0);
                }
                else if (vPair2He[std::make_pair(he.tailVertex().getIndex(), he.tipVertex().getIndex())] -> isWaleOut)
                    currLabel.push_back(2);
            }
            he = he.twin().next();
        } while (he != v.halfedge());

        // Add this dual face
        if (!dualFace.empty()) {
            std::reverse(dualFace.begin(), dualFace.end());
            dualFaces.push_back(dualFace);
        }
        if (!currLabel.empty()){
            //adjust edge labels to account for correct orientation
            //since we reverse the order of vertices in the face
            //need to carefully handle what happens to the edges
            std::reverse(currLabel.begin(), currLabel.end());
            std::vector<int> adjustedLabels(currLabel.size());
            for (int i = 1; i < currLabel.size(); i++){
                adjustedLabels[i - 1] = currLabel[i];
            }
            adjustedLabels[currLabel.size() - 1] = currLabel[0];
            edgeLabels.push_back(adjustedLabels);
        }
    }

    ManifoldSurfaceMesh * dualMesh = new ManifoldSurfaceMesh(dualFaces);
    VertexPositionGeometry * dualGeometry = new VertexPositionGeometry(*dualMesh, dualVertexPositions);

    //std::vector<int> handled(dualMesh -> nEdges(), 0);
    // for (Face f : dualMesh->faces()){
    //     int i = 0;
    //     std::cout << "for dual face " << f.getIndex() << std::endl;
    //     for (Halfedge he : f.adjacentHalfedges()){
    //         //std::cout << "edge between " << he.tailVertex().getIndex() << " and " << he.tipVertex().getIndex() << " is " << edgeLabels[f.getIndex()][i++] << std::endl;
    //         if ((edgeLabels[f.getIndex()][i++] != 1) && (handled[he.edge().getIndex()] == 0)){//this is a course edge
    //             //half it so that Kui can split wales 
    //             scaleEdgeLength(he.edge(), 1.0, dualGeometry);
    //             handled[he.edge().getIndex()] = 1;

    //         }
    //     }
    //     std::cout << "----------------------------" << std::endl;
    // }

    // std::cout << "number of faces in the dual mesh " << dualMesh -> nFaces() << std::endl;
    // std::cout << "size of edge labels = " << edgeLabels.size() << std::endl;
    polyscope::registerSurfaceMesh("dual mesh after adjustment", dualGeometry->inputVertexPositions, dualMesh->getFaceVertexList());

    std::ofstream outfile("stitchMesh.obj");
    for (const auto &v : dualMesh->vertices()){
        Vector3 p = dualGeometry->vertexPositions[v];
        outfile << "v " << p.x << " " << p.y << " " << p.z << std::endl;
    }
    for (auto &f : dualFaces){
        outfile << "f ";
        for (auto &v : f){
            outfile << v + 1 << " ";
        }
        outfile << "\n";
    }
    for (auto &e : edgeLabels){
        outfile << "e ";
        for (auto &l : e){
            outfile << l << " ";
        }
        outfile << "\n";
    } 

    outfile.close(); 


    //print out percentage errors for histograms
    // double period = 0.25;
    // primalGeometry->requireEdgeLengths();
    // std::ofstream error("errors.csv");
    // for (Edge e : primalMesh->edges()){
    //     double percentError = ((primalGeometry->edgeLengths[e] - period)/period) * 100.;
    //     error << percentError << "\n";
    // }
}


//for saving your graph as a line element obj
void KnitGraph::writeLineElementObj(){

    std::ofstream outfile("lineElement.obj");
    for (auto &v : vertices){
        outfile << " v " << v.position.x << " " << v.position.y << " " << v.position.z << std::endl;
    }

    for (auto &v : vertices){
        if (v.row_out != -1) outfile << "l " << v.id + 1 << " " << v.row_out + 1 << std::endl;
        if (v.col_out[0] != -1) outfile << "l " << v.id + 1  << " " << v.col_out[0] + 1 << std::endl;
        if (v.col_out[1] != -1) outfile << "l " << v.id + 1 << " " << v.col_out[1] + 1 << std::endl;
    }

    outfile.close();

}