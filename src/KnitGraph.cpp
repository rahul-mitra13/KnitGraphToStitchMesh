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

    std::vector<Vector3> dups;
    int numDuplicate = 0;
    for (int i = 0; i < edges.size(); i++){
        for (int j = 0; j < edges.size(); j++){
            if (i == j) continue;
            if (edges[i][0] == edges[j][0] && edges[i][1] == edges[j][1]){
                dups.push_back(embeddedVertices[edges[i][0]]);
                dups.push_back(embeddedVertices[edges[i][1]]);
                std::cout << "duplicate edges " << std::endl;
                std::cout << "duplicate edge between " << edges[i][0] << " and " << edges[i][1] << std::endl;
                numDuplicate++;
                //exit(1);
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

    if (!isDisk){
   
        //now trace the faces 
        for (auto &he : halfedges){
            if (he->isVisited) continue;//skip halfedges we've already visited
            //need to come up with a way to trace bottom coure, wale, top course, wale
            std::vector<size_t> currFace;
            KnitGraphHalfedge * currHe = he;
            currHe -> isVisited = true;
            int startVertex = currHe->tail;
            do{
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

        std::cout << "Constructing primal mesh " << std::endl;
        SurfaceMesh * primalMesh = new SurfaceMesh(faces);
        VertexPositionGeometry * primalGeometry = new VertexPositionGeometry(*primalMesh, primalVertexPositions);
        std::cout << "Finished constructing primal mesh " << std::endl;
    
        // std::cout << "Number of faces in the primal mesh = " << primalMesh->nFaces() << std::endl;
        // std::cout << "Number of vertices in the primal mesh = " << primalMesh->nVertices() << std::endl;
        // std::cout << "Number of edges in the primal mesh = " << primalMesh->nEdges() << std::endl;
        // std::cout << "Vertices in the knit graph = " << vertices.size() << std::endl;
        // std::cout << "Is primalMesh manifold? = " << primalMesh -> isManifold() << std::endl;
        // std::cout << "Is primalMesh edge manifold? = " << primalMesh -> isEdgeManifold() << std::endl;
        // //can't be oriented if it's not manifold
        // std::cout << "Is primalMesh oriented? = " << primalMesh -> isOriented() << std::endl;
        auto primalPSMesh = polyscope::registerSurfaceMesh("primal mesh", primalVertexPositions, faces);
        FaceData<int> primalMeshFaceLabels(*primalMesh, 0);

        // for (Face f : primalMesh -> faces()){
        //     int ctr = 0;
        //     for (Halfedge he : f.adjacentHalfedges()){
        //         ctr++;
        //     }
        //     if (ctr == 5) primalMeshFaceLabels[f] = -1.0;
        //     else if (ctr == 3) primalMeshFaceLabels[f] = 1.0;
        // }
        // primalPSMesh->addFaceScalarQuantity("face labels", primalMeshFaceLabels);

    
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
        
        std::vector<std::vector<size_t>> dualFaces;
        std::vector<std::vector<int>> edgeLabels;
        for (Vertex v : primalMesh->vertices()){
            if (v.isBoundary()) continue;
            std::vector<size_t> dualFace;
            std::vector<int> currLabel;
            // Traverse the faces incident on this vertex
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

        std::cout << "Constructing dual mesh " << std::endl;
        SurfaceMesh * dualMesh = new SurfaceMesh(dualFaces);
        VertexPositionGeometry * dualGeometry = new VertexPositionGeometry(*dualMesh, dualVertexPositions);
        std::cout << "Finished constructing dual mesh " << std::endl;

        // std::cout << "number of faces in the dual mesh " << dualMesh -> nFaces() << std::endl;
        // std::cout << "size of edge labels = " << edgeLabels.size() << std::endl;
        polyscope::registerSurfaceMesh("dual mesh", dualGeometry->inputVertexPositions, dualMesh->getFaceVertexList());

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
    }

    else{

        //now trace the faces in the disk setting 
        for (auto &he : halfedges){
            if (he->isVisited) continue;//skip halfedges we've already visited
            //need to come up with a way to trace bottom coure, wale, top course, wale
            std::vector<size_t> currFace;
            KnitGraphHalfedge * currHe = he;
            currHe -> isVisited = true;
            currHe -> twin -> isVisited = true;
            int startVertex = currHe->tail;
            std::vector<Vector3> tracedPoints;
            std::cout << "Size of faces traced = " << faces.size() << std::endl;
            std::cout << "start vertex = " << startVertex << std::endl;
            if (he -> isRowOut){
                do{
                    std::cout << "curr vertex = " << currHe->tail << std::endl;
                    currFace.push_back(currHe->tail);
                    tracedPoints.emplace_back(vertices[currHe->tail].position);
                
                    if (currHe -> isRowOut){
                        std::cout << "In row out he " << std::endl;
                        if (vertices[currHe->tip].col_out[0] != -1){//we're not at the top most row
                            //go to the next wale out
                            currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_out[0])];
                            currHe -> isVisited = true;
                            //currHe -> twin -> isVisited = true;
                            if (currHe -> tail == startVertex) break;
                            else{
                                continue;
                            }
                        }
                    }
                    if (currHe -> isRowIn){
                        std::cout << "In row in he " << std::endl;
                        if (vertices[currHe->tip].col_in[0] != -1){//we're not at the bottom most row
                            //go to the next wale in
                            currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_in[0])];
                            currHe -> isVisited = true;
                            //currHe -> twin -> isVisited = true;
                            if (currHe -> tail == startVertex) break;
                            else{
                                continue;
                            }
                        }     
                    }
                    if (currHe -> isWaleIn){
                        std::cout << "In wale in he " << std::endl;
                        if (vertices[currHe->tip].row_out != -1){//we're not at a short row 
                            //go to the next row out
                            currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].row_out)];
                            currHe -> isVisited = true;
                            //currHe -> twin -> isVisited = true;
                            if (currHe -> tail == startVertex) break;
                            else{
                                continue;
                            }
                        }
                        else{
                            //if the next wale in exists, take it, otherwise break
                            if (vertices[currHe->tip].col_in[0] != -1){
                                //hit a short row, take the next wale in 
                                currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_in[0])];
                                currHe -> isVisited = true;
                                //currHe -> twin -> isVisited = true;
                                if (currHe -> tail == startVertex) break;
                                else{
                                    continue;
                                }
                            }
                            else{
                                break;
                            }
                            if (currHe -> tail == startVertex) break;
                            else{
                                continue;;
                            }
                        }
                    }
                    if (currHe -> isWaleOut){
                        std::cout << "In wale out he " << std::endl;
                        if (vertices[currHe->tip].row_in != -1){//we're not at a short row
                            //go to the next row in
                            currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].row_in)];
                            currHe -> isVisited = true;
                            //currHe -> twin -> isVisited = true;
                            if (currHe -> tail == startVertex) break;
                            else{
                                continue;
                            }
                        }
                        else{
                            //if the next wale out exists, take it, otherwise break
                            if (vertices[currHe->tip].col_out[0] != -1){
                                //hit a short row, take the next wale out 
                                currHe = vPair2He[std::make_pair(currHe->tip, vertices[currHe->tip].col_out[0])];
                                currHe -> isVisited = true;
                                //currHe -> twin -> isVisited = true;
                            }
                            else{
                                break;
                            }
                            if (currHe -> tail == startVertex) break;
                            else{
                                continue;
                            }
                        }
                    }
                }while(currHe -> tail != startVertex);
                if (currFace.size() > 3)
                    faces.push_back(currFace);
            }
        }
    
        Eigen::MatrixXd primalVertexPositions(vertices.size(), 3);
        for (const auto& v : vertices){
            primalVertexPositions(v.id, 0) = v.position[0];
            primalVertexPositions(v.id, 1) = v.position[1];
            primalVertexPositions(v.id, 2) = v.position[2];
        }

        // ============================================================
        // Remove unreferenced vertices and remap faces 
        // ============================================================

        // Collect referenced vertices
        std::unordered_set<size_t> referenced;
        for (const auto& f : faces) {
            for (size_t v : f) {
                referenced.insert(v);
            }
        }

        // Build old → new index map
        std::unordered_map<size_t, size_t> oldToNew;
        std::vector<size_t> newToOld;
        newToOld.reserve(referenced.size());

        for (size_t v : referenced) {
            size_t newIdx = newToOld.size();
            oldToNew[v] = newIdx;
            newToOld.push_back(v);
        }

        // Remap faces
        std::vector<std::vector<size_t>> compactFaces;
        compactFaces.reserve(faces.size());

        for (const auto& f : faces) {
            std::vector<size_t> nf;
            nf.reserve(f.size());
            for (size_t v : f) {
                nf.push_back(oldToNew[v]);
            }
            compactFaces.push_back(nf);
        }

        // Build compact vertex positions
        Eigen::MatrixXd compactPositions(newToOld.size(), 3);
        for (size_t i = 0; i < newToOld.size(); ++i) {
            const auto& v = vertices[newToOld[i]];
            compactPositions(i, 0) = v.position[0];
            compactPositions(i, 1) = v.position[1];
            compactPositions(i, 2) = v.position[2];
        }

        // sanity logging
        std::cout << "Compacted vertices: " << newToOld.size()
          << " / " << vertices.size() << std::endl;
        std::cout << "Faces: " << compactFaces.size() << std::endl;

        std::cout << "Constructing primal mesh " << std::endl;
        SurfaceMesh* primalMesh = new SurfaceMesh(compactFaces);
        VertexPositionGeometry* primalGeometry =
            new VertexPositionGeometry(*primalMesh, compactPositions);
        std::cout << "Finished constructing primal mesh " << std::endl;
        auto primalPSMesh = polyscope::registerSurfaceMesh("primal mesh", compactPositions, compactFaces);
        FaceData<int> primalMeshFaceLabels(*primalMesh, 0);

        std::ofstream outfile("primalMesh.obj");

        // ------------------ vertices ------------------
        for (int i = 0; i < compactPositions.rows(); ++i) {
            outfile << "v "
                << compactPositions(i,0) << " "
                << compactPositions(i,1) << " "
                << compactPositions(i,2) << "\n";
            }

        // ------------------ faces ------------------
        for (const auto& f : compactFaces) {
            outfile << "f ";
            for (size_t v : f) {
            // OBJ is 1-indexed, and we want v/vt
            outfile << (v + 1) << "/" << (v + 1) << " ";
            }
            outfile << "\n";
        }
        outfile.close();
    }

}

//trace the short rows in the graph 
void KnitGraph::traceShortRows(){   

    int ctr = 0;
    for (auto &v : vertices){
        if (v.row_in == -1){
            std::vector<Vector3> pos;
            std::vector<std::array<int, 2>> edges;
            auto walker = v;
            while(walker.row_out != -1){
                pos.push_back(walker.position);
                walker = vertices[walker.row_out];
            }
            for (int i = 0; i < (int)pos.size() - 1; i++){
                edges.push_back(std::array<int, 2>{i, i + 1});
            }
            if (pos.size() > 1000) return;
            polyscope::registerCurveNetwork("traced short row " + std::to_string(ctr), pos, edges)->setRadius(0.003);
            ctr++;
        }
    }
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

//---------------------------------------------//
// Implementing DAG for integer UV coordinates
//---------------------------------------------//

//Collect rows
std::vector<Chain> collectRows(const std::vector<knitGraphVertex>& vertices) {
    std::vector<Chain> rows;
    std::vector<char> seen(vertices.size(), false);

    for (const auto& v : vertices) {
        if (v.row_in == -1 && v.row_out != -1) {          // row start
        Chain row;
        row.start = v.id;

        int cur = v.id;
        std::unordered_set<int> guard;            
        while (cur != -1 && !seen[cur]) {
            if (guard.count(cur)) break;
            guard.insert(cur);

            seen[cur] = true;
            row.verts.push_back(cur);

            int nxt = vertices[cur].row_out;
            cur = (nxt != -1 ? nxt : -1);
        }

        // only keep nontrivial rows
        if (!row.verts.empty()) rows.push_back(std::move(row));
        }
    }

    std::cout << "Number of row chains " << rows.size() << std::endl;
    return rows;
}

//Collect columns 
std::vector<Chain> collectCols(const std::vector<knitGraphVertex>& vertices) {
  
    std::vector<Chain> cols;
    std::vector<char> seen(vertices.size(), false);

    for (const auto& v : vertices) {
        if (v.col_in[0] == -1 && v.col_out[0] != -1) {    // col start
        Chain col;
        col.start = v.id;

        int cur = v.id;
        std::unordered_set<int> guard;
        while (cur != -1 && !seen[cur]) {
            if (guard.count(cur)) break;
            guard.insert(cur);

            seen[cur] = true;
            col.verts.push_back(cur);

            int nxt = vertices[cur].col_out[0];
            cur = (nxt != -1 ? nxt : -1);
        }

        if (!col.verts.empty()) cols.push_back(std::move(col));
        }
    }

    std::cout << "Number of col chains " << cols.size() << std::endl;
    return cols;
}

// Build vertex->chainId map from chains (rows or cols)
std::vector<int> buildVertexToChainMap(
    int nVerts,
    const std::vector<Chain>& chains
) {
  std::vector<int> v2c(nVerts, -1);
  for (int cid = 0; cid < (int)chains.size(); ++cid) {
    for (int v : chains[cid].verts) {
      if (v < 0 || v >= nVerts) continue;
      v2c[v] = cid;
    }
  }
  return v2c;
}

// Row DAG edges: use wale links v -> col_out[0]
DAG buildRowDAG(const std::vector<knitGraphVertex>& vertices,
                const std::vector<Chain>& rows,
                const std::vector<int>& v2row) {

  DAG g;
  g.n = (int)rows.size();
  g.adj.assign(g.n, {});
  g.indeg.assign(g.n, 0);

  std::unordered_set<std::pair<int,int>, PairHash> seen;

  for (int v = 0; v < (int)vertices.size(); ++v) {
    int r0 = v2row[v];
    if (r0 < 0) continue;

    int w = vertices[v].col_out[0];
    if (w == -1) continue;

    if (w < 0 || w >= (int)vertices.size()) continue;
    int r1 = v2row[w];
    if (r1 < 0) continue;

    if (r0 == r1) continue; // ignore self

    std::pair<int,int> e{r0, r1};
    if (seen.insert(e).second) {
      g.adj[r0].push_back(r1);
      g.indeg[r1] += 1;
    }
  }

  return g;
}

// Column DAG edges: use course links v -> row_out
DAG buildColDAG(const std::vector<knitGraphVertex>& vertices,
                const std::vector<Chain>& cols,
                const std::vector<int>& v2col) {

  DAG g;
  g.n = (int)cols.size();
  g.adj.assign(g.n, {});
  g.indeg.assign(g.n, 0);

  std::unordered_set<std::pair<int,int>, PairHash> seen;

  for (int v = 0; v < (int)vertices.size(); ++v) {
    int c0 = v2col[v];
    if (c0 < 0) continue;

    int u = vertices[v].row_out;
    if (u == -1) continue;

    if (u < 0 || u >= (int)vertices.size()) continue;
    int c1 = v2col[u];
    if (c1 < 0) continue;

    if (c0 == c1) continue;

    std::pair<int,int> e{c0, c1};
    if (seen.insert(e).second) {
      g.adj[c0].push_back(c1);
      g.indeg[c1] += 1;
    }
  }

  return g;
}

//Kahn topo sort (will throw if cycle)
std::vector<int> topoSort(const DAG& g) {
  std::queue<int> q;
  std::vector<int> indeg = g.indeg;
  for (int i = 0; i < g.n; ++i) if (indeg[i] == 0) q.push(i);

  std::vector<int> order;
  order.reserve(g.n);

  while (!q.empty()) {
    int u = q.front(); q.pop();
    order.push_back(u);
    for (int v : g.adj[u]) {
      if (--indeg[v] == 0) q.push(v);
    }
  }

  if ((int)order.size() != g.n) {
    throw std::runtime_error("DAG topoSort failed: cycle detected (or missing nodes).");
  }
  return order;
}