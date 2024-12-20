#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "KnitGraph.h"
#include "StitchMesh.h"
#include "imgui.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

int main(int argc, char **argv) {

  // // Initialize polyscope
  polyscope::init();

  
  //read in the knit graph 
  std::string filename = argv[1]; // File path
  std::ifstream file(filename); // Open the file for reading

  if (!file.is_open()) {
    std::cerr << "Error: Could not open the file " << filename << std::endl;
    return 1;
  }

  std::vector<std::vector<double>> matrix; // To store the matrix
  std::string line;

  // Read file line by line
  while (getline(file, line)) {
      std::stringstream ss(line); // Stream to parse line
      std::vector<double> row;    // Temporary row
      double value;

      // Extract values from the line
      while (ss >> value) {
          row.push_back(value);
      }

      matrix.push_back(row); // Add the row to the matrix
  }

  file.close(); // Close the file

  //make a knit graph 
  KnitGraph graph;
  for (int i = 0; i < matrix.size(); i++){
    knitGraphVertex v;
    v.id = matrix[i][0];
    v.position = Vector3{matrix[i][1], matrix[i][2], matrix[i][3]};
    v.row_in = matrix[i][4];
    v.row_out = matrix[i][5];
    v.col_in[0] = matrix[i][6];
    v.col_in[1] = matrix[i][7];
    v.col_out[0] = matrix[i][8];
    v.col_out[1] = matrix[i][9];
    graph.vertices.push_back(v);
  }

  graph.renderGraph();

  graph.traceFaces();

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
