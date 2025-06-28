#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "KnitGraph.h"
#include "imgui.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

float currAngle = 0.0f;

void rotateObjects() {
  currAngle += 0.007f; // radians per frame

  // Set the view rotation (about Y axis)
  polyscope::view::lookAt(
  glm::vec3(cos(currAngle) * 3.0f, 0.5f, sin(currAngle) * 3.0f), // camera position
  glm::vec3(0.0f, 0.0f, 0.0f),                           // target
  glm::vec3(0.0f, 1.0f, 0.0f)                            // up vector
  );
}

int main(int argc, char **argv) {

  // // Initialize polyscope
  polyscope::init();

  
  //read in the knit graph 
  std::string knitGraphName = argv[1]; // File path
  std::ifstream file(knitGraphName); // Open the file for reading

  if (!file.is_open()) {
    std::cerr << "Error: Could not open the file " << knitGraphName << std::endl;
    return 1;
  }

  std::string meshName = argv[2];
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = readManifoldSurfaceMesh(meshName);
  auto psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(meshName), geometry->inputVertexPositions, mesh -> getFaceVertexList());

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

  //graph.traceFaces();

  graph.traceShortRows();

  graph.writeLineElementObj();

  //for rendering figures
  //has to go before the show() call
  // std::string viewerString = R"({"farClipRatio":20.0,"fov":45.0,"nearClipRatio":0.005,"projectionMode":"Perspective","viewMat":[-0.916695833206177,1.74622982740402e-09,0.399592250585556,-0.924567401409149,0.0970502495765686,0.970054030418396,0.222643107175827,26.9699611663818,-0.387625247240067,0.24287411570549,-0.889246106147766,-205.221466064453,0.0,0.0,0.0,1.0],"windowHeight":1440,"windowWidth":2560} )";
  // polyscope::view::setViewFromJson(viewerString, false);

  // Set the callback
  polyscope::state::userCallback = rotateObjects;

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
