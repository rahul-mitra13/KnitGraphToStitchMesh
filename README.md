This code computes the stich mesh (such as those used by yarn-level rendering pipelines) by taking the dual of a knit graph. Also, used to quickly view my graphs

This knitgraph is stored in a txt file with the following format:

`vID vPos.x vPos.y vPos.z v.rowIn v.rowOut v.colIn[0] v.colIn[1] v.colOut[0] v.colOut[1]`

The output is an obj file with appropriate halfedge labels denoted by e. 

# Build instructions and usage

```
git clone --recursive https://github.com/rahul-mitra13/KnitGraphToStitchMesh.git
cd KnitGraphToStitchMesh 
mkdir build; cd build 
cmake ../ 
make -j
```

Then to run the executable use `./bin/stitch_mesh [knitgraph.txt]`

The polysope GUI shows the knit graph, the primal mesh (generated from the knit graph), and its dual. It also writes an obj file where the stitch mesh is saved.

