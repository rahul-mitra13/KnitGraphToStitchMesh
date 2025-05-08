This code computes the stich mesh (such as those used by yarn-level rendering pipelines) by taking the dual of a knit graph. 

This knitgraph is stored in a txt file with the following format:

`vID vPos.x vPos.y vPos.z v.rowIn v.rowOut v.colIn[0] v.colIn[1] v.colOut[0] v.colOut[1]`

The output is an obj file with appropriate halfedge labels denoted by e. 

Usage

`./bin/stitch_mesh [knitgraph.txt]`
