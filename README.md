# Minimal Surfaces

Here is a small program for visualizing and processing triangulated surfaces with the intrinsic discrete Laplace-Beltrami operator introduced by Bobenko and Springborn, see [1,2]. This operator differs from usual combinatorial Laplacians in that it is defined via the intrinsic Delaunay triangulation, which can differ considerably from the immersed triangulation. 

In particular, the program implements a mean curvature flow based on the intrinsic Laplace-Beltrami operator, along with a version of the Pinkall-Polthier algorithm, see [], that allows for changing combinatorics of the surface. For comparison, various other flows and algorithms for mesh processing are also implemented.

## Building and Running
In order to run the program, you will need the following libraries.
* [GLEW](http://glew.sourceforge.net)
* [GLFW](http://www.glfw.org)
* [Eigen](http://eigen.tuxfamily.org/)

Then compile and run minimal_surfaces.cpp.

## Usage
````
Left mouse button drag   : orbit camera
Right mouse button drag  : dolly camera
f                        : compute the intrinsic Delaunay triangulation
g                        : flow to minimize the elastic energy
h                        : a time step of mean curvature flow
j                        : a time step of the Pinkall-Polthier algorithm
k                        : a time step of the Pinkall-Polthier algorithm, with Delaunay retriangulation
````

For building your own initial meshes, the program includes some helper methods to generate random triangulations of a planar disk with holes, immerse them into R^3 with arbitrary functions. See example_surfaces.cpp.

## References
[1] A.I. Bobenko and B. Springborn. "A discrete Laplace-Beltrami operator for simplicial surfaces,"	Discrete Comput. Geom. 38:4 (2007) 740-756
[2] "An algorithm for the construction of intrinsic delaunay triangulations with applications to digital geometry processing"
[3] Pinkall Polthier.
