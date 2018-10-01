# Minimal Surfaces
<img align="right" width="300" src="https://github.com/LittleBadger/MinimalSurfaces/blob/master/PSurface.png">

Here is a small program for visualizing and processing triangulated surfaces with the intrinsic discrete Laplace-Beltrami operator introduced by Bobenko and Springborn, see [1,2]. This operator is defined via the intrinsic Delaunay triangulation of the surface, and consequently enjoys many remarkable properties analogous to smooth Laplacians.

This program focuses in particular on the construction of discrete minimal surfaces. It implements a mean curvature flow based on the intrinsic Laplace-Beltrami operator, along with a version of the Pinkall-Polthier algorithm (see [3]) that allows for the changing of combinatorics. Some other standard techniques in mesh processing are also implemented for comparison.

The program also includes some utility functions for constructing random triangulations and discrete surfaces described by immersion functions.







## Building and Running
In order to build the program, you will first need the following:
* [OpenGL](https://www.opengl.org/)
* [GLEW](http://glew.sourceforge.net)
* [GLFW](http://www.glfw.org)
* [Eigen](http://eigen.tuxfamily.org/)




Once you are sure these are installed, you can build using GCC, for example, by executing:

```
g++ -std=gnu++11 -o minimal_surfaces src\minimal_surfaces.cpp src\mesh.cpp src\example_surfaces.cpp -lglfw3 -lgdi32 -lglew32  -lopengl32 -IC:\eigen-eigen\
```

## Usage


### Interface
````
Left mouse button drag   : orbit camera
Right mouse button drag  : dolly camera
f                        : retriangulate using the combinatorics of the intrinsic Delaunay triangulation
g                        : a time step of flow to minimize the elastic energy
h                        : a time step of mean curvature flow using the standard cotan Laplacian
j                        : a time step of mean curvature flow using the intrinsic Laplacian.
k                        : an iteration of the Pinkall-Polthier algorithm (using the standard cotan Laplacian)
l                        : an iteration of the Pinkall-Polthier algorithm, with Delaunay retriangulation
````

### Constructing Surfaces

Meshes can be specified in many ways. See example_surfaces.cpp for examples and details. For convenience, the program has functions for creating random triangulations of a disk (or a disk with holes), which can then be immersed in R^3 by a specified function. By random triangulation, we mean the Delaunay triangulation of a random set of points in the disk. 

The Schwarz surface above, for example, was created by applying the Pinkall-Polthier algorithm to a surface constructed by first randomly triangulating a disk with 5 holes and then mapping stereographically to sphere. The points are sampled on the disk according to the pushforward of the uniform distribution of the sphere.

<p align = "center">
<img width="600" src="https://github.com/LittleBadger/MinimalSurfaces/blob/master/mapping.png"> 
</p>


## Acknowledgements
This joint work with A.I. Bobenko is supported by the DFG Collaborative Research Center TRR 109 “Discretization in Geometry and Dynamics.”

## References
[1] A.I. Bobenko and B. Springborn. "A discrete Laplace-Beltrami operator for simplicial surfaces,"	Discrete Comput. Geom. 38:4 (2007) 740-756.\
[2] M. Fisher, B. Springborn, P. Schröder, A. I. Bobenko. "An algorithm for the construction of intrinsic delaunay triangulations with applications to digital geometry processing." ACM SIGGRAPH 2006 Courses, 69-74.\
[3] U. Pinkall, K. Polthier. Experiment. Math., Volume 2, Issue 1 (1993), 15-36.
