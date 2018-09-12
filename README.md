# Minimal Surfaces
<img align="right" width="300" src="https://github.com/LittleBadger/MinimalSurfaces/blob/master/PSurface.png">

Here is a small program for visualizing and processing triangulated surfaces with the intrinsic discrete Laplace-Beltrami operator introduced by Bobenko and Springborn, see [1,2]. The operator differs from usual combinatorial Laplacians in that it is defined via the intrinsic Delaunay triangulation, which can differ considerably from the immersed triangulation.

In particular, we focushe program implements a mean curvature flow based on the intrinsic Laplace-Beltrami operator, along with a version of the Pinkall-Polthier algorithm, see [], that allows for changing combinatorics of the surface. For comparison, various other flows and algorithms for mesh processing are also implemented. 







## Building and Running
In order to build the program, you will first need the following:
* OpenGL
* [GLEW](http://glew.sourceforge.net)
* [GLFW](http://www.glfw.org)
* [Eigen](http://eigen.tuxfamily.org/)




Once you are sure these are installed, you can build by running for example using GCC:

```
g++ -std=gnu++11 -o minimal_surfaces src\minimal_surfaces.cpp src\mesh.cpp src\example_surfaces.cpp -lglfw3 -lgdi32 -lglew32  -lopengl32 -IC:\eigen-eigen\
```

## Usage


### Interface
````
Left mouse button drag   : orbit camera
Right mouse button drag  : dolly camera
f                        : compute the intrinsic Delaunay triangulation
g                        : flow to minimize the elastic energy
h                        : a time step of mean curvature flow
j                        : a time step of the Pinkall-Polthier algorithm
k                        : a time step of the Pinkall-Polthier algorithm, with Delaunay retriangulation
````

### Constructing Surfaces




The program includes some utility functions for creating random triangulations of a disk or disk with holes, that can then be immersed in R^3. Please see example_surfaces.cpp for examples and details. The Schwarz surface above, for example, was created by applying the Pinkall-Polthier algorithm to a randomly triangulated disk with 5 holes, mapped to stereographically sphere. 
<img align="center width="300" src="https://github.com/LittleBadger/MinimalSurfaces/blob/master/diskholes.png">
<img align="center width="300" src="https://github.com/LittleBadger/MinimalSurfaces/blob/master/sphere.png">




## References
[1] A.I. Bobenko and B. Springborn. "A discrete Laplace-Beltrami operator for simplicial surfaces,"	Discrete Comput. Geom. 38:4 (2007) 740-756\
[2] "An algorithm for the construction of intrinsic delaunay triangulations with applications to digital geometry processing"\
[3] Pinkall Polthier.
