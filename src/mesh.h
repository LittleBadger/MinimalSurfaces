#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <functional>
using namespace glm;

class mesh {
public:
	// the vertices, faces, and edges are enumerated.
	// the first SizeOfBoundary vertices are designated as boundary vertices.
	std::vector<glm::vec3> V; // an array of vertex positions
	std::vector<glm::vec3> N; // an array of vertex normals
	std::vector<float> VA; // an array to be used for storing vertex areas.
	std::vector<int> bc; // an array giving the boundary component of the vertex. -1 for interior vertices.

	std::vector<int> E; // an array where E[2 i], E[2 i + 1] are the vertices of the ith edge
	std::vector<double> L; // an array where L[i] is the length of the ith edge
	std::vector<int> EF; // an array where E[2 i], E[2 i + 1] are the faces adjacent to the ith edge
	std::vector<int> F; // an array of faces, where F[3i], F[3i+1], F[3i+2] are the vertices of the ith face.
	std::vector<int> FE; // an array where FE[3 i], FE[3 i + 1], FE[3 i + 2] are the edges adjacent to the ith face

	int SizeOfInterior;
	
	
	// Various flow/jump algorithms.
	void PPFlip();
	void PP();
	void ElasticFlow(float dt);
	void MeanCurvatureFlow(float dt);

	// These take as input maps f: R^2 -> R^3 either in polar coordinates
	// or cartesian coordinates, and map the mesh by first projecting to the xy plane
	// and then mapping by f.
	void MapDomainPolar(std::function<vec3(double,double)> c);
	void MapDomain(std::function<vec3(double,double)> c);
	
	// This takes as input a map f: [0,1]-> R^3, and maps the c-th boundary component
	// of the mesh by f.
	void MapBoundaryComponent(std::function<vec3(double)> curve, int c);

	// the only difference between these methods is how they recompute the edge metric.
	// the intrinsic sets the flipped edge length to be the intrinsic distance
	// the extrinsic sets the flipped edge length to be the extrinsic distance
	void ExtrinsicDelEdgeFlipAlg();
	void IntrinsicDelEdgeFlipAlg();
	
	void InitRegularRectTriangulation(int L, int W); // Initializes a regular rectangular triangulation of square with L x W vertices.
	
	// In the following, a random triangulation means the Delaunay triangulation of a randomly chosen set of points.
	
	// For Delaunay triangulating a set of points in the xy plane.
	void PlanarDelaunaryTriangulate();
	
	// Initializes a random triangulation of a disk.
	void InitRandDiskTriangulation(int M, int N, long seed); 
	
	// Initializes a random triangulation of a disk of radius r[0] and center centers[0] with the disks r[i] and centers[i] for i > 0 removed.
	void InitRandDiskTriangulation(int M, int N, std::vector<double> radii, std::vector<vec2> centers, long seed);
	
	// This uses the distribution dist to sample the radius of points.
	void InitRandDiskTriangulation(int M, int N, std::vector<double> radii, std::vector<vec2> centers, std::function<double(double)> dist, long seed);

	double TanHalfAngle(int e, int f); // tan of half the angle opposite to edge e in f
	float CotanWeight(int e); // cotan weight of edge e
	float Area(int f); // area of face f
	float SurfaceArea(); // sum of all face areas.
	
	void ComputeEdgeLengths();
	void translateMesh(vec3 v);
	
	void PrintEdgeLengths();
	void PrintTriangleQNumber();

	// returns the edge adjacent to face f and vertices v1 and v2
	int EdgeFace(int v1, int v2, int f); 
	
	// if e = [v1,v2] and belongs to triangle f = [v1,v2,v3], returns v3. 
	int TheOtherVertex(int e, int f);  
	
	// if e = [v1,v2] and belongs to triangle f = [v1,v2,v3], returns edges [v1,v3], [v2,v3] (in order)
	std::pair<int,int> TheOtherEdges(int e, int f); 
	
	// return the index of the (first ) edge adjacent to v1 and v2.
	int FindEdge(int v1, int v2); 
	
	// either finds the edge, or adds the edge [v1,v2].
	int FindOrInsertEdge(int v1, int v2);
	
	// flips the edge e. make sure e is not a boundary edge.
	void EdgeFlip(int e); 
	
	// add a vertex with the coords, and boundary component.
	void Vertex(float x, float y, float z, int bcv = -1);

	// add a face [v1,v2,v3], and the necessary edges. this function adds the edges if needed.
	void Face(int v1, int v2, int v3);
	
	void PrintEdge( int i );
	void PrintFace( int i );
	void PrintVertex(int i);
	void PrintFaceVertex(int i);
	
	void ComputeBCV(); 	//boundary covariance matrix
	glm::mat2 cv; glm::vec3 sh; // for drawing the shadow
};

#endif



