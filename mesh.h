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
	std::vector<float> VA; // an array to be used vertex areas.
	std::vector<int> bc; // an array giving the boundary component of the vertex

	std::vector<int> E; // an array where E[2 i], E[2 i + 1] are the vertices of the ith edge
	std::vector<double> L; // an array where L[i] is the length of the ith edge
	std::vector<int> EF; // an array where E[2 i], E[2 i + 1] are the faces adjacent to the ith edge
	std::vector<int> F; // an array of faces, where F[3i], F[3i+1], F[3i+2] are the vertices of the ith face.
	std::vector<int> FE; // an array where FE[3 i], FE[3 i + 1], FE[3 i + 2] are the edges adjacent to the ith face

	void ComputeBCV(); 	//boundary covariance matrix
	glm::mat2 cv; glm::vec3 sh; // for drawing the shadow

	int SizeOfInterior;

	void PrintEdgeLengths();
	void PrintTriangleQNumber();
	float Area(int f);
	float SurfaceArea();
	void PPFlip();
	void PP();
	
	void ElasticFlow(float dt);
	void MeanCurvatureFlow(float dt);

	void MapDomainPolar(std::function<vec3(double,double)> c);
	void MapDomain(std::function<vec3(double,double)> c);
	void MapBoundaryComponent(std::function<vec3(double)> curve, int c);

	// the only difference between these methods is how they recompute the edge metric.
	// the intrinsic sets the flipped edge length to be the intrinsic distance
	// the extrinsic sets the flipped edge length to be the extrinsic distance
	void ExtrinsicDelEdgeFlipAlg();
	void IntrinsicDelEdgeFlipAlg();

	void InitRegularRectTriangulation(int L, int W);
	void InitRandDiskTriangulation(int M, int N, long seed);
	void InitRandDiskTriangulation(int M, int N, std::vector<double> radii, std::vector<vec2> centers, long seed);
	void InitRandDiskTriangulation(int M, int N, std::vector<double> radii, std::vector<vec2> centers, std::function<double(double)> dist, long seed);

	double TanHalfAngle(int e, int f);
	float CotanWeight(int i);
	void ComputeEdgeLengths();

	void PlanarDelaunaryTriangulate();

	void translateMesh(vec3 v);

	int EdgeFace(int v1, int v2, int f);
	void EdgeFlip(int e);
	int TheOtherVertex(int e, int f);
	std::pair<int,int> TheOtherEdges(int e, int f);
	int FindEdge(int v1, int v2);
	int FindOrInsertEdge(int v1, int v2);
	void PrintEdge( int i );
	void PrintFace( int i );
	void PrintVertex(int i);
	void PrintFaceVertex(int i);

	void Vertex(float x, float y, float z);
	void Vertex(float x, float y, float z, int bcv);

	void Face(int x, int y, int z);

};

#endif



