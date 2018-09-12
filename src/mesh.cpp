#include "mesh.h"
#include <vector>
#include <climits>
#include <iostream>
#include <cstdlib>
#include <random>
#include <Eigen/Dense>
#include <unordered_map>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;

double eps = .0001;
double eps2 = .1;

// ------------------------- Embedding Functionality ---------------//
void mesh::MapDomainPolar(std::function<vec3(double,double)> c) {

	for (int i = 0; i < V.size(); ++i) {

		double r = glm::length(V[i]);
		double t = glm::atan(V[i].x,V[i].z);

		V[i] = c(r,t);
	}
}

void mesh::MapDomain(std::function<vec3(double,double)> c) {
	for (int i = 0; i < V.size(); ++i) {
		V[i] = c(V[i].x,V[i].z);
	}
}

void mesh::MapBoundaryComponent(std::function<vec3(double)> curve, int c) {

	int i = 0;
	while ( i < V.size() && bc[i] != c) {
		++i;
	}
	int s = i;
	while ( i < V.size() && bc[i] == c)
		++i;

	int e = i;
	int M = e-s;

	for (i = s; i < e; ++i) {
		double t = 1.0*(s-i) / M;
		V[i] = curve(t);
	}
}

// ------------------------- Pinkall-Polthier Jump algorithms ------------------//
void mesh::PP() {
	int nv = SizeOfInterior;

	MatrixXd lap(nv,nv);
	VectorXd ppx(nv); VectorXd ppy(nv); VectorXd ppz(nv);

	for (int i = 0; i < nv; ++i) {
		for (int j = 0; j < nv; ++j) {
			lap(i,j) = 0;
		}
		ppx(i) = 0; ppy(i) = 0; ppz(i) = 0;
	}

	ComputeEdgeLengths();

	for (int i = 0; i < E.size()/2; ++i) {
		int v1 = E[2*i]; int v2 = E[2*i+1];

		double w = CotanWeight(i);

		if ( v1 < SizeOfInterior && v2 < SizeOfInterior ) { lap(v1,v1) += w; lap(v1,v2) -= w; lap(v2,v1) -= w; lap(v2,v2) += w; }
		else if ( v1 < SizeOfInterior ) {
			lap(v1,v1) += w;
			ppx(v1) += w*V[v2].x; ppy(v1) += w*V[v2].y; ppz(v1) += w*V[v2].z;
		}
		else if ( v2 < SizeOfInterior ) {
			lap(v2,v2) += w;
			ppx(v2) += w*V[v1].x; ppy(v2) += w*V[v1].y; ppz(v2) += w*V[v1].z;
		}
	}

	VectorXd hx = lap.llt().solve(ppx); VectorXd hy = lap.llt().solve(ppy); VectorXd hz = lap.llt().solve(ppz);

	for (int i = 0; i < nv; ++i) {
		V[i].x = hx[i]; V[i].y = hy[i]; V[i].z = hz[i];
	}
}

void mesh::PPFlip() {

	int nv = SizeOfInterior;

	MatrixXd lap(nv,nv);

	VectorXd ppx(nv); VectorXd ppy(nv); VectorXd ppz(nv);

	for (int i = 0; i < nv; ++i) {
		for (int j = 0; j < nv; ++j) {
			lap(i,j) = 0;
		}
		ppx(i) = 0; ppy(i) = 0; ppz(i) = 0;
	}

	IntrinsicDelEdgeFlipAlg();

	for (int i = 0; i < E.size()/2; ++i) {
		int v1 = E[2*i];
		int v2 = E[2*i+1];

		double w = CotanWeight(i);

		if ( v1 < SizeOfInterior && v2 < SizeOfInterior ) { lap(v1,v1) += w; lap(v1,v2) -= w; lap(v2,v1) -= w; lap(v2,v2) += w; }
		else if ( v1 < SizeOfInterior ) {
			lap(v1,v1) += w;
			ppx(v1) += w*V[v2].x; ppy(v1) += w*V[v2].y;	ppz(v1) += w*V[v2].z;
		}
		else if ( v2 < SizeOfInterior ) {
			lap(v2,v2) += w;
			ppx(v2) += w*V[v1].x; ppy(v2) += w*V[v1].y; ppz(v2) += w*V[v1].z;
		}
	}

	VectorXd hx = lap.llt().solve(ppx);	VectorXd hy = lap.llt().solve(ppy); VectorXd hz = lap.llt().solve(ppz);

	for (int i = 0; i < nv; ++i) {
		V[i].x = hx[i]; V[i].y = hy[i]; V[i].z = hz[i];
	}
}

// ------------------------------- Flows ---------------------------------------//

// Delaunay triangulates, then minimizes spring energy
void mesh::ElasticFlow(float dt) {
	//ExtrinsicDelEdgeFlipAlg();
	//ComputeEdgeLengths();

	for (int i = 0; i < N.size(); ++i)
		N[i] = vec3(0,0,0);


	for (int i = 0; i < E.size()/2; ++i) {
		vec3 p = V[E[2*i]]-V[E[2*i+1]];
		N[E[2*i]] -= p;
		N[E[2*i+1]] += p;
	}

	for (int i = 0; i < SizeOfInterior; ++i) {
		V[i] += (N[i])*dt;
	}

	for (int i = SizeOfInterior; i < V.size(); ++i)
		N[i] = vec3(0,0,0);
}

// Retriangulates so that the triangulation is extrinsically delaunay.
// Then computes the laplacian and flows.
void mesh::MeanCurvatureFlow(float dt) {
	ExtrinsicDelEdgeFlipAlg();
	ComputeEdgeLengths();
	
	for (int i = 0; i < N.size(); ++i) {
		N[i] = vec3(0,0,0);
		VA[i] = 0;
	}

	for (int i = 0; i < E.size()/2; ++i) {
		float w = CotanWeight(i);
		N[E[2*i]] += (V[E[2*i+1]]  - V[E[2*i]])*w;
		N[E[2*i+1]] -= (V[E[2*i+1]]  - V[E[2*i]])*w;
		float EA = 0;
		if ( EF[2*i] > 0 ) EA += Area(EF[2*i]);
		if ( EF[2*i+1] > 0 ) EA += Area(EF[2*i+1]);
	
		VA[E[2*i]] += EA;
	}

	for (int i = 0; i < SizeOfInterior; ++i) {
		V[i] += (N[i])*dt/VA[i];
	}
}

// ------------------------- Mesh Metrics -----------------------//

float mesh::Area(int i) {
	float l1 = L[FE[3*i]]; float l2 = L[FE[3*i+1]];	float l3 = L[FE[3*i+2]];

	float s = (l1+l2+l3)/2;
	return std::sqrt(s*(s-l1)*(s-l2)*(s-l3));
}

float mesh::SurfaceArea() {

	ComputeEdgeLengths();

	float A = 0;

	for (int i = 0; i < F.size()/3; ++i) {
		A += Area(i);
	}

	return A;
}

void mesh::PrintTriangleQNumber() {
	for (int i = 0; i < F.size()/3; ++i) {

		float l1 = L[FE[3*i]];
		float l2 = L[FE[3*i+1]];
		float l3 = L[FE[3*i+2]];

		float M = std::max(std::max(l1,l2),l3);
		float m = std::min(std::min(l1,l2),l3);

		std::cout<<M/m<<",";
	}

	std::cout<<std::endl;

}

void mesh::PrintEdgeLengths() {
	ComputeEdgeLengths();
	for (int i = 0; i < E.size()/2; ++i) {
		std::cout<<L[i]<<",";
	}

	std::cout<<std::endl;
}

//-------------------- Cotan Weights ---------------------------------- //

// Computes the cotan weights based on the stored length ( which may differ from the
// extrinsic lengths in the case that there have been intrinsic edge flips ).
float mesh::CotanWeight(int e) {
	double w = 0;
	int c = 0;
	
	for (int i = 0; i < 2; ++i) {
		int f = EF[2*e + i];
		if ( f != -1 ) {
			c+= 1;
			double tanHA = TanHalfAngle(e,f);
			w += (1 - tanHA*tanHA)/ ( 2*tanHA );
		}

	}
	return w/c;
}

double mesh::TanHalfAngle(int e, int f) {
	std::pair<int,int> oe = TheOtherEdges(e,f);
	double l1 = L[e]; double l2 = L[oe.first]; double l3 = L[oe.second];
	double da = (l1-l2 + l3)*(l1+l2-l3)/((l1+l2+l3)*(-l1+l2+l3));
	if ( da < 0 ) {
		return 0;
	}
	return std::sqrt(da);
}

void mesh::ComputeEdgeLengths() {
	for (int i = 0; i < E.size()/2; ++i) {
		L[i] = length(V[E[2*i]] - V[E[2*i+1]]);
	}
}

// -------------------- Delaunay Edge Flip ---------------------------- //
void mesh::IntrinsicDelEdgeFlipAlg() {
	bool allDelaunay = false;

	ComputeEdgeLengths();

	
	// the following should be implemented with a queue,
	// but at the moment it is implemented with a brute force
	// loop for flexibility reasons.
	
	while (!allDelaunay) {
		allDelaunay = true;
		for (int i = 0; i < E.size()/2; ++i) {
			if ( EF[2*i] != -1 && EF[2*i+1] != - 1) // not a boundary face
			{
				if ( CotanWeight(i) < -eps ) // not delaunay, problem
				{
					allDelaunay = false;
					EdgeFlip(i);
				}
			}
		}
	}
}

void mesh::ExtrinsicDelEdgeFlipAlg() {
	bool allDelaunay = false;

	ComputeEdgeLengths();

	while (!allDelaunay) {

		allDelaunay = true;

		for (int i = 0; i < E.size()/2; ++i) {
			if ( EF[2*i] != -1 && EF[2*i+1] != - 1) // not a boundary face
			{
				if ( CotanWeight(i) < -eps ) // not delaunay, problem
				{
					allDelaunay = false;
					float sa1 = SurfaceArea();
					EdgeFlip(i);
					float sa2 = SurfaceArea();
					L[i] = length(V[E[2*i]]-V[E[2*i+1]]);
					std::cout<<"Flip! Area diff: "<<sa1 - sa2 <<std::endl;
				}
			}
		}
	}
}

// -------------------- Constructing Initial Meshes ------------------- //
void mesh::InitRandDiskTriangulation(int M, int N, long seed) {
	std::vector<double> radii; radii.push_back(1);
	std::vector<vec2> centers; centers.push_back(vec2(0,0));
	InitRandDiskTriangulation(M,N,radii,centers,seed);
}

void mesh::InitRandDiskTriangulation(int M, int N, std::vector<double> radii, std::vector<vec2> centers, long seed) {
	double r = radii[0];
	std::function<double(double)> dist = [r](double x) { return r*sqrt(x); };
	InitRandDiskTriangulation(M,N,radii,centers,dist,seed);
}

// a reminder: if you want to sample with distribution rho, map the uniform U[0,1] by F(x)
// where F(x) is the inverse of the cumulative of rho.
// for example for distribution 2 r, the cumulative is r^2, so F(x) = sqrt(x)
//
// for the stereographic metric, you should use F[x] = Sqrt[x]/Sqrt[1 - x]
void mesh::InitRandDiskTriangulation(int M, int N, std::vector<double> radii, std::vector<vec2> centers, std::function<double(double)> dist, long seed) {
	srand(76123429);
	
	SizeOfInterior = N;

	int k = 0;
	while ( k < SizeOfInterior ) {

		double r0 = radii[0];

		float th = rand()/(1.0*RAND_MAX)*2*3.14591;
		float r = dist((float)(rand()/(1.0*RAND_MAX)));

		vec2 v( r*cos(th), r*sin(th));

		bool accept = true;

		accept = length(v-centers[0]) < radii[0]-eps2;

		for (int i = 1; i < radii.size(); ++i)
			accept = accept & (length(v-centers[i]) > radii[i]+eps2);
		if ( accept ) {
			k++;
			Vertex(v.x,0,v.y);
		}
	}

	for (int q = 0; q < centers.size(); q++) {
		for (int i = 0; i < M; ++i) {
			float th = 1.0*i/M*2.0*3.14591;
			double r = radii[q];
			vec2 c = centers[q];
			Vertex( r*cos(th)+c.x, 0, r*sin(th)+c.y, q );
		}
	}

	PlanarDelaunaryTriangulate();
}


void mesh::PlanarDelaunaryTriangulate() {
	bool delaunay;

	for (int i = 0; i < V.size(); ++i) {
		for (int j = i+1; j < V.size(); ++j) {
			for (int k = j+1; k < V.size(); ++k) {


				delaunay = true;

				float Ax = V[i].x; float Ay = V[i].z;
				float Bx = V[j].x; float By = V[j].z;
				float Cx = V[k].x; float Cy = V[k].z;

				float s2 = (Bx-Ax)*(Cy-Ay)-(By-Ay)*(Cx-Ax);

				if ( abs(s2) > eps ) {
					for (int m = 0; m < V.size(); m++) {

						float Dx = V[m].x; float Dy = V[m].z;

						float BD = (Bx-Dx)*(Bx-Dx)+(By-Dy)*(By-Dy);
						float AD = (Ax-Dx)*(Ax-Dx)+(Ay-Dy)*(Ay-Dy);
						float CD = (Cx-Dx)*(Cx-Dx)+(Cy-Dy)*(Cy-Dy);
						float s1 = (Ax-Dx)*((By-Dy)*CD- (Cy-Dy)*BD)
																													 - (Bx-Dx)*((Ay-Dy)*CD- (Cy-Dy)*AD)
																													 + (Cx-Dx)*((Ay-Dy)*BD- (By-Dy)*AD);

						if ( s1*s2 > 0) {
							delaunay = false;
							break;
						}
					}
					if (delaunay && !(bc[i] == bc[j] && bc[j] == bc[k] && bc[i] > 0 ) ) {
						if ( true )
						{
							//PrintVertex(i);
							//PrintVertex(j);
							//PrintVertex(k);

							//std::cout<<std::endl;
						}
						if ( s2 < 0 ) Face(i,j,k);
						else Face(j,i,k); 
					}
				}
			}
		}
	}
}


//-------------------- Combinatorial Stuff ---------------------------//
void mesh::InitRegularRectTriangulation(int L, int W) {
	// This mesh is a bit annoying because the boundary vertices should be last
	// in the array of vertices. This makes hooking up the triangles annoying.
	// For simplicity, I just hash the vertex indices.
	std::unordered_map<int,int> m;
	int k = 0;
	
	
	SizeOfInterior = (L-2)*(W-2);
	// add (L-2) x (W-2) vertices, labeling those on the boundary as boundary vertices.
	for (int i = 1; i < L-1; ++i) {
		for (int j = 1; j < W-1; ++j) {
			Vertex(2.0*i/(L-1) - 1, 0, 2.0*j/(W-1) - 1, 0);
			m[i*W+j] = k++;
		}
	}
	
	// add the boundary vertices in clockwise order
	for (int i = 0; i < L-1; ++i) {
		int j = 0;
			Vertex(2.0*i/(L-1) - 1, 0, 2.0*j/(W-1) - 1, 0);
		m[i*W+j] = k++;
	}
	for (int j = 0; j < W-1; ++j) {
		int i = L-1;
			Vertex(2.0*i/(L-1) - 1, 0, 2.0*j/(W-1) - 1, 0);
		m[i*W+j] = k++;
	}
	for (int i = L-1; i > 0; i--) {
		int j = W-1;
			Vertex(2.0*i/(L-1) - 1, 0, 2.0*j/(W-1) - 1, 0);
		m[i*W+j] = k++;
	}
	for (int j = W-1; j > 0; j--) {
		int i = 0;
			Vertex(2.0*i/(L-1) - 1, 0, 2.0*j/(W-1) - 1, 0);
		m[i*W+j] = k++;
	}
	
	// add the faces.
	for (int i = 0; i < L-1; ++i) {
		for (int j = 0; j < W-1; ++j) {
			Face(m[i*W+j], m[i*W+j+1], m[i*W+j+W]);
			Face(m[i*W+j+W], m[i*W+j+1], m[i*W+j+W+1]);
		}
	}
}

int mesh::EdgeFace(int v1, int v2, int f) {
	for (int i = 0; i < 3; ++i) {
		int k = FE[3*f + i];
		if ((E[2*k] == v1 && E[2*k+1] == v2) || (E[2*k] == v2 && E[2*k+1] == v1) ) return k; }

	return 0;
}


int mesh::TheOtherVertex(int e, int f) {
	int v1 = F[3*f]; int v2 = F[3*f+1]; int v3 = F[3*f+2];
	int w1 = E[2*e]; int w2 = E[2*e+1];

	if ( (w1 == v1 && w2 == v2) || (w1 == v2 && w2 == v1) ) return v3;
	if ( (w1 == v1 && w2 == v3) || (w1 == v3 && w2 == v1) ) return v2;
	if ( (w1 == v2 && w2 == v3) || (w1 == v3 && w2 == v2) ) return v1;

	return 0;
}


std::pair<int,int> mesh::TheOtherEdges(int e, int f) {

	std::pair<int,int> oe;

	int e1 = FE[3*f]; int e2 = FE[3*f+1]; int e3 = FE[3*f+2];

	if ( e1 == e ) {
		if ( E[2*e2] == E[2*e] || E[2*e2+1] == E[2*e]) { oe.first = e2; oe.second = e3; return oe; }
		else 										   { oe.first = e3; oe.second = e2; return oe; }
	}
	if ( e2 == e ) {
		if ( E[2*e1] == E[2*e] || E[2*e1+1] == E[2*e]) { oe.first = e1; oe.second = e3; return oe; }
		else 								  		   { oe.first = e3; oe.second = e1; return oe; }
	}
	if ( e3 == e ) {
		if ( E[2*e1] == E[2*e] || E[2*e1+1] == E[2*e]) { oe.first = e1; oe.second = e2; return oe; }
		else 										   { oe.first = e2; oe.second = e1; return oe; }
	}
	
	return oe;
}

int mesh::FindEdge(int v1, int v2) {
	for (size_t e = 0; e < E.size()/2.0; e++) {
		if ( (E[2*e] == v1 && E[2*e+1] == v2)  || (E[2*e] == v2 && E[2*e+1] == v1) )
			return e;
	}
	return -1;
}

void mesh::EdgeFlip(int e) {

	int f1 = EF[2*e];
	int f2 = EF[2*e+1];

	int A = 0; int B = 0; int C = E[2*e]; int D = E[2*e+1];

	// f1 is A C D
	A = TheOtherVertex(e, f1);
	// f2 is B C D
	B = TheOtherVertex(e, f2);

	int AD, BC, AC, BD;
	AD = EdgeFace(A,D,f1);
	BC = EdgeFace(B,C,f2);
	AC = EdgeFace(A,C,f1);
	BD = EdgeFace(B,D,f2);

	//gamma is the angle BD
	double thg = TanHalfAngle(BD,f2);
	//delta is the angle AD
	double thd = TanHalfAngle(AD,f1);

	// tan of (gamma + delta)/2
	double thgd = (thg + thd)/(1-thd*thg);
	double cosgd = (1 - thgd*thgd)/(1 + thgd*thgd);
	double l2 = L[BC]*L[BC] + L[AC]*L[AC] - 2*L[BC]*L[AC]*cosgd;

	E[2*e] = A; E[2*e+1] = B;

	L[e] = std::sqrt(l2);

	if ( EF[2*AD] == f1 ) EF[2*AD] = f2; else if ( EF[2*AD+1] == f1 ) EF[2*AD+1] = f2;
	if ( EF[2*BC] == f2 ) EF[2*BC] = f1; else if ( EF[2*BC+1] == f2 ) EF[2*BC+1] = f1;

	FE[3*f1 + 0] = e; FE[3*f1+1] = AC; FE[3*f1+2] = BC;
	FE[3*f2 + 0] = e; FE[3*f2+1] = AD; FE[3*f2+2] = BD;
	
	// Keep track of orientation. *!! This is should be done better but a quick fix for now !!*
	/*
	F[3*f1 + 0] = A; F[3*f1 + 1] = C; F[3*f1 + 2] = B;
	F[3*f2 + 0] = A; F[3*f2 + 1] = B; F[3*f2 + 2] = D;
	*/
	
	if ( F[3*f1 + 0] == D ) { F[3*f1 + 0] = B; }
	else if ( F[3*f1 + 1] == D ) { F[3*f1 + 1] = B; }
	else if ( F[3*f1 + 2] == D ) { F[3*f1 + 2] = B; }
	
	if ( F[3*f2 + 0] == C ) { F[3*f2 + 0] = A; }
	else if ( F[3*f2 + 1] == C ) { F[3*f2 + 1] = A; }
	else if ( F[3*f2 + 2] == C ) { F[3*f2 + 2] = A; }

}


void mesh::Vertex(float x, float y, float z, int bcv) {
	V.push_back(vec3(x,y,z));
	N.push_back(vec3(0,0,0));
	bc.push_back(-1);
	VA.push_back(0);
	bc[bc.size()-1] = bcv;
}

void mesh::Face(int v1, int  v2, int v3) {
	int n = (F.size() / 3.0);
	F.push_back(v1); F.push_back(v2); F.push_back(v3);

	int e;
	e = FindOrInsertEdge(v1,v2); if ( EF[2*e] == -1 ) EF[2*e] = n; else EF[2*e+1] = n; FE.push_back(e);
	e = FindOrInsertEdge(v2,v3); if ( EF[2*e] == -1 ) EF[2*e] = n; else EF[2*e+1] = n; FE.push_back(e);
	e = FindOrInsertEdge(v3,v1); if ( EF[2*e] == -1 ) EF[2*e] = n; else EF[2*e+1] = n; FE.push_back(e);
}

int mesh::FindOrInsertEdge(int v1, int v2) {
	int e = FindEdge(v1,v2);
	if ( e != -1 ) { return e; }
	else { E.push_back(v1); E.push_back(v2); EF.push_back(-1); EF.push_back(-1); L.push_back(0); return (E.size()/2 - 1); };
	return e;
}

void mesh::translateMesh(vec3 t) {
	for (int i = 0; i < V.size(); ++i)
		V[i] += t;
}


void mesh::ComputeBCV() {
	cv[0][0] = 0; cv[0][1] = 0; cv[1][0] = 0; cv[1][1] = 0;
	sh = glm::vec3(0,V[0].y,0);

	for (int i = 0; i < V.size(); ++i) {
		cv[0][0] += V[i].x*V[i].x;
		cv[1][0] += V[i].x*V[i].z;
		cv[0][1] += V[i].x*V[i].z;
		cv[1][1] += V[i].z*V[i].z;

		sh.x += V[i].x;
		sh.y = std::min(sh.y,V[i].y);
		sh.z += V[i].z;
	}

	cv *= (1.0)/(V.size());
	sh.x *= (1.0)/(V.size());
	sh.z *= (1.0)/(V.size());

	cv[0][0] -= sh.x*sh.x;
	cv[1][0] -= sh.x*sh.z;
	cv[0][1] -= sh.x*sh.z;
	cv[1][1] -= sh.z*sh.z;
}