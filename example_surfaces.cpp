#include "example_surfaces.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <functional>


void graph(mesh &m1) {
	m1.InitRandDiskTriangulation(20,30,423452);

	std::function<vec3(double, double)> randFunc = [](double x, double y) {
		return vec3(x, 1.0*rand()/(1.0*RAND_MAX) , y);
	};

	m1.MapDomain(randFunc);

	std::function<vec3(double)> embedding = [](double t) {

		float th = 2 * pi * t;
		return vec3(cos(th), .9*cos(2*th) , -sin(th));
	};

	m1.MapBoundaryComponent(embedding,0);
}

void enneper(mesh &m1) {
	m1.InitRandDiskTriangulation(100,200, 1234132);

		std::function<vec3(double)> embedding = [](double t) {
			double r = 1.9;
			float th = 2 * pi * t;
			return vec3( cos(th) - 1/3.0f*r*r*cos(3*th), -sin(th)-1/3.0f*r*r*sin(3*th), r*cos(2*th));
		};

		m1.MapBoundaryComponent(embedding,0);
}

void annulusHV(mesh &m1) {
	std::vector<double> radii; std::vector<vec2> centers;
	radii.push_back(1); centers.push_back(vec2(0,0));
	radii.push_back(.5); centers.push_back(vec2(0,0));
	m1.InitRandDiskTriangulation(20,60,radii,centers, 1234234);
	std::function<vec3(double)> embedding = [](double t) {
		return vec3( 0, .8*cos(2*pi*t), -.7*sin(2*pi*t));
	};

	m1.MapBoundaryComponent(embedding,1);

}

void psurface(mesh &m1) {
	std::vector<double> radii;
	std::vector<vec2> centers;

	double r0 = 4.79129; double r1 = 0.436436; double r2 = 0.208712;
	double c = 1.09109;
	radii.push_back(r0); centers.push_back(vec2(0,0));
	radii.push_back(r1); centers.push_back(vec2(c,0));
	radii.push_back(r1); centers.push_back(vec2(-c,0));
	radii.push_back(r1); centers.push_back(vec2(0,-c));
	radii.push_back(r1); centers.push_back(vec2(0,c));
	radii.push_back(r2); centers.push_back(vec2(0,0));

	std::function<double(double)> stdist = [](double x) {
		return sqrt(x)/sqrt(1-x);
	};
	m1.InitRandDiskTriangulation(20,100,radii,centers,stdist,34534);

	std::function<vec3(double,double)> stereographic = [](double x, double y) {
		return vec3(2*x/(1 + x*x+y*y),(-1 + x*x + y*y)/(1+x*x+y*y),2*y/(1+x*x+y*y));
	};

	m1.MapDomain(stereographic);
}

void house(mesh &m1) {

	std::function<vec3(double,double)> embedding = [](double r, double t) {
		if ( t > -pi/4.0 && t < pi/4.0 )
			return vec3(r ,.0, r*tan(t));
		if ( t > pi/4.0 && t < 3*pi/4.0 )
			return vec3(r/tan(t), 0, r);
		if ( t > 3*pi/4.0 || t < -3*pi/4.0 )
			return vec3(-r,0, -r*tan(t));
		if ( t > -3*pi/4.0 || t < -pi/4.0 )
			return vec3(-r/tan(t),0, -r);
	};

	std::function<vec3(double,double)> embedding2 = [](double x, double y) {
		return vec3(x,3-3*abs(y),y);
	};

	m1.InitRandDiskTriangulation(48,200,123492);
	m1.MapDomainPolar(embedding);
	m1.MapDomain(embedding2);
}

void box(mesh &m1) {

	std::function<vec3(double,double)> embedding = [](double r, double t) {
		if ( t > -pi/4.0 && t < pi/4.0 )
			return vec3(r ,.0, r*tan(t));
		if ( t > pi/4.0 && t < 3*pi/4.0 )
			return vec3(r/tan(t), 0, r);
		if ( t > 3*pi/4.0 || t < -3*pi/4.0 )
			return vec3(-r,0, -r*tan(t));
		if ( t > -3*pi/4.0 || t < -pi/4.0 )
			return vec3(-r/tan(t),0, -r);
	};

	std::function<vec3(double,double)> embedding2 = [](double x, double y) {
		if ( y > .4)
			return vec3(.4*x,1-y,.4);
		if ( y > -.4)
			return vec3(.4*x,.6,y);

		return vec3(.4*x,1+y,-.4);
	};

	m1.InitRandDiskTriangulation(48,50,123492);
	m1.MapDomainPolar(embedding);
	m1.MapDomain(embedding2);
}

void drop(mesh &m1) {
	std::function<vec3(double,double)> embedding = [](double r, double t) {
		return vec3(2*(r-.75), .6*cos(t), .6*sin(t));
			};
	std::vector<double> radii; std::vector<vec2> centers;
	radii.push_back(1); centers.push_back(vec2(0,0));
	radii.push_back(.5); centers.push_back(vec2(0,0));
	m1.InitRandDiskTriangulation(10,10,radii,centers,32434354);
	m1.MapDomainPolar(embedding);
}

void pinch(mesh &m1) {
	m1.InitRandDiskTriangulation(60,120,42323452);
	std::function<vec3(double,double)> embedding = [](double r, double t) {
		return vec3(r*cos(t), 1*exp(-8*r*r), r*sin(t));
	};
	m1.MapDomainPolar(embedding);
}

void cylinder(mesh &m1) {
	std::function<vec3(double,double)> embedding = [](double r, double t) {
		return vec3(1.3*(r-.75), .6*cos(t), .6*sin(t));
			};
	std::vector<double> radii; std::vector<vec2> centers;
	radii.push_back(1); centers.push_back(vec2(0,0));
	radii.push_back(.5); centers.push_back(vec2(0,0));
	m1.InitRandDiskTriangulation(50,100,radii,centers,32434354);
	m1.MapDomainPolar(embedding);
}

void doublyWindedCurve(mesh &m1) {
	std::function<vec3(double,double)> doubleEmb = [](double r, double t) {
		return vec3(r*cos(2*t), .5*sin(t),r*sin(2*t));
		};

	m1.InitRandDiskTriangulation(22,40,123424);
	m1.MapDomainPolar(doubleEmb);
}


void RegularBox(mesh &m1) {
	m1.InitRegularRectTriangulation(10, 10);
	
	std::function<vec3(double,double)> embedding2 = [](double x, double y) {
		if ( y > .201)
			return vec3(.4*x,1-y,.4);
		if ( y > -.201)
			return vec3(.4*x,.6,y);

		return vec3(.4*x,1+y,-.4);
	};

	m1.MapDomain(embedding2);
}
