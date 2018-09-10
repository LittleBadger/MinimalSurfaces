#define GLEW_STATIC

#include <windows.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include "Shader.h"
#include <vector>
#include <Eigen/Dense>
#include <functional>
#include "example_surfaces.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include "mesh.h"

#define pi 3.141592654f

using namespace Eigen;
using namespace glm;
using namespace std;
float dt = .01;

float a1 = pi/4.0, a2 = -pi/2.0f +1.1;
float view_rad = -.2, sf = 1;

int mx0, my0;
bool left_mouse; 

void display();

mesh m1;

GLuint VBO, NBO, VAO, EBOF, EBOE;
GLuint shaderProgram;
Shader *mainShader, *lineShader, *normalShader, *groundShader, *skyShader;

GLuint GVAO, GVBO, GEBO;
GLuint SVAO, SVBO, SEBO;

void BufferMesh();

float rot = 0;

// two volumes on minimal surfaces (springer bonn), hildebrandt
// house better triangulation
// for the doubly winded disk, what is known in the smooth setting

void SetupMesh() {

	cylinder(m1);
	BufferMesh();
	return;
}

void display()
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	

	glm::mat4 viewMat = glm::translate( glm::mat4(1.0), glm::vec3(0,0,view_rad))*glm::rotate( glm::rotate( glm::mat4(1.0), a2, glm::vec3(1,0,0)), a1, glm::vec3(0,1,0));
	viewMat = glm::scale( viewMat,  glm::vec3(sf));
	glm::mat4 inviewMat = glm::inverse(viewMat);
	// draw triangles

	
	/* --------- draw the surface --------------------------------------------------- */
	glDisable(GL_CULL_FACE);
	mainShader->Use();
	glUniformMatrix4fv(glGetUniformLocation(mainShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOF);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glPolygonOffset( 1, 1 );
	glDrawElements(GL_TRIANGLES, m1.F.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	/* ------------------------------------------------------------------------------- */

	
	/* -------- draw a sky. turn this off for a black sky, or edit skyshader.frag --- */
	/*
	skyShader->Use();
	glUniformMatrix4fv(glGetUniformLocation(groundShader->Program, "viewTransform"), 1, false, &inviewMat[0][0]);

	glBindVertexArray(SVAO);
	glBindBuffer(GL_ARRAY_BUFFER, SVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SEBO);
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	*/
	/* ------------------------------------------------------------------------------- */
	
	/* --------- draw the ground plane and shadows ---------------------------------- */
	groundShader->Use();
	m1.ComputeBCV();
	glm::mat2 icv = glm::inverse(m1.cv);
	glUniformMatrix4fv(glGetUniformLocation(groundShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glUniformMatrix2fv(glGetUniformLocation(groundShader->Program, "invCov"), 1, false, &icv[0][0]);
	glUniform3fv(glGetUniformLocation(groundShader->Program, "shift"), 1, &m1.sh[0]);

	glEnable(GL_CULL_FACE);
	glBindVertexArray(GVAO);
	glBindBuffer(GL_ARRAY_BUFFER, GVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GEBO);
	glPolygonMode( GL_FRONT, GL_FILL );
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	/* ------------------------------------------------------------------------------- */


	/* ------- draw triangle edges --------------------------------------------------- */
	lineShader->Use();
	glUniformMatrix4fv(glGetUniformLocation(lineShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOE);
	glDrawElements(GL_LINES, m1.E.size(), GL_UNSIGNED_INT, 0);
	/* ------------------------------------------------------------------------------- */
	
	glBindVertexArray(0);;
	glFlush();
}

void BufferMesh() {
	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, m1.V.size()*(sizeof(vec3)), &(m1.V[0].x), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (GLvoid*)0);

	glBindBuffer(GL_ARRAY_BUFFER, NBO);
	glBufferData(GL_ARRAY_BUFFER, m1.N.size()*(sizeof(vec3)), &(m1.N[0].x), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (GLvoid*)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOF);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m1.F.size()*sizeof(int), &(m1.F[0]), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOE);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m1.E.size()*sizeof(int), &(m1.E[0]), GL_STATIC_DRAW);
}

void keyfun(unsigned char key, int x, int y) {
	if ( key == 'p' ) { m1.translateMesh(vec3(0,.1f,0)); BufferMesh(); }
	if ( key == ';' ) { m1.translateMesh(vec3(0,-.1f,0)); BufferMesh(); }
	if ( key == 'l' ) { m1.translateMesh(vec3(-.1f,0,0)); BufferMesh(); }
	if ( key == '\'') { m1.translateMesh(vec3(.1f,0,0)); BufferMesh(); }
	if ( key == '.' ) { m1.translateMesh(vec3(0,0,-.1f)); BufferMesh(); }
	if ( key == '/' ) { m1.translateMesh(vec3(0,0,.1f)); BufferMesh(); }

	
	if ( key == 'f' ) {
		m1.ExtrinsicDelEdgeFlipAlg();
		BufferMesh();
	}
	if ( key == 'h' ) {
		m1.MeanCurvatureFlow(.01);
		BufferMesh();
	}
	if ( key == 'g' ) {
		m1.ElasticFlow(.01);
		BufferMesh();
	}
	if ( key == 'j' ) {
		m1.PP();
		BufferMesh();
	}
	if ( key == 'k' ) {
		m1.PPFlip();
		BufferMesh();
	}
}

void mousedragfun(int mx, int my) {
	if ( left_mouse ) {
		a2 -= .01*(my-my0);
		a1 -= .01*(mx-mx0);
	} else {
		sf *= 1-.001*(my-my0);
	}

	mx0 = mx;
	my0 = my;
}

void mousefun(int button, int state, int x, int y) {
	mx0 = x; my0 = y; 
	left_mouse = (button == GLUT_LEFT_BUTTON);
}

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB | GLUT_MULTISAMPLE);
	glutInitWindowSize(700, 700);
	glutCreateWindow(argv[0]);

	glewExperimental = GL_TRUE;
	glewInit();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH_CLAMP);
	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LEQUAL);

	glDepthRange(0,100);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth(1.5);

	mainShader = new Shader("vert.vs","frag.frag");
	lineShader = new Shader("vert.vs","fragL.frag");
	skyShader = new Shader("vertSky.vs","fragSky.frag");
	groundShader = new Shader("vert.vs","fragGr.frag");
	//normalShader = new Shader("vert.vs","geom.geo","fragLB.frag");

	/* ground */
	glGenVertexArrays(1, &GVAO);
	glBindVertexArray(GVAO);
	glGenBuffers(1, &GVBO);
	glGenBuffers(1, &GEBO);

	vector<float> ground({-5,-1.1,-5,   5,-1.1,-5,    5,-1.1,5,        -5,-1.1,5});
	vector<int> gebo({0,1,2,0,2,3});

	glBindBuffer(GL_ARRAY_BUFFER, GVBO);
	glBufferData(GL_ARRAY_BUFFER, ground.size()*(sizeof(float)), &(ground[0]), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GEBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, gebo.size()*(sizeof(int)), &(gebo[0]), GL_STATIC_DRAW);

	/*sky background*/

	glGenVertexArrays(1, &SVAO);
	glBindVertexArray(SVAO);
	glGenBuffers(1, &SVBO);
	glGenBuffers(1, &SEBO);

	vector<float> sky({-1,-1,5,      -1,1,5,  1,1,5,       1,-1,5});
	vector<int> sebo({0,1,2,0,2,3});

	glBindBuffer(GL_ARRAY_BUFFER, SVBO);
	glBufferData(GL_ARRAY_BUFFER, sky.size()*(sizeof(float)), &(sky[0]), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, SEBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sebo.size()*(sizeof(int)), &(sebo[0]), GL_STATIC_DRAW);

	/* actual mesh */
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &NBO);
	glGenBuffers(1, &EBOF);
	glGenBuffers(1, &EBOE);

	SetupMesh();

	glBindVertexArray(0);

	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutKeyboardFunc(keyfun);
	glutMouseFunc(mousefun);
	
	glutMotionFunc(mousedragfun);
	glutMainLoop();

	return 0;
}

