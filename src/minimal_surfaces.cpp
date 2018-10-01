#define GLEW_STATIC

#include <windows.h>
#include <GL/glew.h>
#include <iostream>
#include <GLFW/glfw3.h>
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

// Camera angles, and scale
float a1 = 3.14/4.0, a2 = -3.14/2.0f +1.1;
float sf = .5;

// GUI stuff
GLFWwindow* window;
double mx, my;
bool lbutton_down, rbutton_down; 

// The surface
mesh m1;

// GPU objects and buffers
GLuint VBO, NBO, VAO, EBOF, EBOE;
GLuint GVAO, GVBO, GEBO;
GLuint SVAO, SVBO, SEBO;

// Shaders
GLuint shaderProgram;
Shader *mainShader, *backShader, *lineShader, *normalShader, *groundShader, *skyShader;

void display()
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 viewMat = glm::translate( glm::mat4(1.0), glm::vec3(0,0,-.2))*glm::rotate( glm::rotate( glm::mat4(1.0), a2, glm::vec3(1,0,0)), a1, glm::vec3(0,1,0));
	viewMat = glm::scale( viewMat,  glm::vec3(sf));
	glm::mat4 inviewMat = glm::inverse(viewMat);
	
	/* --------- draw the surface --------------------------------------------------- */

	mainShader->Use();
	glUniformMatrix4fv(glGetUniformLocation(mainShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOF);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glPolygonOffset( 1, 1 );
	glDrawElements(GL_TRIANGLES, m1.F.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
	
	/* --------- draw the back --------------------------------------------------- */
	glCullFace(GL_FRONT);
	backShader->Use();
	glUniformMatrix4fv(glGetUniformLocation(mainShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOF);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glPolygonOffset( 1, 1 );
	glDrawElements(GL_TRIANGLES, m1.F.size(), GL_UNSIGNED_INT, 0);
	glCullFace(GL_BACK);
	glBindVertexArray(0);
	
	
	
	
	/* --------- draw the ground plane and shadows ---------------------------------- */
	groundShader->Use();
	m1.ComputeBCV();
	glm::mat2 icv = glm::inverse(m1.cv);
	glUniformMatrix4fv(glGetUniformLocation(groundShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glUniformMatrix2fv(glGetUniformLocation(groundShader->Program, "invCov"), 1, false, &icv[0][0]);
	glUniform3fv(glGetUniformLocation(groundShader->Program, "shift"), 1, &m1.sh[0]);


	glBindVertexArray(GVAO);
	glBindBuffer(GL_ARRAY_BUFFER, GVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GEBO);
	glPolygonMode( GL_FRONT, GL_FILL );
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	/* ------- draw triangle edges --------------------------------------------------- */
	lineShader->Use();
	glUniformMatrix4fv(glGetUniformLocation(lineShader->Program, "viewTransform"), 1, false, &viewMat[0][0]);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOE);
	glDrawElements(GL_LINES, m1.E.size(), GL_UNSIGNED_INT, 0);

	glBindVertexArray(0);
	glFlush();
}

void BufferMesh() {
	glBindVertexArray(VAO);

	/* ------------------ Buffer Vertex Data --------------------------------------- */
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, m1.V.size()*(sizeof(vec3)), &(m1.V[0].x), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (GLvoid*)0);
	
	/* ------------------- Buffer Normal Data ---------------------------------------- */
	/*
	glBindBuffer(GL_ARRAY_BUFFER, NBO);
	glBufferData(GL_ARRAY_BUFFER, m1.N.size()*(sizeof(vec3)), &(m1.N[0].x), GL_STATIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (GLvoid*)0);*/

	/* ------------------- Element buffer ojbects for the faces ---------------------- */
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOF);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m1.F.size()*sizeof(int), &(m1.F[0]), GL_STATIC_DRAW);

	/* ------------------- Element buffer ojbects for the edges ---------------------- */
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOE);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m1.E.size()*sizeof(int), &(m1.E[0]), GL_STATIC_DRAW);
}


void init() {
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH_CLAMP);
	glEnable(GL_CULL_FACE);
	glDepthFunc(GL_LEQUAL);

	glDepthRange(0,100);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_LINE_SMOOTH);
	
	mainShader = new Shader("src/shaders/vert.vs","src/shaders/frag.frag");
	backShader = new Shader("src/shaders/vert.vs","src/shaders/fragBack.frag");
	lineShader = new Shader("src/shaders/vert.vs","src/shaders/fragL.frag");
	//skyShader = new Shader("src/shaders/vert.vs","src/shaders/fragSky.frag");
	groundShader = new Shader("src/shaders/vert.vs","src/shaders/fragGr.frag");
	//normalShader = new Shader("vert.vs","geom.geo","fragLB.frag");
	

	// ------------------- Ground ---------------------------------------------------- //
	glGenVertexArrays(1, &GVAO);
	glBindVertexArray(GVAO);
	glGenBuffers(1, &GVBO);
	glGenBuffers(1, &GEBO);

	std::vector<float> ground({-15,-1.1,-15,   15,-1.1,-15,    15,-1.1,15,        -15,-1.1,15});
	std::vector<int> gebo({0,1,2,0,2,3});

	glBindBuffer(GL_ARRAY_BUFFER, GVBO);
	glBufferData(GL_ARRAY_BUFFER, ground.size()*(sizeof(float)), &(ground[0]), GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GEBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, gebo.size()*(sizeof(int)), &(gebo[0]), GL_STATIC_DRAW);

	// ------------------- Buffers for the surface --------- //
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &NBO);
	glGenBuffers(1, &EBOF);
	glGenBuffers(1, &EBOE);

	// Initialize Mesh
	psurface(m1);
	BufferMesh();

	glBindVertexArray(0);
}


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if ( key == 'P' ) { m1.translateMesh(vec3(0,.1f,0)); BufferMesh(); }
	if ( key == ';' ) { m1.translateMesh(vec3(0,-.1f,0)); BufferMesh(); }
	if ( key == 'L' ) { m1.translateMesh(vec3(-.1f,0,0)); BufferMesh(); }
	if ( key == '\'') { m1.translateMesh(vec3(.1f,0,0)); BufferMesh(); }
	if ( key == '.' ) { m1.translateMesh(vec3(0,0,-.1f)); BufferMesh(); }
	if ( key == '/' ) { m1.translateMesh(vec3(0,0,.1f)); BufferMesh(); }

	
	if ( key == 'F' ) {
		m1.ExtrinsicDelEdgeFlipAlg();
		BufferMesh();
	}
	if ( key == 'H' ) {
		m1.MeanCurvatureFlow(.01);
		BufferMesh();
	}
	if ( key == 'G' ) {
		m1.ElasticFlow(.01);
		BufferMesh();
	}
	if ( key == 'J' ) {
		m1.PP();
		BufferMesh();
	}
	if ( key == 'K' ) {
		m1.PPFlip();
		BufferMesh();
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	glfwGetCursorPos(window,&mx,&my);
	
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if(GLFW_PRESS == action)
            lbutton_down = true;
        else if(GLFW_RELEASE == action)
            lbutton_down = false;
    }
	
	if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if(GLFW_PRESS == action)
            rbutton_down = true;
        else if(GLFW_RELEASE == action)
            rbutton_down = false;
    }
}

void mousedragfun() {
	double x = 0; double y = 0;
	
	glfwGetCursorPos(window,&x,&y);
	if ( lbutton_down ) {
		a2 -= .01*(y-my); a1 -= .01*(x-mx);
	} else {
		sf *= 1-.001*(y-my);
	}
	
	mx = x; my = y;
}


int main(int argc, char** argv) {

	glewExperimental = true; // Needed for core profile
	glfwInit();
	
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); 
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


	window = glfwCreateWindow( 700, 700, "Minimal Surfaces", NULL, NULL);
	glfwMakeContextCurrent(window);
	
	glewExperimental = GL_TRUE;
	glewInit();
	
	//glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	init();
	
	do{

		glClear( GL_COLOR_BUFFER_BIT );

		display();
		
		glfwSwapBuffers(window);
		glfwPollEvents();
		
		if ( lbutton_down || rbutton_down ) mousedragfun();

	} while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 );
		   
	return 0;
}

