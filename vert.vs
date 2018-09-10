#version 330 core

uniform mat4 viewTransform;
uniform mat4 modelTransform;

layout (location = 0) in vec4 vPosition;
layout (location = 1) in vec4 vNormal;

out vData {
    vec3 pos;
    vec3 norm;
} vert;

void main()
{   
	vert.pos = vPosition.xyz;
	vert.norm = vNormal.xyz;
    gl_Position = (viewTransform*vPosition);
}
