#version 330 core

uniform mat2 invCov;
uniform vec3 shift;

in vData {
    vec3 pos;
    vec3 norm;
} frag;

out vec4 color;

void main()
{	
	vec2 pos2 = frag.pos.xz - shift.xz;
	float dis = dot(pos2,invCov*pos2);
	color = vec4(1.0f, 1.0f, 1.0f, 1.0f) - exp(-.25*dis-1.2-shift.y)*vec4(1.0f,1.0f,0.9f,0.0f);
}