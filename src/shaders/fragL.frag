#version 330 core



in vData {
    vec3 pos;
    vec3 norm;
} frag;

out vec4 color;

void main()
{
	color = vec4(.3f, 0.6f, .9f, 1.0f);	
}