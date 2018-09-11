#version 330 core


in vData {
    vec3 pos;
    vec3 norm;
} frag;


out vec4 color;

void main()
{
	vec3 sv = normalize( (frag.pos+vec3(0.0f,1.0f,0.0f)) );

//	color = (1-max(1.0f,min(0.0f,sv.y)) )*vec4(0.8f,0.8f,0.8f,1.0f)+ max(0.0f,sv.y)*vec4(0.8f,0.8f,0.9f,1.0f);;
	//color = vec4(0.8f,0.8f,0.9f,1.0f);
	color = vec4(1.0f,1.0f,1.0f,1.0f);
}