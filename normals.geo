#version 330 core

layout (points) in;
layout (line_strip,  max_vertices=2) out;

in vData {
    vec3 pos;
    vec3 norm;
} vertices[];

out vData {
    vec3 pos;
    vec3 norm;
} frag;

void main()
{
    frag.norm = vertices[0].norm;
    frag.pos = vertices[0].pos;
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    frag.norm = vertices[0].norm;
    frag.pos = vertices[0].pos + vertices[0].norm;
    gl_Position = gl_in[0].gl_Position + vec4(.5*vertices[0].norm,0);
    EmitVertex();

    EndPrimitive();
}