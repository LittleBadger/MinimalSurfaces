#version 330 core


in vData {
    vec3 pos;
    vec3 norm;
} frag;


out vec4 color;

void main()
{
	vec3 normal = normalize( -cross(dFdx(frag.pos), dFdy(frag.pos)) );
	
	float col = .3*normal.x + normal.y + .1*normal.z;
	
	if (col > 0) {
		color = vec4(.3*col+.7,.3*col+.7,.5*col+.8,1.0f);
	} else {
		col = abs(col);
		color = vec4(.2*col+.7,.2*col+.7,.5*col+.7,1.0f);
	}
}