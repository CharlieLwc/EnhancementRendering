// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D texEnv; // TEX_BUFF1
uniform sampler2D screenNormal; // TEX_BUFF1

 // local2 = { n_i/n_t, (n_i/n_t)^2, <unused>, <unused> }
 //    n_i = index of refraction of external material (air)
 //    n_t = index of refraction of internal material
uniform vec2 local2; // TEX_BUFF1


// Function to compute a refraction direction.  Evidently, the built-in Cg refraction
//    computes pseudo-refraction vectors, as per previous work in the area.
vec3 refraction(in vec3 incident,in vec3 normal,in float ni_nt,in float ni_nt_sqr )
{
    float IdotN = dot( -incident, normal );
    float cosSqr = 1.0 - ni_nt_sqr*(1.0 - IdotN*IdotN);
    if (cosSqr < 0.0) 
        cosSqr = 0.0;
    else
        cosSqr = sqrt( cosSqr );
    return  normalize( ni_nt * incident + (ni_nt* IdotN - cosSqr) * normal ); 
}


void main(void) {

    const float PI  = 3.14159265;

    vec3 n = texture2D(screenNormal,gl_TexCoord[0].xy).xyz;
    vec3 v;
    vec3 light;
    if(n==vec3(0.0))
    {
        light = -normalize(vec3((gl_TexCoord[0].x-0.5)*-0.49*5.0, (gl_TexCoord[0].y-0.5)*-0.49*5.0, 1.0));
    }
    else
    {
        v = -normalize(vec3((gl_TexCoord[0].x-0.5)*-0.49*5.0, (gl_TexCoord[0].y-0.5)*-0.49*5.0, 1.0));
        light = refraction( v, n, local2.x, local2.y );

    }

    
        float y = 1.0 - acos(light.y)/PI;
        vec2 coordn = normalize(vec2(light.x, light.z));
        coordn.x = acos(coordn.x)/PI;
        coordn.y = acos(coordn.y)/PI;
        float x = coordn.x>0.5 ? 0.5*coordn.y : 1.0-0.5*coordn.y;

        gl_FragData[0] = vec4(texture2D(texEnv,vec2(x,y)).xyz, 1.0);    






        
}
