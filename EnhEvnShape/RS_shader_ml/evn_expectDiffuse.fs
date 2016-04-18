// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D texEvn;
uniform sampler2D screenNormal;         // TEX_BUFF0 (normal,depth)
uniform sampler2D view;         // TEX_BUFF1 (view,1.0)


uniform vec3 mainDirection;
uniform float xSpan;
uniform float ySpan;


void main(void) {
    const float ZERO = 0.00000001;
    const float PI  = 3.14159265;
    vec3 v, n;

    n = vec3(texture2D(screenNormal, gl_TexCoord[0].xy).xyz);


    




    if(n==vec3(0.0))
    {
        vec3 light = -normalize(vec3((gl_TexCoord[0].x-0.5)*xSpan*5.0, (gl_TexCoord[0].y-0.5)*ySpan*5.0, 1.0));
        
        float y = 1.0 - acos(light.y)/PI;
        vec2 coordn = normalize(vec2(light.x, light.z));
        coordn.x = acos(coordn.x)/PI;
        coordn.y = acos(coordn.y)/PI;
        float x = coordn.x>0.5 ? 0.5*coordn.y : 1.0-0.5*coordn.y;

        gl_FragData[0] = vec4(texture2D(texEvn,vec2(x,y)).xyz, 1.0);           
        
    }
    else
    {


  
        n = normalize(n);
        v = normalize(texture2D(view, gl_TexCoord[0].xy).xyz);
        v = normalize(mainDirection);
        
        gl_FragData[0] =  vec4(v,1.0);
        gl_FragData[0] =  vec4(dot(n,v)*0.6);
    }
    
    
    
}
