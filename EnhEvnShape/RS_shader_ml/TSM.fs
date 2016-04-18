// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D texEnv; // TEX_BUFF1
uniform sampler2D backScreenNormal; // TEX_BUFF1
uniform sampler2D texView; // TEX_BUFF1

uniform float enhanceScale;

void main(void) {

    const float PI  = 3.14159265;

    vec3 texColor = vec3(1.0, 1.0, 1.0);

    float r0 = 0.01701323251417769376181474480151;
    vec3 n = texture2D(backScreenNormal,gl_TexCoord[0].xy).xyz;
    if(n == vec3(0.0))
        gl_FragData[0] = vec4(0.0); 
    else
    {   
        

//        vec3 reverseLight= -normalize(vec3((gl_TexCoord[0].x-0.5)*-0.49*5.0, (gl_TexCoord[0].y-0.5)*-0.49*5.0, 1.0));
//    
//    
//        float y = 1.0 - acos(reverseLight.y)/PI;
//        vec2 coordn = normalize(vec2(reverseLight.x, reverseLight.z));
//        coordn.x = acos(coordn.x)/PI;
//        coordn.y = acos(coordn.y)/PI;
//        float x = coordn.x>0.5 ? 0.5*coordn.y : 1.0-0.5*coordn.y;
//
//        vec3 lightColor  = texture2D(texEnv,vec2(x,y)).xyz;      
//
//        float NDotL = clamp(dot(n,reverseLight),0.0,1.0);
//        float k = pow(1.0-NDotL,5.0);
//        float Ft = 1.0 - (r0 + k*(1.0-r0));
//        vec3 E = Ft*NDotL*texColor*lightColor;

        float curv = texture2D(texView,gl_TexCoord[0].xy).w;
        float intensity = 1.0;

        if(curv > 0.1)
              intensity -= enhanceScale;
        else if(curv < -0.1)
              intensity += enhanceScale;
        gl_FragData[0] = vec4(intensity);     
    }

 
        
}
