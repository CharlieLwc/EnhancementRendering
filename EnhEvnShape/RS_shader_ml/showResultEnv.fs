// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D tex; // TEX_BUFF1


uniform sampler2D dir;         // 
uniform sampler2D col;         //
uniform sampler2D scale;         //

uniform sampler2D screenNormal;         // TEX_BUFF0 (normal,depth)
uniform sampler2D view;         // TEX_BUFF1 (view,1.0)/ 

uniform float increasement;
uniform float enhance;


uniform int bufferSize;


float scaled(in float impfunc,in float beta) {
  float alpha = 0.5;  // RS - controls the scaling invariant point ]0,1[
  float expbeta = exp(-beta);
  return (alpha*expbeta+impfunc*(1.0-alpha-alpha*expbeta)) / (alpha+impfunc*(expbeta-alpha-alpha*expbeta)); 
}


void main(void) {

    const float PI  = 3.14159265;
    
    vec3 n = texture2D(screenNormal,gl_TexCoord[0].xy).xyz;
    if(true||n==vec3(0.0))
    {
        float y = -cos(gl_TexCoord[0].y*PI);
        float theta = 2.0*PI*(1.0 - gl_TexCoord[0].x);
        float length = sqrt(1.0-y*y);
        float z = cos(theta)*length;
        float x = sin(theta)*length;
        vec3 v  = normalize(vec3(x, y, z));
        vec3  l, p;
    
        vec3 m;
        m =  vec3(0.0);
        float angle = -100.0;
    
        float tx = 0.0;
        float ty = 0.0;
        for(int i=0; i<bufferSize; i++)
        {
            p = texture2D(col,vec2(tx,ty)).xyz;
            l = texture2D(dir,vec2(tx,ty)).xyz;
            float s = texture2D(scale,vec2(tx,ty)).x;
            tx+=increasement;
            if(tx>=1.0)
            {
                tx = 0.0;
                ty+=increasement;
            }
            float tempangle = dot(l,v);
            if(tempangle > angle)
            {
                angle = tempangle;
                s = s>0.5?3000.0:0.0;
                m = vec3(p*s);
            }
    
        }
        gl_FragData[0] =  vec4(m, 1.0);
//        gl_FragData[0] =  vec4(1.0);
    }
    else
    {

        
        float curv = texture2D(view,gl_TexCoord[0].xy).w;
        
        vec3 v = normalize(vec3((gl_TexCoord[0].x-0.5)*-0.49*5.0, (gl_TexCoord[0].y-0.5)*-0.49*5.0, 1.0));
        v = texture2D(view,gl_TexCoord[0].xy).xyz;
        vec3  l, p;
        vec3 m;
        m =  vec3(0.0);
    
        float tx = 0.0;
        float ty = 0.0;
        for(int i=0; i<bufferSize; i++)
        {
            p = texture2D(col,vec2(tx,ty)).xyz;
            l = texture2D(dir,vec2(tx,ty)).xyz;
            float s = texture2D(scale,vec2(tx,ty)).x;

            tx+=increasement;
            if(tx>=1.0)
            {
                tx = 0.0;
                ty+=increasement;
            }
            
            vec3 halfvector = normalize(l + v);
//            if(s > 0.5)
//                m += vec3(pow(clamp(dot(l,n),0.0,1.0),20.0)*p)*0.02;
//            else
                m += vec3(pow(clamp(dot(l,n),0.0,1.0),20.0)*p)*0.015 * scaled(max(dot(n,l),0.0),curv*1.0*1.0*enhance);
    
                


        }
        gl_FragData[0] =  vec4(m, 1.0);



    }




    
    

        
}
