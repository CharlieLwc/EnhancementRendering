// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

varying vec3 screenNormal;
varying vec3 view;
uniform float curv_thre_pos;
uniform float curv_thre_neg;
uniform bool needDepth;

void main(void) {

    float a = gl_Color.x>=curv_thre_pos ? 1.0:0.0;    
    a = gl_Color.y>=curv_thre_neg ? -1.0:a;    
   
//    a = gl_Color.x>0.0 ? gl_Color.x:-gl_Color.y;

    if(needDepth)
        gl_FragData[0] = vec4(normalize(screenNormal), -view.z);
    else
        gl_FragData[0] = vec4(normalize(screenNormal), 1.0);
            
    gl_FragData[1] = vec4(view,a);

      
}
