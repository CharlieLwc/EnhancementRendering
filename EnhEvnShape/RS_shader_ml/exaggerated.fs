// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

varying float density;

uniform vec3 color;

void main(void) {
//
//    vec3 color = vec3(0.286, 0.415, 0.568);
//    color = vec3(0.12, 0.43, 0.50);
//    
//    //map
//    color = vec3(0.28627, 0.41568, 0.56863); 
//    //gargo
//    color = vec3(0.855, 0.34, 0.1333);
//    //vase
//    color = vec3(0.486, 0.4196, 0.0745);
//    //persidon
//    color = vec3(0.5, 0.5, 0.5);
//    //adma
//    color = vec3(0.46745, 0.27705, 0.1601);
//    //ball
//    color = vec3(0.345, 0.7255, 0.749);

    gl_FragData[0] = vec4(color * density, 1.0);
    
    
      
}
