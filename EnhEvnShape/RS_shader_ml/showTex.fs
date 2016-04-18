// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D tex; // TEX_BUFF1

void main(void) {

    vec3 n = texture2DLod(tex,gl_TexCoord[0].xy, 0.0).xyz;
    gl_FragColor = vec4(gl_TexCoord[0].xy, 0.0, 1.0);
           
    gl_FragColor = vec4(texture2D(tex,gl_TexCoord[0].xy).w);
    gl_FragColor = vec4(texture2D(tex,gl_TexCoord[0].xy).xyz,1.0);



}
