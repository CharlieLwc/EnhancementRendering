// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable 

void main(void) {
  

    vec4 col = gl_Color;

    if(col.y >= 0.7)
        col = vec4(0.0);
    else if(col.y > 0.3)
        col = vec4(1.0, 0.0, 0.0, 0.0);
    else
        col = vec4(0.0, 1.0, 0.0, 0.0);
    



    gl_FragData[0] = col;

}
