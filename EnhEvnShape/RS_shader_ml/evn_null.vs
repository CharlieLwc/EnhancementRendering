// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120

void main(void)
{
	gl_FrontColor = gl_Color;
  gl_TexCoord[0] = gl_MultiTexCoord0;
  gl_Position    = ftransform();
}
