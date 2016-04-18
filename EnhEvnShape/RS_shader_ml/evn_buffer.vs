// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120


varying vec3  screenNormal;
varying vec3  wordNormal;
varying vec3  view;


void main() {

	//calculate everything in screen space
	screenNormal = gl_NormalMatrix*gl_Normal;
	view = -vec3(gl_ModelViewMatrix*gl_Vertex);

	gl_Position = ftransform();
	gl_FrontColor =gl_Color;


//  const float PI   = 3.14159265;
//  vec3 light = normalize(reflect(view, normal));
//  float y = 1.0 - acos(light.y)/PI;
//  vec2 coordn = normalize(vec2(light.x, light.z));
//  coordn.x = acos(coordn.x)/PI;
//  coordn.y = acos(coordn.y)/PI;
//  float x = coordn.x>0.5 ? 0.5*coordn.y : 1.0-0.5*coordn.y;
//  evnColor = vec3(texture2D(evn,vec2(x,y)).xyz);

}
