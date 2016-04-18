// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120

attribute vec3 smoothedNormal0;
attribute vec3 smoothedNormal1;
attribute vec3 smoothedNormal2;
attribute vec3 smoothedNormal3;
attribute vec3 smoothedNormal4;
attribute vec3 smoothedNormal5;
attribute vec3 smoothedNormal6;
attribute vec3 smoothedNormal7;
attribute vec3 smoothedNormal8;
attribute vec3 smoothedNormal9;
attribute vec3 smoothedNormal10;

uniform float weights[12];    // scaling intensity: varies between 0 and 1 (may be more)

uniform float enhanceScale;
uniform float averageIlluminate;
uniform int layerNum;


varying float density;
uniform vec3  lGlobal;

void main() {

	vec3 n[12];
	
	n[0] = gl_Normal;
	n[1] = smoothedNormal0;
	n[2] = smoothedNormal1;
	n[3] = smoothedNormal2;
	n[4] = smoothedNormal3;
	n[5] = smoothedNormal4;
	n[6] = smoothedNormal5;
	n[7] = smoothedNormal6;
	n[8] = smoothedNormal7;
	n[9] = smoothedNormal8;
	n[10] = smoothedNormal9;
	n[11] = smoothedNormal10;



	density = 0.0;


	for(int i=0; i<layerNum; i++)
	{
		vec3 lDirection = normalize(lGlobal - n[i+1]*dot(n[i+1], lGlobal));
		density += weights[i] * clamp(enhanceScale*dot(n[i],lDirection),-0.45,1.0);
	}
	

	density += weights[layerNum] * clamp(enhanceScale*dot(n[layerNum],lGlobal),-0.45,1.0);

//	density *= 0.5;

//	if(layerNum != 0)
//		density += 0.5;

density += averageIlluminate;


  gl_Position = ftransform();

}
