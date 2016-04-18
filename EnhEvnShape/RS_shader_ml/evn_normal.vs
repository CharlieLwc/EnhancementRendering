// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120



varying vec3 screenNormal;
varying vec3 worldNormal;
varying vec3 view;
varying vec3 worldView;
varying vec3 worldPos;
varying vec3 smoothedNormalScreen0;
varying vec3 smoothedNormalScreen1;
varying float curvatureObject;
varying float depth;
varying float feature;
varying vec3 curvDirection;

attribute vec3 smoothedNormal0;
attribute vec3 smoothedNormal1;
attribute float weightedCurv;
attribute vec3 curDir;
attribute float weightedFeature;

uniform vec3  campos; // camera position
uniform float zmin;   // distance to the closest point 
uniform float zmax;   // distance to the farest point 

uniform bool needWorldNormal;
uniform bool needView;
uniform bool needWorldView;
uniform bool needWorldPos;
uniform bool needLayer0;
uniform bool needLayer1;
uniform bool needCurv;
uniform bool needDepth;
uniform bool needFeature;
uniform bool needCurvDir; 

void main() {

	
	screenNormal = normalize(gl_NormalMatrix*gl_Normal);
	if(needWorldNormal)
		worldNormal = gl_Normal;
	if(needWorldView)
		worldView = vec3(gl_Vertex.xyz) - campos;
	if(needView)
		view = vec3(gl_ModelViewMatrix*gl_Vertex);
	if(needWorldPos)
		worldPos = gl_Vertex.xyz;
	if(needLayer0)
		smoothedNormalScreen0 = gl_NormalMatrix*smoothedNormal0;
	if(needLayer1)
		smoothedNormalScreen1 = gl_NormalMatrix*smoothedNormal1;
	if(needCurv)
		curvatureObject = weightedCurv;
    if(needDepth)
		depth  = log(clamp(-(gl_ModelViewMatrix*gl_Vertex).z,zmin,zmax)/zmin)/log(zmax/zmin);
    if(needFeature)
		feature  = weightedFeature;
    if(needCurvDir)
		curvDirection = curDir;


		
	gl_FrontColor =gl_Color;
	gl_Position = ftransform();
}
