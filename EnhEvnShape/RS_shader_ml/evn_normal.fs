// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable 


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

uniform bool needWorldNormal;
uniform bool needView;
uniform bool needWorldView;
uniform bool needWorldPos;
uniform bool needColor;
uniform bool needLayer0;
uniform bool needLayer1;
uniform bool needCurv;
uniform bool needDepth;
uniform bool needFeature;
uniform bool needCurvDir;


void main(void) {

//    float a = gl_Color.x>=curv_thre_pos ? 1.0:0.0;    
//    a = gl_Color.y>=curv_thre_neg ? -1.0:a;    
   
//    a = gl_Color.x>0.0 ? gl_Color.x:-gl_Color.y;




//	vec3 tempN = normalize(worldNormal);
//	float radiance1 = max(dot(tempN,vec3(1.0, 0.0, 0.0)), 0.0);
//	float radiance2 = max(dot(tempN,vec3(0.0, -1.0, 0.0)), 0.0);
//	float radiance3 = max(dot(tempN,vec3(0.0, 0.0, -1.0)), 0.0);
//	
//	vec3 color1 = vec3(0.0, 0.66, 1.0);
//	vec3 color2 = vec3(1.0, 0.5478, 0.0);
//	vec3 color3 = vec3(1.0, 0.0, 0.8795);
//
//	
//	color1 = vec3(0.0, 1.0, 1.0);
//	color2 = vec3(0.7588, 1.0, 0.0);
//	color3 = vec3(1.0, 0.4243, 1.0);

//	color1 = vec3(1.0);
//	color2 = vec3(0.4663);
//	color3 = vec3(0.0);

//	vec3 tempColor = color1*radiance1 + 
//				color2*radiance2 + 
//				color3*radiance3;



    int increasement = 0;
      
    gl_FragData[increasement++] = vec4(normalize(screenNormal), 1.0);
    if(needWorldNormal)
        gl_FragData[increasement++] = vec4(normalize(worldNormal), 1.0);
    if(needView)
        gl_FragData[increasement++] = vec4(normalize(view), 1.0);
    if(needWorldView)
        gl_FragData[increasement++] = vec4(normalize(worldView), 1.0);
    if(needWorldPos)
        gl_FragData[increasement++] = vec4(worldPos, 1.0);
    if(needLayer0)
        gl_FragData[increasement++] = vec4(normalize(smoothedNormalScreen0), 1.0);
    if(needLayer1)
        gl_FragData[increasement++] = vec4(normalize(smoothedNormalScreen1), 1.0);

    if(needCurv || needDepth || needFeature)
    {
        vec3 CD = vec3(0.0);
        if(needCurv)
            CD.x = curvatureObject;
        else if(needFeature)
            CD.x = feature;

        if(needDepth)
            CD.y = depth;

        if(abs(worldNormal.x) > 0.1 || abs(worldNormal.y) > 0.1 || abs(worldNormal.z) > 0.1)
           CD.z = 1.0;
           
        gl_FragData[increasement++] = vec4(CD, 1.0);
    }

    if(needCurvDir)
        gl_FragData[increasement++] = vec4(normalize(curvDirection), 1.0);

    if(needColor)
        gl_FragData[increasement++] = gl_Color;

}
