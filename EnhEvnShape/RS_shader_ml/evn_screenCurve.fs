// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120

uniform sampler2D norms[3];     // 0,1:snoothedNormal. 2screenNormal;


uniform float weights[3];    // scaling intensity: varies between 0 and 1 (may be more)

uniform sampler2D silhouette;

uniform float sw;     // 1/width
uniform float sh;     // 1/height
uniform float enhanceScale;     // 

vec2 coord;
vec3 v;



float curvature() {
  // |A|B|C|
  // |D|X|E|
  // |F|G|H|

  float curvature = 0.0;
  float ms = texture2D(silhouette,vec2(coord.s,coord.t)).x;
  float ms1 = texture2D(silhouette,vec2(coord.s+sw,coord.t)).x;
  float ms2 = texture2D(silhouette,vec2(coord.s,coord.t-sh)).x;

  for(int lIndex=0; lIndex<3; lIndex++)
  {
    if(weights[lIndex]>0.0)
    {
        vec3 E = texture2D(norms[lIndex],vec2(coord.s+sw,coord.t)).xyz;
        vec3 D = texture2D(norms[lIndex],vec2(coord.s-sw,coord.t)).xyz;
        vec3 B = texture2D(norms[lIndex],vec2(coord.s,coord.t+sh)).xyz;
        vec3 G = texture2D(norms[lIndex],vec2(coord.s,coord.t-sh)).xyz;
        vec3 X = texture2D(norms[lIndex],vec2(coord.s,coord.t)).xyz;
        



       if(!(ms < 0.8 || ms < 0.8 || ms2<0.8 ||E.xyz==vec3(0.0) || D.xyz==vec3(0.0) || B.xyz==vec3(0.0) || G.xyz==vec3(0.0)))
         //curvature += 0.5*(E.x-D.x+B.y-G.y) * weights[lIndex] * smoothstep(0.1, 0.5, dot(-v,X));
         //curvature += 0.5*(E.x-X.x+X.y-G.y) * weights[lIndex] * smoothstep(0.1, 0.5, dot(-v,X));
         //curvature += 0.5*(E.x-D.x+B.y-G.y) * weights[lIndex];
         curvature += 0.5*(E.x-X.x+X.y-G.y) * weights[lIndex];

    } 
  }
  return curvature;
}


void main(void) {

    coord = vec2(gl_TexCoord[0].x, gl_TexCoord[0].y);
          
    float curvature = curvature();
	


    gl_FragData[0] = vec4(curvature*enhanceScale); 

}
