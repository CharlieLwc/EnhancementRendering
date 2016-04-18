// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120

uniform sampler2D depth;

uniform float sw;     // 1/width
uniform float sh;     // 1/height

uniform float threshold;     // 1/height
uniform float hardness;     // 1/height

float silhouetteWeight(in float s)
{
    float t2 = 0.9+threshold;
    float t1 = t2-hardness;

    return smoothstep(t1,t2,max(1.0-s,0.9));
}

float  meanSilhouette() {

// |A|B|C|
// |D|X|E|
// |F|G|H|
    //depth
    float x,a,b,c,d,e,f,g,h;
    // silhouette weights 
    float Asw,Bsw,Csw,Dsw,Esw,Fsw,Gsw,Hsw;
    
    float xc = gl_TexCoord[0].s;
    float yc = gl_TexCoord[0].t;

    
    if(texture2D(depth,vec2(xc,yc)).z < 0.5)
        return 0.0;


    float xm = gl_TexCoord[0].s-sw;
    float xp = gl_TexCoord[0].s+sw;
    
    float ym = gl_TexCoord[0].t-sh;
    float yp = gl_TexCoord[0].t+sh;
    
    
        
    x = texture2D(depth,vec2(xc,yc)).y;
    a = texture2D(depth,vec2(xm,yp)).y;
    b = texture2D(depth,vec2(xc,yp)).y;
    c = texture2D(depth,vec2(xp,yp)).y;
    d = texture2D(depth,vec2(xm,yc)).y;
    e = texture2D(depth,vec2(xp,yc)).y;
    f = texture2D(depth,vec2(xm,ym)).y;
    g = texture2D(depth,vec2(xc,ym)).y;
    h = texture2D(depth,vec2(xp,ym)).y;

    Asw = silhouetteWeight(abs(a-x));
    Bsw = silhouetteWeight(abs(b-x));
    Csw = silhouetteWeight(abs(c-x));
    Dsw = silhouetteWeight(abs(d-x));
    Esw = silhouetteWeight(abs(e-x));
    Fsw = silhouetteWeight(abs(f-x));
    Gsw = silhouetteWeight(abs(g-x));
    Hsw = silhouetteWeight(abs(h-x));

    
    return (Asw+Bsw+Csw+Dsw+Esw+Fsw+Gsw+Hsw)/8.0;

}



void main(void) {
          

    float ms = meanSilhouette();

    gl_FragData[0] = vec4(ms); 
    
    
 
}
