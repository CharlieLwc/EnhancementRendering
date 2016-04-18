// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D texEnv; // TEX_BUFF1
uniform sampler2D screenNormal; // TEX_BUFF1
uniform sampler2D texView; // TEX_BUFF1
uniform sampler2D backScreenNormal; // TEX_BUFF1
uniform sampler2D texBackView; // TEX_BUFF1
uniform sampler2D texTSM; // TEX_BUFF1

uniform vec3 zr3;
uniform vec3 zv3;
uniform vec3 sigma_t;
uniform vec3 sigma_tr;
uniform vec3 Rd0;





vec4 Rd(in vec3 xout, in vec3 texCoords) // texCoords.z = miplevel
{


    vec4 tsmnormal = texture2DLod(backScreenNormal, texCoords.xy, texCoords.z).xyzw;
   
    if(tsmnormal.xyz == vec3(0.0)) return vec4(0.0);

    float w = tsmnormal.w;


    vec3 lightNormal = normalize(tsmnormal.xyz);

    
    vec3 irradiance = texture2DLod(texTSM, texCoords.xy, texCoords.z).xyz/w;



    vec3 xin = texture2DLod(texBackView, texCoords.xy, texCoords.z).xyz/w;
    lightNormal = normalize(xin);
   

    float d = distance(xin, xout);
//    if(d < 0.00001)
//        return vec4(irradiance*1.0/(1.0+d), 1.0);

    vec3 xr = xin - zr3*lightNormal;
    vec3 xv = xin + zv3*lightNormal;



    float dr = length(xr-xout);
    float dv = length(xv-xout);

    vec3 Rd1r = zr3*(sigma_tr*dr+1.0);
    vec3 Rd2r =exp(-sigma_tr*dr) / (sigma_t);
    vec3 Rd1v =zv3*(sigma_tr+1.0);
    vec3 Rd2v =exp(-sigma_tr*dv) / (sigma_t);
    //    vec3 Rd1r = zr3*(sigma_tr*dr+1.0);
//    vec3 Rd2r =exp(-sigma_tr*dr) / (sigma_t*pow(dr,3.0));
//    vec3 Rd1v =zv3*(sigma_tr*dv+1.0);
//    vec3 Rd2v =exp(-sigma_tr*dv) / (sigma_t*pow(dv,3.0));


    vec3 finalColor = Rd0 * (Rd1r * Rd2r + Rd1v * Rd2v);

    return vec4(irradiance * finalColor, 1.0);




//
//    vec4 tsmnormal = texture2DLod(backScreenNormal, texCoords.xy, texCoords.z).xyzw;
//   
//    if(tsmnormal.xyz == vec3(0.0)) return vec4(0.0);
//
//    float w = tsmnormal.w;
//
//
//    vec3 lightNormal = normalize(tsmnormal.xyz);
//    vec3 irradiance = texture2DLod(texTSM, texCoords.xy, texCoords.z).xyz / w;
//    vec3 xin = texture2DLod(texBackView, texCoords.xy, texCoords.z).xyz / w;
//    lightNormal = normalize(xin);
//    irradiance = vec3(2.0);
//
//    float d = distance(xin, xout);
//    if(d < 0.001)
//        return vec4(irradiance*1.0/(1.0+d), 1.0);
//
//    vec3 xr = xin - zr3*lightNormal;
//    vec3 xv = xin + zv3*lightNormal;
//    float dr = length(xr-xout);
//    float dv = length(xv-xout);
//
//    vec3 Rd1r = zr3*(sigma_tr*dr+1.0);
//    vec3 Rd2r =exp(-sigma_tr*dr) / (sigma_t*pow(dr,3.0));
//    vec3 Rd1v =zv3*(sigma_tr*dv+1.0);
//    vec3 Rd2v =exp(-sigma_tr*dv) / (sigma_t*pow(dv,3.0));
//
//    vec3 finalColor = Rd0 * (Rd1r * Rd2r + Rd1v * Rd2v);
//    return vec4(min(finalColor, 0.8)*irradiance, 1.0);
}


vec3 sampleTSM(void)
{

      vec2 delta[4];
      delta[0] = 1.0 / vec2(512.0, 512.0);
      delta[1] = 1.0 / (2.0 * vec2(512.0, 512.0));
      delta[2] = 1.0 / (8.0 * vec2(512.0, 512.0));
      delta[3] = 1.0 / (16.0 *vec2(512.0, 512.0));
  
      float k[4];
      k[0] = 1.0 / (4.0 * 9.0);
      k[1] = 1.0 / (4.0 * 4.0);
      k[2] = 1.0 / (4.0 * 4.0);
      k[3] = 1.0 / (4.0 * 4.0);

      vec3 xout = texture2D(texView,gl_TexCoord[0].xy).xyz;
      vec3 adjTex = vec3(gl_TexCoord[0].xy,0.0);
  
      vec4 color = k[0]*Rd(xout, adjTex);    
      
      
      return color.xyz;


      adjTex.x = gl_TexCoord[0].x+delta[0].x;
      color += k[0]*Rd(xout, adjTex);
      adjTex.x += delta[0].x;
      color += k[0]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x-delta[0].x;
      color += k[0]*Rd(xout, adjTex);
      adjTex.x -= delta[0].x;
      color += k[0]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x;
      adjTex.y = gl_TexCoord[0].y+delta[0].y;
      color += k[0]*Rd(xout, adjTex);
      adjTex.y += delta[0].y;
      color += k[0]*Rd(xout, adjTex);
      adjTex.y = gl_TexCoord[0].y-delta[0].y;
      color += k[0]*Rd(xout, adjTex);
      adjTex.y -= delta[0].y;
      color += k[0]*Rd(xout, adjTex);
      
      adjTex.x = gl_TexCoord[0].x+delta[1].x;
      adjTex.y = gl_TexCoord[0].y+delta[1].y;
      adjTex.z = 1.0;
      color += k[1]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x-delta[1].x;
      color += k[1]*Rd(xout, adjTex);
      adjTex.y = gl_TexCoord[0].y-delta[1].y;
      color += k[1]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x+delta[1].x;
      color += k[1]*Rd(xout, adjTex);
  
      adjTex.x = gl_TexCoord[0].x+delta[2].x;
      adjTex.y = gl_TexCoord[0].y+delta[2].y;
      adjTex.z = 3.0;
      color += k[2]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x-delta[2].x;
      color += k[2]*Rd(xout, adjTex);
      adjTex.y = gl_TexCoord[0].y-delta[2].y;
      color += k[2]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x+delta[2].x;
      color += k[2]*Rd(xout, adjTex);
  
      adjTex.x = gl_TexCoord[0].x+delta[3].x;
      adjTex.y = gl_TexCoord[0].y;
      adjTex.z = 4.0;
      color += k[3]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x-delta[3].x;
      color += k[3]*Rd(xout, adjTex);
      adjTex.x = gl_TexCoord[0].x;
      adjTex.y = gl_TexCoord[0].y+delta[3].y;
      color += k[3]*Rd(xout, adjTex);
      adjTex.y = gl_TexCoord[0].y-delta[3].y;
      color += k[3]*Rd(xout, adjTex);
      
      color /= color.w;
      return color.xyz;
return vec3(0.0);
}

void main(void) {

    const float PI  = 3.14159265;

    vec3 texColor = vec3(1.0, 1.0, 1.0);

    float r0 = 0.01701323251417769376181474480151;
    vec3 n = texture2D(screenNormal,gl_TexCoord[0].xy).xyz;
    vec3 reverseLight= -normalize(vec3((gl_TexCoord[0].x-0.5)*-0.49, (gl_TexCoord[0].y-0.5)*-0.49, 1.0));

    if(n == vec3(0.0))
    {
    
    
        float y = 1.0 - acos(reverseLight.y)/PI;
        vec2 coordn = normalize(vec2(reverseLight.x, reverseLight.z));
        coordn.x = acos(coordn.x)/PI;
        coordn.y = acos(coordn.y)/PI;
        float x = coordn.x>0.5 ? 0.5*coordn.y : 1.0-0.5*coordn.y;

        gl_FragData[0] = vec4(texture2D(texEnv,vec2(x,y)).xyz, 1.0);  
        gl_FragData[0] = vec4(1.0);         
    }
    else
    {   
        
        vec3 outputColor =  sampleTSM();
        
        float NDotL = clamp(dot(n,-reverseLight),0.0,1.0);
        float k = pow(1.0-NDotL,5.0);
        float Ft = 1.0 - (r0 + k*(1.0-r0));
      //  outputColor *= (Ft * 1.0/PI);

      outputColor += vec3(NDotL*0.15);

        gl_FragData[0] = vec4(outputColor, 1.0);     
    }
//    n = texture2DLod(backScreenNormal,gl_TexCoord[0].xy,4.0).xyz;
//    gl_FragData[0] = vec4(n, 1.0);     
        
}
