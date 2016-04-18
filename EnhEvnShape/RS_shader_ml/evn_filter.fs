// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D environmentMap;
uniform sampler2D normal;         // TEX_BUFF0 (normal,depth)
uniform sampler2D view;         // TEX_BUFF1 (view,1.0)
uniform sampler2D evnSampleDir;         // 
uniform sampler2D evnSampleCol;         // 
uniform sampler2D curvature;

uniform float increasement;
uniform int bufferSize;

uniform float xSpan;
uniform float ySpan;
uniform float envsize;           // TONE MAPPING - environement size (number of lights)
uniform float enhanceScale;           // TONE MAPPING - environement size (number of lights)

uniform float diffuseValue;
uniform float specularValue;
uniform float shiness;
uniform vec3 diffColor;
uniform vec3 specColor;
uniform float highlightTher;


// BRDF coeficients
vec3 diffTerm;
vec3 specTerm;


float c = 0.f;


// **** ASNGAN FUNCTIONS ****
void computeASNganCoefs(in vec3 n, in vec3 v, in vec3 l) {
  const float PI_23  = 72.256631032;
  const float PI_8   = 25.132741228;
  


  float f = 0.5;
  float r = (diffColor.x+diffColor.y+diffColor.z)*0.333333;
  vec3  h = normalize(l-v);
  float normalization = (shiness+1.0)/PI_8;

  float VDotH = max(dot(-v,h),0.0);
  float NDotH = max(dot(n,h),0.0);
  

  float NDotV = max(dot(n,-v),0.0);
  float NDotL = max(dot(n,l),0.0);
  
  //specTerm = specColor*max(((pow(NDotH,shiness)/(VDotH*max(NDotL,NDotV)))*normalization*(f+(1.0-f)*pow(VDotH,5.0))),0.0)/envsize;
  specTerm = specColor*max(((pow(NDotH,shiness)/max(VDotH*max(NDotL,NDotV), highlightTher))*normalization*(f+(1.0-f)*pow(VDotH,5.0))),0.0)/envsize;
  diffTerm = diffColor*max((((28.0*r)/PI_23)*(1.0-f)*(1.0-pow(1.0-(NDotL/2.0),5.0))*(1.0-pow(1.0-(NDotV/2.0),5.0))),0.0)/envsize;

   


}


// **** SCALING FUNCTION ****
float scale(in float impfunc,in float beta) {
  float alpha = 0.5;  // RS - controls the scaling invariant point ]0,1[
  float expbeta = exp(-beta);
  return (alpha*expbeta+impfunc*(1.0-alpha-alpha*expbeta)) / (alpha+impfunc*(expbeta-alpha-alpha*expbeta)); 
}





vec3 computeLighting(in vec3 n,in vec3 v) {

    vec3 m =  vec3(0.0);
        
    vec3  l, p;



    float tx = 0.0;
    float ty = 0.0;
    for(int i=0; i<bufferSize; i++)
    {
        p = texture2D(evnSampleCol,vec2(tx,ty)).xyz;
        l = texture2D(evnSampleDir,vec2(tx,ty)).xyz;
   
        tx+=increasement;
        if(tx>=1.0)
        {
            tx = 0.0;
            ty+=increasement;
        }
   
        computeASNganCoefs(n, v, l);


        
        m += ((diffTerm*scale(max(dot(n,l),0.0),c*enhanceScale)*diffuseValue) + (specTerm*scale(max(dot(normalize(l-v),l),0.0),c*enhanceScale)*specularValue))*p;

    }
   
    return m;
}

void main(void) {
    const float ZERO = 0.00000001;
    const float PI  = 3.14159265;
    vec3 v, n;

    
        gl_FragData[0] = vec4(texture2D(evnSampleDir,gl_TexCoord[0].xy).xyz,1.0);


            n = vec3(texture2D(normal, gl_TexCoord[0].xy).xyz);

    if(n==vec3(0.0))
    {

        vec3 light = -normalize(vec3((gl_TexCoord[0].x-0.5)*xSpan*5.0, (gl_TexCoord[0].y-0.5)*ySpan*5.0, 1.0));
        
        float y = 1.0 - acos(light.y)/PI;
        vec2 coordn = normalize(vec2(light.x, light.z));
        coordn.x = acos(coordn.x)/PI;
        coordn.y = acos(coordn.y)/PI;
        float x = coordn.x>0.5 ? 0.5*coordn.y : 1.0-0.5*coordn.y;

        gl_FragData[0] = vec4(texture2D(environmentMap,vec2(x,y)).xyz, 1.0);        
        gl_FragData[0] = vec4(0.0);        
        
    }
    else
    {
  
        n = normalize(n);
        v = normalize(texture2D(view, gl_TexCoord[0].xy).xyz);
        
          if(enhanceScale>0.0)
              c = texture2D(curvature, gl_TexCoord[0].xy).x;

        gl_FragData[0] =  vec4(computeLighting(n,v), 1.0);
        
        


    }
    
    if(false)
    {
        float y = -cos(gl_TexCoord[0].y*PI);
        float theta = 2.0*PI*(1.0 - gl_TexCoord[0].x);
        float length = sqrt(1.0-y*y);
        float z = cos(theta)*length;
        float x = sin(theta)*length;
        vec3 v  = normalize(vec3(x, y, z));
        vec3  l, p;
        
        vec3 m;
        m =  vec3(0.0);

        float tx = 0.0;
        float ty = 0.0;
        for(int i=0; i<bufferSize; i++)
        {
            p = texture2D(evnSampleCol,vec2(tx,ty)).xyz;
            l = texture2D(evnSampleDir,vec2(tx,ty)).xyz;


            tx+=increasement;
            if(tx>=1.0)
            {
                tx = 0.0;
                ty+=increasement;
            }

            //m+=clamp(pow(dot(l,v),500.0), 0.0, 1.0)*vec3(p);
       //     m+=pow(clamp(dot(l,v), 0.0, 1.0),500.0)*vec3(p);
       //    m+=vec3(1000000);
            m += pow(clamp(dot(l,v), 0.0, 1.0),500.0)*p;

        }
        gl_FragData[0] =  vec4(m,1.0);

    }

    

    



}
