// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D screenNormal;         // TEX_BUFF0 (normal,depth)
uniform sampler2D view;         // TEX_BUFF1 (view,1.0)
uniform sampler2D texEnvDir;

uniform float sampleScale; //the size of the environment map
uniform float sampleScaleReverse; //the size of the environment map

uniform float increasement;


void main(void) {
    const float ZERO = 0.00000001;
    const float PI  = 3.14159265;
    

    if(gl_TexCoord[0].x>sampleScale || gl_TexCoord[0].y>sampleScale)
        gl_FragData[0] = vec4(0.0);
    else
    {
        vec3 l = vec3(texture2D(texEnvDir, gl_TexCoord[0].xy*sampleScaleReverse).xyz);
        float intensity = 0.0;
        for(float x=increasement*0.5; x<1.0; x+=increasement)
            for(float y=increasement*0.5; y<1.0; y+=increasement)
            {
                float curv = texture2D(view,vec2(x,y)).w;
               if(curv != 1.0 && curv != -1.0)
                   continue;

                vec3 v = normalize(vec3((x-0.5)*-0.49*5.0, (y-0.5)*-0.49*5.0, 1.0));
                vec3 halfvector = normalize(l + v);

                vec3 n = texture2D(screenNormal,vec2(x,y)).xyz;
                float weight = clamp(pow(dot(n, l),20.0),0.0, 1.0);
                if(weight < 0.99)
                    continue;

                if(curv > 0.0)
                    intensity += weight;
                else
                    intensity -= weight*100.0;
            }
        if(intensity>0.0)
            intensity = 1.0;
        else
            intensity = 0.0;
        gl_FragData[0] = vec4(intensity);
    }
}
