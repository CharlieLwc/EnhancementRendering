// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable

uniform sampler2D tempResult;   
uniform sampler2D targetTex;     

uniform int phase;        // abs of the normal

uniform float xShift;
uniform float yShift;
uniform vec3 averageNormal;



void main(void) {

    if(phase == 0)
    {
    //calculate sum
        if(gl_TexCoord[0].x<xShift && gl_TexCoord[0].y<yShift)
        {
            vec4 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy);
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift));
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, 0.0));
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, yShift));

            gl_FragData[0]  = element[0]+element[1]+element[2]+element[3];
        }
        else
            gl_FragData[0] =  vec4(0.0);

    }
    else if(phase == 1)
    {
    //normal - average normal
    //the first loop of calculate average

//        vec3 element = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
//
//        if(element != vec3(0.0))
//            gl_FragData[0] = vec4(element - averageNormal, 1.0);
//        else
//            gl_FragData[0] = vec4(0.0);



        if(gl_TexCoord[0].x<xShift && gl_TexCoord[0].y<yShift)
        {
            vec3 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift)).xyz;
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, 0.0)).xyz;
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, yShift)).xyz;
            for(int i=0; i<4; i++)
                if(element[i] != vec3(0.0))
                {
                    element[i] -= averageNormal;
                    element[i] *= element[i].x;
                }
            gl_FragData[0] = vec4(element[0]+element[1]+element[2]+element[3], 1.0);
        }
        else if(gl_TexCoord[0].x > 1.0-xShift && gl_TexCoord[0].y<yShift)
        {
            vec3 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift)).xyz;
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(-xShift, 0.0)).xyz;
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(-xShift, yShift)).xyz;
            for(int i=0; i<4; i++)
                if(element[i] != vec3(0.0))
                {
                    element[i] -= averageNormal;
                    element[i] *= element[i].y;
                }
            gl_FragData[0] = vec4(element[0]+element[1]+element[2]+element[3], 1.0);
        }
        else if(gl_TexCoord[0].x<xShift && gl_TexCoord[0].y > 1.0-yShift)
        {
            vec3 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, -yShift)).xyz;
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, 0.0)).xyz;
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, -yShift)).xyz;
            for(int i=0; i<4; i++)
                if(element[i] != vec3(0.0))
                {
                    element[i] -= averageNormal;
                    element[i] *= element[i].z;
                }
            gl_FragData[0] = vec4(element[0]+element[1]+element[2]+element[3], 1.0);
        }
        else 
            gl_FragData[0] =  vec4(0.0);

    }
    else if(phase == 2)
    {
    //calculate average

        if(gl_TexCoord[0].x<xShift && gl_TexCoord[0].y<yShift)
        {
            vec3 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift)).xyz;
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, 0.0)).xyz;
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, yShift)).xyz;

            gl_FragData[0] = vec4(element[0]+element[1]+element[2]+element[3], 1.0);
        }
        else if(gl_TexCoord[0].x > 1.0-xShift && gl_TexCoord[0].y<yShift)
        {
            vec3 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift)).xyz;
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(-xShift, 0.0)).xyz;
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(-xShift, yShift)).xyz;

            gl_FragData[0] = vec4(element[0]+element[1]+element[2]+element[3], 1.0);
        }
        else if(gl_TexCoord[0].x<xShift && gl_TexCoord[0].y > 1.0-yShift)
        {
            vec3 element[4];
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy).xyz;
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, -yShift)).xyz;
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, 0.0)).xyz;
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, -yShift)).xyz;

            gl_FragData[0] = vec4(element[0]+element[1]+element[2]+element[3], 1.0);
        }
        else 
            gl_FragData[0] =  vec4(0.0);
    }




}

    
    
