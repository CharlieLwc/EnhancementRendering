// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable
 
uniform sampler2D normal;   
uniform sampler2D view;   
uniform sampler2D tempResult;  



uniform int phase;        // abs of the normal

uniform float xShift;
uniform float yShift;
uniform vec3 glossyDirection;



void main(void) {

    if(phase == 0)
    {     
        if ((gl_TexCoord[0].x<xShift && gl_TexCoord[0].y<yShift) || (gl_TexCoord[0].x > 1.0-xShift && gl_TexCoord[0].y<yShift))
        {
            float xTemp = gl_TexCoord[0].x<xShift? xShift : -xShift;
            float yTemp = gl_TexCoord[0].y<yShift? yShift : -yShift;
            
            vec4 element[4];
            float weight[4];
            vec3 pointNormals;
            
            
            pointNormals = texture2D(normal, gl_TexCoord[0].xy).xyz;
            weight[0] = max(dot(pointNormals, glossyDirection), 0.0);

            pointNormals = texture2D(normal, gl_TexCoord[0].xy+vec2(0.0, yTemp)).xyz;
            weight[1] = max(dot(pointNormals, glossyDirection), 0.0);

            pointNormals = texture2D(normal, gl_TexCoord[0].xy+vec2(xTemp, 0.0)).xyz;
            weight[2] = max(dot(pointNormals, glossyDirection), 0.0);

            pointNormals = texture2D(normal, gl_TexCoord[0].xy+vec2(xTemp, yTemp)).xyz;
            weight[3] = max(dot(pointNormals, glossyDirection), 0.0);
            
            

            if(gl_TexCoord[0].x<xShift)
            {
                element[0] = texture2D(normal, gl_TexCoord[0].xy);
                element[1] = texture2D(normal, gl_TexCoord[0].xy+vec2(0.0, yTemp));
                element[2] = texture2D(normal, gl_TexCoord[0].xy+vec2(xTemp, 0.0));
                element[3] = texture2D(normal, gl_TexCoord[0].xy+vec2(xTemp, yTemp));
            }
            else
            {
                element[0] = -texture2D(view, gl_TexCoord[0].xy);
                element[1] = -texture2D(view, gl_TexCoord[0].xy+vec2(0.0, yTemp));
                element[2] = -texture2D(view, gl_TexCoord[0].xy+vec2(xTemp, 0.0));
                element[3] = -texture2D(view, gl_TexCoord[0].xy+vec2(xTemp, yTemp));
            }


            gl_FragData[0]  = element[0]*weight[0]+element[1]*weight[1]+element[2]*weight[2]+element[3]*weight[3];

        }
        else 
            gl_FragData[0] =  vec4(0.0);

            


    }
    else if(phase == 1)
    {
    //calculate sum
        vec4 element[4];
        
        

        if(gl_TexCoord[0].x<xShift && gl_TexCoord[0].y<yShift)
        {
            
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy);
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift));
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, 0.0));
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(xShift, yShift));
            gl_FragData[0]  = element[0]+element[1]+element[2]+element[3];

        }
        else if(gl_TexCoord[0].x > 1.0-xShift && gl_TexCoord[0].y<yShift)
        {
            element[0] = texture2D(tempResult, gl_TexCoord[0].xy);
            element[1] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(0.0, yShift));
            element[2] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(-xShift, 0.0));
            element[3] = texture2D(tempResult, gl_TexCoord[0].xy+vec2(-xShift, yShift));
            gl_FragData[0]  = element[0]+element[1]+element[2]+element[3];

        }
        else
        {
            gl_FragData[0] =  vec4(0.0);
        
        }


    }

}

    
    
