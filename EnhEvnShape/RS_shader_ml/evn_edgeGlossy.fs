// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable 


uniform sampler2D ridgeLines;
uniform sampler2D silhouette;
uniform sampler2D normal;


uniform float shineness;     // 1/height
uniform float gamma;     // 1/height

uniform int glossySize;     // 1/width
uniform float sw;     // 1/width
uniform float sh;     // 1/height



int calStart(float start, float center, float step)
{
    int temp = int(center/step);
    if(start < 0.0)
        return glossySize - temp;
    else
        return 0;
}

bool checkBolck(vec2 start, float x, float y)
{
    float startx = start.x;
    float starty = start.y;
    
    float tempx = x - startx;
    float tempy = y - starty;
    
    int xnum = int(abs(tempx/sw));
    int ynum = int(abs(tempy/sh));

    int num;

    float xstep = sw;
    float ystep = sh;

    if(xnum < ynum)
    {
        num = ynum;
        if(tempy < 0) ystep = -ystep;
        xstep = tempx * ystep/tempy;
    }
    else
    {
        num = xnum;
        if(tempx < 0) xstep = -xstep;
        ystep = tempy * xstep/tempx;
    }

    for(int i=0; i<num; i++)
    {
        starty += ystep;
        if(texture2D(ridgeLines,vec2(startx, starty)).y > 0.5 || texture2D(silhouette,vec2(startx, starty)).x<0.9)
            return false;

        startx += xstep;
        if(texture2D(ridgeLines,vec2(startx, starty)).y > 0.5 || texture2D(silhouette,vec2(startx, starty)).x<0.9)
            return false;
    }


    return true;
}


vec3 calColor(vec3 col)
{

    float halfSW = sw*0.5;
    float halfSH = sh*0.5;

    vec2 coord = vec2(gl_TexCoord[0].s, gl_TexCoord[0].t);
    float length2 = glossySize * glossySize * sw * sw; 

    float startx = coord.x - sw * glossySize;
    float starty = coord.y - sh * glossySize; 
    float endx = coord.x + sw * glossySize;
    float endy = coord.y + sh * glossySize;

    if(endx > 1.0) endx = 1.0;
    if(endy > 1.0) endy = 1.0;
    if(startx < halfSW) startx = halfSW;
    if(starty < halfSH) starty = halfSH;

    float illuminate = 0.0;
    

    for(float x = startx; x<endx; x+= sw)
    {
        float tempx = abs(coord.x - x);
        float xx = tempx*tempx;
        for(float y = starty; y<endy; y+=sh)
        {
            float tempy = abs(coord.y - y);
            if(xx + tempy* tempy <= length2)
            {
                vec4 ridge = texture2D(ridgeLines,vec2(x, y));

                if(ridge.x > 0.5)
                {
                    if(checkBolck(coord, x, y))
                    {
                        vec3 nor2 = texture2D(normal,vec2(x,y)).xyz;
                        illuminate += pow(dot(nor2, col),shineness);  
                    
                    }
                
                }
            }      
        }
    }

    if(illuminate == 0)
        return col;
    else
        return vec3(illuminate*gamma);

}




void main(void) {
  
    vec2 coord = vec2(gl_TexCoord[0].s, gl_TexCoord[0].t);


    vec3 nor = texture2D(normal,coord).xyz;
    vec3 col = nor;
    if(nor != vec3(0.0))
    {
    
        vec3 ridge = texture2D(ridgeLines,coord).xyz;
        if(texture2D(silhouette,coord).x < 0.9)
            col = vec3(0.0);
        else if(ridge.y > 0.5)
            col = ridge;
        else if(ridge.x > 0.5)
            col = ridge;
        else
            col = calColor(col);
    }
    
    gl_FragData[0] = vec4(col, 0.0);

}
