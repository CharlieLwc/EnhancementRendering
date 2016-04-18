// Initial software: Radiance Scaling Version 1.0
// Co-authors: Romain VERGNE, Romain PACANOWSKI, Pascal BARLA, Xavier GRANIER and Christophe SCHLICK.
// Owners: INRIA, University of Bordeaux 1 and University of Bordeaux 2.
// Copyright © 2009-2010, spread under the terms and conditions of the license CeCILL B Version 2.0.

#version 120
#extension GL_ARB_draw_buffers : enable 

uniform sampler2D worldNormal;
uniform sampler2D worldView;
uniform sampler2D worldPos;
uniform sampler2D curvDir;
uniform sampler2D modelColor;



uniform bool needColor;


uniform sampler2D areaLight;
uniform float xIncrement;
uniform float yIncrement;

uniform int materialType;     // 1/height
uniform bool isPointLight;


uniform float roughness;     // 1/height
uniform float fresnel;     // 1/height


uniform float diffuseValue;
uniform float specularValue;
uniform vec3 diffColor;
uniform vec3 specColor;


uniform float alphaX;
uniform float alphaY;
uniform vec3 iosDir;

uniform float matelness;
uniform float edginess;
uniform float backscatter;

uniform float fTransparency;













vec3 ATI(vec3 normal, vec3 view, vec3 light)
{


//    cheating
    float LdotN =  dot( light, normal ); 

    float subLamb = smoothstep(-fresnel,1.0,LdotN) - smoothstep(0.0,1.0,LdotN);
    subLamb = max(0.0,subLamb);
    vec3 subContrib = vec3(subLamb);

    
    float diffComp = max(0.0,LdotN);
    float NdotV = 1.0-dot(view,normal);

    vec3 result = subContrib + diffComp * diffColor*diffuseValue + NdotV*specColor*specularValue;



    return result;




}

vec3 SHW(vec3 normal, vec3 view, vec3 light)
{

//	Translated from HLSL shader for VVVV
//	vNoiseVelvety.fx
//	by Desaxismundi 2008
// 
//	Velvet shader originally by ATI
// 
//	GLSL conversion toneburst 2008

    // Retroreflective lobe



  float diffuseT = clamp(dot(light, normal), 0.0, 1.0);

  float cosine = clamp(dot(light, view), 0.0, 1.0); 
  float shiny = pow(cosine, 1.0 / roughness) * backscatter;

  cosine = clamp(dot(normal, view), 0.0, 1.0);
  float sine = sqrt(1.0 - cosine);
  shiny = shiny + pow(sine, edginess) * diffuseT;

  return (diffColor*diffuseT* diffuseValue) + (specColor*shiny* specularValue);

}

vec3 Ashikhmin_Shirley(vec3 normal, vec3 view, vec3 light)
{

    float alphaXTrem = alphaX*128;
    float alphaYTrem = alphaY*128;
 
    // Define the coordinate frame 
    vec3 halfwayVector = normalize(light + view);
 
    vec3 tangent = normalize(iosDir-normal * dot(normal, iosDir));

    vec3 bitangent = normalize(cross(tangent, normal));

    // Generate any useful aliases
    float VdotN = dot( view, normal );
    float LdotN =  dot( light, normal ); 
    
    float Pd = 0.0;
    float Ps = 0.0;


        float HdotN = max(dot( halfwayVector, normal ), 0.f);
        float HdotL = dot( halfwayVector, light );
        float HdotT = dot( halfwayVector, tangent );
        float HdotB = dot( halfwayVector, bitangent );
  
        // Compute the diffuse term
        Pd = (28.0f * diffuseValue) / ( 23.0f * 3.14159f );
        Pd *= (1.0f - specularValue);
        Pd *= (1.0f - pow(1.0f - (LdotN / 2.0f), 5.0f));
        Pd *= (1.0f - pow(1.0f - (VdotN / 2.0f), 5.0f));
 

        // Compute the specular term

        float fr = ( specularValue + (1.0f - specularValue) * pow( 1.0f - HdotL, 5.0f ) );

        float ps_num_exp = alphaXTrem * HdotT * HdotT + alphaYTrem * HdotB * HdotB;
        ps_num_exp /= (1.0f - HdotN * HdotN);
 
        float Ps_num = sqrt( (alphaXTrem + 1) * (alphaYTrem+ 1) );
        Ps_num *= pow( HdotN, ps_num_exp );
 
        float Ps_den = 8.0f * 3.14159f * HdotL + 0.000001;
        Ps_den *= max( LdotN, VdotN ) + 0.000001;
 
        Ps = specularValue * (Ps_num / Ps_den);
        Ps *= fr;
 
        // Composite the final value:
    
    
    

    Pd = clamp(Pd, 0.f, 1.f);

    Ps = clamp(Ps, 0.f, 1.f);
    vec3 finalValue =  diffColor * Pd + specColor * Ps;


    return finalValue;


}


vec3 Ward(vec3 normal, vec3 view, vec3 light)
{
    
    const float PI = 3.14159;


    // renormalize
    vec3 tangent = normalize(cross(iosDir, normal));

    vec3 bitangent = normalize(cross(tangent, normal));
    
        // calculate intermediary values
    vec3 halfwayVector = normalize(light + view);
    float HdotN = dot(halfwayVector, normal);
    float VdotN = max(dot(view, normal), 0.f);
    float HdotTAlphaX = dot(halfwayVector, tangent) / alphaX;
    float HdotBAlphaY = dot(halfwayVector, bitangent) / alphaY;


    // calculate N dot L
    float NdotL = max(dot(normal, light), 0.f);
    

    float diffuse = diffuseValue * NdotL;
    
        
    // calculate the specular value
    float exponent = -2.0 * (HdotTAlphaX * HdotTAlphaX + HdotBAlphaY * HdotBAlphaY) / (1.0 + HdotN);
    // Evaluate the specular denominator
    float s_den  = 4.0f * 3.14159f; 
    s_den       *= alphaX;
    s_den       *= alphaY;
    s_den       *= sqrt( NdotL * VdotN )+0.00001;

    vec3 specular = specColor * ( exp( exponent ) / s_den ) * specularValue;

    
    vec3 finalValue = (diffColor * diffuseValue + specular) * NdotL;
	


    return finalValue;

}



vec3 Strauss(vec3 normal, vec3 view, vec3 light)
{

    // Declare any aliases:
    vec3 halfVector = reflect( light, normal );

    float NdotL   = dot( normal, light );
    float NdotV   = dot( normal, view );
    float HdotV   = dot( halfVector, view );
    
        
    // Schlick approximation
    float fNdotL = 0.014567258064516/(NdotL - 1.12)/(NdotL - 1.12) - 0.0116129032258065;
    
    
    
    float s_cubed = roughness * roughness * roughness;
 
    // Evaluate the diffuse term
    float d  = ( 1.0f - matelness * roughness );
    float Rd = ( 1.0f - s_cubed ) * ( 1.0f - fTransparency );
    vec3 diffuse = diffColor * NdotL * d * Rd;
 
    // Compute the inputs into the specular term
    float r = ( 1.0f - fTransparency ) - Rd;
 
    float shadowNdotL = 0.0001*( 10000.0 - 1.0/(NdotL-1.01)/(NdotL-1.01))*1.000098039215686;
    float shadowNdotV = 0.0001*( 10000.0 - 1.0/(NdotV-1.01)/(NdotV-1.01))*1.000098039215686;

    float j = fNdotL * shadowNdotL * shadowNdotV;
 
    // 'k' is used to provide small off-specular
    // peak for very rough surfaces. Can be changed
    // to suit desired results...
        const float k = 0.1;
    float reflect = min( 1.0, r + j * ( r + k ) );
 
    vec3 C1 = vec3( 1.0, 1.0, 1.0 );
    vec3 Cs = C1 + (diffColor - C1) * matelness * (1.0 - fNdotL);
 
    // Evaluate the specular term
    vec3 specular = Cs * reflect;
    specular *= pow( -HdotV, 3.0 / (1.0 - roughness) );
     // Composite the final result, ensuring
    // the values are >= 0.0f yields better results. Some
    // combinations of inputs generate negative values which
    // looks wrong when rendered...
    diffuse  = max( vec3(0.0), diffuse );
    specular = max( vec3(0.0), specular );

    return diffuse + specular;
}

vec3 Oren_Nayar(vec3 normal, vec3 view, vec3 light)
{

·
    const float PI = 3.14159;
    
    float NdotL = dot(normal, light);
    float NdotV = dot(normal, view); 

    float angleVN = acos(NdotV);
    float angleLN = acos(NdotL);
    



    float alpha = max(angleVN, angleLN);
    float beta = min(angleVN, angleLN);
    float gamma = dot(view - normal * NdotV, light - normal * NdotL);
    
    float roughnessSquared = roughness * roughness;   

    
    // calculate A and B
    float A = 1.0 - 0.5 * (roughnessSquared / (roughnessSquared + 0.57));

    float B = 0.45 * (roughnessSquared / (roughnessSquared + 0.09));
 
    float C = sin(alpha) * tan(beta);
    
    // put it all together
    float L1 = max(0.0, NdotL) * (A + B * max(0.0, gamma) * C);
    
    // get the final color 
    vec3 finalValue = diffColor * diffuseValue * L1;
    
    return finalValue;


}

vec3 Cook_Torrance(vec3 normal, vec3 view, vec3 light)
{
    // fraction of diffuse reflection (specular reflection = 1 - k)
    
    // do the lighting calculation for each fragment.
    float NdotL = max(dot(normal, light), 0.0);
    
    float specular = 0.0;
    if(NdotL > 0.0 && roughness > 0.0)
    {
        // calculate intermediary values
        vec3 halfVector = normalize(light + view);
        float NdotH = max(dot(normal, halfVector), 0.0); 
        float NdotV = max(dot(normal, view), 0.0); // note: this could also be NdotL, which is the same value
        float VdotH = max(dot(view, halfVector), 0.0);
        float mSquared = roughness * roughness;
        
        // geometric attenuation
        float geo_numerator = 2.0 * NdotH;
        float g1 = (geo_numerator * NdotV) / VdotH;
        float g2 = (geo_numerator * NdotL) / VdotH;
        float geoAtt = min(1.0, min(g1, g2));
     
        // roughness (or: microfacet distribution function)
        // beckmann distribution function
        float r1 = 1.0 / ( 4.0 * mSquared * pow(NdotH, 4.0));
        float r2 = (NdotH * NdotH - 1.0) / (mSquared * NdotH * NdotH);
        float roughTerm = r1 * exp(r2);
        
        // fresnel
        // Schlick approximation
        float fresnelTerm = pow(1.0 - VdotH, 5.0);
        fresnelTerm *= (1.0 - fresnel);
        fresnelTerm += fresnel;
        
        specular = (fresnelTerm * geoAtt * roughTerm) / (NdotV * NdotL * 3.14);
    }
    
    vec3 finalValue =  NdotL * (diffColor * diffuseValue + specular * specularValue * specColor);

    return finalValue;


}
vec3 Blinn_Phong(vec3 normal, vec3 view, vec3 light)
{

    // Compute the half vector
    vec3 half_vector = normalize(light + view);
 
    // Compute the angle between the half vector and normal
    float  HdotN = max( 0.0f, dot( half_vector, normal ) );
 
    // Compute the specular colour
    vec3 specular = specColor * pow( HdotN, 1.0/roughness ) * specularValue;
 
    // Compute the diffuse term
    vec3 diffuse = diffColor * max( 0.0f, dot( normal, light ) ) * diffuseValue;
 
    // Determine the final colour    
    return diffuse + specular;
}

void main(void) {
  

    vec2 coord = vec2(gl_TexCoord[0].s, gl_TexCoord[0].t);


    vec3 normal = texture2D(worldNormal,coord).xyz;
    vec3 view = -texture2D(worldView,coord).xyz;
    vec3 pos = texture2D(worldPos,coord).xyz;





    vec3 light;


    
    vec3 col = vec3(0.0); 
    if(needColor)
        col = texture2D(modelColor,coord).xyz;

    if(normal != vec3(0.0))
    {  
        int lightNum = 0;
        for(float tx = xIncrement*0.5; tx<1.0; tx+=xIncrement)
        {
            for(float ty = yIncrement*0.5; ty<1.0; ty+=yIncrement)
            {
            
                light = texture2D(areaLight,vec2(tx,ty)).xyz;
        
                if(isPointLight)
                    light = normalize(light - pos);
                else
                    light = light;
        
                switch(materialType)
                {
                    case 0: col += Blinn_Phong(normal, view, light); break;
                    case 1: col += Cook_Torrance(normal, view, light);break;
                    case 2: col += Oren_Nayar(normal, view, light); break;
                    case 3: col += Strauss(normal, view, light);break;
                    case 4: col += Ward(normal, view, light);break;
                    case 5: col += Ashikhmin_Shirley(normal, view, light);break;
                    case 6: col += SHW(normal, view, light);break;
                    case 7: col += ATI(normal, view, light);break;
                }
                lightNum++;
            }
        
        }
        col = col/lightNum;
    }



    gl_FragData[0] = vec4(col, 0.0);

}
