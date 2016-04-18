#include "stdafx.h"

#include "barItems.h"


#define PI 3.141592654f




/*************************************************************************/
/****************************    Parameter    ******************************/
/*************************************************************************/

int PARA_vertexIndex = 0;
int PARA_sampleNum = 11;
int PARA_size = 10;
float PARA_ridgeThreshold = 2.4f;
float PARA_angle = 0.f;
float PARA_transX = 0.f;
float PARA_transY = 0.f;
float PARA_scale = 1.f;

bool PARA_saveMatrix;
bool PARA_loadMatrix;

bool PARA_useFaceNormal = false;

/*************************************************************************/
/****************************    Draw_Ori    ******************************/
/*************************************************************************/

bool draw_shiny = true;
bool draw_light = true;
bool draw_2side = false;
bool draw_white_bg = false;
bool draw_points = false;
bool draw_edges = false;
bool draw_ridge = false;
bool draw_falsecolor = false;
bool draw_colorSpace = false;


enum    DrawMod { Faces, Edges, Points };
DrawMod draw_drawMod;

int draw_point_size = 1;
int draw_point_size_space = 1;
int draw_line_width = 1;


/*************************************************************************/
/****************************   ModelSeg   ******************************/
/*************************************************************************/

bool MS = false;
bool MS_needColor = false;
bool ms_needRefesh = false;
bool ms_needRefreshParmeters = false;
int ms_layerIndex = 0;
bool ms_showPatch = false;
bool ms_showNormal = false;
bool ms_showNormalColor = false;
bool ms_showSmoothedNormal = false;
bool ms_showSmoothedNormalColor = false;
bool ms_showCurvetureColor = false;
int ms_smoothRadius = 1.f;
float ms_clusterAngle = 11.f;
float ms_maxClusterSize = 0.f;
float ms_minClusterSize = 0.f;





bool ms_useCluster = true;
float ms_smoothness = 1.0;

bool ori = false;
bool clu = false;
bool edg = false;
bool our = false;



/*************************************************************************/
/****************************   Curvature   ******************************/
/*************************************************************************/
bool CURV = false;
bool draw_curvature = false;
bool draw_dcurvature = false;
bool curv_needRefresh = false;

float curv_thre_pos = 1.0;
float curv_thre_neg = 1.0;
float curv_scale = 1.0;


/*************************************************************************/
/****************************      PRT      ******************************/
/*************************************************************************/
bool PRT = false;

bool prt_isTexInit = false;
bool prt_needRefresh = false;
bool prt_needDesignReflect = false;
bool prt_needDesignRefract = false;
bool prt_needRefreshTexture = false;

bool prt_enhSingleFLE = false;
bool prt_enhVertexFLE = false;
bool prt_enhUsing = false;
bool prt_enhSingleFRA = false;
bool prt_enhAllFRA = false;

float prt_sharpnessFLE = 1.0;



bool prt_showLAB = false;
vec3 prt_eyeDirection = vec3(0.f);
vec3 prt_left = vec3(0.f);


bool prt_UsingRawReflect = false;
bool prt_UsingRawRefract = false;
float prt_scaleRadiance = 1.f;

int prt_envPicIndex = 0;

GLuint texCubeMap[6];

//enum    ColorSpaceMod { CS_NO, CS_LAB, CS_RGB, CS_HSV, CS_SRGB, CS_XYZ };



Color::Colorspace PRT_colSpaMod = Color::Colorspace::YCBCR;;
Color::Colorspace PRT_sourceColSpaMod = Color::Colorspace::YCBCR;;

EyeType prt_eyeDir;
MaterialType prt_material;

bool prt_GLSL_addMaterial = false;

bool prt_GLSL_CustomLightColor = false;


enum EvnMap { EM_CUR, EM_DES, EM_ORI, EM_AVG, EM_NONE };
EvnMap PRT_EvnMapType = EM_CUR;
EvnMap PRT_ColSpaSourEvnMapType = EM_NONE;


/*************************************************************************/
/****************************      GLSL     ******************************/
/*************************************************************************/
bool GLSL = false;
bool GLSL_reflection = false;
bool GLSL_refraction = false;
bool GLSL_radianceScaling = false;
bool GLSL_translucent = false;


bool GLSL_CurvNeedRefrash = true;

bool GLSL_showCurv = false;

bool GLSL_showNormal = false;



float GLSL_layerWeight[5] = { 0.f, 0.f, 0.f, 0.f, 1.f };

bool GLSL_needRidgeInfo = false;
bool translucent = false;

float GLSL_enhanceScale = 1.f;

float GLSL_diffuseValue = 0.5f;
float GLSL_specularValue = 0.5f;
float GLSL_hightlightThre = 0.015;
float GLSL_shiness = 10.f;
float GLSL_transparency = 0.5;

float GLSL_silThre = 0.05;
float GLSL_silHardness = 0.001;

vec3 GLSL_diffColor = vec3(1.0, 1.0, 1.0);
vec3 GLSL_specColor = vec3(1.0, 1.0, 1.0);

bool GLSL_enhancePRT = false;
bool GLSL_silhouette = false;

enum    CurMod { ScreenSpace, ObjectSpace, SurfaceReliefFeature };
CurMod GLSL_curvMode;
enum    MoveMod { MoveObject, MoveEye };
CurMod GLSL_moveMode;

//Blinn - Phong		plastics
//Cook_Torrance		metals
//Oren_Nayar		Diffuse - only rough surfaces. brick, concrete and also clothing materials such as velvet.
//Strauss			Easy to control
//Ward				Clothing, for example velvet. 
//					Machined metals such as brushed steel.
//					Wood grain. Although the colour and pattern can be expressed fairly easily it is difficult to express the lighting characteristics without considering the effect the grain has.
//					Paint and varnish applied using a paintbrush will usually leave a directional pattern matching the brush - stroke.
//Ashikhmin_Shirley metals and plastics




enum    MaterialMod { Blinn_Phong, Cook_Torrance, Oren_Nayar, Strauss, Ward, Ashikhmin_Shirley, SHW, ATI };
MaterialMod GLSL_materialMod;
enum    LightMod { Light_Direct, Light_Point, Light_GeoDirect, Light_GeoPoint, Light_CamDirect, Light_CamPoint };
LightMod GLSL_lightMod;

vec3 GLSL_lightDirection = vec3(1.f, 0.f, 0.f);
float GLSL_Roughness = 0.5;
float GLSL_Metalness = 0.5;
float GLSL_Edginess = 4.0;
float GLSL_Backscatter = 0.25;

float GLSL_fresnel = 0.8;

float GLSL_alphaX = 0.5;
float GLSL_alphaY = 0.13;
vec3 GLSL_iosDirection = vec3(1.f, 0.f, 0.f);

float GLSL_lightSizeX = 0.0;
float GLSL_lightSizeY = 0.0;
float GLSL_lightSolution = 0.1;


bool GLSL_ExaggeratedShading = false;
bool GLSL_ExaggeratedShadingOur = false;




/*************************************************************************/
/************************      SurfaceRelief     *************************/
/*************************************************************************/

bool SR = false;
bool SR_showBaseNormal = false;
int SR_normalSmoothRadiance = 0;
bool SR_Height = false;

bool SR_showVertexHeight = false;
bool SR_needRefesh = false;
bool SR_surfaceRelif = false;



/*************************************************************************/
/************************      CurvAndPlan     *************************/
/*************************************************************************/

bool CP = false;
bool CP_showCPsingle = false;
bool CP_showCPall = false;
bool CP_needRefesh = false;

int CP_ridgeIndex;
float CP_curvSize = 0.1;
float CP_planSize = 0.1;



/*************************************************************************/



bool refraction = false;
bool draw_index = false;
bool UsingSH = false;
bool UsingWavelet = false;



bool UsingGlossReflect = false;
bool UsingGlossRefract = false;


bool isIndirect = false;


bool draw_faded = true;
bool test_rv = true;

float sharpnessRefract = 1.0;
float enhScaleReflect = 1.0;
float enhScaleRefract = 1.0;
float gammaAverage = 1.f;
float gammaScaleReflect = 1.f;
float gammaScaleRefract = 1.f;

float feature_size;
float rv_thresh = 0.1;


float m_refindex = 1.3f;
float m_sigma_s[3] = { 0.3f, 0.3f, 0.3f };
float m_sigma_a[3] = { 0.001f, 0.001f, 0.001f };



float cubeSize = 1.5f;

float CubeMapCoord[6][4][3] = {

	-cubeSize, cubeSize, -cubeSize, cubeSize, cubeSize, -cubeSize, cubeSize, cubeSize, cubeSize, -cubeSize, cubeSize, cubeSize,			//T


	-cubeSize, cubeSize, -cubeSize, -cubeSize, cubeSize, cubeSize, -cubeSize, -cubeSize, cubeSize, -cubeSize, -cubeSize, -cubeSize,		//L


	-cubeSize, cubeSize, cubeSize, cubeSize, cubeSize, cubeSize, cubeSize, -cubeSize, cubeSize, -cubeSize, -cubeSize, cubeSize,			//F


	cubeSize, cubeSize, cubeSize, cubeSize, cubeSize, -cubeSize, cubeSize, -cubeSize, -cubeSize, cubeSize, -cubeSize, cubeSize, 		//R


	-cubeSize, -cubeSize, cubeSize, cubeSize, -cubeSize, cubeSize, cubeSize, -cubeSize, -cubeSize, -cubeSize, -cubeSize, -cubeSize,		//D


	-cubeSize, -cubeSize, -cubeSize, cubeSize, -cubeSize, -cubeSize, cubeSize, cubeSize, -cubeSize, -cubeSize, cubeSize, -cubeSize		//B

};

float texCoord[4][2] = { 0.f, 0.f, 1.f, 0.f, 1.f, 1.f, 0.f, 1.f };

vec3 LABmin(0.f, -85.f, -106.f);
vec3 LABmax(98.f, 96.6f, 92.9f);

void drawEyedirOnScreen(TriMesh::BSphere &global_bsph, xform &global_xf)
{
	vec3 eye(0.f, 0.f, 0.f);


	eye = inv(global_xf) * eye;
	prt_eyeDirection = eye - global_bsph.center;
	normalize(prt_eyeDirection);


	vec3 cen = global_xf*global_bsph.center;

	vec3 upRight(0.f, global_bsph.r, 0.f);
	upRight = cen + upRight;
	upRight = inv(global_xf) * upRight;


	vec3 left(-global_bsph.r, 0.f, 0.f);
	left = cen + left;
	left = inv(global_xf) * left;
	prt_left = left - global_bsph.center;
	normalize(prt_left);

	
	glLineWidth(5.0);

	/*
	glBegin(GL_LINES);

	glVertex3fv(global_bsph.center);
	glVertex3fv(eye);

	glVertex3fv(global_bsph.center);
	glVertex3fv(upRight);

	glVertex3fv(global_bsph.center);
	glVertex3fv(left);

	glEnd();

	

	vec3 tempAixeZ(prt_eyeDirection[0], 0.f, prt_eyeDirection[2]);
	float xLength = len(tempAixeZ);
	float yLength = prt_eyeDirection[1];
	vec3 tempAixeY = tempAixeZ * (-yLength / xLength);
	tempAixeY[1] = xLength;
	vec3 tempAixeX = tempAixeY.cross(prt_eyeDirection);
	
	glBegin(GL_LINES);
	glVertex3fv(global_bsph.center);
	glVertex3fv(global_bsph.center + tempAixeY * global_bsph.r);
	glVertex3fv(global_bsph.center);
	glVertex3fv(global_bsph.center + tempAixeX * global_bsph.r);
	glEnd();
	*/
	
}

void SLmatrix(xform &global_xf)
{
	if (PARA_saveMatrix)
	{
		PARA_saveMatrix = false;
		
		global_xf.write(meshes[0]->dFName.getDN(DFNameType::BestPosition));

	}

	if (PARA_loadMatrix)
	{
		PARA_loadMatrix = false;

		global_xf.read(meshes[0]->dFName.getDN(DFNameType::BestPosition));
	}






}

void initeCubeMapTexture()
{
	if (prt_isTexInit)
		return;
	prt_isTexInit = true;
	//generate texture 
	//space
	glGenTextures(6, texCubeMap);
	for (int i = 0; i < 6; i++)
	{
		glBindTexture(GL_TEXTURE_2D, texCubeMap[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	}
}

void meshSegPreparetion()
{

	if (ms_needRefesh)
	{



		ms_needRefesh = false;
		if (ms_showPatch)
			meshes[0]->setColorCluster(ms_layerIndex);

		if (ms_showNormal)
			meshes[0]->setNormalClusterNormal(ms_layerIndex);

		if (ms_showNormalColor)
			meshes[0]->setColorClusterNormal(ms_layerIndex);


		if (ms_showSmoothedNormalColor)
			meshes[0]->setColorSmoothedNormal(ms_layerIndex, ms_smoothRadius, ms_useCluster, ms_smoothness);

		if (ms_showCurvetureColor)
			meshes[0]->getColorSmoothedCurvature(ms_layerIndex, ms_smoothRadius, curv_thre_pos, curv_thre_neg, curv_scale);

	}

	if (ms_needRefreshParmeters)
	{
		int downLayer = ms_layerIndex > 0 ? ms_layerIndex - 1 : 0;
		ms_needRefreshParmeters = false;
		meshes[0]->need_cluster(ms_layerIndex);
		ms_smoothRadius = meshes[0]->smoothRadius[ms_layerIndex];
		ms_clusterAngle = meshes[0]->clusterAngle[ms_layerIndex];
		ms_minClusterSize = meshes[0]->lowClusterScale[downLayer];
		ms_maxClusterSize = meshes[0]->upClusterScale[downLayer];
	}



}

void normalPreparetion()
{
	// Normals

	if ((!meshes[0]->normals.empty() || ms_showSmoothedNormal || ms_showNormal) && !draw_index) {

		glEnableClientState(GL_NORMAL_ARRAY);
		if (ms_showSmoothedNormal)
		{
			meshes[0]->smoothNormal(ms_layerIndex, ms_smoothRadius, ms_useCluster, ms_smoothness);

			glNormalPointer(GL_FLOAT,
				sizeof(meshes[0]->smoothedNormal[ms_layerIndex][0]),
				&meshes[0]->smoothedNormal[ms_layerIndex][0][0]);
		}
		else if (ms_showNormal)
		{

			glNormalPointer(GL_FLOAT,
				sizeof(meshes[0]->clusterNormalForVertex[0]),
				&meshes[0]->clusterNormalForVertex[0][0]);
		}
		else
		{
			glNormalPointer(GL_FLOAT,
				sizeof(meshes[0]->normals[0]),
				&meshes[0]->normals[0][0]);
		}
	}
	else
	{
		glDisableClientState(GL_NORMAL_ARRAY);
	}
}

void curvPreparetion()
{
	if (curv_needRefresh)
	{
		curv_needRefresh = false;
		if (draw_curvature)
			meshes[0]->setColorCurvature(curv_thre_pos, curv_thre_neg, curv_scale);
		if (draw_dcurvature)
			meshes[0]->setColorDcurvature();
	}
}

void surfaceReliefPreparetion(TriMesh::BSphere &global_bsph, xform &global_xf)
{
	if (SR_needRefesh)
	{
		SR_needRefesh = false;

		if (SR_showVertexHeight)
			meshes[0]->setColorVertexSamplePoly(PARA_vertexIndex, PARA_sampleNum, PARA_angle, false);
		if (SR_showBaseNormal)
			meshes[0]->setColorBaseNormal(SR_normalSmoothRadiance);
		if (SR_Height)
			meshes[0]->setColorHeight();
		if (SR_showVertexHeight)
			meshes[0]->setColorVertexSamplePoly(PARA_vertexIndex, PARA_sampleNum, PARA_angle, true);
		if (SR_surfaceRelif)
			meshes[0]->setColorSurfaceRelif(PARA_sampleNum, curv_scale, curv_thre_pos, curv_thre_neg);

		


	}


	//height curve
	if (!SR_needRefesh && SR_showVertexHeight)
	{
		glPushMatrix();
		glLoadIdentity();



		vec3 newCenter = global_xf * global_bsph.center;

		//				glMultMatrixd(global_xf);
		glTranslatef(newCenter[0], newCenter[1], newCenter[2]);
		glColor3f(1.0, 1.0, 1.0);

		float xStart = 0.0;
		float xEnd = 1.2f*global_bsph.r;
		float yStart = -xEnd;
		float yEnd = -0.6*global_bsph.r;
		float xStep = xEnd / (float)(PARA_sampleNum - 1);
		float yStep = yEnd - yStart;

		glRectf(xStart, yStart, xEnd, yEnd);


		glTranslatef(0, 0, global_bsph.r * 0.01);
		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i<PARA_sampleNum; i++)
		{
			float a = meshes[0]->halfEdge.sampleHeight[i];
			//glVertex2f(-90.0+2.0*i,-80.0+1.0*((sampleHeight[i]-0.5)*50+0.5)); ami
			//glVertex2f(-50.0+1.0*i,-60.0+0.7*((sampleHeight[i]-0.5)*50+0.5));//cir
			glVertex2f(xStart + xStep*i, yStart + meshes[0]->halfEdge.sampleHeight[i] * yStep);//gar
		}
		glEnd();
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j<PARA_sampleNum; j++)
		{
			int i = j + PARA_sampleNum;
			//glVertex2f(-90.0+2.0*j,-80.0+1.0*((sampleHeight[i]-0.5)*50+0.5));
			//glVertex2f(-50.0+1.0*j,-60.0+0.7*((sampleHeight[i]-0.5)*50+0.5));
			glVertex2f(xStart + xStep*j, yStart + meshes[0]->halfEdge.sampleHeight[i] * yStep);//gar
		}
		glEnd();

		glPopMatrix();
	}

}

void curveAndPlanePreparetion(TriMesh::BSphere &global_bsph, xform &global_xf)
{

	if (CP_needRefesh && (CP_showCPsingle || CP_showCPall))
	{
		//		meshes[0]->need_curvatures();
		CP_needRefesh = false;
		if (CP_showCPsingle)
			CP_ridgeIndex = meshes[0]->setColorRidgeVertexCurvPlane(PARA_vertexIndex, PARA_sampleNum, curv_scale, PARA_ridgeThreshold, CP_curvSize, 1.f - CP_planSize);
		if (CP_showCPall)
			meshes[0]->setColorRidgeCurvPlane(PARA_sampleNum, curv_scale, PARA_ridgeThreshold, CP_curvSize, 1.f - CP_planSize);

	}


	//height curve
	if (!CP_needRefesh && CP_showCPsingle)
	{
		glPushMatrix();
		glLoadIdentity();

		std::vector<vec2> &leftScrPosition = meshes[0]->halfEdge.leftScrPosition;
		std::vector<vec2> &rightScrPosition = meshes[0]->halfEdge.rightScrPosition;
		std::vector<vec2> &leftScrNormal = meshes[0]->halfEdge.leftScrNormal;
		std::vector<vec2> &rightScrNormal = meshes[0]->halfEdge.rightScrNormal;
		std::vector<bool> &leftSampleIsCurv = meshes[0]->halfEdge.leftIsSampleCurv[CP_ridgeIndex];
		std::vector<bool> &rightSampleIsCurv = meshes[0]->halfEdge.rightIsSampleCurv[CP_ridgeIndex];

		HalfEdge &ccc = meshes[0]->halfEdge;



		vec3 newCenter = global_xf * global_bsph.center;

		//				glMultMatrixd(global_xf);
		glTranslatef(newCenter[0], newCenter[1], newCenter[2]);
		glColor3f(1.0, 1.0, 1.0);

		float xStart = 0.0;
		float xEnd = 1.2f*global_bsph.r;
		float yStart = -xEnd;
		float yEnd = -0.6*global_bsph.r;
		float xStep = xEnd - xStart;
		float yStep = yEnd - yStart;

		glRectf(xStart, yStart, xEnd, yEnd);


		glTranslatef(0, 0, global_bsph.r * 0.01);
		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < leftScrPosition.size(); i++)
		{
			if (i >= leftSampleIsCurv.size() || leftSampleIsCurv[i])
				glColor3f(1.0, 0.0, 0.0);
			else
				glColor3f(0.0, 0.0, 0.0);


			glVertex2f(xStart + leftScrPosition[i][0] * xStep, yStart + leftScrPosition[i][1] * yStep);//gar
		}
		glEnd();


		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < rightScrPosition.size(); i++)
		{
			if (i >= rightSampleIsCurv.size() || rightSampleIsCurv[i])
				glColor3f(1.0, 0.0, 0.0);
			else
				glColor3f(0.0, 0.0, 0.0);

			glVertex2f(xStart + rightScrPosition[i][0] * xStep, yStart + rightScrPosition[i][1] * yStep);//gar
		}
		glEnd();



		glBegin(GL_LINES);
		{

			for (int i = 0; i < leftScrPosition.size(); i++)
			{

				if (i >= leftSampleIsCurv.size() || leftSampleIsCurv[i])
					glColor3f(1.0, 0.0, 0.0);
				else
					glColor3f(0.0, 0.0, 0.0);
				glVertex2f(xStart + leftScrPosition[i][0] * xStep, yStart + leftScrPosition[i][1] * yStep);//gar
				glVertex2f(xStart + leftScrNormal[i][0] * xStep, yStart + leftScrNormal[i][1] * yStep);//gar


			}



			for (int i = 0; i < rightScrPosition.size(); i++)
			{
				if (i >= rightSampleIsCurv.size() || rightSampleIsCurv[i])
					glColor3f(1.0, 0.0, 0.0);
				else
					glColor3f(0.0, 0.0, 0.0);
				glVertex2f(xStart + rightScrPosition[i][0] * xStep, yStart + rightScrPosition[i][1] * yStep);//gar
				glVertex2f(xStart + rightScrNormal[i][0] * xStep, yStart + rightScrNormal[i][1] * yStep);//gar

			}


		}


		glEnd();



		glPopMatrix();
	}
}

void PRTPreparetion(TriMesh::BSphere &global_bsph)
{

	// Cube Map
	if (PRT)
	{
		initeCubeMapTexture();

		if (prt_needDesignReflect)
		{

			prt_needDesignReflect = false;
			if (prt_enhVertexFLE)
				transport.designLightReflectRidgeVertexColor(prt_eyeDirection, PARA_vertexIndex, prt_sharpnessFLE, enhScaleReflect, curv_scale, gammaAverage, gammaScaleReflect, UsingGlossReflect);
			if (PRT_EvnMapType == EM_DES)
				transport.readingExistingLight();

			prt_needRefresh = true;
			prt_needRefreshTexture = true;

		}

		if (prt_needDesignRefract)
		{
			prt_needDesignRefract = false;
			transport.designLightRefract(meshes[0]->curv1, isIndirect, false, PARA_vertexIndex, sharpnessRefract, enhScaleRefract, curv_scale, gammaScaleRefract, UsingGlossRefract);
			prt_needRefresh = true;
			prt_needRefreshTexture = true;

		}


		// calculate radiance using precomputed lighting transport
		if (prt_needRefresh)
		{
			prt_needRefresh = false;

			if (UsingWavelet)
				transport.setColorWav(meshes[0]->getColor(), isIndirect, false, false, false, UsingGlossReflect);
			else if (prt_UsingRawReflect || prt_UsingRawRefract)
				transport.setColorRaw(meshes[0]->getColor(), isIndirect, false, prt_UsingRawReflect, prt_UsingRawRefract, (PRT_EvnMapType == EM_DES || prt_enhVertexFLE),
				prt_enhAllFRA, PARA_vertexIndex, prt_scaleRadiance, UsingGlossRefract, prt_material, prt_eyeDir);
			else if (UsingSH)
				shLight.setColor(meshes[0]->getColor());
		}

		if (prt_needRefreshTexture)
		{
			prt_needRefreshTexture = false;

			if (UsingSH)
				shLight.initTexture(texCubeMap, PRT_EvnMapType == EM_ORI);
			else
				transport.initTexture(texCubeMap, PRT_EvnMapType == EM_ORI, PRT_EvnMapType == EM_AVG, UsingWavelet, (PRT_EvnMapType == EM_DES || prt_enhVertexFLE), prt_enhAllFRA);
		}



		if (PRT_EvnMapType != EM_NONE)
		{
			glEnable(GL_TEXTURE_2D);

			for (int face = 0; face < 6; face++)
			{
				glBindTexture(GL_TEXTURE_2D, texCubeMap[face]);
				glBegin(GL_QUADS);
				for (int i = 0; i < 4; i++)
				{
					glTexCoord2fv(texCoord[i]);
					glVertex3f(CubeMapCoord[face][i][0] * global_bsph.r + global_bsph.center[0], CubeMapCoord[face][i][1] * global_bsph.r + global_bsph.center[1], CubeMapCoord[face][i][2] * global_bsph.r + global_bsph.center[2]);		//LB


				}
				glEnd();
			}

			glDisable(GL_TEXTURE_2D);
		}






	}

}

bool GLSLPreparetion(TriMesh::BSphere &global_bsph, GLCamera &camera, void(*draw_tstrips)(const TriMesh *, const bool), xform &global_xf)
{


	meshes[0]->need_faceNormals();

	if (GLSL)
	{

		bool GLSL_needNormal = false;
		bool GLSL_needObjCurv = false;
		bool GLSL_needScreenCurve = false;

		bool GLSL_needWorldView = false;
		bool GLSL_needWorldNormal = false;
		bool GLSL_needView = false;

		bool GLSL_needWorldPos = false;

		bool GLSL_needDepth = false;
		bool GLSL_needBack = false;
		bool GLSL_needSilhouette = false;

		bool GLSL_needCurDir = false;

		bool GLSL_needFeature = false;

		bool GLSL_needColor = false;

		vec3 lightPosition = GLSL_lightDirection;

		if (GLSL_showCurv || (GLSL_radianceScaling && GLSL_enhanceScale>0.f))
		{
			if (GLSL_curvMode == ScreenSpace)
				GLSL_needScreenCurve = true;
			else if ((GLSL_curvMode == ObjectSpace))
				GLSL_needObjCurv = true;
		}
		

		 


		GLSL_needFeature = GLSL_radianceScaling && GLSL_curvMode == SurfaceReliefFeature;

		GLSL_needView = (GLSL_radianceScaling && (GLSL_moveMode == MoveObject));
		GLSL_needWorldView = (GLSL_radianceScaling && (GLSL_moveMode == MoveEye)) || prt_GLSL_addMaterial;
		GLSL_needWorldNormal = (GLSL_radianceScaling && (GLSL_moveMode == MoveEye)) || prt_GLSL_addMaterial || GLSL_showNormal;

		GLSL_needWorldPos = prt_GLSL_addMaterial;
		
		GLSL_needRidgeInfo = GLSL_enhancePRT;

		GLSL_needSilhouette = (GLSL_silhouette || GLSL_enhancePRT || GLSL_needScreenCurve);

		GLSL_needWorldNormal = GLSL_needWorldNormal || (GLSL_needSilhouette || GLSL_enhancePRT);

		GLSL_needDepth = GLSL_needScreenCurve || GLSL_needSilhouette;

		GLSL_needNormal = GLSL_showNormal || GLSL_radianceScaling || GLSL_needObjCurv || GLSL_needScreenCurve || GLSL_needWorldNormal || GLSL_needDepth || GLSL_needWorldView || GLSL_needWorldPos || GLSL_needCurDir || GLSL_needFeature;
		
		GLSL_needColor = prt_GLSL_addMaterial && PRT;


		if (prt_GLSL_addMaterial)
		{
			normalize(lightPosition);
			if (GLSL_lightMod == Light_CamPoint || GLSL_lightMod == Light_CamDirect || GLSL_lightMod == Light_GeoDirect || GLSL_lightMod == Light_GeoPoint)
			{
				vec3 cen = global_xf*global_bsph.center;
				lightPosition = cen + lightPosition;
				lightPosition = inv(global_xf) * lightPosition;
				lightPosition = lightPosition - global_bsph.center;
				normalize(lightPosition);
			}
		}
		

		float tempLayerweight[5] = { -1.f, -1.f, -1.f, -1.f, -1.f};
		if (GLSL_needNormal)
		{

			int tempLayerIndex[3] = { -1, -1, 4 };
			int layerNum = 0;
			float sumWeight = GLSL_layerWeight[4];
			for (int lIndex = 0; lIndex < 5; lIndex++)
			{
				if (GLSL_layerWeight[lIndex] <= 0.f || (layerNum >= 2 && lIndex != 4))
					GLSL_layerWeight[lIndex] = 0.f;
				else if (lIndex != 4)
				{
					sumWeight += GLSL_layerWeight[lIndex];
					tempLayerIndex[layerNum++] = lIndex;
				}
			}



			sumWeight = 1.0 / sumWeight;
			for (int lIndex = 0; lIndex < 3; lIndex++)
			{
				if (tempLayerIndex[lIndex] >= 0)
					tempLayerweight[lIndex] = GLSL_layerWeight[tempLayerIndex[lIndex]] * sumWeight;
			}



			GLint sn[10];

			if (GLSL_CurvNeedRefrash && GLSL_needObjCurv)
			{
				meshes[0]->needWeightedCurvature(tempLayerweight, curv_scale);
				GLSL_CurvNeedRefrash = false;
			}

			if (GLSL_needCurDir)
				meshes[0]->need_curvatures();

			radianceScaling.getNormal(GLSL_needWorldNormal, GLSL_needView, GLSL_needWorldView, tempLayerIndex, GLSL_needObjCurv, GLSL_needDepth, 
				GLSL_needFeature, GLSL_needBack, GLSL_needWorldPos, GLSL_needCurDir, GLSL_needColor, 
				prt_GLSL_addMaterial, lightPosition, global_bsph,
				sn, camera.neardist, camera.fardist, draw_tstrips, PARA_useFaceNormal);

		}


		if (GLSL_needSilhouette)
		{
			radianceScaling.calSilhouette(GLSL_silThre, GLSL_silHardness);
		}


		if (GLSL_needScreenCurve)
			radianceScaling.calScreenCurv(tempLayerweight, curv_scale);

		if (GLSL_radianceScaling)
			radianceScaling.EvnRS(camera.field_of_view, (GLSL_curvMode == ObjectSpace || GLSL_curvMode == SurfaceReliefFeature), (GLSL_moveMode == MoveObject),
				GLSL_enhanceScale, GLSL_diffuseValue, GLSL_specularValue, GLSL_shiness, GLSL_diffColor, GLSL_specColor, GLSL_hightlightThre);

		if (GLSL_needRidgeInfo)
		{
			void(*func)(void);
			func = DrawRidgeLines;

			draw_ridge = true;
			radianceScaling.getRidgeInfo(draw_tstrips, func, PARA_useFaceNormal);
			draw_ridge = false;

		}


		if (GLSL_enhancePRT)
		{
			radianceScaling.addEdgeGlossy(PARA_size, GLSL_shiness, GLSL_enhanceScale);

		}
		
		if (prt_GLSL_addMaterial)
		{

			if (GLSL_lightMod == Light_GeoDirect || GLSL_lightMod == Light_GeoPoint)
			{
				vec3 averageNormal, averageView;
				radianceScaling.anlyNormal(lightPosition, averageNormal, averageView);

				lightPosition = averageNormal;

				switch (GLSL_materialMod)
				{
				case Blinn_Phong:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;
					break;

				case Cook_Torrance:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;


					break;

				case Oren_Nayar:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;
					break;

				case Strauss:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;

					break;

				case Ward:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;

					break;

				case Ashikhmin_Shirley:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;

					break;

				case SHW:
					lightPosition = averageNormal * averageView.dot(averageNormal)*2.f - averageView;

					break;

				default:
					break;
				}

			}



			radianceScaling.getAdjacentLightDir(lightPosition, GLSL_lightSizeX / GLSL_lightSolution, GLSL_lightSizeY / GLSL_lightSolution, GLSL_iosDirection, GLSL_lightSolution, 
				(GLSL_lightMod == Light_Point || GLSL_lightMod == Light_CamPoint || GLSL_lightMod == Light_GeoPoint), global_bsph);


			Color averageColor = transport.wavTranEnvironment.getAverageColor(lightPosition);

			radianceScaling.renderMaterial(GLSL_materialMod, (GLSL_lightMod == Light_Point || GLSL_lightMod == Light_CamPoint || GLSL_lightMod == Light_GeoPoint), lightPosition, 
				GLSL_Roughness, GLSL_Metalness, GLSL_Edginess, GLSL_Backscatter, 
				GLSL_fresnel, GLSL_alphaX, GLSL_alphaY, GLSL_iosDirection, GLSL_transparency,
				GLSL_lightSizeX / GLSL_lightSolution, GLSL_lightSizeY / GLSL_lightSolution,
				GLSL_diffuseValue, GLSL_specularValue, GLSL_diffColor, GLSL_specColor, GLSL_needColor, averageColor, prt_GLSL_CustomLightColor);
		}

		if (GLSL_ExaggeratedShading || GLSL_ExaggeratedShadingOur)
		{



				vec3 cen = global_xf*global_bsph.center;

				vec3 upRight(-global_bsph.r, global_bsph.r, 0.f);
				upRight = cen + upRight;
				upRight = inv(global_xf) * upRight;

				vec3 global = upRight - global_bsph.center;

				normalize(global);

			float tempLayerweightEx[5] = { GLSL_layerWeight[4], GLSL_layerWeight[0], GLSL_layerWeight[1], GLSL_layerWeight[2], GLSL_layerWeight[3]};


			int tempLayerIndex[5] = { -1, -1, -1, -1, -1};
			int layerNum = 0;
			float sumWeight = 0.f;
			for (int lIndex = 0; lIndex < 5; lIndex++)
			{
				if (tempLayerweightEx[lIndex] <= 0.f)
					tempLayerweightEx[lIndex] = 0.f;
				else 
				{
					sumWeight += tempLayerweightEx[lIndex];
					tempLayerIndex[layerNum++] = lIndex;
				}
			}

			sumWeight = 1.0 / sumWeight;
			for (int lIndex = 0; lIndex < layerNum; lIndex++)
			{
				tempLayerweightEx[lIndex] = tempLayerweightEx[tempLayerIndex[lIndex]] * sumWeight;
			}

			radianceScaling.exaggeratedshading(layerNum - 1, tempLayerweightEx, GLSL_ExaggeratedShadingOur, GLSL_enhanceScale, gammaAverage, draw_tstrips, GLSL_diffColor, tempLayerIndex, global, PARA_useFaceNormal);
		}

		radianceScaling.showResult(GLSL_showNormal, GLSL_showCurv, (GLSL_curvMode == ObjectSpace), GLSL_radianceScaling, GLSL_enhancePRT, GLSL_silhouette, prt_GLSL_addMaterial, GLSL_ExaggeratedShading || GLSL_ExaggeratedShadingOur);

		/*

		radianceScaling.getFrontNormal(curv_thre_pos, curv_thre_neg, false);
		draw_tstrips(themesh);
		glFinish();

		glDepthFunc(GL_GREATER);
		glClearDepth(0);
		glDisable(GL_CULL_FACE);
		radianceScaling.getBackNormal(curv_thre_pos, curv_thre_neg, false);
		draw_tstrips(themesh);
		glFinish();
		glCullFace(GL_BACK);

		glEnable(GL_CULL_FACE);
		glClearDepth(1);
		glDepthFunc(GL_LESS);

		radianceScaling.afterWriteMesh(camera.field_of_view, draw_curvature, refraction, GLSL_translucent, m_refindex, m_sigma_s, m_sigma_a, enhanceScale);
		*/

		return true;
	}
	return false;
}

bool DrawLabSpace(TriMesh::BSphere &global_bsph)
{


	if (draw_colorSpace && PRT_colSpaMod != Color::Colorspace::YCBCR)
	{
		glPushMatrix();
		float scale;
				
		switch (PRT_colSpaMod)
		{
		case Color::Colorspace::HSV:
			glTranslatef(global_bsph.center[0], global_bsph.center[1], global_bsph.center[2]);
			scale = global_bsph.r * 0.8;
			break;
		case  Color::Colorspace::CIELAB:
			glTranslatef(global_bsph.center[0], global_bsph.center[1], global_bsph.center[2]);
			scale = global_bsph.r / LABmax[0] * 0.8;
			break;
		case Color::Colorspace::RGB:
			glTranslatef(global_bsph.center[0] - global_bsph.r*0.8, global_bsph.center[1] - global_bsph.r*0.8, global_bsph.center[2] - global_bsph.r*0.8);
			scale = global_bsph.r * 1.6;
			break;

		case Color::Colorspace::XYZ:
			glTranslatef(global_bsph.center[0] - global_bsph.r*0.8, global_bsph.center[1] - global_bsph.r*0.8, global_bsph.center[2] - global_bsph.r*0.8);
			scale = global_bsph.r * 1.6;
			break;

		case Color::Colorspace::SRGB:
			glTranslatef(global_bsph.center[0] - global_bsph.r*0.8, global_bsph.center[1] - global_bsph.r*0.8, global_bsph.center[2] - global_bsph.r*0.8);
			scale = global_bsph.r * 1.6;
			break;

		default:
			break;
		}
		glScalef(scale, scale, scale);
		glPointSize(draw_point_size);
		glBegin(GL_POINTS);



		Color start;
		vec3 step = (LABmax - LABmin)*0.05f;
		if (PRT_colSpaMod == Color::Colorspace::HSV)
		{
			switch (PRT_sourceColSpaMod)
			{
			case Color::Colorspace::RGB:
			case Color::Colorspace::XYZ:
			case Color::Colorspace::SRGB:
				for (float x = 0; x < 1.f; x += 0.05f)
				for (float y = 0; y < 1.f; y += 0.05f)
				for (float z = 0; z < 1.f; z += 0.05f)
				{
					Color rgb(x, y, z);
					glColor3fv(rgb); // Used iff unlit
					rgb = rgb.convert(PRT_sourceColSpaMod, PRT_colSpaMod);


					float x1, y1;

					x1 = rgb[1] * cos(rgb[0]);
					y1 = rgb[1] * sin(rgb[0]);
					glVertex3f(x1, y1, rgb[2]);
				}
				break;
				
			case Color::Colorspace::CIELAB:
				for (float x = LABmin[0]; x < LABmax[0]; x += step[0])
				{
					start[0] = 0.f;
					for (float y = LABmin[1]; y < LABmax[1]; y += step[1])
					{
						start[1] = 0.f;
						for (float z = LABmin[2]; z < LABmax[2]; z += step[2])
						{
							Color rgb(x, y, z);
							glColor3f(start[0], start[1], 0); // Used iff unlit
							rgb = rgb.convert(rgb.CIELAB, PRT_colSpaMod);


							float x1, y1;

							x1 = rgb[1] * cos(rgb[0]);
							y1 = rgb[1] * sin(rgb[0]);

							glVertex3f(x1, y1, rgb[2]);
							start[1] += 0.05;
						}
						start[0] += 0.05;
					}
				}
				break;
			case Color::Colorspace::HSV:
				for (float v = 0; v < 1.f; v += 0.05f)
				for (float s = 0.1; s < 1.f; s += 0.1f)
				{
					float start = 0.f;
					float step = 0.1f / s;
					float step2 = step / (PI*2.f);
					for (float h = 0; h < PI*2.f; h += step)
					{
						float h2 = h + PARA_angle;
						glColor3f(s, start, 0.f); // Used iff unlit

						float x, y, z;

						x = s*cos(h2);
						y = s*sin(h2);
						z = v;

						glVertex3f(x,y,z);


						start += step2;
					}
				}
				break;

			default:
				break;
			}

		}
		else
		{
			switch (PRT_sourceColSpaMod)
			{
			case Color::Colorspace::RGB:
			case Color::Colorspace::XYZ:
			case Color::Colorspace::SRGB:
				for (float x = 0; x < 1.f; x += 0.05f)
				for (float y = 0; y < 1.f; y += 0.05f)
				for (float z = 0; z < 1.f; z += 0.05f)
				{
					Color rgb(x, y, z);
					glColor3f(x, y, z); // Used iff unlit
					

					rgb = rgb.convert(PRT_sourceColSpaMod, PRT_colSpaMod);
					glVertex3fv(rgb);
				}
				break;     
			case Color::Colorspace::CIELAB:
				for (float x = LABmin[0]; x < LABmax[0]; x += step[0])
				{
					start[0] = 0.f;
					for (float y = LABmin[1]; y < LABmax[1]; y += step[1])
					{
						start[1] = 0.f;
						for (float z = LABmin[2]; z < LABmax[2]; z += step[2])
						{
							Color rgb(x, y, z);
							glColor3f(start[0], start[1], 0); // Used iff unlit
							rgb = rgb.convert(rgb.CIELAB, PRT_colSpaMod);

							glVertex3fv(rgb);
							start[1] += 0.05;
						}
						start[0] += 0.05;
					}
				}
				break;

			case Color::Colorspace::HSV:
				for (float v = 0.1f; v < 1.f; v += 0.05f)
				for (float s = 0.1; s < 1.f; s += 0.1f)
				{
					float start = 0.f;
					float step = 0.1f / s;
					float step2 = step / (PI*2.f);
					for (float h = 0; h < PI*2.f; h += step)
					{
						float h2 = h + PARA_angle;

						glColor3f(s, start, 0.f); // Used iff unlit
						
						Color rgb(h2, s, v);
						rgb = rgb.convert(rgb.HSV, PRT_colSpaMod);
						glVertex3fv(rgb);


						start += step2;
					}
					/*
					for (float h = PI; h < PI*1.5f; h += step)
					{
						float h2 = h + PARA_angle;

						glColor3f(s, start, 0.f); // Used iff unlit

						Color rgb(h2, s, v);
						rgb = rgb.convert(rgb.HSV, spaceType);
						glVertex3fv(rgb);


						start += step2;
					}
					*/
				}
				break;

			default:
				break;
			}
		}

		glEnd();
		glPointSize(draw_point_size_space);

		glBegin(GL_POINTS);
		switch (PRT_colSpaMod)
		{
		case Color::Colorspace::RGB:
		case Color::Colorspace::XYZ:
		case Color::Colorspace::SRGB:
			for (float x = 0; x < 1.f; x += 0.05f)
			for (float y = 0; y < 1.f; y += 0.05f)
			for (float z = 0; z < 1.f; z += 0.05f)
			{
				Color position(x, y, z);
				Color rgb = position.convert(PRT_colSpaMod, position.RGB);
				glColor3fv(rgb); // Used iff unlit
				glVertex3fv(position);
			}
			break;

		case Color::Colorspace::HSV:
			for (float v = 0; v < 1.f; v += 0.05f)
			for (float s = 0.1; s < 1.f; s += 0.1f)
			{
				float step = 0.1f / s;
				for (float h = 0; h < PI*2.f; h += step)
				{
					float x, y, z;

					x = s*cos(h);
					y = s*sin(h);
					z = v;

					Color rgb(h, s, v);
					rgb = rgb.convert(rgb.HSV, rgb.RGB);
					glColor3fv(rgb); // Used iff unlit
					glVertex3f(x, y, z);

				}
			}
			break;
		case Color::Colorspace::CIELAB:
			
			for (float x = LABmin[0]; x < LABmax[0]; x += step[0])
			for (float y = LABmin[1]; y < LABmax[1]; y += step[1])
			for (float z = LABmin[2]; z < LABmax[2]; z += step[2])
			{
				Color position(x, y, z);
				Color rgb = position.convert(position.CIELAB, position.RGB);


				glColor3fv(rgb); // Used iff unlit
				glVertex3fv(position);
			}
			break;


		default:
			break;
		}







		glEnd();
		glPointSize(1.f);
		glPopMatrix();
	}




	if (false && PRT_colSpaMod != Color::Colorspace::YCBCR)
	{
		glPushMatrix();
		glTranslatef(global_bsph.center[0], global_bsph.center[1], global_bsph.center[2]);
		float scale = global_bsph.r / LABmax[0] * 0.8;
		glScalef(scale, scale, scale);

		glBegin(GL_LINES);

		glColor3f(1.0f, 0.0f, 0.f);
		glVertex3f(LABmax[0], LABmax[1], LABmax[2]);
		glVertex3f(LABmin[0], LABmax[1], LABmax[2]);
		glVertex3f(LABmax[0], LABmax[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmax[1], LABmin[2]);
		glVertex3f(LABmax[0], LABmin[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmin[1], LABmin[2]);
		glVertex3f(LABmax[0], LABmin[1], LABmax[2]);
		glVertex3f(LABmin[0], LABmin[1], LABmax[2]);

		glColor3f(0.0f, 1.0f, 0.f);
		glVertex3f(LABmax[0], LABmax[1], LABmax[2]);
		glVertex3f(LABmax[0], LABmin[1], LABmax[2]);
		glVertex3f(LABmax[0], LABmax[1], LABmin[2]);
		glVertex3f(LABmax[0], LABmin[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmax[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmin[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmax[1], LABmax[2]);
		glVertex3f(LABmin[0], LABmin[1], LABmax[2]);

		glColor3f(0.0f, 0.0f, 1.f);
		glVertex3f(LABmax[0], LABmax[1], LABmax[2]);
		glVertex3f(LABmax[0], LABmax[1], LABmin[2]);
		glVertex3f(LABmax[0], LABmin[1], LABmax[2]);
		glVertex3f(LABmax[0], LABmin[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmin[1], LABmax[2]);
		glVertex3f(LABmin[0], LABmin[1], LABmin[2]);
		glVertex3f(LABmin[0], LABmax[1], LABmax[2]);
		glVertex3f(LABmin[0], LABmax[1], LABmin[2]);
		glEnd();

		glPopMatrix();
	}



	if (prt_showLAB)
	{
			vector<Color> &col = meshes[0]->getColor();

			glPushMatrix();

			float scale = 0.f;
			
			switch (PRT_colSpaMod)
			{
			case Color::Colorspace::HSV:
				glTranslatef(global_bsph.center[0], global_bsph.center[1], global_bsph.center[2]);
				scale = global_bsph.r * 0.8;
				break;
			case Color::Colorspace::CIELAB:
				glTranslatef(global_bsph.center[0], global_bsph.center[1], global_bsph.center[2]);
				scale = global_bsph.r / LABmax[0] * 0.8;
				break;
			case Color::Colorspace::RGB:
				glTranslatef(global_bsph.center[0] - global_bsph.r*0.8, global_bsph.center[1] - global_bsph.r*0.8, global_bsph.center[2] - global_bsph.r*0.8);
				scale = global_bsph.r * 1.6;
				break;

			case Color::Colorspace::XYZ:
				glTranslatef(global_bsph.center[0] - global_bsph.r*0.8, global_bsph.center[1] - global_bsph.r*0.8, global_bsph.center[2] - global_bsph.r*0.8);
				scale = global_bsph.r * 1.6;
				break;

			case Color::Colorspace::SRGB:
				glTranslatef(global_bsph.center[0] - global_bsph.r*0.8, global_bsph.center[1] - global_bsph.r*0.8, global_bsph.center[2] - global_bsph.r*0.8);
				scale = global_bsph.r * 1.6;
				break;

			default:
				break;
			}
			glScalef(scale, scale, scale);
			glPointSize(draw_point_size);

			glBegin(GL_POINTS);


			float gamma;
			int number;
			VectorXf* lightVec = transport.wavTranEnvironment.getRawLightVec();



			switch (PRT_ColSpaSourEvnMapType)
			{
			case EM_CUR:
				lightVec = transport.wavTranEnvironment.getRawLightVec();
				break;
			case EM_DES:
				lightVec = transport.wavTranEnvironment.getRawScaledLightVec();
				break;
			case EM_AVG:
				lightVec = transport.wavTranEnvironment.getRawAverageLightVec();
				break;
			default:
				break;
			}




			switch (PRT_ColSpaSourEvnMapType)
			{
			case EM_CUR:
			case EM_DES:
			case EM_AVG:

				gamma = transport.wavTranEnvironment.getLightDirNum();
				number = transport.wavTranEnvironment.getLightDirNum();

				switch (PRT_colSpaMod)
				{
				case Color::Colorspace::RGB:
				case Color::Colorspace::XYZ:
				case Color::Colorspace::SRGB:
				case Color::Colorspace::CIELAB:

					for (int vecInd = 0; vecInd < number; vecInd++)
					{

						Color rgb(lightVec[0][vecInd] * gamma, lightVec[1][vecInd] * gamma, lightVec[2][vecInd] * gamma);
						glColor3fv(rgb); // Used iff unlit
						Color position = rgb.convert(rgb.RGB, PRT_colSpaMod);
						glVertex3fv(position);
					}
					break;

				case Color::Colorspace::HSV:
					for (int vecInd = 0; vecInd < number; vecInd++)
					{
						Color rgb(lightVec[0][vecInd] * gamma, lightVec[1][vecInd] * gamma, lightVec[2][vecInd] * gamma);
						glColor3fv(rgb); // Used iff unlit
						Color position = rgb.convert(rgb.RGB, PRT_colSpaMod);

						float x1, y1;

						x1 = position[1] * cos(position[0]);
						y1 = position[1] * sin(position[0]);
						glVertex3f(x1, y1, position[2]);

					}
					break;

				default:
					break;
				}
				break;

			case EM_NONE:
			case EM_ORI:

				switch (PRT_colSpaMod)
				{
				case Color::Colorspace::RGB:
				case Color::Colorspace::XYZ:
				case Color::Colorspace::SRGB:
				case Color::Colorspace::CIELAB:

					for (int i = 0; i < col.size(); i++)
					{
						Color rgb(col[i][0], col[i][1], col[i][2]);
						Color position = rgb.convert(rgb.RGB, PRT_colSpaMod);
						glColor3fv(col[i]); // Used iff unlit
						glVertex3fv(position);
					}
					break;

				case Color::Colorspace::HSV:
					for (int i = 0; i < col.size(); i++)
					{
						Color rgb(col[i][0], col[i][1], col[i][2]);
						Color position = rgb.convert(rgb.RGB, PRT_colSpaMod);
						glColor3fv(col[i]); // Used iff unlit

						float x1, y1;

						x1 = position[1] * cos(position[0]);
						y1 = position[1] * sin(position[0]);
						glVertex3f(x1, y1, position[2]);

					}
					break;

				default:
					break;
				}





			}

			glEnd();


			if (PRT_ColSpaSourEvnMapType == EM_AVG && PRT_colSpaMod == Color::Colorspace::RGB)
			{
				vec3 *envScale = transport.wavTranEnvironment.getEnvScale();

				glBegin(GL_LINES);
				for (int i = 0; i < 3; i++)
				{
					glColor3fv(envScale[i]);
					glVertex3fv(envScale[i]);

					glColor3f(0.f, 0.f, 0.f);
					glVertex3f(0.f, 0.f, 0.f);

					glColor3fv(envScale[i]);
					glVertex3fv(envScale[i]);

					glColor3f(1.f, 1.f, 1.f);
					glVertex3f(1.f, 1.f, 1.f);
				}

				glColor3f(0.f, 0.f, 0.f);
				glVertex3f(0.f, 0.f, 0.f);
				glColor3f(1.f, 1.f, 1.f);
				glVertex3f(1.f, 1.f, 1.f);

				glColor3fv(envScale[0]);
				glVertex3fv(envScale[0]);
				glColor3fv(envScale[1]);
				glVertex3fv(envScale[1]);

				glColor3fv(envScale[1]);
				glVertex3fv(envScale[1]);
				glColor3fv(envScale[2]);
				glVertex3fv(envScale[2]);

				glColor3fv(envScale[2]);
				glVertex3fv(envScale[2]);
				glColor3fv(envScale[0]);
				glVertex3fv(envScale[0]);

				glEnd();

			}


			glPointSize(1.f);
			glPopMatrix();

	}
	return draw_colorSpace;
}

void DrawRidgeLines()
{
	// Ridge drawing pass
	if (draw_ridge) {

		meshes[0]->need_curvatures();
		float tempTh = PARA_ridgeThreshold * meshes[0]->averageCurv1;

		glLineWidth(float(draw_line_width));


		glDisableClientState(GL_COLOR_ARRAY);

		/*
		glDisable(GL_COLOR_MATERIAL);
		GLfloat global_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat light0_diffuse[] = { 0.8f, 0.8f, 0.8f, 0.0f };
		GLfloat light1_diffuse[] = { -0.2f, -0.2f, -0.2f, 0.0f };
		GLfloat light0_specular[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
		GLfloat mat_diffuse[4] = { 1.0f, 0.0f, 0.0f, 1.0f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glColor3f(0, 0, 1); // Used iff unlit
		*/


		meshes[0]->need_RidgeLines(true);
		meshes[0]->need_RidgeLines(false);


		glBegin(GL_LINE_STRIP);


		int totalSampleIndex = 0;
		for (int lineIndex = 0; lineIndex < meshes[0]->ridgeLines.ridgeLinesPoint.size(); lineIndex++)
		{
			vector<vec3> &currentLine = meshes[0]->ridgeLines.ridgeLinesPoint[lineIndex];
			vec3 colorStart = vec3(1.0, 0.5, 0.0);
			vec3 step = vec3(-1.f / meshes[0]->ridgeLines.ridgeLinesPoint[lineIndex].size(), 0.f, 1.f / meshes[0]->ridgeLines.ridgeLinesPoint[lineIndex].size());
			for (int sIndex = 0; sIndex < meshes[0]->ridgeLines.ridgeLinesPoint[lineIndex].size(); sIndex++)
			{
				colorStart += step;
				glColor3fv(colorStart);
				
				if (abs(meshes[0]->ridgeLines.ridgeLinesPointPcurv[totalSampleIndex]) >= tempTh)
				{
					glVertex3fv(meshes[0]->ridgeLines.ridgeLinesPoint[lineIndex][sIndex]);
				}
				else
				{
					glEnd();
					glBegin(GL_LINE_STRIP);
				}

				
				totalSampleIndex++;
			}
			glEnd();
			glBegin(GL_LINE_STRIP);
		}
		totalSampleIndex = 0;
		for (int lineIndex = 0; lineIndex < meshes[0]->ridgeLines.valleyLinesPoint.size(); lineIndex++)
		{

			vec3 colorStart = vec3(1.0, 0.0, 0.0);
			vec3 step = vec3(-1.f / meshes[0]->ridgeLines.valleyLinesPoint[lineIndex].size(), 0.f, 1.f / meshes[0]->ridgeLines.valleyLinesPoint[lineIndex].size());
			for (int sIndex = 0; sIndex < meshes[0]->ridgeLines.valleyLinesPoint[lineIndex].size(); sIndex++)
			{
				colorStart += step;
				glColor3fv(colorStart);

				if (abs(meshes[0]->ridgeLines.valleyLinesPointPcurv[totalSampleIndex]) >= tempTh)
				{
					glVertex3fv(meshes[0]->ridgeLines.valleyLinesPoint[lineIndex][sIndex]);
				}
				else
				{
					glEnd();
					glBegin(GL_LINE_STRIP);
				}
				totalSampleIndex++;
			}
			glEnd();
			glBegin(GL_LINE_STRIP);
		}

		glEnd();
		
		/*
		glEnableClientState(GL_NORMAL_ARRAY);
		for (unsigned int i=0; i<themeshes[0]->ridgeLinesPoint.size(); i++)
		{
		glVertexPointer(3, GL_FLOAT,
		sizeof(themeshes[0]->ridgeLinesPoint[i][0]),
		&themeshes[0]->ridgeLinesPoint[i][0][0]);
		glNormalPointer(GL_FLOAT,
		sizeof(themeshes[0]->ridgeLinesNormal[i][0]),
		&themeshes[0]->ridgeLinesNormal[i][0][0]);
		glDrawArrays(GL_LINE_STRIP, 0, themeshes[0]->ridgeLinesPoint[i].size());

		}
		*/
	}
}

bool IsFaceNormal() { return PARA_useFaceNormal; }

bool ISdraw_white_bg(){return draw_white_bg;}

bool ISdraw_falsecolor(){return draw_falsecolor;}

bool ISdraw_shiny(){return draw_shiny;}

bool ISdraw_2side(){return draw_2side;}

bool ISdraw_points(){return draw_points;}

bool ISdraw_edges(){return draw_edges;}

bool ISdraw_ridge(){return draw_ridge;}

bool ISUsingSH(){return UsingSH;}

int getdraw_point_size(){return draw_point_size;}

int getdraw_line_width(){return draw_line_width;}

bool ISdisableLighting(){ return meshes[0]->normals.empty() || CP_showCPsingle || !draw_light || draw_colorSpace || PRT || CURV || SR || PRT_colSpaMod != Color::Colorspace::YCBCR; }

bool useColor(){return !draw_falsecolor && !draw_index && (CURV || PRT || MS_needColor || SR || CP);}

void TW_CALL getInt(void *value, void *clientData){*((int*)value) = *((int*)clientData);}

void TW_CALL getBool(void *value, void *clientData){ *((bool*)value) = *((bool*)clientData); }

void TW_CALL getFloat(void *value, void *clientData){ *((float*)value) = *((float*)clientData); }

void TW_CALL getMod(void *value, void *clientData){ *((float*)value) = *((float*)clientData); }

void TW_CALL savePlaneNormal(void *clientData){ meshes[0]->writePlaneNormal(PARA_sampleNum, PARA_ridgeThreshold, CP_curvSize, 1.f - CP_planSize); }

void TW_CALL saveTransferMattix(void *clientData){ transport.writeTransferMatrixForLRNormals(); }

void TW_CALL saveEvnScale(void *clientData) { transport.wavTranEnvironment.outputEnvScale(); }

void TW_CALL saveCluster(void *clientData){ meshes[0]->save_cluster(ms_layerIndex); }

void TW_CALL setInt(const void *value, void *clientData)
{
	*((int*)clientData) = *((int*)value);


	if (clientData == &prt_envPicIndex)
	{

		transport.changEnvironmentMap(picName[prt_envPicIndex]);
		prt_needRefreshTexture = true;
		prt_needRefresh = true;
	}

	

	if (clientData == &ms_smoothRadius)
	{
		GLSL_CurvNeedRefrash = true;
		ms_needRefesh = true;
	}

	if (clientData == &ms_layerIndex)
	{
		ms_needRefesh = true;
		ms_needRefreshParmeters = true;
	}


	if (clientData == &PARA_sampleNum)
	{
		if (CP_showCPall || CP_showCPsingle)
			CP_needRefesh = true;
		if (SR_showVertexHeight)
			SR_needRefesh = true;
	}



	if (clientData == &PARA_vertexIndex)
	{
		if (CP_showCPsingle)
			CP_needRefesh = true;
		else if (SR_showVertexHeight)
			SR_needRefesh = true;
		else if (prt_enhVertexFLE)
			prt_needDesignReflect = true;
	}

}

void TW_CALL setBool(const void *value, void *clientData)
{
	*((bool*)clientData) = *((bool*)value);


	if (clientData == &draw_colorSpace)
	{
		setBarVisible();
		return;
	}



	if (ms_showSmoothedNormal)
	{
		if (clientData == &ori || clientData == &clu || clientData == &edg || clientData == &our)
		{
			ms_needRefesh = true;
			meshes[0]->smoothedNormal[ms_layerIndex].clear();
			if (clientData == &ori)meshes[0]->smoothNormalWhithOutCluster(ms_layerIndex, ms_smoothRadius);
			if (clientData == &clu)meshes[0]->smoothNormalWhithCluster(ms_layerIndex, ms_smoothRadius);
			if (clientData == &edg)meshes[0]->smoothNormalWithEdge(ms_layerIndex, ms_smoothRadius);
			if (clientData == &our)meshes[0]->smoothNormal(ms_layerIndex, ms_smoothRadius, ms_useCluster, ms_smoothness);
		}
	}

	MS = PRT = CURV = GLSL = SR = CP = false;

	/***************   ModelSeg  *******************/
	if (clientData == &ms_showPatch || clientData == &ms_showNormalColor || clientData == &ms_showSmoothedNormalColor ||
		clientData == &ms_showNormal || clientData == &ms_showSmoothedNormal || clientData == &ms_showCurvetureColor)
	{
		if (clientData == &ms_showNormal || clientData == &ms_showSmoothedNormal)
			ms_showNormal = ms_showSmoothedNormal = false;


		if (clientData == &ms_showPatch || clientData == &ms_showNormalColor || clientData == &ms_showSmoothedNormalColor || clientData == &ms_showCurvetureColor)
			ms_showPatch = ms_showNormalColor = ms_showSmoothedNormalColor = ms_showCurvetureColor = false;

		*((bool*)clientData) = *((bool*)value);

		MS = (ms_showPatch || ms_showNormalColor || ms_showSmoothedNormalColor || ms_showNormal || ms_showSmoothedNormal || ms_showCurvetureColor);
		if (MS)
			ms_needRefesh = true;
		MS_needColor = (ms_showPatch || ms_showNormalColor || ms_showSmoothedNormalColor || ms_showCurvetureColor);
	}
	else
		ms_showNormal = ms_showSmoothedNormal = ms_showPatch = ms_showNormalColor = ms_showSmoothedNormalColor = ms_showCurvetureColor = false;



	/***************   PRT  *******************/

	if (clientData == &UsingSH || clientData == &UsingWavelet || clientData == &prt_UsingRawReflect || clientData == &prt_UsingRawRefract ||
		clientData == &prt_enhVertexFLE || clientData == &UsingGlossReflect ||
		clientData == &prt_enhAllFRA || clientData == &UsingGlossRefract ||
		clientData == &prt_showLAB || clientData == &prt_GLSL_addMaterial)
	{

		if (clientData == &UsingSH || clientData == &UsingWavelet || clientData == &prt_UsingRawReflect || clientData == &prt_UsingRawRefract)
		{
			UsingWavelet = UsingSH = false;
			if (clientData == &UsingSH || clientData == &UsingWavelet)
				prt_UsingRawReflect = prt_UsingRawRefract = false;

			*((bool*)clientData) = *((bool*)value);

			prt_needRefreshTexture = true;
		}

		PRT = (prt_UsingRawReflect || prt_UsingRawRefract || UsingSH || UsingWavelet);
		if (PRT)
			prt_needRefresh = true;
		else
			UsingWavelet = UsingSH = prt_UsingRawReflect = prt_UsingRawRefract =
			prt_enhVertexFLE = UsingGlossReflect =
			prt_enhAllFRA = UsingGlossRefract =
			prt_showLAB = false;



		if (clientData == &prt_enhVertexFLE || clientData == &UsingGlossReflect)
		{
			if (!prt_UsingRawReflect)
				*((bool*)clientData) = false;
			if (prt_enhVertexFLE)
				prt_needDesignReflect = true;
		}

		if (clientData == &prt_enhAllFRA || clientData == &UsingGlossRefract)
		{
			if (!prt_UsingRawRefract)
				*((bool*)clientData) = false;
			if (prt_enhAllFRA)
				prt_needDesignRefract = true;
		}
	}
	else
		UsingWavelet = UsingSH = prt_UsingRawReflect = prt_UsingRawRefract =
		prt_enhVertexFLE = UsingGlossReflect =
		prt_enhAllFRA = UsingGlossRefract =
		prt_showLAB = prt_GLSL_addMaterial = false;

	/***************   CURV  *******************/

	if (clientData == &draw_curvature || clientData == &draw_dcurvature)
	{
		draw_curvature = draw_dcurvature = false;
		*((bool*)clientData) = *((bool*)value);
		CURV = (draw_curvature || draw_dcurvature);
		if (CURV)
			curv_needRefresh = true;
	}
	else
		draw_curvature = draw_dcurvature = false;



	/***************   GLSL  *******************/

	
	if (clientData == &GLSL_reflection || clientData == &GLSL_refraction || clientData == &GLSL_radianceScaling || 
		clientData == &GLSL_translucent || clientData == &GLSL_showNormal || clientData == &GLSL_showCurv || clientData == &GLSL_enhancePRT ||
		clientData == &GLSL_silhouette || clientData == &GLSL_ExaggeratedShading || 
		clientData == &GLSL_ExaggeratedShadingOur || clientData == &prt_GLSL_addMaterial)
	{
		GLSL_reflection = GLSL_refraction = GLSL_radianceScaling = GLSL_enhancePRT = GLSL_translucent = GLSL_showNormal = GLSL_showCurv = GLSL_silhouette =
			GLSL_ExaggeratedShading = GLSL_ExaggeratedShadingOur = prt_GLSL_addMaterial = false;
		*((bool*)clientData) = *((bool*)value);
		GLSL = (GLSL_reflection || GLSL_refraction || GLSL_radianceScaling || GLSL_enhancePRT || GLSL_translucent || GLSL_showNormal || GLSL_showCurv || GLSL_silhouette ||
			GLSL_ExaggeratedShading || GLSL_ExaggeratedShadingOur || prt_GLSL_addMaterial);
	}
	else
		GLSL_reflection = GLSL_refraction = GLSL_radianceScaling = GLSL_enhancePRT = GLSL_translucent = GLSL_showNormal = GLSL_showCurv = GLSL_silhouette =
			GLSL_ExaggeratedShading = GLSL_ExaggeratedShadingOur = prt_GLSL_addMaterial = false;


	/***************   SR  *******************/

	if (clientData == &SR_showBaseNormal || clientData == &SR_Height || clientData == &SR_showVertexHeight || clientData == &SR_surfaceRelif)
	{
		SR_showBaseNormal = SR_Height = SR_showVertexHeight = SR_surfaceRelif = false;
		*((bool*)clientData) = *((bool*)value);
		SR = (SR_showBaseNormal || SR_Height || SR_showVertexHeight || SR_surfaceRelif);

		SR_needRefesh = true;

	}
	else
		SR_showBaseNormal = SR_Height = SR_showVertexHeight = SR_surfaceRelif = false;


	/***************   CP  *******************/

	if (clientData == &CP_showCPsingle || clientData == &CP_showCPall)
	{
		CP_showCPsingle = CP_showCPall = false;
		*((bool*)clientData) = *((bool*)value);
		CP = (CP_showCPsingle || CP_showCPall);

		if (CP_showCPsingle || CP_showCPall)
			CP_needRefesh = true;
	}
	else
		CP_showCPsingle = CP_showCPall = false;

	setBarVisible();
}

void TW_CALL setFloat(const void *value, void *clientData)
{

	*((float*)clientData) = *((float*)value);

	if (clientData == &curv_scale || clientData == &curv_thre_pos || clientData == &curv_thre_neg)
	{
		if (GLSL_showCurv || GLSL_radianceScaling)
			GLSL_CurvNeedRefrash = true;
		if (ms_showCurvetureColor)
			ms_needRefesh = true;
		if (draw_curvature)
			curv_needRefresh = true;
		if (SR_surfaceRelif)
			SR_needRefesh = true;
	}

	if (prt_enhVertexFLE && (clientData == &enhScaleReflect || clientData == &prt_sharpnessFLE || clientData == &curv_scale || clientData == &gammaAverage || clientData == &gammaScaleReflect))
		prt_needDesignReflect = true;

	if (prt_enhAllFRA && (clientData == &enhScaleRefract || clientData == &sharpnessRefract || clientData == &curv_scale || clientData == &gammaScaleRefract))
		prt_needDesignRefract = true;


	if ((prt_UsingRawRefract || prt_UsingRawReflect) && clientData == &prt_scaleRadiance)
		prt_needRefresh = true;


	if (MS && (clientData == &ms_clusterAngle || clientData == &ms_minClusterSize || clientData == &ms_maxClusterSize))
	{
		if (clientData == &ms_clusterAngle || clientData == &ms_minClusterSize || clientData == &ms_maxClusterSize)
			meshes[0]->creatLayer(ms_layerIndex, ms_clusterAngle, ms_minClusterSize, ms_maxClusterSize);
		ms_needRefesh = true;									  
		GLSL_CurvNeedRefrash = true;
	}



	if (clientData == &GLSL_layerWeight[0] || clientData == &GLSL_layerWeight[1] || clientData == &GLSL_layerWeight[2] || clientData == &GLSL_layerWeight[3] || clientData == &GLSL_layerWeight[4])
		GLSL_CurvNeedRefrash = true;


	if (clientData == &PARA_angle)
	{
		if (SR_showVertexHeight)
			SR_needRefesh = true;
	}

	if (clientData == &PARA_transX || clientData == &PARA_transY || clientData == &PARA_scale || clientData == &PARA_angle)
	{
		if (prt_showLAB)
		{
			transport.transHVSToRGB(meshes[0]->getColor(), transport.rawScaledReflectResult, PARA_transX, PARA_transY, PARA_angle, PARA_scale);
		}
	}

	if (clientData == &ms_smoothness)
	{
		GLSL_CurvNeedRefrash = true;
		ms_needRefesh = true;
	}


	if (clientData == &CP_planSize || clientData == &CP_curvSize)
		CP_needRefesh = true;
}

void TW_CALL setMod(const void *value, void *clientData)
{
	*((float*)clientData) = *((float*)value);

	if (clientData == &prt_eyeDir || clientData == &prt_material)
		prt_needRefresh = true;


	if (clientData == &PRT_EvnMapType)
	{
		prt_needRefreshTexture = true;
		if (PRT_EvnMapType == EM_DES || PRT_EvnMapType == EM_CUR)
		{
			prt_needRefresh = true;
			prt_needDesignReflect = true;
		}

	}


	setBarVisible();
}

void setBarVisible()
{
	int paraVis;

	/***************   ModelSeg  *******************/
	paraVis = (ms_showPatch || ms_showNormalColor || ms_showSmoothedNormalColor || ms_showCurvetureColor || ms_showNormal || ms_showSmoothedNormal) ? 1 : 0;
	TwSetParam(bar, "useCluster", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "smoothness", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "layer", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "Angle", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "minSize", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "maxSize", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "Sradius", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "saveCluster", "visible", TW_PARAM_INT32, 1, &paraVis);
	/***************   ModelSeg  *******************/



	/***************   Plane&Normal  *******************/
	paraVis = CP_showCPall ? 1 : 0;
	TwSetParam(bar, "savePlaNor", "visible", TW_PARAM_INT32, 1, &paraVis);
	/***************   Plane&Normal  *******************/




	/***************   PRT  *******************/
	paraVis = (prt_UsingRawReflect || prt_UsingRawRefract) ? 1 : 0;
	TwSetParam(bar, "scalRad", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "EvnType", "visible", TW_PARAM_INT32, 1, &paraVis);
	
	paraVis = (prt_enhVertexFLE) ? 1 : 0;
	TwSetParam(bar, "sharpFLE", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = prt_enhVertexFLE ? 1 : 0;
	TwSetParam(bar, "enhScaleFLE", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "scaleFLE", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_enhVertexFLE || GLSL_ExaggeratedShading || GLSL_ExaggeratedShadingOur) ? 1 : 0;
	TwSetParam(bar, "avgFLE", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = prt_UsingRawReflect ? 1 : 0;
	TwSetParam(bar, "saveTmat", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "saveEvnScale", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "showLAB", "visible", TW_PARAM_INT32, 1, &paraVis);
	/***************   PRT  *******************/



	/***************   GLSL  *******************/
	paraVis = prt_GLSL_addMaterial ? 1 : 0;
	TwSetParam(bar, "CosLightCol", "visible", TW_PARAM_INT32, 1, &paraVis);
	/***************   GLSL  *******************/



	/***************   ORI  *******************/


	paraVis = (draw_colorSpace || prt_showLAB) ? 1 : 0;
	TwSetParam(bar, "ColSpaMod", "visible", TW_PARAM_INT32, 1, &paraVis);


	paraVis = draw_colorSpace ? 1 : 0;
	TwSetParam(bar, "SurColMod", "visible", TW_PARAM_INT32, 1, &paraVis);
	/***************   ORI  *******************/




	/***************   PRT  *******************/



	/***************   PRT  *******************/




	paraVis = (CP_showCPsingle || SR_showVertexHeight || prt_enhVertexFLE) ? 1 : 0;
	TwSetParam(bar, "verInd", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (CP_showCPsingle || CP_showCPall || SR_showVertexHeight || SR_surfaceRelif) ? 1 : 0;
	TwSetParam(bar, "sampleNum", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (SR_showVertexHeight || prt_showLAB || draw_colorSpace) ? 1 : 0;
	TwSetParam(bar, "angle", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = 0;
	TwSetParam(bar, "tranX", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "tranY", "visible", TW_PARAM_INT32, 1, &paraVis);
	paraVis = (prt_showLAB) ? 1 : 0;
	TwSetParam(bar, "scale", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "SurEvnType", "visible", TW_PARAM_INT32, 1, &paraVis);
	

	paraVis = (CP_showCPsingle || CP_showCPall) ? 1 : 0;
	TwSetParam(bar, "curvSize", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "planSize", "visible", TW_PARAM_INT32, 1, &paraVis);


	paraVis = (draw_colorSpace || draw_points || prt_showLAB) ? 1 : 0;
	TwSetParam(bar, "Psize", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "PSsize", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (draw_edges || draw_ridge || GLSL_enhancePRT || GLSL_silhouette) ? 1 : 0;
	TwSetParam(bar, "Ewidth", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (draw_ridge || GLSL_enhancePRT) ? 1 : 0;
	TwSetParam(bar, "RidThr", "visible", TW_PARAM_INT32, 1, &paraVis);


	paraVis = (GLSL_enhancePRT || GLSL_silhouette) ? 1 : 0;
	TwSetParam(bar, "silThre", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "silHard", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_enhancePRT) ? 1 : 0;
	TwSetParam(bar, "size", "visible", TW_PARAM_INT32, 1, &paraVis);
	

	paraVis = (SR_surfaceRelif || draw_curvature || ms_showCurvetureColor || GLSL_showCurv || GLSL_radianceScaling) ? 1 : 0;
	TwSetParam(bar, "CurScale", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (SR_surfaceRelif || draw_curvature || ms_showCurvetureColor) ? 1 : 0;
	TwSetParam(bar, "CurThreP", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "CurThreN", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_showCurv || GLSL_radianceScaling) ? 1 : 0;
	TwSetParam(bar, "CurvMode", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_showCurv || GLSL_radianceScaling || GLSL_ExaggeratedShadingOur) ? 1 : 0;
	TwSetParam(bar, "oriWei", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "lay0Wei", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "lay1Wei", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "lay2Wei", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "lay3Wei", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_radianceScaling) ? 1 : 0;
	TwSetParam(bar, "spot", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "MoveMod", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_GLSL_addMaterial) ? 1 : 0;
	TwSetParam(bar, "MaterialGS", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "LightType", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "LightDir", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "lightXSize", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "lightYSize", "visible", TW_PARAM_INT32, 1, &paraVis);

	TwSetParam(bar, "lightFine", "visible", TW_PARAM_INT32, 1, &paraVis);
	
	paraVis = (GLSL_radianceScaling || GLSL_enhancePRT ) ? 1 : 0;
	TwSetParam(bar, "shiness", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_radianceScaling || prt_GLSL_addMaterial || GLSL_ExaggeratedShadingOur || GLSL_ExaggeratedShading) ? 1 : 0;
	TwSetParam(bar, "difCol", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_radianceScaling || (prt_GLSL_addMaterial && (GLSL_materialMod != Strauss))) ? 1 : 0;
	TwSetParam(bar, "difTerm", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (GLSL_radianceScaling || (prt_GLSL_addMaterial && (GLSL_materialMod != Strauss || GLSL_materialMod != Oren_Nayar))) ? 1 : 0;
	TwSetParam(bar, "SpeTerm", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "SpeCol", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_GLSL_addMaterial && (GLSL_materialMod == Cook_Torrance || GLSL_materialMod == ATI)) ? 1 : 0;
	TwSetParam(bar, "Fresnel", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_GLSL_addMaterial && (GLSL_materialMod == Oren_Nayar || GLSL_materialMod == Cook_Torrance || GLSL_materialMod == Strauss || GLSL_materialMod == Blinn_Phong)) ? 1 : 0;
	TwSetParam(bar, "Roughness", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_GLSL_addMaterial && GLSL_materialMod == Strauss) ? 1 : 0;
	TwSetParam(bar, "Metalness", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_GLSL_addMaterial && GLSL_materialMod == SHW) ? 1 : 0;
	TwSetParam(bar, "Edginess", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "backscatter", "visible", TW_PARAM_INT32, 1, &paraVis);
	
	



	paraVis = (prt_GLSL_addMaterial && (GLSL_materialMod == Ward || GLSL_materialMod == Ashikhmin_Shirley)) ? 1 : 0;
	TwSetParam(bar, "alphaX", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "alphaY", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "IosDir", "visible", TW_PARAM_INT32, 1, &paraVis);

	paraVis = (prt_GLSL_addMaterial && (GLSL_materialMod == Strauss)) ? 1 : 0;
	TwSetParam(bar, "transparency", "visible", TW_PARAM_INT32, 1, &paraVis);
	
	paraVis = (GLSL_radianceScaling || GLSL_enhancePRT || GLSL_ExaggeratedShading || GLSL_ExaggeratedShadingOur) ? 1 : 0;
	TwSetParam(bar, "enhScale", "visible", TW_PARAM_INT32, 1, &paraVis);




	



}

void initTwBar()
{
	TwInit(TW_OPENGL, NULL);
	TwGLUTModifiersFunc(glutGetModifiers);

	bar = TwNewBar("Parameters");
	TwDefine("Parameters size='170 750'");


	int paraVis = 0;

	

	TwAddVarCB(bar, "sigleC", TW_TYPE_BOOL8, setBool, getBool, &CP_showCPsingle, "group=CurePlane");
	TwAddVarCB(bar, "allC", TW_TYPE_BOOL8, setBool, getBool, &CP_showCPall, "group=CurePlane");
	TwAddButton(bar, "savePlaNor", savePlaneNormal, NULL, "group=CurePlane");


	TwAddVarCB(bar, "Reflect", TW_TYPE_BOOL8, setBool, getBool, &prt_UsingRawReflect, "group=PRT");



	TwAddVarCB(bar, "envIndex", TW_TYPE_INT32, setInt, getInt, &prt_envPicIndex, "min=0 max=7 step=1 group=PRT");
	
	// Current; designed; Orignal environment map; average; none.
	TwEnumVal EvnMapModEV[] = { { EM_CUR, "CUR" }, { EM_DES, "DES" }, { EM_ORI, "ORI" }, { EM_AVG, "AVG" }, { EM_NONE, "None" } };
	TwType EvnModType = TwDefineEnum("EvnMapModEV", EvnMapModEV, 5);
	TwAddVarCB(bar, "EvnType", EvnModType, setMod, getMod, &PRT_EvnMapType, " group=PRT");




	TwEnumVal eyeModEV[] = { { current, "current" }, { around, "around" }};
	TwType eyeModType = TwDefineEnum("EyeMod", eyeModEV, 2);


	TwAddVarCB(bar, "EyeDir", eyeModType, setMod, getMod, &prt_eyeDir, " group=PRT");
	TwEnumVal materialModEV[] = { { diffuse, "diffuse" }, { ward, "ward" } };
	TwType materialModType = TwDefineEnum("MaterialMod", materialModEV, 2);
	TwAddVarCB(bar, "MaterialP", materialModType, setMod, getMod, &prt_material, " group=PRT");

	TwAddVarCB(bar, "showLAB", TW_TYPE_BOOL8, setBool, getBool, &prt_showLAB, "group=PRT");
	TwAddVarRW(bar, "SurEvnType", EvnModType, &PRT_ColSpaSourEvnMapType, " group=PRT");
	TwAddVarCB(bar, "Gloss", TW_TYPE_BOOL8, setBool, getBool, &UsingGlossReflect, "group=PRT");
	TwAddVarCB(bar, "SH", TW_TYPE_BOOL8, setBool, getBool, &UsingSH, "group=PRT");
	TwAddVarCB(bar, "Wavelet", TW_TYPE_BOOL8, setBool, getBool, &UsingWavelet, "group=PRT");
	TwSetParam(bar, "Gloss", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "Wavelet", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "SH", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwAddVarCB(bar, "enhVertex", TW_TYPE_BOOL8, setBool, getBool, &prt_enhVertexFLE, "group=PRT");
	TwSetParam(bar, "enhVertex", "visible", TW_PARAM_INT32, 1, &paraVis);


	//output transfer matrix for Left and Right Normals
	TwAddButton(bar, "saveTmat", saveTransferMattix, NULL, "group=PRT");
	//output the scale of environment light. 
	TwAddButton(bar, "saveEvnScale", saveEvnScale, NULL, "group=PRT");


	TwAddVarCB(bar, "showNorl", TW_TYPE_BOOL8, setBool, getBool, &GLSL_showNormal, "group=GLSL");
	TwAddVarCB(bar, "showCurv", TW_TYPE_BOOL8, setBool, getBool, &GLSL_showCurv, "group=GLSL");
	
	TwAddVarCB(bar, "RS", TW_TYPE_BOOL8, setBool, getBool, &GLSL_radianceScaling, "group=GLSL");
	TwAddVarCB(bar, "showMaterial", TW_TYPE_BOOL8, setBool, getBool, &prt_GLSL_addMaterial, "group=GLSL");
	TwAddVarRW(bar, "CosLightCol", TW_TYPE_BOOL8, &prt_GLSL_CustomLightColor, "group=GLSL");
	TwAddVarCB(bar, "silhouette", TW_TYPE_BOOL8, setBool, getBool, &GLSL_silhouette, "group=GLSL");
	TwAddVarCB(bar, "EX", TW_TYPE_BOOL8, setBool, getBool, &GLSL_ExaggeratedShading, "group=GLSL");
	TwAddVarCB(bar, "EXOur", TW_TYPE_BOOL8, setBool, getBool, &GLSL_ExaggeratedShadingOur, "group=GLSL");
	TwAddVarRW(bar, "reflect", TW_TYPE_BOOL8, &GLSL_reflection, "group=GLSL");
	TwAddVarCB(bar, "Glass", TW_TYPE_BOOL8, setBool, getBool, &GLSL_refraction, "group=GLSL");
	TwAddVarCB(bar, "enhPRT", TW_TYPE_BOOL8, setBool, getBool, &GLSL_enhancePRT, "group=GLSL");
	TwSetParam(bar, "reflect", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "Glass", "visible", TW_PARAM_INT32, 1, &paraVis);
	TwSetParam(bar, "enhPRT", "visible", TW_PARAM_INT32, 1, &paraVis);




	TwAddVarRW(bar, "Psize", TW_TYPE_INT32, &draw_point_size, "min=1 max=20 step=1 group=Parameter");
	TwAddVarRW(bar, "PSsize", TW_TYPE_INT32, &draw_point_size_space, "min=1 max=20 step=1 group=Parameter");
	TwAddVarRW(bar, "Ewidth", TW_TYPE_INT32, &draw_line_width, "min=1 max=6 step=1 group=Parameter");
	TwAddVarCB(bar, "verInd", TW_TYPE_INT32, setInt, getInt, &PARA_vertexIndex, "min=0 max=10000000 step=1 group=Parameter");
	TwAddVarCB(bar, "sampleNum", TW_TYPE_INT32, setInt, getInt, &PARA_sampleNum, "min=3 max=101 step=2 group=Parameter");
	TwAddVarCB(bar, "angle", TW_TYPE_FLOAT, setFloat, getFloat, &PARA_angle, "min=-3.1415926 max=3.1415926 step=0.1f group=Parameter");
	TwAddVarCB(bar, "tranX", TW_TYPE_FLOAT, setFloat, getFloat, &PARA_transX, "min=-100.0 max=100.0 step=1.f group=Parameter");
	TwAddVarCB(bar, "tranY", TW_TYPE_FLOAT, setFloat, getFloat, &PARA_transY, "min=-100.0 max=100.0 step=1.f group=Parameter");
	TwAddVarCB(bar, "scale", TW_TYPE_FLOAT, setFloat, getFloat, &PARA_scale, "min=-1000 max=1000.0 step=0.1f group=Parameter");
	TwAddVarRW(bar, "size", TW_TYPE_INT32, &PARA_size, "min=0 max=1000 step=1f group=Parameter");
	TwAddVarCB(bar, "RidThr", TW_TYPE_FLOAT, setFloat, getFloat, &PARA_ridgeThreshold, "min=0.01 max=10.0 step=0.1 group=Parameter");
	int vertexNum = meshes[0]->vertices.size();
	TwSetParam(bar, "verInd", "max", TW_PARAM_INT32, 1, &vertexNum);
	TwAddVarCB(bar, "curvSize", TW_TYPE_FLOAT, setFloat, getFloat, &CP_curvSize, "min=0 max=1.0 step=0.001f group=Parameter");
	TwAddVarCB(bar, "planSize", TW_TYPE_FLOAT, setFloat, getFloat, &CP_planSize, "min=0 max=1.0 step=0.001f group=Parameter");
	TwEnumVal materialModeEV[] = { { Blinn_Phong, "Blinn_Phong" }, { Cook_Torrance, "Cook_Torrance" }, { Oren_Nayar, "Oren_Nayar" }, 
	{ Strauss, "Strauss" }, { Ward, "Ward" }, { Ashikhmin_Shirley, "Ashikhmin_Shirley" }, { SHW, "SHW" }, { ATI, "ATI" } };
	TwType materialModeType = TwDefineEnum("Material", materialModeEV, 8);
	TwAddVarCB(bar, "MaterialGS", materialModeType, setMod, getMod, &GLSL_materialMod, " group=Parameter");
	TwEnumVal lightModeEV[] = { { Light_Direct, "Direct" }, { Light_Point, "Point" }, { Light_GeoDirect, "DirectGeo" }, { Light_GeoPoint, "PointGeo" }, { Light_CamDirect, "DirectCam" }, { Light_CamPoint, "PointCam" } };
	TwType lightModeType = TwDefineEnum("Light", lightModeEV, 6);
	TwAddVarRW(bar, "LightType", lightModeType, &GLSL_lightMod, " group=Parameter");
	TwAddVarRW(bar, "LightDir", TW_TYPE_DIR3F, &GLSL_lightDirection,
		" label='Light direction' opened=true help='Change the light direction.' ");
	TwAddVarRW(bar, "Roughness", TW_TYPE_FLOAT, &GLSL_Roughness, "min=0.0 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "Metalness", TW_TYPE_FLOAT, &GLSL_Metalness, "min=0.0 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "Edginess", TW_TYPE_FLOAT, &GLSL_Edginess, "min=0.0 max=10.0 step=0.1f group=Parameter");
	TwAddVarRW(bar, "backscatter", TW_TYPE_FLOAT, &GLSL_Backscatter, "min=0.0 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "Fresnel", TW_TYPE_FLOAT, &GLSL_fresnel, "min=0.0 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "alphaX", TW_TYPE_FLOAT, &GLSL_alphaX, "min=0.0 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "alphaY", TW_TYPE_FLOAT, &GLSL_alphaY, "min=0.0 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "IosDir", TW_TYPE_DIR3F, &GLSL_iosDirection,
		" label='Light direction' opened=true help='Change the light direction.' ");
	TwAddVarRW(bar, "lightXSize", TW_TYPE_FLOAT, &GLSL_lightSizeX, "min=0.0 max=1.5 step=0.1f group=Parameter");
	TwAddVarRW(bar, "lightYSize", TW_TYPE_FLOAT, &GLSL_lightSizeY, "min=0.0 max=1.5 step=0.1f group=Parameter");
	TwAddVarRW(bar, "lightFine", TW_TYPE_FLOAT, &GLSL_lightSolution, "min=0.01 max=0.1 step=0.01f group=Parameter");
	TwAddVarRW(bar, "difTerm", TW_TYPE_FLOAT, &GLSL_diffuseValue, "min=0.0 max=100.0 step=0.1 group=Parameter");
	TwAddVarRW(bar, "SpeTerm", TW_TYPE_FLOAT, &GLSL_specularValue, "min=0.0 max=10.0 step=0.1 group=Parameter");
	TwAddVarRW(bar, "shiness", TW_TYPE_FLOAT, &GLSL_shiness, "min=0.0 max=1000.0 step=1 group=Parameter");
	TwAddVarRW(bar, "transparency", TW_TYPE_FLOAT, &GLSL_transparency, "min=0.0 max=1.0 step=0.01 group=Parameter");
	TwAddVarRW(bar, "spot", TW_TYPE_FLOAT, &GLSL_hightlightThre, "min=0.0 max=0.2 step=0.001 group=Parameter");
	TwAddVarRW(bar, "difCol", TW_TYPE_COLOR3F, &GLSL_diffColor, "group=Parameter");
	TwAddVarRW(bar, "SpeCol", TW_TYPE_COLOR3F, &GLSL_specColor, "group=Parameter");
	TwAddVarRW(bar, "enhScale", TW_TYPE_FLOAT, &GLSL_enhanceScale, "min=-100.0 max=100.0 step=0.1 group=Parameter");
	TwAddVarCB(bar, "oriWei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[4], "min=0.0 max=100.0 step=1.0 group=Parameter");
	TwAddVarCB(bar, "lay0Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[0], "min=0.0 max=100.0 step=1.0 group=Parameter");
	TwAddVarCB(bar, "lay1Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[1], "min=0.0 max=100.0 step=1.0 group=Parameter");
	TwAddVarCB(bar, "lay2Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[2], "min=0.0 max=100.0 step=1.0 group=Parameter");
	TwAddVarCB(bar, "lay3Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[3], "min=0.0 max=100.0 step=1.0 group=Parameter");
	TwEnumVal moveModeEV[] = { { MoveObject, "object" }, { MoveEye, "eye" } };
	TwType moveModeType = TwDefineEnum("MoveMod", moveModeEV, 2);
	TwAddVarRW(bar, "MoveMod", moveModeType, &GLSL_moveMode, " group=Parameter");
	TwEnumVal curvModeEV[] = { { ScreenSpace, "sreen" }, { ObjectSpace, "object" }, { SurfaceReliefFeature, "SR" } };
	TwType curvModeType = TwDefineEnum("CurvMode", curvModeEV, 3);
	TwAddVarRW(bar, "CurvMode", curvModeType, &GLSL_curvMode, " group=Parameter");
	TwAddVarRW(bar, "silThre", TW_TYPE_FLOAT, &GLSL_silThre, "min=0.0 max=0.1 step=0.001 group=Parameter");
	TwAddVarRW(bar, "silHard", TW_TYPE_FLOAT, &GLSL_silHardness, "min=0.0 max=0.1 step=0.001 group=Parameter");
	TwAddVarCB(bar, "CurScale", TW_TYPE_FLOAT, setFloat, getFloat, &curv_scale, "min=-10000.0 max=10000.0 step=1.0 group=Parameter");
	TwAddVarCB(bar, "CurThreP", TW_TYPE_FLOAT, setFloat, getFloat, &curv_thre_pos, "min=0.02 max=1.0 step=0.02 group=Parameter");
	TwAddVarCB(bar, "CurThreN", TW_TYPE_FLOAT, setFloat, getFloat, &curv_thre_neg, "min=0.02 max=1.0 step=0.02 group=Parameter");
	TwAddVarCB(bar, "scalRad", TW_TYPE_FLOAT, setFloat, getFloat, &prt_scaleRadiance, "min=-1000 max=10000.0 step=0.1 group=Parameter");
	TwAddVarCB(bar, "sharpFLE", TW_TYPE_FLOAT, setFloat, getFloat, &prt_sharpnessFLE, "min=1.0 max=100.0 step=1.0 group=Parameter");
	TwAddVarCB(bar, "enhScaleFLE", TW_TYPE_FLOAT, setFloat, getFloat, &enhScaleReflect, "min=0.1 max=10.0 step=0.1 group=Parameter");
	TwAddVarCB(bar, "avgFLE", TW_TYPE_FLOAT, setFloat, getFloat, &gammaAverage, "min=0.1 max=4.0 step=0.1 group=Parameter");
	TwAddVarCB(bar, "scaleFLE", TW_TYPE_FLOAT, setFloat, getFloat, &gammaScaleReflect, "min=0.1 max=100.0 step=0.1 group=Parameter");

	TwAddVarRW(bar, "saveCam", TW_TYPE_BOOL8, &PARA_saveMatrix, "group=Parameter");
	TwAddVarRW(bar, "loadCam", TW_TYPE_BOOL8, &PARA_loadMatrix, "group=Parameter");

/*
	TwAddVarRW(bar, "weight1", TW_TYPE_FLOAT, &GLSL_EXlayerWeight[0], "min=0.01 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "weight2", TW_TYPE_FLOAT, &GLSL_EXlayerWeight[1], "min=0.01 max=1.0 step=0.01f group=Parameter");
	TwAddVarRW(bar, "weight3", TW_TYPE_FLOAT, &GLSL_EXlayerWeight[2], "min=0.01 max=1.0 step=0.01f group=Parameter");

	TwAddVarRW(bar, "smRa1", TW_TYPE_FLOAT, &GLSL_lightSolution, "min=0 max=100 step=1f group=Parameter");
	TwAddVarRW(bar, "smRa2", TW_TYPE_FLOAT, &GLSL_lightSolution, "min=0 max=100 step=1f group=Parameter");
	TwAddVarRW(bar, "smRa3", TW_TYPE_FLOAT, &GLSL_lightSolution, "min=0 max=100 step=1f group=Parameter");*/
	/*
	TwAddVarRW(bar, "LightDir", TW_TYPE_DIR3F, &prt_eyeDirection,
	" label='Light direction' opened=true help='Change the light direction.' ");
	*/





	/*

	TwAddVarRW(bar, "reflect2", TW_TYPE_BOOL8, &CP_reflection, "group=CurePlane");
	TwAddVarCB(bar, "Glass", TW_TYPE_BOOL8, setBool, getBool, &GLSL_refraction, "group=CurePlane");
	TwAddVarCB(bar, "RS", TW_TYPE_BOOL8, setBool, getBool, &GLSL_radianceScaling, "group=CurePlane");
	TwAddVarCB(bar, "showCurv", TW_TYPE_BOOL8, setBool, getBool, &GLSL_showCurv, "group=CurePlane");
	TwAddVarRW(bar, "MoveObj/Eye", TW_TYPE_BOOL8, &GLSL_moveObject, "group=CurePlane");
	TwAddVarRW(bar, "ObjCurv/Scr", TW_TYPE_BOOL8, &GLSL_useObjectCurv, "group=CurePlane");
	TwAddVarRW(bar, "enhScale", TW_TYPE_FLOAT, &GLSL_enhanceScale, "min=0.0 max=100.0 step=1.0 group=CurePlane");
	TwAddVarCB(bar, "lay0Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[0], "min=0.0 max=100.0 step=1.0 group=CurePlane");
	TwAddVarCB(bar, "lay1Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[1], "min=0.0 max=100.0 step=1.0 group=CurePlane");
	TwAddVarCB(bar, "lay2Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[2], "min=0.0 max=100.0 step=1.0 group=CurePlane");
	TwAddVarCB(bar, "lay3Wei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[3], "min=0.0 max=100.0 step=1.0 group=CurePlane");
	TwAddVarCB(bar, "oriWei", TW_TYPE_FLOAT, setFloat, getFloat, &GLSL_layerWeight[4], "min=0.0 max=100.0 step=1.0 group=CurePlane");
	*/









	/*

	TwAddVarCB(bar, "m_refindex", TW_TYPE_FLOAT, setFloat, getFloat, &m_refindex, "min=0.1 max=100.0 step=0.1 group=Waves");
	TwAddVarRW(bar, "m_sigma_s", TW_TYPE_COLOR3F, &m_sigma_s, "group=Waves");
	TwAddVarRW(bar, "m_sigma_a", TW_TYPE_COLOR3F, &m_sigma_a, "group=Waves");
	TwAddVarRW(bar, "enhanceScale", TW_TYPE_FLOAT, &curv_scale, "min=0.001 max=1.0 step=0.001 group=Waves");
	*/





	//	TwAddVarCB(bar, "inDirect", TW_TYPE_BOOL8, setBool, getBool, &isIndirect, "group=PRT");


	//	TwAddVarCB(bar, "LightFace", TW_TYPE_INT32, setInt, getInt, &lightFace, "min=0 max=5 step=1 group=PRT");
	//	TwAddVarCB(bar, "LightX", TW_TYPE_INT32, setInt, getInt, &lightX, "min=0 max=15 step=1 group=PRT");
	//	TwAddVarCB(bar, "LightY", TW_TYPE_INT32, setInt, getInt, &lightY, "min=0 max=15 step=1 group=PRT");
	//	TwAddVarCB(bar, "verInd", TW_TYPE_INT32, setInt, getInt, &vertexInd, "min=0 max=695 step=1 group=PRT");



	TwAddVarRW(bar, "Shiny", TW_TYPE_BOOL8, &draw_shiny, "group=ORI");
	TwAddVarRW(bar, "light", TW_TYPE_BOOL8, &draw_light, "group=ORI");
	TwAddVarRW(bar, "2Sides", TW_TYPE_BOOL8, &draw_2side, "group=ORI");
	TwAddVarRW(bar, "White", TW_TYPE_BOOL8, &draw_white_bg, "group=ORI");
	TwAddVarCB(bar, "Point", TW_TYPE_BOOL8, setBool, getBool, &draw_points, "group=ORI");
	TwAddVarCB(bar, "Edge", TW_TYPE_BOOL8, setBool, getBool, &draw_edges, "group=ORI");
	TwAddVarRW(bar, "FColor", TW_TYPE_BOOL8, &draw_falsecolor, "group=ORI");
	TwEnumVal colSpaModEV[] = { { Color::Colorspace::YCBCR, "NONE" }, { Color::Colorspace::CIELAB, "LAB" }, { Color::Colorspace::RGB, "RGB" },
	{ Color::Colorspace::HSV, "HSV" }, { Color::Colorspace::SRGB, "SRGB" }, { Color::Colorspace::XYZ, "XYZ" } };
	TwType colSpaModType = TwDefineEnum("Color::Colorspace", colSpaModEV, 6);
	TwAddVarCB(bar, "colorSpace", TW_TYPE_BOOL8, setBool, getBool, &draw_colorSpace, "group=ORI");
	TwAddVarRW(bar, "SurColMod", colSpaModType, &PRT_sourceColSpaMod, " group=ORI");
	TwAddVarRW(bar, "ColSpaMod", colSpaModType, &PRT_colSpaMod, " group=ORI");
	TwAddVarCB(bar, "Ridge", TW_TYPE_BOOL8, setBool, getBool, &draw_ridge, " group=ORI");
	TwAddVarCB(bar, "Cur", TW_TYPE_BOOL8, setBool, getBool, &draw_curvature, "group=ORI");
	TwAddVarCB(bar, "DCur", TW_TYPE_BOOL8, setBool, getBool, &draw_dcurvature, "group=ORI");
	TwAddVarCB(bar, "FaceNor", TW_TYPE_BOOL8, setBool, getBool, &PARA_useFaceNormal, "group=ORI");
	




	TwAddVarCB(bar, "patch", TW_TYPE_BOOL8, setBool, getBool, &ms_showPatch, "group=ModelSeg");
	TwAddVarCB(bar, "norCol", TW_TYPE_BOOL8, setBool, getBool, &ms_showNormalColor, "group=ModelSeg");
	TwAddVarCB(bar, "SNorCol", TW_TYPE_BOOL8, setBool, getBool, &ms_showSmoothedNormalColor, "group=ModelSeg");
	TwAddVarCB(bar, "CurvCol", TW_TYPE_BOOL8, setBool, getBool, &ms_showCurvetureColor, "group=ModelSeg");
	TwAddVarCB(bar, "normal", TW_TYPE_BOOL8, setBool, getBool, &ms_showNormal, "group=ModelSeg");
	TwAddVarCB(bar, "SNor", TW_TYPE_BOOL8, setBool, getBool, &ms_showSmoothedNormal, "group=ModelSeg");
	TwAddVarRW(bar, "useCluster", TW_TYPE_BOOL8, &ms_useCluster, "group=ModelSeg");
	TwAddVarCB(bar, "smoothness", TW_TYPE_FLOAT, setFloat, getFloat, &ms_smoothness, "min=0.0 max=1.0 step=0.1 group=ModelSeg");
	//	TwAddVarCB(bar, "ori", TW_TYPE_BOOL8, setBool, getBool, &ori, "group=ModelSeg");
	//	TwAddVarCB(bar, "clu", TW_TYPE_BOOL8, setBool, getBool, &clu, "group=ModelSeg");
	//	TwAddVarCB(bar, "edg", TW_TYPE_BOOL8, setBool, getBool, &edg, "group=ModelSeg");
	//	TwAddVarCB(bar, "our", TW_TYPE_BOOL8, setBool, getBool, &our, "group=ModelSeg");
	TwAddVarCB(bar, "layer", TW_TYPE_INT32, setInt, getInt, &ms_layerIndex, "min=0 max=3 step=1 group=ModelSeg");
	TwAddVarCB(bar, "Angle", TW_TYPE_FLOAT, setFloat, getFloat, &ms_clusterAngle, "min=1.0 max=100.0 step=1.0 group=ModelSeg");
	TwAddVarCB(bar, "minSize", TW_TYPE_FLOAT, setFloat, getFloat, &ms_minClusterSize, "min=0.0 max=100.0 step=0.01f group=ModelSeg");
	TwAddVarCB(bar, "maxSize", TW_TYPE_FLOAT, setFloat, getFloat, &ms_maxClusterSize, "min=0.0 max=100.0 step=0.01f group=ModelSeg");
	TwAddVarCB(bar, "Sradius", TW_TYPE_INT32, setInt, getInt, &ms_smoothRadius, "min=1.0 max=1000.0 step=1.0 group=ModelSeg");
	TwAddButton(bar, "saveCluster", saveCluster, NULL, "group=ModelSeg");



	TwAddVarCB(bar, "base", TW_TYPE_BOOL8, setBool, getBool, &SR_showBaseNormal, "group=SRelief");
	TwAddVarRW(bar, "radiu", TW_TYPE_INT32, &SR_normalSmoothRadiance, "min=0 max=100 step=1 group=SRelief");
	TwAddVarCB(bar, "Height", TW_TYPE_BOOL8, setBool, getBool, &SR_Height, "group=SRelief");
	TwAddVarCB(bar, "sigle", TW_TYPE_BOOL8, setBool, getBool, &SR_showVertexHeight, "group=SRelief");
	TwAddVarCB(bar, "relief", TW_TYPE_BOOL8, setBool, getBool, &SR_surfaceRelif, "group=SRelief");

	/*

	TwAddVarCB(bar, "RawRefract", TW_TYPE_BOOL8, setBool, getBool, &prt_UsingRawRefract, "group=PRT_fract");
	TwAddVarCB(bar, "enh2", TW_TYPE_BOOL8, setBool, getBool, &prt_enhAllFRA, "group=PRT_fract");
	TwAddVarCB(bar, "Gloss2", TW_TYPE_BOOL8, setBool, getBool, &UsingGlossRefract, "group=PRT_fract");
	TwAddVarCB(bar, "sharp2", TW_TYPE_FLOAT, setFloat, getFloat, &sharpnessRefract, "min=1.0 max=100.0 step=2.0 group=PRT_fract");
	TwAddVarCB(bar, "enhScale2", TW_TYPE_FLOAT, setFloat, getFloat, &enhScaleRefract, "min=0.0 max=10.0 step=1.0 group=PRT_fract");
	TwAddVarCB(bar, "scale2", TW_TYPE_FLOAT, setFloat, getFloat, &gammaScaleRefract, "min=0.1 max=100.0 step=0.5 group=PRT_fract");
	*/

	// Show or hide the light variable in the Lights tweak bar



	TwWindowSize(800, 800);


	setBarVisible();



}

