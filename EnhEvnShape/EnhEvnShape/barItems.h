#ifndef INCLUDE_G_VALUE
#define INCLUDE_G_VALUE

#include "stdafx.h"


#include "AntTweakBar.h"

#include "../TriMesh2/XForm.h"



#include "radianceScaling.h"
#include"../GL/glut.h"
#include "../TriMesh2/GLCamera.h"

extern TwBar *bar;
extern vector<TriMesh *> meshes;
extern SHLighting shLight;
extern Transport transport;
extern RadianceScaling radianceScaling;

extern vector<char*> picName;



	void setBarVisible();
	void TW_CALL getFloat(void *value, void *clientData);
	void TW_CALL setFloat(const void *value, void *clientData);
	void TW_CALL setFloatAll(const void *value, void *clientData);
	void TW_CALL savePlaneNormal(void *clientData);
	void TW_CALL getInt(void *value, void *clientData);
	void TW_CALL setInt(const void *value, void *clientData);
	void TW_CALL getBool(void *value, void *clientData);
	void TW_CALL setBool(const void *value, void *clientData);


	void SLmatrix(xform &global_xf);

	void initTwBar();
	void initeCubeMapTexture();
	void drawEyedirOnScreen(TriMesh::BSphere &global_bsph, xform &global_xf);


	void meshSegPreparetion();
	void normalPreparetion();
	void curvPreparetion();
	void surfaceReliefPreparetion(TriMesh::BSphere &global_bsph, xform &global_xf);
	void curveAndPlanePreparetion(TriMesh::BSphere &global_bsph, xform &global_xf);
	void PRTPreparetion(TriMesh::BSphere &global_bsph);
	bool GLSLPreparetion(TriMesh::BSphere &global_bsph, GLCamera &camera, void(*draw_tstrips)(const TriMesh *, const bool), xform &global_xf);
	bool DrawLabSpace(TriMesh::BSphere &global_bsph);
	void DrawRidgeLines();

	bool useColor();

	bool ISdraw_white_bg();
	bool ISdraw_falsecolor();
	bool ISdraw_shiny();
	bool ISdraw_2side();
	bool ISdraw_points();
	bool ISdraw_edges();
	bool ISdraw_ridge();
	bool ISUsingSH();
	int getdraw_point_size();
	int getdraw_line_width();
	bool ISdisableLighting();
	bool IsFaceNormal();

#endif