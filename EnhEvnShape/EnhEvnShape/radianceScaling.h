#include "stdafx.h"


#include"../GL/glew.h"
#include"../GL/glut.h"

#pragma comment(lib, "glew32s.lib")


//model normal for exaggreted shading;
class exNormal
{
public:
	exNormal(){ smoothRadius = -1; }
	~exNormal(){}

	void getNormal(float tempSR, TriMesh* mesh, bool our, int ms_layerIndex)
	{
		if (smoothRadius != tempSR || our != isOur || layerIndex != ms_layerIndex)
		{
			smoothRadius = tempSR;
			isOur = our;
			layerIndex = ms_layerIndex;
			if (isOur)
				mesh[0].getSmoothNormalCluster(sommthedNormal, ms_layerIndex);
			else
				mesh[0].getSmoothNormal(sommthedNormal, smoothRadius);
		}

	}



	vector<vec3> sommthedNormal;
	float smoothRadius;
	bool isOur;
	int layerIndex;
};




class ShadingProgram
{
public:
	ShadingProgram(){inited = false;}
	~ShadingProgram();
	void inite(char* v, char* f);
	bool unInited(){return !inited;}
	GLuint getID(){return programId;}
private:
	//rutine;
	void printShaderInfoLog(GLuint obj);
	void printProgramInfoLog(GLuint obj);

	GLuint vsId, fsId;
	GLuint programId;
	bool inited;
};

//environment light;
class EvnLight
{
public:
	EvnLight(){inited = false;}
	~EvnLight(){};

	//sample Number;
	int getDrectionNum(){return directionNum;}

	
	void inite(GLuint texTempEvn, GLuint texEnvCol, GLuint texEnvDir, char* filename, int interNum);
	bool unInited(){return !inited;}




	GLint ImageWidth;//图像宽度;
	GLint ImageHeight;//图像高度;

	void getSampleColor(GLuint texEnvCol, ::std::vector<Color> &PixelDataf, ::std::vector<vec> &directions);

	//creat samples;
	void getSamplePoints(GLuint texEnvDir, ::std::vector<vec> &directions);
	void readBmp(GLuint texTempEvn, char* fileName, ::std::vector<Color> &PixelDataf);



	::std::vector<Color> PixelDataf; //像素数据;

	//the direction of samples;
	::std::vector<vec> directions;

	int directionNum;
	//use to control the number of the samples;
	int interNum;

	::std::vector<ivec3> facePoints;

private:
	bool inited;

};


//store all the GLSL files
//and shading using GLSL
class RadianceScaling
{
public:

	RadianceScaling(){}
	~RadianceScaling(){};
	//radiance scaling;
	//get model normal and so on

	void getNormal(bool needWorldNormal, bool needView, bool needWorldView, int* layerIndex, bool needCurv, bool needDepth, 
		bool needFeature, bool needBack, bool needPos, bool needCurvDir, bool needColor,
		bool needLightDirection, vec3 lightDirection, TriMesh::BSphere &global_bsph,
		GLint* sn, int zmin, int zmax, void(*draw_tstrips)(const TriMesh *, bool), bool useFaceNormal);

	void showResult(bool GLSL_showNormal, bool GLSL_showCurv, bool GLSL_useObjectCurv, bool GLSL_radianceScaling, bool GLSL_enhancePRT, bool GLSL_silhouette, bool showMaterial, bool GLSL_ExaggeratedShading);


	void getFrontNormal(float curv_thre_pos, float curv_thre_neg, bool needDepth);
	//get model normal on the back side
	void getBackNormal(float curv_thre_pos, float curv_thre_neg, bool needDepth);
	//do what you need
	void afterWriteMesh(float field_of_view, bool isGlossy, bool isRefract, bool isTranslucent, float m_refindex, float m_sigma_s[], float m_sigma_a[], float);

	void inite(TriMesh* m);

	void resizeTexture(void);

	void renderRS(float *);

	void calScreenCurv(float *, float enhanceScale);
	void getRidgeInfo(void(*draw_tstrips)(const TriMesh *, const bool), void(*DrawRidgeLines)(void), bool useFaceNormal);

	void calSilhouette(float silThre, float silHard);

	void addEdgeGlossy(int lightScale, float shiness, float enhanceScale);


	void renderMaterial(int materialType, bool isPoint, vec3 lightSource, float roughness, float metalness, float edginess, float backscatter, 
		float fresnel, float alphaX, float alphaY, vec3 IosDir, float transparency,
		int lightSizeX, int lightSizeY, 
		float diffuseValue, float specularValue, vec3 diffColor, vec3 specColor, bool needColor, Color averageColor, bool useEvnColor);

	void EvnRS(float field_of_view, bool useObjectCurv, bool movingObject, float enhanceScale, float diffuseValue, float specularValue, float shiness, vec3 diffColor, vec3 specColor, float highlightTher);

	EvnLight evnlight;

	//analysis the sum of the normal near vec3.
	void anlyNormal(vec3, vec3&, vec3&);
	void getAdjacentLightDir(vec3 lightDir, int sizeX, int sizeY, vec3 xDirection, float resolution, bool isPointLight, TriMesh::BSphere &global_bsph);

	void exaggeratedshading(int layerNum, float layerWeight[], bool isOur, float enhanceScale, float averageIlluminate, void(*draw_tstrips)(const TriMesh *, bool), vec3 color, int layerIndex[], vec3 lGlobal, bool useFaceNormal);
private:


	void initeFrameBuffer(void);

	//show the shading result or texture
	void showText(GLuint);

	//show the result environment map
	void showResultEnv(GLuint);

	//filter environment light get shading result
	void filtEvn(float field_of_view);
	//show expectation results for diffuse material
	void expectDiffuse(float field_of_view);
	//show expectation results for glossy material
	void expectGlossy(float field_of_view);
	//show expectation results for transparent material
	void transparent(float field_of_view);
	//show expectation results for translucent material
	void translucent(float field_of_view, float m_refindex, float m_sigma_s[], float m_sigma_a[]);
	//show expectation results for transparent material
	void createTSM(float field_of_view, float);

	//calculate the sum of pixels on a texture
	void PCATex(GLuint, float, float []);



	//design environmap for glossy material
	void remapGlossy();



	GLfloat* normalTexture;
	GLuint farameBuffers;


	GLuint texScreenNormal;
	GLuint* texSmoothedScreenNormal;
	GLuint texCurvAndDepth;
	GLuint texCurvDir;
	GLuint texScreenCurv;
	GLuint texRidgeLines;
	GLuint texSilhouette;
	GLuint texResult;
	GLuint texColor;
	


	GLuint texBackScreenNormal;
	GLuint texWordNormal;
	GLuint texWordView;
	GLuint texWordPos;
	
	GLuint texView;
	GLuint texBackView;
	GLuint texEnvCol;
	GLuint texEnvDir;
	GLuint texAeraLightDir;

	GLuint texMapedEvn;
	GLuint texTempEvn;
	GLuint texEvn;
	GLuint texResultEvn;
	GLuint texTSM;

	GLuint texSum;
	GLuint texTempSum;

	GLint creatShaders(char* v, char* f);
	GLint getUniLoc(GLuint program, const GLchar *name);

	ShadingProgram evnNormal;
	ShadingProgram evnScreenCurv;
	ShadingProgram evnRadianScaling;
	ShadingProgram evnRidgeLines;
	ShadingProgram evnSilhouette;
	ShadingProgram evnEdgeGlossy;
	ShadingProgram showMaterial;
	
	ShadingProgram showTexture;
	ShadingProgram refraction;
	ShadingProgram tralTSM;
	ShadingProgram subScattering;
	ShadingProgram showTexResultEnv;
	ShadingProgram analyNormal;
	ShadingProgram PCATexture;
	ShadingProgram buffProgramEvn;
	ShadingProgram evnFilter;
	ShadingProgram evnExpectation;
	ShadingProgram evnRemapGlossy;
	ShadingProgram tralucentTSM;
	ShadingProgram exaggeratedShading;

	GLuint depthbuffer;


	GLint viewPort[4];

	float enhanceScale;

	TriMesh* mesh;
	int layerNum;

	float *oldCurvature; 


	exNormal exSmoothedNormals[12];
};



