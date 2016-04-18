#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <string>
#include "Environment.h"
#include "Eigen/Sparse"
#include "BVH.h"

#include <cmath>
#include "Eigen/Core"

USING_PART_OF_NAMESPACE_EIGEN
using namespace std;



class RefractionParemeter
{
public:
	RefractionParemeter(float Teta, Color Tsigma_a, Color Tsigma_s) :eta(Teta), sigma_a(Tsigma_a), sigma_s(Tsigma_s)
	{

		if (eta >= 1)
			Fdr = -1.4399f / (eta*eta) + 0.7099f / eta + 0.6681f + 0.0636f*eta;
		else
			Fdr = -0.4399f + 0.710f / eta - 0.3319f / (eta*eta) + 0.0636f / (eta*eta*eta);

		A = (1.0f + Fdr) / (1.0f - Fdr);

		sigma_t = sigma_a + sigma_s;
		sigma_tr = Color(sqrt(3.0f*sigma_a[0] * sigma_t[0]), sqrt(3.0f*sigma_a[1] * sigma_t[1]), sqrt(3.0f*sigma_a[2] * sigma_t[2]));

		alphap = sigma_s / sigma_t;
		zpos = Color(1.0) / sigma_t;
		zneg = -zpos *(1.f + (4.f / 3.f)*A);
		Rd0 = alphap / (4.0f * 3.141592653589f);
		zpos2 = zpos*zpos;
		zneg2 = zneg*zneg;

	};
	~RefractionParemeter(){};


	float eta;
	Color sigma_a;
	Color sigma_s;

	float Fdr;
	float A;
	Color sigma_t;
	Color sigma_tr;
	Color alphap;
	Color zpos;
	Color zneg;
	Color zpos2;
	Color zneg2;
	Color Rd0;
	
};

// Norm of an RGB color
inline float ColorNorm(float red, float green, float blue)
{
	//return sqrt(red*red+green*green+blue*blue);
	return max(fabs(red), max(fabs(green), fabs(blue)));
	//return min(fabs(red), max(fabs(green),fabs(blue)));
}

/* Class for transforming in and out of the Haar
wavelet basis. */
class WaveletCompressor
{
public:
	WaveletCompressor(int dim = 0);
	void Compress(const MatrixXf& dataIn, MatrixXf& dataOut);
	void Decompress(const MatrixXf& dataIn, MatrixXf& dataOut);

	MatrixXf areaWeights;

private:
	int dimension;
	MatrixXf xform;
	MatrixXf inverse;

	// Computing the wavelet transform
	void Iterate(int iter, MatrixXf& mat);

	// Computing wavelet area weights
	void Diagonal(int imin, int jmin, int imax, int jmax, float weight);
	void Horizontal(int imin, int jmin, int imax, int jmax, float weight);
	void Vertical(int imin, int jmin, int imax, int jmax, float weight);
};



// Render black-and-white images using just the red channel
//#define RED_ONLY

class CubeDirections
{
public:
	CubeDirections(int cSL):cubeSideLength(cSL){isInit = false;}
	~CubeDirections(void);
	vec3 getDirection(int f, int x, int y);
	
private:
	int cubeSideLength;
	float **a;
	float **b;
	float **c;
	bool isInit;


};


/*  class for lighting WavTranEnvironments */
class WavTranEnvironment : public Environment
{
public:
	WavTranEnvironment(char* fileDN, int cSL = 16, float envThreshold = 0.05,
		int samplePP = 16, bool tWS = false, int nWL = 100, bool uAW = true)
		: Environment(fileDN, cSL, samplePP),
		threshold(envThreshold), comp(cSL), transWeightSelect(tWS), numWaveletLights(nWL), useAreaWeights(uAW),
		isRawInit(false), isWavInit(false), rawHasRotted(true), wavHasRotted(true), isScaledLightInit(false), isScaledLightRefractInit(false), isAverageLightInit(false), 
		rawAverageLight(0)
	{

	}

	~WavTranEnvironment(){};

	void initRaw();
	void initWav();
	void RotateRaw(Eigen::Quaternionf& q);
	void RotateWav(Eigen::Quaternionf& q);

	void initTextureAvg(unsigned int texCubeMap[]);
	void initTextureRaw(unsigned int texCubeMap[]);
	void initTextureWav(unsigned int texCubeMap[]);
	void initTextureScaledRaw(unsigned int texCubeMap[]);
	void initTextureScaledWav(unsigned int texCubeMap[]);

	float getRawLightAvg(void);
	void initAverageLightVec();
	void initScaledLights();
	void initScaledLightsRefract(int vertexNum);

	void outputEnvScale();
	//return the raw light vector
	VectorXf* getRawLightVec()
	{
		initRaw();
		return rawLightVec;
	}

	//return the vawelate light vector
	Eigen::SparseVector<float>* getWavLightVec()
	{
		initWav();
		return wavLightVec;
	}


	//return the raw light vector
	VectorXf* getRawScaledLightVec()
	{
		initScaledLights();
		return rawScaledLightVec;
	}


	//return the raw light vector
	VectorXf* getRawAverageLightVec()
	{
		initAverageLightVec();
		return rawAverageLightVec;
	}

	
	//return the raw light vector
	VectorXf* getRawScaledLightVecRefract(int vertexNum)
	{
		initScaledLightsRefract(vertexNum);
		return rawScaledRefractLightVecRefract;
	}

	//return the raw light vector
	VectorXf& getRawMidLightVec()
	{
		if (rawMidLightVec.size() == 0)
			rawMidLightVec.resize(6 * cubeSideLength*cubeSideLength);
		return rawMidLightVec;
	}
	


	VectorXf* getLabLightVecColor()
	{
		if (rawLabLightVecColor[0].size() == 0)
		{
			int csl2 = cubeSideLength*cubeSideLength;
			rawLabLightVecColor[0].resize(6 * csl2);
			rawLabLightVecColor[1].resize(6 * csl2);
			rawLabLightVecColor[2].resize(6 * csl2);
		}

		return rawLabLightVecColor;
	}

	//return the raw light vector
	VectorXf& getRawMidRefractLightVec(int vertexNum)
	{
		initScaledLightsRefract(vertexNum);
		return rawMidLightVecRefract;
	}

	//return the raw light vector
	VectorXf& getRawRefractNearDistance(int vertexNum)
	{
		initScaledLightsRefract(vertexNum);
		return lightNearDistance;
	}

	//return the vawelate light vector
	Eigen::SparseVector<float>* getWavScaledLightVec()
	{
		initScaledLights();
		return wavScaledLightVec;
	}



	bool rawNeelRelight()
	{
		bool temp = rawHasRotted;
		rawHasRotted = false;
		return temp;
	};

	bool wavNeelRelight()
	{
		bool temp = wavHasRotted;
		wavHasRotted = false;
		return temp;
	};

	vec3 *getEnvScale()
	{
		initAverageLightVec();
		return envScale;
	}


	float** transColEnergy;

	void changeWaveEnvironment(char*);


	EIGEN_MAKE_ALIGNED_OPERATOR_NEW




private:

	VectorXf rawLightVec[3];
	bool isRawInit, isWavInit, isScaledLightInit, isScaledLightRefractInit, isAverageLightInit;

	bool rawHasRotted, wavHasRotted;


	//the original lights for enhancement rendering
	//raw and wavelete
	Eigen::SparseVector<float> wavLightVec[3];



	//the designed lights for enhancement rendering
	//raw and wavelete
	VectorXf rawAverageLightVec[3];
	VectorXf rawScaledLightVec[3];
	VectorXf rawScaledLightVecRefract[3];
	VectorXf rawScaledRefractLightVecRefract[3];
	VectorXf rawMidLightVec;
	VectorXf rawLabLightVecColor[3];
	VectorXf rawMidLightVecRefract;
	VectorXf lightNearDistance;

	vec3 envScale[3];


	Eigen::SparseVector<float> wavScaledLightVec[3];

	float threshold;
	bool transWeightSelect;
	bool useAreaWeights;
	int numWaveletLights;
	WaveletCompressor comp;

	MatrixXf waveletFaces[6][3];


	float rawAverageLight;

	inline bool IsBetterWaveletLight(MatrixXf* face, int i, int j, int x, int y, float* val)
	{
		return ColorNorm(transColEnergy[i][0], transColEnergy[i][1], transColEnergy[i][2])*ColorNorm(face[0](y, x), face[1](y, x), face[2](y, x))
			> ColorNorm(transColEnergy[j][0], transColEnergy[j][1], transColEnergy[j][2])*ColorNorm(val[0], val[1], val[2]);
	}

};





/* Base class for light transport matrices */
class Transport
{
public:
	Transport(vector<TriMesh *> & m, char* FileND, int cSL = 16, float thre = 0.0075)
		:cubeSideLength(cSL), meshes(m), cubeDirections(cSL), threshold(thre),
		isRawRefractInit(false), isRawReflectInit(false), isUnOclInit(false),
		rawNeedInitRGB(true), isWavInit(false), wavNeedInitRGB(true),
		rawUsingRGB(false), rawUsingIndirect(false), wavUsingRGB(false), wavUsingIndirect(false),
		calResultColor(false), needRelight(true), oldSharpness(-1.f), oldEnhScale(-1.f), oldAverage(-1.f), oldScale(-1.f),
		needRefract(false), isVertexScaledLightInited(false), 
		RGBnum(rad), wavTranEnvironment(FileND, cSL)
	{

		colEnergy = new float*[3];
		wavTranEnvironment.transColEnergy = colEnergy;
		vertexNum = 0;
	}
	~Transport();

	//creat memory for scaled lights
	void initScaledLights();

	//calculate the environment map texture
	void initTexture(unsigned int texCubeMap[], bool isOri, bool isAvg, bool isWavelet, bool isEnhReflect, bool isEnhRefract);

	//Set the color of each point on the mesh. using Raw transport matrix.
	void setColorRaw(vector<Color>& color, bool isIndirect, bool isRGB, bool isReflect, bool isRefract, bool enhReflect, bool enhRefract, int evnLightIndex, float radianceScale, bool isRefractGlossy, MaterialType useMtaerial, EyeType useEye);


	void setColorRrefract(vector<Color>& color, bool isIndirect, bool isRefraction, bool enhCur, bool isRGB, int evnLightIndex, float radianceScale);

	//Set the color of each point on the mesh. using vawelate transport matrix.
	void setColorWav(vector<Color>& color, bool isIndirect, bool isRefraction, bool enhCur, bool isRGB, bool isGlossy);

	//design environment light based curvature. and scale the lights. For both raw and wav.
	void designLightReflectAllvertexRadiance(vector<float> &cur, bool isIndirect, bool isRGB, int evnLightIndex, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy);
	void designLightReflectRidgeRadiance(vector<float> &cur, bool isIndirect, bool isRGB, int evnLightIndex, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy);
	void designLightReflectRidgeVertexColor(vec3 eyeDirection, int vertexIndex, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy);
	void designLightReflectRidgeColor(vec3 eyeDirection, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy);

	void writeTransferMatrixForLRNormals();

	void readingExistingLight(void);

	void transHVSToRGB(vector<Color>& color, VectorXf *source, float transX, float transY, float rotAngle, float scale);

	void designLightRefract(vector<float> &cur, bool isIndirect, bool isRGB, int evnLightIndex, float sharpness, float enhScale, float curv_scale, float scale, bool isGlossy);//, vec3 eye);

	WavTranEnvironment wavTranEnvironment;
	

	float designReflectLightTransA;
	float designReflectLightTransB;


	VectorXf rawScaledReflectResult[3];

	void changEnvironmentMap(char* fileName) { wavTranEnvironment.changeWaveEnvironment(fileName);
	needRelight = true;
	}


protected:
	void getVertexLightOcl(bool*** vis, int face, BVH& meshBVH);
	void getVertexVertexInsideOcl(bool** vis, int face, BVH& meshBVH);
	void calDirectMat(MatrixXf(*matrices)[3], int face, BVH& meshBVH, bool*** vertLitOcl,  bool isUnOclReflect);
	void calIndirectMat(MatrixXf(*indirectMatrices)[3], int face, BVH& meshBVH);
	void calRefractionMat(MatrixXf(*refractionMatrices)[3], int face, BVH& meshBVH, bool isGlossy);
	void getRefractionVertexAsLightMat(MatrixXf *refractionVertexMatrices, bool isGlossy);
	void getTranMat(MatrixXf(*matrices)[3], int face, bool isDirect, bool isRefraction, bool isUnOclReflect);
	void compressTranMat(bool isGlossy);

	int cubeSideLength;
	int vertexNum;
	vector<TriMesh *> &meshes;
	CubeDirections cubeDirections;
	

	float oldSharpness, oldEnhScale, oldAverage, oldScale;
	float oldRefractSharpness, oldRefractEnhScale, oldRefractAverage, oldRefractScale;
	int oldEvnLightIndex;
	bool oldReflectGlossy, oldRefractGlossy;

	bool isWavInit, rawNeedInitRGB, wavNeedInitRGB, needRelight;
	bool isRawRefractInit, isRawReflectInit, isUnOclInit;
	bool rawUsingRGB, rawUsingIndirect, rawUsingRefLight, rawUsingRefVertex, wavUsingRGB, wavUsingIndirect, wavUsingRefraction, rawRefractGlossy;
	bool calResultColor;
	bool needRefract;
	bool isVertexScaledLightInited;

	enum LightType { rad = 1, col = 3 };

	LightType RGBnum;
	MaterialType curMaterial;
	EyeType curEye;

	void setColor(vector<Color>& color, VectorXf *result, VectorXf* lightVec, MatrixXf *transportMat, float radianceScale);

	void setAB(VectorXf *result, VectorXf* lightVec, MatrixXf *transportMat);
	void setL(VectorXf *result, VectorXf* lightVec, MatrixXf *transportMat);


	void LoadRaw(bool isRefracetVertex, bool isRefracetLight, bool isReflect, bool isIndirect, bool isRGB, bool isRefractGlossy, bool isUnOclReflect, MaterialType useMtaerial, EyeType useEye);



	void LoadWav(bool isIndirect, bool isRefraction, bool isRGB, bool isGlossy);



	float threshold;
	float **colEnergy;
	Eigen::SparseMatrix<float, Eigen::ColMajor> wavTransportMat[3];
	Eigen::SparseMatrix<float> wavResult[3];
	MatrixXf rawRefractTransportMat[3];
	MatrixXf rawReflectTransportMat[3];
	MatrixXf rawReflectTransportMatNoOcl;
	VectorXf rawRefractResult[3];
	VectorXf rawReflectResult[3];




	//the result of enhancement rendering
	VectorXf rawScaledRefractResult[3];

	MatrixXf designReflectLightMid;
	MatrixXf designRefractLightMid;


	Color* resultColor;

	Eigen::SparseMatrix<float> wavScaledResult[3];


};

#endif		// TRANSPORT_H