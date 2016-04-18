#include "stdafx.h"

#include "BVH.h"
#include "Environment.h"


//Evaluate an Associated Legendre Polynomial P(l, m) at x
double P(int l, int m, double x);

//Calculate the normalisation constant for an SH function
double K(int l, int m);

//Sample a spherical harmonic basis function Y(l, m) at a point on the unit sphere
double SH(int l, int m, double theta, double phi);


//Calculate n!
int Factorial(int n);



class SAMPLE
{
public:
	//Spherical polar coords
	double theta;
	double phi;

	//Cartesian direction
	vec3 direction;
	Color color;
	//Values of each SH function at this point
	double *shValues;

	SAMPLE():shValues(NULL)
	{}
	~SAMPLE()
	{
		if(shValues)
			delete [] shValues;
	}
};

class SampleInfo
{
public:
	SampleInfo(Environment &environment, int sN, int numB);
	~SampleInfo();
	vector<SAMPLE>& getSamples(){return samples;}

private:
	bool GenerateSamples(Environment &environment);

	int cubeSideLength;
	int numBands;
	vector<SAMPLE> samples;
	int numSamples;
};



class SHVertex
{
public:

	//SH coefficients for transfer function
	double * unshadowedCoeffs;
	double * shadowedCoeffs;


	SHVertex()	: unshadowedCoeffs(NULL), shadowedCoeffs(NULL)
	{}
	~SHVertex()
	{
		if(unshadowedCoeffs)
			delete [] unshadowedCoeffs;
		unshadowedCoeffs=NULL;

		if(shadowedCoeffs)
			delete [] shadowedCoeffs;
		shadowedCoeffs=NULL;
	}
};


class  SHDirectCoeffs
{
public:
	SHDirectCoeffs(vector<TriMesh *> &objects, int numB);
	void calculateSHDirectCoeffs(SampleInfo &sampleInfo);
	~SHDirectCoeffs();


	vector<vector<SHVertex>> vertexSHCoffs;

private:

	void calculateSHDirectCoeffsSingle(vector<SAMPLE> &samples, TriMesh* object, int objectIndex);

	int numBands;
	vector<TriMesh *> &meshes;
};


class SHLighting
{
public:

	SHLighting(vector<TriMesh *> &ms, char* FileND, int cSL)
		:environment(FileND, cSL),
		meshes(ms), numBands(4), cubeSideLength(cSL),
		sampleNum(6 * cSL * cSL),
		shDiredctCoeffs(ms, 4), envLightImage(FileND)
	{
		isInit = false;
		theta = 0.0;
		phi = 0.0;
		needRote = true;
	}
	~SHLighting();

	void rote(int thetaIncrement, int phiIncrement);
	void init();

	void setColor(vector<Color>& color);

	void initTexture(unsigned int texCubeMap[], bool isOri);

private:

	void RotateSHCoefficients();
	void getSHCoefficients(SampleInfo& sampleInfo);
	void GetZRotationMatrix(int band, double * entries, double angle);
	void GetX90DegreeRotationMatrix(int band, double *entries);
	void ApplyMatrix(	int size, double * matrix, bool transpose,
		double * inVector, double * outVector);

	double light(double theta, double phi);


	bool needRote;
	int cubeSideLength;
	int numBands;
	int sampleNum;
	double theta;
	double phi;
	double (*lightCoeffs)[3];
	double (*rotatedLightCoeffs)[3];
	bool isInit;
	vector<TriMesh *> &meshes;
	SHDirectCoeffs shDiredctCoeffs;
	char* envLightImage;
	vector<Color> PixelDataf;
	Environment environment;

};
