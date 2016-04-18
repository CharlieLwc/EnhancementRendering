#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "stdafx.h"

#include <string>

#include "Eigen/Geometry"
USING_PART_OF_NAMESPACE_EIGEN
using namespace std;

class Environment
{
public:
	Environment(char* fileDN, int cSL = 16, int samplePP = 16)
		: cubeSideLength(cSL), samplesPerPixel(samplePP), rot(Eigen::Quaternionf::Identity()),		
		isBaseInit(false), isRadiance(true)
	{



		char* filepath = "..\\Image\\";
		char* fileAfter = ".pfm";

		strcpy_s(fileDirectionName, 80, filepath);
		strcat_s(fileDirectionName, 80, fileDN);
		strcat_s(fileDirectionName, 80, fileAfter);
		strcpy_s(fileName, 80, fileDN);
	}

	~Environment();

	//something used to load the image
	void initBase();

	void Export(string filename);
	void initTexture(unsigned int texCubeMap[]);
	void getDirectionColor(int f, int x, int y, vec3& outdir, Color& color);
	Color getAverageColor(vec3 dir);
	bool getSourceColor(int f, int x, int y, Color& sourceColor, float weight);

	float maxGamma;
	int getLightDirNum(void) { return 6 * cubeSideLength *cubeSideLength; }
	int getSamplePP(void){ return samplesPerPixel; }
	Vector3f**** getSampleDirs(void){ return sampleDirs; }

	void initTextureHole(unsigned int texCubeMap[], vec3 direction);
	char* getFileName() { return fileName; }


	EIGEN_MAKE_ALIGNED_OPERATOR_NEW



protected:
	int cubeSideLength;
	int sourceSideLength;
	int samplesPerPixel;
	Vector3f*** sampleDirs[6];
	Eigen::Quaternionf rot;
	MatrixXf faces[6][3];

	bool isRadiance, isBaseInit;

	void CubeCoordToDir(int f, float u, float v, Vector3f& outdir);
	void DirToCubeCoord(int& f, float& u, float& v, Vector3f& indir);
	char fileDirectionName[80];

	void changeEnvironmentMap(char *newFile);
	char fileName[80];

	void LoadImage();
};


#endif		// ENVIRONMENT_H

