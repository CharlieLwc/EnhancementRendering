#ifndef DFNAME_H
#define DFNAME_H

enum DFNameType
{
	Obj,
	Normal,
	FaceNormal,
	Curvature,
	Dcurvature,
	Tstrips,
	Bsphere,
	Pointareas,
	Neighbors,
	BVHDN,
	VertexVertexOutsideOcl,
	VertexVertexInsideOcl,
	VertexVertexInsideDis,
	AvrageEdgeLength,
	Halfedge,
	EdgeLength,
	VertexHeight,
	RidgeSample,
	PlaneNormals,
	PlaNorLines,
	SurfaceRelief,
	BestPosition

};

enum MaterialType { diffuse, ward };
enum EyeType { current, around };


namespace trimesh {


	class DFName {

	public:

		DFName() : fileDirection(NULL), fileName(NULL)
		{
			strcpy_s(fileAfter, 10, ".obj");
			str = new char[180];
			std::cout << "errorInDFName";
		}

		DFName(const char* d, const char *n) : fileDirection(d), fileName(n)
		{
			strcpy_s(fileDirectionTemp, 80, d);
			strcat_s(fileDirectionTemp, 80, "temp\\");

			strcpy_s(fileAfter, 10, ".obj");
			str = new char[180];
		}

		~DFName()
		{
			delete[] str;
		}

		const char *fileDirection;
		char fileDirectionTemp[80];
		const char *fileName;
		char fileAfter[10];
		char* str;

		char *getDN(DFNameType dfnType)
		{
			strcpy_s(str, 180, fileDirection);
			strcat_s(str, 180, fileName);

			switch (dfnType)
			{
			case Obj:
				break;
			case VertexHeight:
				strcat_s(str, 180, "_VertexHeight");
				break;
			case RidgeSample:
				strcat_s(str, 180, "_RidgeSample");
				break;
			case PlaneNormals:
				strcat_s(str, 180, "_PlaneNormals");
				break;
			case PlaNorLines:
				strcat_s(str, 180, "_PlaNorLines");
				break;
			case BestPosition:
				strcat_s(str, 180, "_BestPosition");
				break;
				


			default:
				break;
			}

			strcat_s(str, 180, fileAfter);
			return str;
		}

		char *getDNT(DFNameType dfnType)
		{
			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);

			switch (dfnType)
			{
			case Obj:
				break;
			case Normal:
				strcat_s(str, 180, "_normal");
				break;
			case FaceNormal:
				strcat_s(str, 180, "_faceNormal");
				break;
			case Curvature:
				strcat_s(str, 180, "_curvature");
				break;
			case Dcurvature:
				strcat_s(str, 180, "_dcurvature");
				break;
			case Tstrips:
				strcat_s(str, 180, "_tstrips");
				break;
			case Bsphere:
				strcat_s(str, 180, "_bsphere");
				break;
			case Pointareas:
				strcat_s(str, 180, "_pointareas");
				break;
			case Neighbors:
				strcat_s(str, 180, "_neighbors");
				break;
			case BVHDN:
				strcat_s(str, 180, "_BVH");
				break;
			case VertexVertexOutsideOcl:
				strcat_s(str, 180, "_VertVertOutsOcl");
				break;
			case VertexVertexInsideOcl:
				strcat_s(str, 180, "_VertVertInsOcl");
				break;
			case VertexVertexInsideDis:
				strcat_s(str, 180, "_VertVertInsDis");
				break;
			case AvrageEdgeLength:
				strcat_s(str, 180, "_AvrageEdgeLength");
				break;
			case Halfedge:
				strcat_s(str, 180, "_HalfEdge");
				break;
			case EdgeLength:
				strcat_s(str, 180, "_EdgeLength");
				break;
			case VertexHeight:
				strcat_s(str, 180, "_VertexHeight");
				break;
			case RidgeSample:
				strcat_s(str, 180, "_RidgeSample");
				break;
			case PlaneNormals:
				strcat_s(str, 180, "_PlaneNormals");
				break;
			case PlaNorLines:
				strcat_s(str, 180, "_PlaNorLines");
				break;
				
			default:
				break;
			}

			strcat_s(str, 180, fileAfter);
			return str;
		}

		char *getObjDN(){ return getDN(Obj); }
		char *getVertexHeightDN(){ return getDN(VertexHeight); }
		char *getPlaNorLinesDN(){ return getDN(PlaNorLines); }
		char *getPlaneNormalsDN(){ return getDN(PlaneNormals); }

		char *getBVHDN(){ return getDNT(BVHDN); }
		char *getCurvatureDN(){ return getDNT(Curvature); }
		char *getDcurvatureDN(){ return getDNT(Dcurvature); }
		char *getTstripsDN(){ return getDNT(Tstrips); }
		char *getBsphereDN(){ return getDNT(Bsphere); }
		char *getPointareasDN(){ return getDNT(Pointareas); }
		char *getNormalDN(){ return getDNT(Normal); }
		char *getFaceNormalDN() { return getDNT(FaceNormal); }
		char *getVertexVertexOutsideOclDN(void){ return getDNT(VertexVertexOutsideOcl); }
		char *getVertexVertexInsideOclDN(void){ return getDNT(VertexVertexInsideOcl); }
		char *getVertexVertexInsideDisDN(void){ return getDNT(VertexVertexInsideDis); }
		char *getNeighborsDN(){ return getDNT(Neighbors); }
		char *getAvrageEdgeLengthDN(){ return getDNT(AvrageEdgeLength); }
		char *getHalfEdgeDN(){ return getDNT(Halfedge); }
		char *getEdgeLengthDN(){ return getDNT(EdgeLength); }
		char *getRidgeSampleDN(){ return getDNT(RidgeSample); }

		char *getOptLightDN(int cubeSideLength, MaterialType curMaterial, EyeType curEye, char* envName){
			char str2[25];
			sprintf_s(str2, "%d_", cubeSideLength);

			strcpy_s(str, 180, fileDirection);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_OptLight_");
			
			switch (curMaterial)
			{
			case diffuse:
				strcat_s(str, 180, "Diffuse_");
				break;
			case ward:
				strcat_s(str, 180, "Ward_");
				break;
			default:
				break;
			}

			switch (curEye)
			{
			case current:
				strcat_s(str, 180, "CurEye_");
				break;
			case around:
				strcat_s(str, 180, "Around_");
				break;
			default:
				break;
			}

			strcat_s(str, 180, str2);
			strcat_s(str, 180, envName);
			strcat_s(str, 180, fileAfter);
			return str;
		
		}

		char *getClusterDN(int layerIndex)
		{
			char str2[25];
			sprintf_s(str2, "%d", layerIndex);
			strcpy_s(str, 180, fileDirection);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_cluster");
			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}

		char *getRidgeLinesDN(bool isRidge)
		{
			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			if (isRidge)
				strcat_s(str, 180, "_ridgeLines");
			else
				strcat_s(str, 180, "_valleyLines");
			strcat_s(str, 180, fileAfter);
			return str;
		}
		char *getSHDirectCoeffsDN(int numSamples)
		{
			char str2[25];
			sprintf_s(str2, "%d", numSamples);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_SHDirectCoeffs_");
			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}


		char *getVertexLightOclDN(int cubeSideLength, int face)
		{
			char str2[25];
			sprintf_s(str2, "%d_face%d", cubeSideLength, face + 1);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_VertLitOcl_");
			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}

		char *getTranMatDN(int cubeSideLength, int face, bool isIndirect, bool isRefraction, int RGBnum, bool isUnOclReflect, MaterialType curMaterial, EyeType curEye)
		{
			char str2[25];
			sprintf_s(str2, "%d_face%d", cubeSideLength, face + 1);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			if (isIndirect)
				strcat_s(str, 180, "_TranMat_Indirect_");
			else if (isRefraction)
				strcat_s(str, 180, "_TranMat_Refraction_");
			else
				strcat_s(str, 180, "_TranMat_Direct_");


			switch (curMaterial)
			{
			case diffuse:
				strcat_s(str, 180, "Diffuse_");
				break;
			case ward:
				strcat_s(str, 180, "Ward_");
				break;
			default:
				break;
			}

			switch (curEye)
			{
			case current:
				strcat_s(str, 180, "CurEye_");
				break;
			case around:
				strcat_s(str, 180, "Around_");
				break;
			default:
				break;
			}
			if (isUnOclReflect)
				strcat_s(str, 180, "UnCol_");

			if (RGBnum == 3)
				strcat_s(str, 180, "RGB_");

			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;

		}

		char *getTranMatRefractVertexAsLightDN(int RGBnum, bool isGlossy)
		{

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_RefractVertexAsLight_");

			if (isGlossy)
				strcat_s(str, 180, "Glossy_");
			else
				strcat_s(str, 180, "Diffuse_");

			if (RGBnum == 3)
				strcat_s(str, 180, "_RGB_");
			strcat_s(str, 180, fileAfter);
			return str;
		}

		char *getWaveCofDN(int cubeSideLength, bool isIndirect, int RGBnum, bool isGlossy)
		{
			char str2[25];
			sprintf_s(str2, "%d", cubeSideLength);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			if (isIndirect)
				strcat_s(str, 180, "_WaveCof_Indirect_");
			else
				strcat_s(str, 180, "_WaveCof_Direct_");

			if (isGlossy)
				strcat_s(str, 180, "Glossy_");
			else
				strcat_s(str, 180, "Diffuse_");

			if (RGBnum == 3)
				strcat_s(str, 180, "_RGB_");
			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}
		char *getWaveCofThreDN(int cubeSideLength, float thre, bool isIndirect, int RGBnum, bool isGlossy)
		{
			char str2[25];
			sprintf_s(str2, "%d_%f", cubeSideLength, thre);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			if (isIndirect)
				strcat_s(str, 180, "_WaveCof_Indirect_");
			else
				strcat_s(str, 180, "_WaveCof_Direct_");

			if (isGlossy)
				strcat_s(str, 180, "Glossy_");
			else
				strcat_s(str, 180, "Diffuse_");

			if (RGBnum == 3)
				strcat_s(str, 180, "_RGB_");
			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}
		char *getDesignedLightDN(int cubeSideLength, float sharpness, bool isGlossy)
		{
			char str2[25];
			sprintf_s(str2, "%d_%f", cubeSideLength, sharpness);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_DesignedLight_");
			if (isGlossy)
				strcat_s(str, 180, "Glossy_");
			else
				strcat_s(str, 180, "Diffuse_");

			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}
		char *getDesignedRefractLightDN(float sharpness, bool isGlossy)
		{
			char str2[25];
			sprintf_s(str2, "%f", sharpness);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_DesignedRefractLight_");

			if (isGlossy)
				strcat_s(str, 180, "Glossy_");
			else
				strcat_s(str, 180, "Diffuse_");
			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}



		char *getFeatureToLightDN(int cubeSideLength, int featureNum)
		{
			char str2[25];
			sprintf_s(str2, "%d_%d", cubeSideLength, featureNum);

			strcpy_s(str, 180, "..\\..\\matlabSoveNonlinearEquations\\");
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_FeatureToLight_");

			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;
		}

		char *getReliefDN(int sampleNum)
		{

			char str2[25];
			sprintf_s(str2, "%d", sampleNum);

			strcpy_s(str, 180, fileDirectionTemp);
			strcat_s(str, 180, fileName);
			strcat_s(str, 180, "_SurfaceRelief_");

			strcat_s(str, 180, str2);
			strcat_s(str, 180, fileAfter);
			return str;

		}


	};


















}

#endif