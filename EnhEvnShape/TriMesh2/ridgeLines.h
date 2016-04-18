//半边数据结构

#ifndef RIDGELINES_H
#define RIDGELINES_H

#include "math.h"
#include "vec.h"
#include <iostream>

#include <vector>

using namespace std;
namespace trimesh {



	struct RidgeLines
	{


		//points on the lines
		::std::vector<::std::vector<vec3>> ridgeLinesPoint;
		//the normal of the points.
		::std::vector<vec> ridgeLinesPointNormal;
		//the curv of the points.
		::std::vector<float> ridgeLinesPointPcurv;
		//the normal of the points.
		::std::vector<vec> ridgeLinesPointPDir;
		::std::vector<int> ridgeLinesPointEdge;
		::std::vector<bool> ridgeLinesPointIsVertex;



		//points on the lines
		::std::vector<::std::vector<vec3>> valleyLinesPoint;
		//the normal of the points.
		::std::vector<vec> valleyLinesPointNormal;
		//the curv of the points.
		::std::vector<float> valleyLinesPointPcurv;
		//the normal of the points.
		::std::vector<vec> valleyLinesPointPDir;

		::std::vector<int> valleyLinesPointEdge;
		::std::vector<bool> valleyLinesPointIsVertex;

		int changeToVaildIndex(int sampleIndex, float curvThreshold)
		{
			int tempSampleNum = -1;



			for (unsigned int i = 0; i < ridgeLinesPointPcurv.size(); i++)
			{

				if (abs(ridgeLinesPointPcurv[i]) >= curvThreshold)
				{
					tempSampleNum++;
					if (sampleIndex == tempSampleNum)
					{
						return i;
					}
				}
			}

			for (unsigned int i = 0; i < valleyLinesPointPcurv.size(); i++)
			{

				if (abs(valleyLinesPointPcurv[i]) >= curvThreshold)
				{
					tempSampleNum++;
					if (sampleIndex == tempSampleNum)
					{
						return i + ridgeLinesPointPcurv.size();
					}
				}
			}
			return -1;
		}

		float getCurv(int sampleIndex)
		{
			if (sampleIndex < (int)ridgeLinesPointNormal.size())
				return ridgeLinesPointPcurv[sampleIndex];

			sampleIndex -= ridgeLinesPointNormal.size();

			if (sampleIndex >(int)valleyLinesPointNormal.size())
				cout << "error at ridgelines.h" << endl;


			return valleyLinesPointPcurv[sampleIndex];

		}

		point getPointPosition(int sampleIndex)
		{
			if (sampleIndex < (int)ridgeLinesPointNormal.size())
			{
				for (int i = 0; i < (int)ridgeLinesPoint.size(); i++)
				{
					if (sampleIndex < (int)ridgeLinesPoint[i].size())
						return ridgeLinesPoint[i][sampleIndex];
					sampleIndex -= ridgeLinesPoint[i].size();
				}
				cout << "error at ridgelines.h" << endl;
			}

			sampleIndex -= ridgeLinesPointNormal.size();

			if (sampleIndex >(int)valleyLinesPointNormal.size())
				cout << "error at ridgelines.h" << endl;

			for (int i = 0; i < (int)valleyLinesPoint.size(); i++)
			{
				if (sampleIndex < (int)valleyLinesPoint[i].size())
					return valleyLinesPoint[i][sampleIndex];
				sampleIndex -= valleyLinesPoint[i].size();
			}
			cout << "error at ridgelines.h" << endl;
			return vec3(-1.f);
		}

		
		vec3 getLineDirection(int sampleIndex)
		{
			vec3 direction;
			int lineIndex;
			if (sampleIndex < (int)ridgeLinesPointNormal.size())
			{
				for (lineIndex = 0; lineIndex < (int)ridgeLinesPoint.size(); lineIndex++)
				{
					if (sampleIndex < (int)ridgeLinesPoint[lineIndex].size())
						break;
					sampleIndex -= ridgeLinesPoint[lineIndex].size();
				}
				if (lineIndex == (int)ridgeLinesPoint.size())
					cout << "error at ridgelines.h" << endl;
				

				if (sampleIndex == 0)
					direction = ridgeLinesPoint[lineIndex][sampleIndex + 1] - ridgeLinesPoint[lineIndex][sampleIndex];
				else
					direction = ridgeLinesPoint[lineIndex][sampleIndex] - ridgeLinesPoint[lineIndex][sampleIndex - 1];

				return direction;
			}

			sampleIndex -= ridgeLinesPointNormal.size();
			if (sampleIndex > (int)valleyLinesPointNormal.size())
				cout << "error at ridgelines.h" << endl;

			for (lineIndex = 0; lineIndex < (int)valleyLinesPoint.size(); lineIndex++)
			{
				if (sampleIndex < (int)valleyLinesPoint[lineIndex].size())
					break;
				sampleIndex -= valleyLinesPoint[lineIndex].size();
			}

			if (lineIndex == valleyLinesPoint.size())
				cout << "error at ridgelines.h" << endl;

			if (sampleIndex == 0)
				direction = valleyLinesPoint[lineIndex][sampleIndex + 1] - valleyLinesPoint[lineIndex][sampleIndex];
			else
				direction = valleyLinesPoint[lineIndex][sampleIndex] - valleyLinesPoint[lineIndex][sampleIndex - 1];

			normalize(direction);
			return direction;

		}

		void getPointInfo(int sampleIndex, point& position, vec3& normal, float& curv, vec3& pDir, int& edgeIndex, bool& isVertex)
		{
			position = getPointPosition(sampleIndex);

			if (sampleIndex < (int)ridgeLinesPointNormal.size())
			{
				normal = ridgeLinesPointNormal[sampleIndex];
				curv = ridgeLinesPointPcurv[sampleIndex];
				pDir = ridgeLinesPointPDir[sampleIndex];
				edgeIndex = ridgeLinesPointEdge[sampleIndex];
				isVertex = ridgeLinesPointIsVertex[sampleIndex];


				return;
			}

			sampleIndex -= ridgeLinesPointNormal.size();

			if (sampleIndex > (int)valleyLinesPointNormal.size())
				cout << "error at ridgelines.h" << endl;

			normal = valleyLinesPointNormal[sampleIndex];
			curv = valleyLinesPointPcurv[sampleIndex];
			pDir = valleyLinesPointPDir[sampleIndex];
			edgeIndex = valleyLinesPointEdge[sampleIndex];
			isVertex = valleyLinesPointIsVertex[sampleIndex];
		}




	};


};

#endif

