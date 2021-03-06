// ShowResult.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include <iostream>;
#include <fstream>;
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{

	int fileNum = 13;
	FILE *fin;
	

	char *fileName[13] = { "..\\Debug\\result_DSJ_mouse.dat", 
		"..\\Debug\\result_Ismail_Mouse.dat", 
		"..\\Debug\\result_LS_Mouse.dat", 
		"..\\Debug\\result_SD_Mouse.dat",
		"..\\Debug\\result_ZAZ_mouse.dat",
		"..\\Debug\\result1.dat",
		"..\\Debug\\result2.dat",
		"..\\Debug\\result3.dat",
		"..\\Debug\\result4.dat",
		"..\\Debug\\result5.dat",
		"..\\Debug\\result6.dat",
		"..\\Debug\\result7.dat",
		"..\\Debug\\result8.dat"
	};
	
	int resultMatrix1[13][4][3][3][3] = { 0 };
	int resultMatrix2[13][4][3][3][3] = { 0 };
	double timeMatrix[13][4][3][3][3] = {0};

	//map po line ball
	int model[36] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
		2, 2, 2,
		-1, -1, -1, -1, -1, -1,
		2, 2, 2,
		1, 1, 1, 1, 1, 1, 1, 1, 1,
		3, 3, 3, 3, 3, 3 };

	int method[3] = { 0, 2, 1 };
	int distance[36] = {2, 2, 2, 1, 1, 1, 0, 0, 0, 
		1, 1, 1,
		-1, -1, -1, -1, -1, -1,
		2, 2, 2,
		2, 2, 2, 1, 1, 1, 0, 0, 0,
		1, 1, 1, 2, 2, 2};



	int timeNum[4][3] = { 0 };

	double times[4][3] = { 0 };


	int win1[4][3][3] = {0};

	int win2[4][3][3] = { 0 };


	int timeNumIn[4][3][3] = { 0 };

	double timeIn[4][3][3] = { 0 };


	for (int f = 0; f < fileNum; f++)
	{
		int dateNum = 30;
		if (f < 5)
			dateNum = 15;


		fopen_s(&fin, fileName[f], "rb");
		if (fin)
		{

			int picIndices[30][2];
			int result[30];
			int result2[30];
			double timeLast[30];

			fread(&picIndices[0], sizeof(int), dateNum*2, fin);

			fread(&result[0], sizeof(int), dateNum, fin);
			fread(&result2[0], sizeof(int), dateNum, fin);
			fread(&timeLast[0], sizeof(double), dateNum, fin);



			for (int i = 0; i < dateNum; i++)
			{

				int curIndex = picIndices[i][0] - 133;
				int curIndex2 = picIndices[i][1] - 133;
				int mo = model[curIndex];
				int dis = distance[curIndex];

				int me = method[curIndex % 3];
				int me2 = method[curIndex2 % 3];

				resultMatrix1[f][mo][dis][me][me2] = 1 - result[i];
				resultMatrix1[f][mo][dis][me2][me] = result[i];

				resultMatrix2[f][mo][dis][me][me2] = 1 - result2[i];
				resultMatrix2[f][mo][dis][me2][me] = result2[i];


				if (result[i] == 0)
					win1[mo][dis][me]++;
				else
					win1[mo][dis][me2]++;

				if (result2[i] == 0)
					win2[mo][dis][me]++;
				else
					win2[mo][dis][me2]++;

				if (timeLast[i] == 0)
					continue;

				if (abs(timeLast[i]) > 2000 || abs(timeLast[i]) < 1)
					continue;

				times[mo][dis] += timeLast[i];
				timeNum[mo][dis]++;

				timeIn[mo][dis][me] += timeLast[i];
				timeIn[mo][dis][me2] += timeLast[i];
				timeNumIn[mo][dis][me]++;
				timeNumIn[mo][dis][me2]++;

				timeMatrix[f][mo][dis][me][me2] = timeLast[i];
				timeMatrix[f][mo][dis][me2][me] = timeLast[i];

			}


		}
	}



	ofstream result_of1("resultMatrix1.txt");
	ofstream result_of2("resultMatrix2.txt");
	ofstream time_of("timeMatrix.txt");

	for (int p = 0; p < 13; p++)
	{
		for (int f = 0; f < 4; f++)
		{
			for (int d = 0; d < 3; d++)
			{
				for (int m1 = 0; m1 < 3; m1++)
				{
					for (int m2 = 0; m2 < 3; m2++)
					{
						result_of1 << resultMatrix1[p][f][d][m1][m2] << " ";
						result_of2 << resultMatrix2[p][f][d][m1][m2] << " ";
						time_of << timeMatrix[p][f][d][m1][m2] << " ";
					}
					result_of1 << endl;
					result_of2 << endl;
					time_of << endl;

				}
				result_of1 << endl;
				result_of2 << endl;
				time_of << endl;
			}
			result_of1 << endl;
			result_of2 << endl;
			time_of << endl;
		}
		result_of1 << endl;
		result_of2 << endl;
		time_of << endl;
	}



	



	for (int m = 0; m < 4; m++)
	for (int d = 0; d < 3; d++)
	{
		times[m][d] /= timeNum[m][d];
		for (int me = 0; me < 3; me++)
			timeIn[m][d][me] /= timeNumIn[m][d][me];

	}


	return 0;
}

