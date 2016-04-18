/*
Szymon Rusinkiewicz
Princeton University

TriMesh_normals.cc
Compute per-vertex normals for TriMeshes

For meshes, uses average of per-face normals, weighted according to:
  Max, N.
  "Weights for Computing Vertex Normals from Facet Normals,"
  Journal of Graphics Tools, Vol. 4, No. 2, 1999.

For raw point clouds, fits plane to k nearest neighbors.
*/

#include "TriMesh.h"
#include "KDtree.h"
#include "lineqn.h"
using namespace std;

void quickShort(int start, int end, float *valuee, int *pointer)
{
	if (start >= end)
		return;

	int middle = start + end;
	middle /= 2;

	int StandardP = pointer[middle];
	float StandardPV = valuee[StandardP];

	pointer[middle] = pointer[start];

	int TStart = start + 1;
	int TEnd = end;
	int TEmpty = start;


	while (TStart <= TEnd)
	{

		while (valuee[pointer[TEnd]] <= StandardPV && TStart <= TEnd)
			TEnd--;
		if (TStart > TEnd)
			break;
		pointer[TEmpty] = pointer[TEnd];
		TEmpty = TEnd;
		TEnd--;

		while (valuee[pointer[TStart]] >= StandardPV && TStart <= TEnd)
			TStart++;
		if (TStart > TEnd)
			break;
		pointer[TEmpty] = pointer[TStart];
		TEmpty = TStart;
		TStart++;
	}

	pointer[TEmpty] = StandardP;

	middle = TEmpty - 1;

	for (int i = middle; i >= start; i--)
	{
		if (valuee[pointer[i]] == StandardPV)
		{
			StandardP = pointer[i];
			pointer[i] = pointer[middle];
			pointer[middle] = StandardP;
			middle--;
		}
	}
	quickShort(start, middle, valuee, pointer);

	middle = TEmpty + 1;
	for (int i = middle; i <= end; i++)
	{
		if (valuee[pointer[i]] == StandardPV)
		{
			StandardP = pointer[i];
			pointer[i] = pointer[middle];
			pointer[middle] = StandardP;
			middle++;
		}
	}
	quickShort(middle, end, valuee, pointer);
	return;

}

void bubbleSort(int clusterNum, float* clusterScale, int *clusterPointer)
{
	int tempPointer;
	for (int i = 1; i<clusterNum; i++)
	{
		for (int j = i; j>0; j--)
		{
			if (clusterScale[clusterPointer[j]] > clusterScale[clusterPointer[j - 1]])
			{
				tempPointer = clusterPointer[j];
				clusterPointer[j] = clusterPointer[j - 1];
				clusterPointer[j - 1] = tempPointer;
			}
		}
	}
}

void bubbleSortVec(int clusterNum, vector<float> &clusterScale, vector<int> &clusterPointer)
{
	int tempPointer;
	for (int i = 1; i<clusterNum; i++)
	{
		for (int j = i; j>0; j--)
		{
			if (clusterScale[clusterPointer[j]] > clusterScale[clusterPointer[j - 1]])
			{
				tempPointer = clusterPointer[j];
				clusterPointer[j] = clusterPointer[j - 1];
				clusterPointer[j - 1] = tempPointer;
			}
		}
	}
}




namespace trimesh {



	// Compute per-vertex normals
	void TriMesh::need_cluster(int layerIndex)
	{
		clamp(layerIndex, 0, 3);
		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];
		if (!currentVertexCluster.empty())
			return;

		smoothRadius[layerIndex] = 0;
		ifstream cluster_if(dFName.getClusterDN(layerIndex));

		if (cluster_if.is_open())
		{
			currentVertexCluster.resize(vertices.size());
			cout << "Reading CLuster.....";
			cluster_if >> clusterNum[layerIndex];
			cluster_if >> lowClusterScale[layerIndex];
			cluster_if >> upClusterScale[layerIndex];
			averageClusterScale[layerIndex].resize(1);
			cluster_if >> averageClusterScale[layerIndex][0];
			cluster_if >> clusterAngle[layerIndex];

			for (unsigned int i = 0; i<vertices.size(); i++)
				cluster_if >> currentVertexCluster[i];

			clusterScale[layerIndex].resize(clusterNum[layerIndex]);
			for (unsigned int i = 0; i < clusterNum[layerIndex]; i++)
				cluster_if >> clusterScale[layerIndex][i];


			clusterNormal[layerIndex].resize(clusterNum[layerIndex]);
			for (unsigned int i = 0; i < clusterNum[layerIndex]; i++)
			{
				cluster_if >> clusterNormal[layerIndex][i][0];
				cluster_if >> clusterNormal[layerIndex][i][1];
				cluster_if >> clusterNormal[layerIndex][i][2];
			}


			cout << "Done" << endl;
			cluster_if.close();

			return;
		}

		clusterAngle[layerIndex] = 0.f;
		creatLayer(layerIndex, -1.f, -1.f, -1.f);

	}

	void TriMesh::save_cluster(int layerIndex)
	{
		clamp(layerIndex, 0, 3);
		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];
		if (currentVertexCluster.empty())
			return;


		ofstream cluster_of(dFName.getClusterDN(layerIndex));

		if (cluster_of.is_open())
		{
			cout << "Writing CLuster.....";
			cluster_of << clusterNum[layerIndex] << endl;
			cluster_of << lowClusterScale[layerIndex]<<" ";
			cluster_of << upClusterScale[layerIndex] << " ";
			cluster_of << averageClusterScale[layerIndex][0] << " ";
			cluster_of << clusterAngle[layerIndex] << endl;

			for (unsigned int i = 0; i<vertices.size(); i++)
				cluster_of << currentVertexCluster[i] << " ";

			cluster_of << endl;


			for (unsigned int i = 0; i < clusterNum[layerIndex]; i++)
				cluster_of << clusterScale[layerIndex][i] << " ";

			cluster_of << endl;


			for (unsigned int i = 0; i < clusterNum[layerIndex]; i++)
			{
				cluster_of << clusterNormal[layerIndex][i][0] << " ";
				cluster_of << clusterNormal[layerIndex][i][1] << " ";
				cluster_of << clusterNormal[layerIndex][i][2] << " ";
			}





			cout << "Done" << endl;
			cluster_of.close();

			return;
		}


	}

	void TriMesh::creatLayer(int layerIndex, float angle, float minCLusterSize, float maxClusterSize)
	{

		clamp(layerIndex, 0, 3);
		if (angle <= 0)
		{
			if (layerIndex == 0)
				angle = 10;
			else
				angle = 15.f;
		}

		minCLusterSize = minCLusterSize < 0.f ? 0.f : minCLusterSize;
		maxClusterSize = maxClusterSize < 0.f ? 0.f : maxClusterSize;

		if (layerIndex == 0 && clusterAngle[layerIndex] == angle)
			return;

		if (layerIndex != 0 && clusterAngle[layerIndex] == angle && lowClusterScale[layerIndex - 1] == minCLusterSize && upClusterScale[layerIndex-1] == maxClusterSize)
			return;


		clusterAngle[layerIndex] = angle;

		if (layerIndex != 0 && minCLusterSize > 0.f)
			lowClusterScale[layerIndex - 1] = minCLusterSize;

		if (layerIndex != 0 && maxClusterSize > 0.f)
			upClusterScale[layerIndex - 1] = maxClusterSize;

		if (layerIndex == 0)
		{
			need_neighbors();
			if (vertexCluster[layerIndex].empty())
				clusterVertex(angle);
			reclusterVertex(angle);
		}
		else
			clusterCluster(layerIndex, angle);



		lowClusterScale[layerIndex] = -1;
		upClusterScale[layerIndex] = -1;
		clusterSize[layerIndex].clear();
		clusterTopology[layerIndex].clear();
		clusterHalfArc[layerIndex].clear();
		clusterNormal[layerIndex].clear();
		smoothedNormal[layerIndex].clear();
		averageClusterScale[layerIndex].clear();

		need_upAndLowScale(layerIndex, 0, clusterNum[layerIndex] - 1);

	}

	void TriMesh::clusterVertex(float angle)
	{
		float angleTher = cos(M_PIf*angle / 180.f);

		cout << "clustering vertex. Threhould: " << angle << " degree. " << endl;
		int vertexNum = vertices.size();

		int *tempVertexClusterIndex = new int[vertexNum];

		/////////////////////////////////////////////
		//分出平缓小块。;
		////////////////////////////////////////////////


		//分块缓存;

		vertexCluster[0].resize(vertices.size());

		for (int i = 0; i < vertexNum; i++)
			tempVertexClusterIndex[i] = -1;

		int tempClusterNum = 0;
		for (int vIndex = 0; vIndex < vertexNum; vIndex++)
		{
			//建立新块;
			if (tempVertexClusterIndex[vIndex] >= 0)
				continue;

			tempVertexClusterIndex[vIndex] = tempClusterNum;
			vec curNormal = normals[vIndex];

			vector<int> unspreadVertexes;
			unspreadVertexes.push_back(vIndex);
			//行走;
			while (!unspreadVertexes.empty())
			{
				//处理末尾节点;
				int curVertex = unspreadVertexes.back();
				unspreadVertexes.pop_back();





				for (int vcIndex = 0; vcIndex < neighbors[curVertex].size(); vcIndex++)
				{
					int nextVertexIndex = neighbors[curVertex][vcIndex];

					if (tempVertexClusterIndex[nextVertexIndex] >= 0 ||
						tempVertexClusterIndex[nextVertexIndex] == -2 - tempClusterNum)
						continue;
					float cosAngle = curNormal.dot(normals[nextVertexIndex]);

					if (cosAngle>angleTher)
					{
						tempVertexClusterIndex[nextVertexIndex] = tempClusterNum;
						unspreadVertexes.push_back(nextVertexIndex);
					}
					else
						tempVertexClusterIndex[nextVertexIndex] = -2 - tempClusterNum;
					//下一个相邻点;
				}
			}
			tempClusterNum++;
		}

		//存储结果;
		cout << "There are " << tempClusterNum << " seed clusters. " << endl;


		clusterNum[0] = tempClusterNum;

		clusterScale[0].clear();
		clusterScale[0].resize(tempClusterNum);

		vector<float> &currentClusterScale = clusterScale[0];

		vector<int> &vertexClusterIndex = vertexCluster[0];




		float *tempClusterScale = new float[tempClusterNum];
		calScale(tempClusterScale, tempVertexClusterIndex, tempClusterNum);

		int *clusterPointer = new int[tempClusterNum];
		for (int cIndex = 0; cIndex < tempClusterNum; cIndex++)
			clusterPointer[cIndex] = cIndex;

		quickShort(0, tempClusterNum - 1, tempClusterScale, clusterPointer);

		int *clusterSequence = new int[tempClusterNum];
		for (int i = 0; i < tempClusterNum; i++)
		{
			clusterSequence[clusterPointer[i]] = i;
			currentClusterScale[i] = tempClusterScale[clusterPointer[i]];
		}
		delete[] clusterPointer;
		delete[] tempClusterScale;

		for (int i = 0; i < vertexNum; i++)
			vertexClusterIndex[i] = clusterSequence[tempVertexClusterIndex[i]];

		delete[] tempVertexClusterIndex;
		delete[] clusterSequence;





	}

	void TriMesh::reclusterVertex(float angle)
	{


		clusterAngle[0] = angle;

		float angleTher = cos(M_PIf*angle / 180.f);

		cout << "Reclustering vertex. Threhould: " << angle << " degree. " << endl;

		need_cluster(0);
		need_clusterNormal(0);

		vector<vec> &seedClusterNormal = clusterNormal[0];

		int seedClusterNum = clusterNum[0];


		int *startVertex = new int[seedClusterNum];
		float *biggestCos = new float[seedClusterNum];

		for (int cIndex = 0; cIndex<seedClusterNum; cIndex++)
			biggestCos[cIndex] = 0.f;

		for (int cIndex = 0; cIndex<seedClusterNum; cIndex++)
			startVertex[cIndex] = -1;
		int vertexNum = vertices.size();
		vector<int> &vertexClusterIndex = vertexCluster[0];
		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
		{
			int cIndex = vertexClusterIndex[vIndex];
			float cos = seedClusterNormal[cIndex].dot(normals[vIndex]);
			if (cos > biggestCos[cIndex])
			{
				biggestCos[cIndex] = cos;
				startVertex[cIndex] = vIndex;
			}
		}
		delete[] biggestCos;


		/////////////////////////////////////////////
		//二次聚类，得到平缓的标准区域以及小区域。;
		////////////////////////////////////////////////


		int *tempVertexClusterIndex = new int[vertexNum];
		for (int i = 0; i<vertexNum; i++)
			tempVertexClusterIndex[i] = -1;

		int *seedClusterVertexNum = new int[seedClusterNum];
		for (int cIndex = 0; cIndex<seedClusterNum; cIndex++)
			seedClusterVertexNum[cIndex] = 0;
		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			seedClusterVertexNum[vertexClusterIndex[vIndex]]++;

		int tempClusterNum = 0;
		for (int cIndex = 0; cIndex<seedClusterNum; cIndex++)
		{
			//寻找块;
			if (seedClusterVertexNum[cIndex] <= 0)
				continue;

			//寻找起点;
			int curIndex = startVertex[cIndex];
			vec curNormal = seedClusterNormal[cIndex];

			//找，草;
			if (tempVertexClusterIndex[curIndex] >= 0)
			{
				curIndex = -1;
				float biggestCos = 0.f;

				for (int i = 0; i<vertexNum; i++)
				if (vertexClusterIndex[i] == cIndex && tempVertexClusterIndex[curIndex] < 0)
				{
					float cos = curNormal.dot(normals[i]);
					if (cos > biggestCos)
					{
						biggestCos = cos;
						curIndex = i;
					}
				}
			}
			if (curIndex == -1)
				continue;

			//建立新块;
			tempVertexClusterIndex[curIndex] = tempClusterNum;
			curNormal = normals[curIndex];
			//压栈;
			vector<int> unspreadVertexes;
			unspreadVertexes.push_back(curIndex);
			while (!unspreadVertexes.empty())
			{
				//处理末尾节点;
				int curVertex = unspreadVertexes.back();
				unspreadVertexes.pop_back();



				for (int vcIndex = 0; vcIndex < neighbors[curVertex].size(); vcIndex++)
				{
					int nextVertexIndex = neighbors[curVertex][vcIndex];

					if (tempVertexClusterIndex[nextVertexIndex] >= 0 ||
						tempVertexClusterIndex[nextVertexIndex] == -2 - tempClusterNum)
						continue;
					float cosAngle = curNormal.dot(normals[nextVertexIndex]);

					if (cosAngle>angleTher)
					{
						tempVertexClusterIndex[nextVertexIndex] = tempClusterNum;
						unspreadVertexes.push_back(nextVertexIndex);

						if (vertexClusterIndex[nextVertexIndex] != cIndex)
							seedClusterVertexNum[vertexClusterIndex[nextVertexIndex]]--;
					}
					else
						tempVertexClusterIndex[nextVertexIndex] = -2 - tempClusterNum;
				}
			}
			tempClusterNum++;
		}


		delete seedClusterVertexNum;
		delete[] startVertex;


		//对剩下的点聚类一下;
		for (int curIndex = 0; curIndex<vertexNum; curIndex++)
		{
			if (tempVertexClusterIndex[curIndex] >= 0)
				continue;

			//建立新块;

			tempVertexClusterIndex[curIndex] = tempClusterNum;
			vec curNormal = normals[curIndex];

			//压栈;
			vector<int> unspreadVertexes;
			unspreadVertexes.push_back(curIndex);
			while (!unspreadVertexes.empty())
			{
				//处理末尾节点;
				int curVertex = unspreadVertexes.back();
				unspreadVertexes.pop_back();



				for (int vcIndex = 0; vcIndex < neighbors[curVertex].size(); vcIndex++)
				{
					int nextVertexIndex = neighbors[curVertex][vcIndex];
					if (tempVertexClusterIndex[nextVertexIndex] >= 0 ||
						tempVertexClusterIndex[nextVertexIndex] == -2 - tempClusterNum)
						continue;
					float cosAngle = curNormal.dot(normals[nextVertexIndex]);

					if (cosAngle>angleTher)
					{
						tempVertexClusterIndex[nextVertexIndex] = tempClusterNum;
						unspreadVertexes.push_back(nextVertexIndex);
					}
					else
						tempVertexClusterIndex[nextVertexIndex] = -2 - tempClusterNum;
				}
			}
			tempClusterNum++;
		}

		cout << "Layer0: There are " << tempClusterNum << " clusters." << endl;

		//冒泡继续排序;

		vector<float> &currentClusterScale = clusterScale[0];


		currentClusterScale.resize(tempClusterNum);


		float *tempClusterScale = new float[tempClusterNum];
		calScale(tempClusterScale, tempVertexClusterIndex, tempClusterNum);


		int *clusterPointer = new int[tempClusterNum];

		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
			clusterPointer[cIndex] = cIndex;

		bubbleSort(tempClusterNum, tempClusterScale, clusterPointer);



		int *clusterSequence = new int[tempClusterNum];
		for (int i = 0; i<tempClusterNum; i++)
		{
			clusterSequence[clusterPointer[i]] = i;
			currentClusterScale[i] = tempClusterScale[clusterPointer[i]];
		}

		delete[] tempClusterScale;
		delete[] clusterPointer;

		for (int i = 0; i<vertexNum; i++)
			vertexClusterIndex[i] = clusterSequence[tempVertexClusterIndex[i]];

		delete[] tempVertexClusterIndex;
		delete[] clusterSequence;
		
		clusterNum[0] = tempClusterNum;

	}

	void TriMesh::clusterCluster(int layerIndex, float angle)
	{

		clusterAngle[layerIndex] = angle;

		float angleTher = cos(M_PIf*angle / 180.f);

		int beforClusterNum = clusterNum[layerIndex];

		int downLayer = layerIndex - 1;
		need_topology(downLayer);
		vector<vector<int>> &downTopology = clusterTopology[downLayer];
		
		float maxClusterSize = upClusterScale[downLayer];
		float minClusterSize = lowClusterScale[downLayer];


		int downClusterNum = clusterNum[downLayer];

		need_clusterScale(downLayer);
		vector<float> &downClusterScale = clusterScale[downLayer];

		need_cluster(downLayer);
		vector<int> &downVertexClusterIndex = vertexCluster[downLayer];

		need_clusterNormal(downLayer);
		vector<vec> &downClusterNormal = clusterNormal[downLayer];


		vec *downTotalClusterNormal = new vec[downClusterNum];
		for (int cIndex = 0; cIndex < downClusterNum; cIndex++)
			downTotalClusterNormal[cIndex] += downClusterNormal[cIndex] * downClusterScale[cIndex];

		int* clusterCLusterIndex = new int[downClusterNum];
		int *tempTopology = new int[downClusterNum];
		int seedClusterNum = 0;

		cout << "Clustering clusters: [" << minClusterSize << "-" << maxClusterSize << "]~[" <<
			downClusterScale[downClusterNum - 1] << "-" << downClusterScale[0] << "]. " << endl;
		
		//计算每个聚类周围一圈以及两圈的法向;
		vec *suroundClusterNormal = new vec[downClusterNum];	//法向;

		for (int cIndex = 0; cIndex<downClusterNum; cIndex++)
		{
			if (downClusterScale[cIndex] >= maxClusterSize)
			{
				clusterCLusterIndex[cIndex] = seedClusterNum;
				seedClusterNum++;
				continue;
			}
			clusterCLusterIndex[cIndex] = -1;

			suroundClusterNormal[cIndex] = vec(0.f);

			for (int i = 0; i<downClusterNum; i++)
				tempTopology[i] = 0;

			//第一圈;

			for (int vcIndex = 0; vcIndex < downTopology[cIndex].size(); vcIndex++)
			{
				int nextCluster = downTopology[cIndex][vcIndex];

				if (downClusterScale[nextCluster] >= maxClusterSize)
					continue;
				tempTopology[nextCluster] = 1;
				//mark second

				for (int vcIndexS = 0; vcIndexS < downTopology[nextCluster].size(); vcIndexS++)
				{
					int nextClusterS = downTopology[nextCluster][vcIndexS];
					tempTopology[nextClusterS] = 1;
				}


			}

			//第二圈;
			for (int i = 0; i<downClusterNum; i++)
			if (tempTopology[i] == 1 && downClusterScale[i] < maxClusterSize)
				suroundClusterNormal[cIndex] += downTotalClusterNormal[i];
			normalize(suroundClusterNormal[cIndex]);
		}
		delete[] tempTopology;
		//找种子1种子 -1已分块 0未分块;


		cout << seedClusterNum << " oversize. ";

		int oversizeclusterNum = seedClusterNum;
		//将法向近似的种子块合并成标准块;
		for (int cIndex = oversizeclusterNum; cIndex<downClusterNum; cIndex++)
		{
			//建立新块;

			if (clusterCLusterIndex[cIndex] >= 0)
				continue;

			clusterCLusterIndex[cIndex] = seedClusterNum;

			float* seedNormal = suroundClusterNormal[cIndex];

			vector<int> undoneClusters;
			undoneClusters.push_back(cIndex);

			while (!undoneClusters.empty())
			{
				int curCluster = undoneClusters.back();
				undoneClusters.pop_back();


				for (int vcIndex = 0; vcIndex < downTopology[curCluster].size(); vcIndex++)
				{
					int nextCluster = downTopology[curCluster][vcIndex];

					if (clusterCLusterIndex[nextCluster] >= 0 ||
						clusterCLusterIndex[nextCluster] == -2 - seedClusterNum)
						continue;


					float length = suroundClusterNormal[cIndex].dot(suroundClusterNormal[nextCluster]);
//					float length = downClusterNormal[cIndex].dot(downClusterNormal[nextCluster]);
					

					if (length >= angleTher)
					{
						clusterCLusterIndex[nextCluster] = seedClusterNum;
						undoneClusters.push_back(nextCluster);
					}
					else
						clusterCLusterIndex[nextCluster] = -2 - seedClusterNum;
				}
			}
			seedClusterNum++;
		}

		cout << seedClusterNum << " seed. ";

		delete[] suroundClusterNormal;
		//去除小于阈值的块;
		int vertexNum = vertices.size();



		vector<float> tempClusterScale;
		tempClusterScale.resize(seedClusterNum);

		vector<int> vertexClusterIndex;
		vertexClusterIndex.resize(vertexNum);




		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			vertexClusterIndex[vIndex] = clusterCLusterIndex[downVertexClusterIndex[vIndex]];


		calScaleVec(tempClusterScale, vertexClusterIndex, seedClusterNum);

		int unCLustered = 0;
		for (int cIndex = 0; cIndex<downClusterNum; cIndex++)
		if (clusterCLusterIndex[cIndex] < 0)
			unCLustered++;



		for (int cIndex = 0; cIndex<downClusterNum; cIndex++)
		{
			int fIndex = clusterCLusterIndex[cIndex];
			if (fIndex >= 0 && tempClusterScale[fIndex] < minClusterSize)
			{
				clusterCLusterIndex[cIndex] = -1;
				unCLustered++;
			}
		}

		int tempClusterNum = oversizeclusterNum;

		int *clusterSequence = new int[seedClusterNum];
		for (int i = 0; i<seedClusterNum; i++)
			clusterSequence[i] = i;

		for (int cIndex = oversizeclusterNum; cIndex<seedClusterNum; cIndex++)
		{
			if (tempClusterScale[cIndex] >= minClusterSize)
			{
				clusterSequence[cIndex] = tempClusterNum;
				tempClusterNum++;
			}
		}

		clusterNum[layerIndex] = tempClusterNum;
		cout << tempClusterNum << " clusters." << endl;


		vector<vec> &currentClusterNormal = clusterNormal[layerIndex];
		currentClusterNormal.resize(tempClusterNum);
		vector<float> &currentClusterScale = clusterScale[layerIndex];
		currentClusterScale.resize(tempClusterNum);


		for (int cIndex = 0; cIndex < tempClusterNum; cIndex++)
			currentClusterNormal[cIndex] = vec(0.f);

		for (int cIndex = 0; cIndex<downClusterNum; cIndex++)
		{
			int fIndex = clusterCLusterIndex[cIndex];
			if (fIndex >= 0)
			{
				fIndex = clusterSequence[clusterCLusterIndex[cIndex]];
				clusterCLusterIndex[cIndex] = fIndex;
				currentClusterNormal[fIndex] += downTotalClusterNormal[cIndex];
			}
		}
		delete[] downTotalClusterNormal;
		delete[] clusterSequence;
		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
			normalize(currentClusterNormal[cIndex]);

		/////////////////////////////////////////////
		//将未分类的点吸收入标准块;
		////////////////////////////////////////////////


		//直到聚类个数稳定 逐渐扩大误差阈值;
		while (unCLustered >0)
		{

			//对每个聚类继续聚类;

			bool nothing = true;
			bool noTher = false;
			if (angleTher < 0)
			{
				noTher = true;
				cout << "no thre!";
			}


			for (int cIndex = 0; cIndex<downClusterNum; cIndex++)
			if (clusterCLusterIndex[cIndex] < 0)
				clusterCLusterIndex[cIndex] = -1;

			for (int cIndex = 0; cIndex<downClusterNum; cIndex++)
			{
				//寻找未分类的点;
				if (clusterCLusterIndex[cIndex] >= 0)
					continue;

				vector<int> undoneVertexes;
				//将其归类;
				int curCluster;


				for (int vcIndex = 0; vcIndex < downTopology[cIndex].size(); vcIndex++)
				{
					int nextCluster = downTopology[cIndex][vcIndex];

					curCluster = clusterCLusterIndex[nextCluster];

					if (curCluster >= 0)
					{
						float angle;
						if (!noTher)
							angle = downClusterNormal[cIndex].dot(currentClusterNormal[curCluster]);

						if (noTher || angle>angleTher)
						{
							clusterCLusterIndex[cIndex] = curCluster;
							undoneVertexes.push_back(cIndex);
							nothing = false;
							unCLustered--;
							break;
						}
					}
				}

				//是否归类;
				if (clusterCLusterIndex[cIndex] < 0)
					continue;

				//行走;
				while (!undoneVertexes.empty())
				{
					//处理末尾节点;
					int curVertex = undoneVertexes.back();

					undoneVertexes.pop_back();

					for (int vcIndex = 0; vcIndex < downTopology[cIndex].size(); vcIndex++)
					{
						int nextCluster = downTopology[cIndex][vcIndex];
						if (clusterCLusterIndex[nextCluster] >= 0 || clusterCLusterIndex[nextCluster] == -2 - cIndex)
							continue;

						float angle;
						if (!noTher)
							angle = downClusterNormal[cIndex].dot(currentClusterNormal[curCluster]);

						if (noTher || angle>angleTher)
						{
							clusterCLusterIndex[nextCluster] = curCluster;
							undoneVertexes.push_back(nextCluster);
							unCLustered--;
						}
						else
							clusterCLusterIndex[nextCluster] = -2 - cIndex;
					}

				}
			}
			//存储结果;
			cout << unCLustered << " ";

			if (nothing)
				angleTher -= 0.05f;
		}
		cout << endl;






		vector<int> clusterPointer;
		clusterPointer.resize(tempClusterNum);


		tempClusterScale.resize(tempClusterNum);
		

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			vertexClusterIndex[vIndex] = clusterCLusterIndex[downVertexClusterIndex[vIndex]];
		calScaleVec(tempClusterScale, vertexClusterIndex, tempClusterNum);

		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
			clusterPointer[cIndex] = cIndex;

		bubbleSortVec(tempClusterNum, tempClusterScale, clusterPointer);

		clusterSequence = new int[tempClusterNum];
		for (int i = 0; i<tempClusterNum; i++)
		{
			clusterSequence[clusterPointer[i]] = i;
			currentClusterScale[i] = tempClusterScale[clusterPointer[i]];
		}

		
		vector<int> &currentVertexCluster = vertexCluster[layerIndex];
		currentVertexCluster.resize(vertexNum);

		for (int i = 0; i<vertexNum; i++)
			currentVertexCluster[i] = clusterSequence[clusterCLusterIndex[downVertexClusterIndex[i]]];


		delete[] clusterSequence;
		delete[] clusterCLusterIndex;



		cout << "Layer" << layerIndex << ": There are " << tempClusterNum << " clusters." << endl;




	}
	
	void TriMesh::calScale(float* tempClusterScale, int* tempVertexClusterIndex, int tempClusterNum)
	{
		int* tempClusterLength = new int[tempClusterNum];


		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
		{
			tempClusterScale[cIndex] = 0.f;
			tempClusterLength[cIndex] = 0;
		}
		int vertexNum = vertices.size();

		bool* isEdge = new bool[vertexNum];

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			isEdge[vIndex] = false;

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
		{
			int cIndex = tempVertexClusterIndex[vIndex];
			if (cIndex<0)
				continue;
			if (isEdge[vIndex])
			{
				tempClusterLength[cIndex]++;
				continue;
			}


			tempClusterScale[cIndex] += 1.f;



			for (int vcIndex = 0; vcIndex < neighbors[vIndex].size(); vcIndex++)
			{
				int nextVertexIndex = neighbors[vIndex][vcIndex];
				if (tempVertexClusterIndex[nextVertexIndex] != cIndex)
				{
					isEdge[vIndex] = true;
					tempClusterLength[cIndex]++;
					isEdge[nextVertexIndex] = true;
					break;
				}
			}
		}

		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
			tempClusterScale[cIndex] /= (float)tempClusterLength[cIndex];

		delete[] tempClusterLength;
		delete[] isEdge;
	}


	void TriMesh::calScaleVec(vector<float> &tempClusterScale, vector<int> &tempVertexClusterIndex, int tempClusterNum)
	{



		vector<int> tempClusterLength;
		tempClusterLength.resize(tempClusterNum);





		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
		{
			tempClusterScale[cIndex] = 0.f;
			tempClusterLength[cIndex] = 1;
		}
		int vertexNum = vertices.size();

		bool* isEdge = new bool[vertexNum];

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			isEdge[vIndex] = false;

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
		{




			int cIndex = tempVertexClusterIndex[vIndex];


			if (cIndex<0)
				continue;
			if (isEdge[vIndex])
			{
				tempClusterLength[cIndex]++;
				continue;
			}


			tempClusterScale[cIndex] += 1.f;


			for (int vcIndex = 0; vcIndex < neighbors[vIndex].size(); vcIndex++)
			{
				int nextVertexIndex = neighbors[vIndex][vcIndex];
				if (tempVertexClusterIndex[nextVertexIndex] != cIndex)
				{
					isEdge[vIndex] = true;
					tempClusterLength[cIndex]++;
					isEdge[nextVertexIndex] = true;
					break;
				}
			}
		}

		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
			tempClusterScale[cIndex] /= (float)tempClusterLength[cIndex];

		delete[] isEdge;
	}




	void TriMesh::need_clusterScale(int layerIndex)
	{

		clamp(layerIndex, 0, 3);


		vector<float> &tempClusterScale = clusterScale[layerIndex];

		if (!tempClusterScale.empty())
			return;


		need_cluster(layerIndex);
		vector<int> &tempVertexClusterIndex = vertexCluster[layerIndex];


		int tempClusterNum = clusterNum[layerIndex];
		tempClusterScale.resize(tempClusterNum);

		int* tempClusterLength = new int[tempClusterNum];


		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
		{
			tempClusterScale[cIndex] = 0.f;
			tempClusterLength[cIndex] = 0;
		}
		int vertexNum = vertices.size();

		bool* isEdge = new bool[vertexNum];

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			isEdge[vIndex] = false;

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
		{
			int cIndex = tempVertexClusterIndex[vIndex];
			if (cIndex<0)
				continue;
			if (isEdge[vIndex])
			{
				tempClusterLength[cIndex]++;
				continue;
			}


			tempClusterScale[cIndex] += 1.f;



			for (int vcIndex = 0; vcIndex < neighbors[vIndex].size(); vcIndex++)
			{
				int nextVertexIndex = neighbors[vIndex][vcIndex];
				if (tempVertexClusterIndex[nextVertexIndex] != cIndex)
				{
					isEdge[vIndex] = true;
					tempClusterLength[cIndex]++;
					isEdge[nextVertexIndex] = true;
					break;
				}
			}
		}

		for (int cIndex = 0; cIndex<tempClusterNum; cIndex++)
			tempClusterScale[cIndex] /= (float)tempClusterLength[cIndex];

		delete[] tempClusterLength;
		delete[] isEdge;



	}

	void TriMesh::need_clusterSize(int layerIndex)
	{
		clamp(layerIndex, 0, 3);
		need_cluster(layerIndex);








		::std::vector<int> &currentVertexClusterSize = clusterSize[layerIndex];
		if (!currentVertexClusterSize.empty())
			return;
		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];

		int currentClusterNum = clusterNum[layerIndex];
		currentVertexClusterSize.resize(currentClusterNum);


		for (int i = 0; i < currentClusterNum; i++)
			currentVertexClusterSize[i] = 0;

		for (unsigned int i = 0; i < vertices.size(); i++)
			currentVertexClusterSize[currentVertexCluster[i]]++;





	}

	// Compute per-vertex normals
	void TriMesh::setColorCluster(int layerIndex)
	{
		clamp(layerIndex, 0, 3);

		need_cluster(layerIndex);
			

		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];

		need_clusterSize(layerIndex);
		::std::vector<int> &currentVertexClusterSize = clusterSize[layerIndex];

		int currentClusterNum = clusterNum[layerIndex];
		vector<vec3> clusterColor;
		clusterColor.resize(currentClusterNum);
	
		for (int i = 0; i < currentClusterNum; i++)
		{
			clusterColor[i][0] = (float)rand() / RAND_MAX;
			clusterColor[i][1] = (float)rand() / RAND_MAX;
			clusterColor[i][2] = (float)rand() / RAND_MAX;
		}
		::std::vector<int> &currentClusterSize = clusterSize[layerIndex];
		

//		float maxCLusterSize = 1.f / currentVertexClusterSize[1];


		if (colors.size() != vertices.size())
			colors.clear();
		if (colors.empty())
			colors.resize(vertices.size());

		for (unsigned int i = 0; i<vertices.size(); i++)
		for (int d = 0; d<3; d++)
			colors[i][d] = clusterColor[currentVertexCluster[i]][d];

	}

	void TriMesh::setNormalClusterNormal(int layerIndex)
	{
		int vertexNum = vertices.size();
		if (clusterNormalForVertex.empty())
			clusterNormalForVertex.resize(vertexNum);
		clamp(layerIndex, 0, 3);
		need_cluster(layerIndex);
		need_clusterNormal(layerIndex);

		vector<int> &curretVertexCluster = vertexCluster[layerIndex];
		vector<vec> &curretClusterNormal = clusterNormal[layerIndex];
		for (int verIndex = 0; verIndex < vertexNum; verIndex++)
			clusterNormalForVertex[verIndex] = curretClusterNormal[curretVertexCluster[verIndex]];
	}

	void TriMesh::need_topology(int layerIndex)
	{
		clamp(layerIndex, 0, 3);
		need_cluster(layerIndex);

		::std::vector<vector<int>> &currentClusterTopology = clusterTopology[layerIndex];
		if (!currentClusterTopology.empty())
			return;

		cout << "Building Layer Topology.....";
		int	currentClusterNum = clusterNum[layerIndex];
		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];
		int vertexNum = vertices.size();
		currentClusterTopology.resize(currentClusterNum);
		need_neighbors();

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
		{
			int currentClusterIndex = currentVertexCluster[vIndex];

			for (unsigned int verIndex = 0; verIndex < neighbors[vIndex].size(); verIndex++)
			{
				vector<int> &me = currentClusterTopology[currentClusterIndex];
				int nextClusterIndex = currentVertexCluster[neighbors[vIndex][verIndex]];
				if (find(me.begin(), me.end(), nextClusterIndex) == me.end())
					currentClusterTopology[currentClusterIndex].push_back(nextClusterIndex);
			}
		}
		
		cout << "Done!" << endl;
	}

	void TriMesh::need_averageEdgeLength()
	{
		if (!averageEdgeLength.empty())
			return;


		float currentAverageEdgeLength = 0.f;

		FILE *fin;
		fopen_s(&fin, dFName.getAvrageEdgeLengthDN(), "rb");
		if (fin)
		{
			fread(&currentAverageEdgeLength, sizeof(int), 1, fin);
			fclose(fin);
			averageEdgeLength.push_back(currentAverageEdgeLength);
		}


		cout << "Calculating averageEdgeLength.....";

		int faceNum = faces.size();
		for (int fIndex = 0; fIndex<faceNum; fIndex++)
		{
			int start[3];
			int end[3];
			for (int d = 0; d<3; d++)
				start[d] = faces[fIndex][d];
			end[0] = faces[fIndex][1];
			end[1] = faces[fIndex][2];
			end[2] = faces[fIndex][0];

			for (int d = 0; d<3; d++)
			if (start[d] > end[d])
				currentAverageEdgeLength += len(vertices[start[d]] - vertices[end[d]]);
		}
		currentAverageEdgeLength = currentAverageEdgeLength / 3.f / (float)faceNum;

		averageEdgeLength.push_back(currentAverageEdgeLength);

		FILE *fout;
		fopen_s(&fout, dFName.getAvrageEdgeLengthDN(), "wb");

		fwrite(&averageEdgeLength[0], sizeof(int), 1, fout);

		fclose(fout);



	}
	

	void TriMesh::need_clusterNormal(int layerIndex)
	{
		clamp(layerIndex, 0, 3);

		::std::vector<vec> &currentClusterNormal = clusterNormal[layerIndex];
		if (!currentClusterNormal.empty())
			return;
		int vertexNum = vertices.size();

		need_cluster(layerIndex);
		currentClusterNormal.resize(clusterNum[layerIndex]);
		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];


		for (int cIndex = 0; cIndex < clusterNum[layerIndex]; cIndex++)
			currentClusterNormal[cIndex] = vec3(0.f);

		for (int vIndex = 0; vIndex < vertexNum; vIndex++)
			currentClusterNormal[currentVertexCluster[vIndex]] += normals[vIndex];

		for (int cIndex = 0; cIndex < clusterNum[layerIndex]; cIndex++)
			normalize(currentClusterNormal[cIndex]);
	}

	void TriMesh::need_upAndLowScale(int layerIndex, int startIndex, int endIndex)
	{

		if (startIndex == endIndex)
			return;

		clamp(layerIndex, 0, 3);

		if (!averageClusterScale[layerIndex].empty())
			return;


		//确定上限和下限;

		float tempAverageClusterScale = -1.f;

		int total = 64;

		int sixtyfour[64] = { 0 };

		vector<float> currentClusterScale = clusterScale[layerIndex];

		float max = currentClusterScale[startIndex];
		float min = currentClusterScale[endIndex];

		int highest = 0;
		int curLeft = 0;
		int nextLeft = 0;
		int right = 0;
		bool hasHighest = false;

		float length = max - min + 0.0001f;
		length = 1.f / length;

		for (int i = startIndex; i<endIndex; i++)
		{
			float index = (float)currentClusterScale[i] - min;
			index *= length;
			index *= total;
			if (index >= total || index <0)
				cout << "error";
			sixtyfour[(int)index]++;
		}

		//反的；
		for (int i = 1; i<total; i++)
		{
			if (sixtyfour[i] > sixtyfour[i - 1])
				hasHighest = false;

			if (sixtyfour[i] < sixtyfour[i - 1])
			{
				if (hasHighest)
					right = i;

				if (sixtyfour[i - 1] > highest)
				{
					highest = sixtyfour[i - 1];
					curLeft = nextLeft;
					hasHighest = true;
				}
				nextLeft = i;

			}
		}

		for (int i = 0; i<curLeft; i++)
			endIndex -= sixtyfour[i];


		for (int i = total - 1; i>right; i--)
			startIndex += sixtyfour[i];

		//	if((clusterNum - endIndex + startIndex)*10 > clusterNum)
		if (startIndex * 10 > clusterNum[layerIndex])
		{
			upClusterScale[layerIndex] = max;
			lowClusterScale[layerIndex] = min;

			tempAverageClusterScale = 0.f;
			int temp = 0;
			for (int cIndex = 0; cIndex<clusterNum[layerIndex]; cIndex++)
			{
				if (currentClusterScale[cIndex] < max && currentClusterScale[cIndex]>min)
				{
					tempAverageClusterScale += currentClusterScale[cIndex];
					temp++;
				}
			}
			tempAverageClusterScale /= (float)temp;
		}
		else
			need_upAndLowScale(layerIndex, startIndex, endIndex);


		averageClusterScale[layerIndex].push_back(tempAverageClusterScale);

		return;

	}
	
	void TriMesh::need_calHalfArc(int layerIndex)
	{
		clamp(layerIndex, 0, 3);

		::std::vector<vector<float>> &currentClusterHalfArc = clusterHalfArc[layerIndex];

		if (!currentClusterHalfArc.empty())
			return;

		need_topology(layerIndex);
		need_averageEdgeLength();
		need_clusterNormal(layerIndex);

		float currentAverageEdgeLength = averageEdgeLength[0]*2.f;
		currentAverageEdgeLength = 0.5f / currentAverageEdgeLength*2.f;

		::std::vector<vector<int>> &currentClusterTopology = clusterTopology[layerIndex];
		::std::vector<vec> &currentClusterNormal = clusterNormal[layerIndex];
		currentClusterHalfArc.resize(currentClusterTopology.size());

		for (unsigned int currentClusterIndex = 0; currentClusterIndex < currentClusterTopology.size(); currentClusterIndex++)
		{
			currentClusterHalfArc[currentClusterIndex].resize(currentClusterTopology[currentClusterIndex].size());
			for (unsigned int clusterIndex = 0; clusterIndex < currentClusterTopology[currentClusterIndex].size(); clusterIndex++)
			{
				int nextClusterIndex = currentClusterTopology[currentClusterIndex][clusterIndex];
				currentClusterHalfArc[currentClusterIndex][clusterIndex] = fabs(acos(currentClusterNormal[currentClusterIndex].dot(currentClusterNormal[nextClusterIndex]))*currentAverageEdgeLength);


			}
		}
	}

	void TriMesh::smoothNormal(int layerIndex, int radius, bool useCluster, float smoothness)
	{

		useCluster = !useCluster;
		need_averageEdgeLength();
		clamp(layerIndex, 0, 3);

		::std::vector<vec> &currentSmoothedNormal = smoothedNormal[layerIndex];

		need_clusterNormal(layerIndex);

		int vertexNum = vertices.size();

		if (currentSmoothedNormal.empty())
		{
			smoothRadius[layerIndex] = -1;
			currentSmoothedNormal.resize(vertices.size());
		}
		if (radius <= 0)
		{
			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				currentSmoothedNormal[vIndex] = normals[vIndex];
			radius = -1;
		}
		
		if (radius == smoothRadius[layerIndex])
			return;

		smoothRadius[layerIndex] = radius;

		if (radius <= 0)
			return;


		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];


		need_curvatures();

		int tempR = 2 * radius;

		int *waitVertex = new int[vertexNum];		//待处理的点;
		float* curvatureWeight = new float[vertexNum];
		int* loopEnd = new int[tempR];	//点的距离;
		float* distantWeight = new float[tempR];//点的权重;



		float *reverseCuvature = new float[vertexNum];
		for (int vIndex = 0; vIndex < vertexNum; vIndex++)
		{
			if (_isnan(curv1[vIndex]))
				reverseCuvature[vIndex] = 0.1;
			else
				reverseCuvature[vIndex] = fabs(smoothness / curv1[vIndex]);
		}


		for (int dIndex = 0; dIndex<tempR; dIndex++)
			distantWeight[dIndex] = (float)exp(-(float)dIndex*(float)dIndex / (2.f*(float)tempR*(float)tempR));

		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
			curvatureWeight[vertexIndex] = -1;

		need_calHalfArc(layerIndex);
		need_neighbors();
		::std::vector<vector<int>> &currentClusterTopology = clusterTopology[layerIndex];
		::std::vector<vector<float>> &currentClusterHalfArc = clusterHalfArc[layerIndex];


		//每个顶点滤波;
		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		{
			int waitPush = 0;									//待处理点 进;
			curvatureWeight[vertexIndex] = 1.f;

			waitVertex[waitPush++] = vertexIndex;			//待处理点;
			loopEnd[0] = 1;
			int disStartIndex = 0;
			//宽度遍历radius圈;
			for (int curDistance = 0; curDistance<tempR - 1;)
			{
				//遍历一圈;
				for (int waitPop = disStartIndex; waitPop<loopEnd[curDistance]; waitPop++)
				{
					int curVertex = waitVertex[waitPop];
					float curWeight = curvatureWeight[curVertex];
					int curCluster = currentVertexCluster[curVertex];

					for (unsigned int tvIndex = 0; tvIndex < neighbors[curVertex].size(); tvIndex++)
					{
						int nextVertex = neighbors[curVertex][tvIndex];
						//如果没有处理过;
						if (curvatureWeight[nextVertex] >= 0)
							continue;
						//同一块或在边缘;
						int nextCLuster = currentVertexCluster[nextVertex];

						if (useCluster || nextCLuster == curCluster)// || curDistance == 0)
						{
							curvatureWeight[nextVertex] = curWeight;
							waitVertex[waitPush++] = nextVertex;	//加入待处理队列;
							continue;
						}

						//计算弧长;
						unsigned int nextClusterIndex;		//待处理kai的相邻kuai开始;
						for (nextClusterIndex = 0; nextClusterIndex < currentClusterTopology[curCluster].size(); nextClusterIndex++)
						if (currentClusterTopology[curCluster][nextClusterIndex] == nextCLuster)
							break;

						//在弧长范围之外;
						float arcLength = currentClusterHalfArc[curCluster][nextClusterIndex] * reverseCuvature[curVertex];



						if (curDistance >= arcLength)
							continue;
						//在弧长范围之内;
						curvatureWeight[nextVertex] = curWeight*(1.f - curDistance / arcLength);
						
						waitVertex[waitPush++] = nextVertex;	//加入待处理队列;

					}
				}
				disStartIndex = loopEnd[curDistance];
				loopEnd[++curDistance] = waitPush;
			}


			vec3 tempNormal;
			currentSmoothedNormal[vertexIndex] = normals[vertexIndex];

			curvatureWeight[vertexIndex] = -1.f;
			for (int curDistance = 1; curDistance<tempR; curDistance++)
			{
				tempNormal = vec(0.f);
				for (int topStaIndex = loopEnd[curDistance - 1]; topStaIndex<loopEnd[curDistance]; topStaIndex++)
				{
					int curVertex = waitVertex[topStaIndex];
					
					tempNormal += normals[curVertex] * curvatureWeight[curVertex];

					curvatureWeight[curVertex] = -1.f;
				}
				currentSmoothedNormal[vertexIndex] += tempNormal *distantWeight[curDistance];
			}
			normalize(currentSmoothedNormal[vertexIndex]);
		}

		delete[] curvatureWeight;
		delete[] waitVertex;
		delete[] reverseCuvature;
	}

	void TriMesh::getSmoothNormalCluster(vector<vec3>& resultNormal, int layerIndex)
	{

		clamp(layerIndex, 0, 3);

		::std::vector<vec> &currentSmoothedNormal = smoothedNormal[layerIndex];
		resultNormal.resize(vertices.size());

		int vertexNum = vertices.size();
		if (currentSmoothedNormal.empty())
		{

			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				resultNormal[vIndex] = normals[vIndex];


		}
		else
		{

			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				resultNormal[vIndex] = currentSmoothedNormal[vIndex];
		}

	}


	void TriMesh::needSmoothNormal(int radius)
	{

		need_averageEdgeLength();

		int vertexNum = vertices.size();

		if (smoothedNormalOrigin.empty())
		{
			smoothedNormalOrigin.resize(vertices.size());
		}
		if (radius <= 0)
		{
			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				smoothedNormalOrigin[vIndex] = normals[vIndex];
			radius = -1;
		}
		
		if (radius <= 0)
			return;



		int tempR = 2 * radius;

		int *waitVertex = new int[vertexNum];		//待处理的点;
		float* curvatureWeight = new float[vertexNum];
		int* loopEnd = new int[tempR];	//点的距离;
		float* distantWeight = new float[tempR];//点的权重;

		for (int dIndex = 0; dIndex<tempR; dIndex++)
			distantWeight[dIndex] = (float)exp(-(float)dIndex*(float)dIndex / (2.f*(float)tempR*(float)tempR));

		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
			curvatureWeight[vertexIndex] = -1;

		need_neighbors();

		//每个顶点滤波;
		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		{
			int waitPush = 0;									//待处理点 进;
			curvatureWeight[vertexIndex] = 1.f;

			waitVertex[waitPush++] = vertexIndex;			//待处理点;
			loopEnd[0] = 1;
			int disStartIndex = 0;
			//宽度遍历radius圈;
			for (int curDistance = 0; curDistance<tempR - 1;)
			{
				//遍历一圈;
				for (int waitPop = disStartIndex; waitPop<loopEnd[curDistance]; waitPop++)
				{
					int curVertex = waitVertex[waitPop];
					float curWeight = curvatureWeight[curVertex];

					for (unsigned int tvIndex = 0; tvIndex < neighbors[curVertex].size(); tvIndex++)
					{
						int nextVertex = neighbors[curVertex][tvIndex];
						//如果没有处理过;
						if (curvatureWeight[nextVertex] >= 0)
							continue;
						//同一块或在边缘;

						curvatureWeight[nextVertex] = curWeight;
						waitVertex[waitPush++] = nextVertex;	//加入待处理队列;
						continue;
					}
				}
				disStartIndex = loopEnd[curDistance];
				loopEnd[++curDistance] = waitPush;
			}


			vec3 tempNormal;
			smoothedNormalOrigin[vertexIndex] = normals[vertexIndex];

			curvatureWeight[vertexIndex] = -1.f;
			for (int curDistance = 1; curDistance<tempR; curDistance++)
			{
				tempNormal = vec(0.f);
				for (int topStaIndex = loopEnd[curDistance - 1]; topStaIndex<loopEnd[curDistance]; topStaIndex++)
				{
					int curVertex = waitVertex[topStaIndex];

					tempNormal += normals[curVertex];

					curvatureWeight[curVertex] = -1.f;
				}
				smoothedNormalOrigin[vertexIndex] += tempNormal *distantWeight[curDistance];
			}
			normalize(smoothedNormalOrigin[vertexIndex]);
		}

		delete[] curvatureWeight;
		delete[] waitVertex;
	}

	void TriMesh::getSmoothNormal(vector<vec3>& resultNormal, float radius)
	{

		need_averageEdgeLength();

		int vertexNum = vertices.size();

		if (resultNormal.empty())
		{
			resultNormal.resize(vertices.size());
		}
		if (radius <= 0.5)
		{
			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				resultNormal[vIndex] = normals[vIndex];
			radius = -1;
		}

		if (radius <= 0)
			return;



		int tempR = 2.f * radius;

		int *waitVertex = new int[vertexNum];		//待处理的点;
		float* curvatureWeight = new float[vertexNum];
		int* loopEnd = new int[tempR];	//点的距离;
		float* distantWeight = new float[tempR];//点的权重;

		for (int dIndex = 0; dIndex<tempR; dIndex++)
			distantWeight[dIndex] = (float)exp(-(float)dIndex*(float)dIndex / (2.f*(float)tempR*(float)tempR));

		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
			curvatureWeight[vertexIndex] = -1;

		need_neighbors();

		//每个顶点滤波;
		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		{
			int waitPush = 0;									//待处理点 进;
			curvatureWeight[vertexIndex] = 1.f;

			waitVertex[waitPush++] = vertexIndex;			//待处理点;
			loopEnd[0] = 1;
			int disStartIndex = 0;
			//宽度遍历radius圈;
			for (int curDistance = 0; curDistance<tempR - 1;)
			{
				//遍历一圈;
				for (int waitPop = disStartIndex; waitPop<loopEnd[curDistance]; waitPop++)
				{
					int curVertex = waitVertex[waitPop];
					float curWeight = curvatureWeight[curVertex];

					for (unsigned int tvIndex = 0; tvIndex < neighbors[curVertex].size(); tvIndex++)
					{
						int nextVertex = neighbors[curVertex][tvIndex];
						//如果没有处理过;
						if (curvatureWeight[nextVertex] >= 0)
							continue;
						//同一块或在边缘;

						curvatureWeight[nextVertex] = curWeight;
						waitVertex[waitPush++] = nextVertex;	//加入待处理队列;
						continue;
					}
				}
				disStartIndex = loopEnd[curDistance];
				loopEnd[++curDistance] = waitPush;
			}


			vec3 tempNormal;
			resultNormal[vertexIndex] = normals[vertexIndex];

			curvatureWeight[vertexIndex] = -1.f;
			for (int curDistance = 1; curDistance<tempR; curDistance++)
			{
				tempNormal = vec(0.f);
				for (int topStaIndex = loopEnd[curDistance - 1]; topStaIndex<loopEnd[curDistance]; topStaIndex++)
				{
					int curVertex = waitVertex[topStaIndex];

					tempNormal += normals[curVertex];

					curvatureWeight[curVertex] = -1.f;
				}
				resultNormal[vertexIndex] += tempNormal *distantWeight[curDistance];
			}
			normalize(resultNormal[vertexIndex]);
		}

		delete[] curvatureWeight;
		delete[] waitVertex;
	}

	void TriMesh::setColorClusterNormal(int layerIndex)

	{
		if (colors.size() != vertices.size())
			colors.resize(vertices.size());

		clamp(layerIndex, 0, 3);
		
		need_clusterNormal(layerIndex);
		vector<int> &currentVertexCluster = vertexCluster[layerIndex];
		vector<vec> &currentClusterNormal = clusterNormal[layerIndex];

		for (unsigned int vIndex = 0; vIndex<vertices.size(); vIndex++)
		for (int d = 0; d<3; d++)
			colors[vIndex][d] = abs(currentClusterNormal[currentVertexCluster[vIndex]][d] + 1.f) / 2.f;
	}

	void TriMesh::setColorSmoothedNormal(int layerIndex, int radius, bool useCluster, float smoothness)
	{

		clamp(layerIndex, 0, 3);

		if (colors.size() != vertices.size())
			colors.resize(vertices.size());

		smoothNormal(layerIndex, radius, useCluster, smoothness);

		for (unsigned int vIndex = 0; vIndex<vertices.size(); vIndex++)
		for (int d = 0; d<3; d++)
			colors[vIndex][d] = abs(smoothedNormal[layerIndex][vIndex][d] + 1.f) / 2.f;
	}
	
	void TriMesh::smoothNormalWhithOutCluster(int layerIndex, int radius)
	{
		clamp(layerIndex, 0, 3);

		::std::vector<vec> &currentSmoothedNormal = smoothedNormal[layerIndex];

		need_clusterNormal(layerIndex);

		int vertexNum = vertices.size();

		if (currentSmoothedNormal.empty())
		{
			smoothRadius[layerIndex] = -1;
			currentSmoothedNormal.resize(vertices.size());
		}
		if (radius <= 0)
		{
			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				currentSmoothedNormal[vIndex] = normals[vIndex];
		}

		if (radius == smoothRadius[layerIndex])
			return;

		smoothRadius[layerIndex] = radius;

		if (radius <= 0)
			return;


		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];


		need_curvatures();

		int tempR = 2 * radius;

		int *waitVertex = new int[vertexNum];		//待处理的点;
		float* curvatureWeight = new float[vertexNum];
		int* loopEnd = new int[tempR];	//点的距离;
		float* distantWeight = new float[tempR];//点的权重;



		float *reverseCuvature = new float[vertexNum];
		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			reverseCuvature[vIndex] = fabs(1.f / curv1[vIndex]);

		for (int dIndex = 0; dIndex<tempR; dIndex++)
			distantWeight[dIndex] = (float)exp(-(float)dIndex*(float)dIndex / (2.f*(float)tempR*(float)tempR));

		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
			curvatureWeight[vertexIndex] = -1;

		need_calHalfArc(layerIndex);
		need_neighbors();
		::std::vector<vector<int>> &currentClusterTopology = clusterTopology[layerIndex];
		::std::vector<vector<float>> &currentClusterHalfArc = clusterHalfArc[layerIndex];


		//每个顶点滤波;
		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		{
			int waitPush = 0;									//待处理点 进;
			curvatureWeight[vertexIndex] = 1.f;

			waitVertex[waitPush++] = vertexIndex;			//待处理点;
			loopEnd[0] = 1;
			int disStartIndex = 0;
			//宽度遍历radius圈;
			for (int curDistance = 0; curDistance<tempR - 1;)
			{
				//遍历一圈;
				for (int waitPop = disStartIndex; waitPop<loopEnd[curDistance]; waitPop++)
				{
					int curVertex = waitVertex[waitPop];
					float curWeight = curvatureWeight[curVertex];
					int curCluster = currentVertexCluster[curVertex];

					for (unsigned int tvIndex = 0; tvIndex < neighbors[curVertex].size(); tvIndex++)
					{
						int nextVertex = neighbors[curVertex][tvIndex];
						//如果没有处理过;
						if (curvatureWeight[nextVertex] >= 0)
							continue;
						//同一块或在边缘;
						int nextCLuster = currentVertexCluster[nextVertex];


						curvatureWeight[nextVertex] = curWeight;
						waitVertex[waitPush++] = nextVertex;	//加入待处理队列;

					}
				}
				disStartIndex = loopEnd[curDistance];
				loopEnd[++curDistance] = waitPush;
			}


			vec3 tempNormal;
			currentSmoothedNormal[vertexIndex] = normals[vertexIndex];

			curvatureWeight[vertexIndex] = -1.f;
			for (int curDistance = 1; curDistance<tempR; curDistance++)
			{
				tempNormal = vec(0.f);
				for (int topStaIndex = loopEnd[curDistance - 1]; topStaIndex<loopEnd[curDistance]; topStaIndex++)
				{
					int curVertex = waitVertex[topStaIndex];

					tempNormal += normals[curVertex] * curvatureWeight[curVertex];

					curvatureWeight[curVertex] = -1.f;
				}
				currentSmoothedNormal[vertexIndex] += tempNormal * distantWeight[curDistance];
			}
			normalize(currentSmoothedNormal[vertexIndex]);
		}

		delete[] curvatureWeight;
		delete[] waitVertex;
	}

	void TriMesh::smoothNormalWithEdge(int layerIndex, int radius)
	{
		clamp(layerIndex, 0, 3);

		::std::vector<vec> &currentSmoothedNormal = smoothedNormal[layerIndex];

		need_clusterNormal(layerIndex);

		int vertexNum = vertices.size();

		if (currentSmoothedNormal.empty())
		{
			smoothRadius[layerIndex] = -1;
			currentSmoothedNormal.resize(vertices.size());
		}
		if (radius <= 0)
		{
			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				currentSmoothedNormal[vIndex] = normals[vIndex];
		}

		if (radius == smoothRadius[layerIndex])
			return;

		smoothRadius[layerIndex] = radius;

		if (radius <= 0)
			return;


		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];


		need_curvatures();

		int tempR = 2 * radius;

		int *waitVertex = new int[vertexNum];		//待处理的点;
		float* curvatureWeight = new float[vertexNum];
		int* loopEnd = new int[tempR];	//点的距离;
		float* distantWeight = new float[tempR];//点的权重;



		float *reverseCuvature = new float[vertexNum];
		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			reverseCuvature[vIndex] = fabs(1.f / curv1[vIndex]);

		for (int dIndex = 0; dIndex<tempR; dIndex++)
			distantWeight[dIndex] = (float)exp(-(float)dIndex*(float)dIndex / (2.f*(float)tempR*(float)tempR));

		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
			curvatureWeight[vertexIndex] = -1;

		need_calHalfArc(layerIndex);
		need_neighbors();
		::std::vector<vector<int>> &currentClusterTopology = clusterTopology[layerIndex];
		::std::vector<vector<float>> &currentClusterHalfArc = clusterHalfArc[layerIndex];


		//每个顶点滤波;
		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		{

			bool isEdge = false;
			for (int nextVertex = 0; nextVertex < neighbors[vertexIndex].size(); nextVertex++)
				if (currentVertexCluster[neighbors[vertexIndex][nextVertex]] > currentVertexCluster[vertexIndex])
					isEdge = true;


			int waitPush = 0;									//待处理点 进;
			curvatureWeight[vertexIndex] = 1.f;

			waitVertex[waitPush++] = vertexIndex;			//待处理点;
			loopEnd[0] = 1;
			int disStartIndex = 0;
			//宽度遍历radius圈;
			for (int curDistance = 0; curDistance<tempR - 1;)
			{
				//遍历一圈;
				for (int waitPop = disStartIndex; waitPop<loopEnd[curDistance]; waitPop++)
				{
					int curVertex = waitVertex[waitPop];
					float curWeight = curvatureWeight[curVertex];
					int curCluster = currentVertexCluster[curVertex];

					for (unsigned int tvIndex = 0; tvIndex < neighbors[curVertex].size(); tvIndex++)
					{
						int nextVertex = neighbors[curVertex][tvIndex];
						//如果没有处理过;
						if (curvatureWeight[nextVertex] >= 0)
							continue;
						//同一块或在边缘;
						int nextCLuster = currentVertexCluster[nextVertex];

						if (nextCLuster == curCluster || isEdge)
						{
							curvatureWeight[nextVertex] = curWeight;
							waitVertex[waitPush++] = nextVertex;	//加入待处理队列;
							continue;
						}

					}
				}
				disStartIndex = loopEnd[curDistance];
				loopEnd[++curDistance] = waitPush;
			}


			vec3 tempNormal;
			currentSmoothedNormal[vertexIndex] = normals[vertexIndex];

			curvatureWeight[vertexIndex] = -1.f;
			for (int curDistance = 1; curDistance<tempR; curDistance++)
			{
				tempNormal = vec(0.f);
				for (int topStaIndex = loopEnd[curDistance - 1]; topStaIndex<loopEnd[curDistance]; topStaIndex++)
				{
					int curVertex = waitVertex[topStaIndex];

					tempNormal += normals[curVertex] * curvatureWeight[curVertex];

					curvatureWeight[curVertex] = -1.f;
				}
				currentSmoothedNormal[vertexIndex] += tempNormal * distantWeight[curDistance];
			}
			normalize(currentSmoothedNormal[vertexIndex]);
		}

		delete[] curvatureWeight;
		delete[] waitVertex;
	}

	void TriMesh::smoothNormalWhithCluster(int layerIndex, int radius)
	{
		clamp(layerIndex, 0, 3);

		::std::vector<vec> &currentSmoothedNormal = smoothedNormal[layerIndex];

		need_clusterNormal(layerIndex);

		int vertexNum = vertices.size();

		if (currentSmoothedNormal.empty())
		{
			smoothRadius[layerIndex] = -1;
			currentSmoothedNormal.resize(vertices.size());
		}
		if (radius <= 0)
		{
			for (int vIndex = 0; vIndex<vertexNum; vIndex++)
				currentSmoothedNormal[vIndex] = normals[vIndex];
		}

		if (radius == smoothRadius[layerIndex])
			return;

		smoothRadius[layerIndex] = radius;

		if (radius <= 0)
			return;


		::std::vector<int> &currentVertexCluster = vertexCluster[layerIndex];


		need_curvatures();

		int tempR = 2 * radius;

		int *waitVertex = new int[vertexNum];		//待处理的点;
		float* curvatureWeight = new float[vertexNum];
		int* loopEnd = new int[tempR];	//点的距离;
		float* distantWeight = new float[tempR];//点的权重;



		float *reverseCuvature = new float[vertexNum];
		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
			reverseCuvature[vIndex] = fabs(1.f / curv1[vIndex]);

		for (int dIndex = 0; dIndex<tempR; dIndex++)
			distantWeight[dIndex] = (float)exp(-(float)dIndex*(float)dIndex / (2.f*(float)tempR*(float)tempR));

		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
			curvatureWeight[vertexIndex] = -1;

		need_calHalfArc(layerIndex);
		need_neighbors();
		::std::vector<vector<int>> &currentClusterTopology = clusterTopology[layerIndex];
		::std::vector<vector<float>> &currentClusterHalfArc = clusterHalfArc[layerIndex];


		//每个顶点滤波;
		for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		{
			int waitPush = 0;									//待处理点 进;
			curvatureWeight[vertexIndex] = 1.f;

			waitVertex[waitPush++] = vertexIndex;			//待处理点;
			loopEnd[0] = 1;
			int disStartIndex = 0;
			//宽度遍历radius圈;
			for (int curDistance = 0; curDistance<tempR - 1;)
			{
				//遍历一圈;
				for (int waitPop = disStartIndex; waitPop<loopEnd[curDistance]; waitPop++)
				{
					int curVertex = waitVertex[waitPop];
					float curWeight = curvatureWeight[curVertex];
					int curCluster = currentVertexCluster[curVertex];

					for (unsigned int tvIndex = 0; tvIndex < neighbors[curVertex].size(); tvIndex++)
					{
						int nextVertex = neighbors[curVertex][tvIndex];
						//如果没有处理过;
						if (curvatureWeight[nextVertex] >= 0)
							continue;
						//同一块或在边缘;
						int nextCLuster = currentVertexCluster[nextVertex];

						if (nextCLuster == curCluster)
						{
							curvatureWeight[nextVertex] = curWeight;
							waitVertex[waitPush++] = nextVertex;	//加入待处理队列;
							continue;
						}
					}
				}
				disStartIndex = loopEnd[curDistance];
				loopEnd[++curDistance] = waitPush;
			}


			vec3 tempNormal;
			currentSmoothedNormal[vertexIndex] = normals[vertexIndex];

			curvatureWeight[vertexIndex] = -1.f;
			for (int curDistance = 1; curDistance<tempR; curDistance++)
			{
				tempNormal = vec(0.f);
				for (int topStaIndex = loopEnd[curDistance - 1]; topStaIndex<loopEnd[curDistance]; topStaIndex++)
				{
					int curVertex = waitVertex[topStaIndex];

					tempNormal += normals[curVertex] * curvatureWeight[curVertex];

					curvatureWeight[curVertex] = -1.f;
				}
				currentSmoothedNormal[vertexIndex] += tempNormal * distantWeight[curDistance];
			}
			normalize(currentSmoothedNormal[vertexIndex]);
		}

		delete[] curvatureWeight;
		delete[] waitVertex;
	}

	void TriMesh::need_smoothedCurvature(int layerIndex, int radius)
	{

		clamp(layerIndex, 0, 3);

		if (radius < 0)
			radius = 0;
		smoothNormal(layerIndex, radius, true, 1.0);
		

		::std::vector<float>  &currentsmoothedCurvature = smoothedCurvature[layerIndex];
		if (smoothRadiusForCurv[layerIndex] != radius || currentsmoothedCurvature.empty())
		{
			smoothRadiusForCurv[layerIndex] = radius;
			normals = smoothedNormal[layerIndex];
			reCalCurvatures();
			currentsmoothedCurvature = curv1;
			curv1.clear();
			normals.clear();
			need_normals();
		}

	}

	void TriMesh::getColorSmoothedCurvature(int layerIndex, int radius, float curv_thre_pos, float curv_thre_neg, float scaleCur)
	{

		clamp(layerIndex, 0, 3);

		need_smoothedCurvature(layerIndex, radius);


		if (colors.size() != vertices.size())
			colors.resize(vertices.size());

		smoothNormal(layerIndex, radius, true, 1.0);

		::std::vector<float>  &currentsmoothedCurvature = smoothedCurvature[layerIndex];

		float scale;
		if (scaleCur < 0.f)
			scale = -1.0 / scaleCur;
		else 
			scale = scaleCur;


		for (unsigned int i = 0; i<vertices.size(); i++)
		{
			colors[i][0] = currentsmoothedCurvature[i] * scale;
			colors[i][1] = -currentsmoothedCurvature[i] * scale;
			colors[i][2] = 0.0;

			float cur = colors[i][0] + colors[i][1];
			if (colors[i][0] >= curv_thre_pos)
				colors[i][2] = 1.0;
			if (colors[i][1] >= curv_thre_neg)
				colors[i][2] = 1.0;

		}

	}

	void TriMesh::needWeightedCurvature(float* weights, float scaleCur)
	{

		float weightSum = 0.f;
		for (int lIndex = 0; lIndex < 5; lIndex++)
		{
			if (weights[lIndex] <= 0)
				continue;
			weightSum += weights[lIndex];
		}


		weightedCurvature.resize(vertices.size());

		for (int vIndex = 0; vIndex < vertices.size(); vIndex++)
			weightedCurvature[vIndex] = 0;

		for (int lIndex = 0; lIndex < 4; lIndex++)
		{
			if (weights[lIndex] <= 0)
				continue;

			need_smoothedCurvature(lIndex, smoothRadius[lIndex]);

			float scale;

			if (scaleCur < 0.f)
				scale = -1.0/scaleCur * weights[lIndex] / weightSum;
			else
				scale = scaleCur * weights[lIndex] / weightSum;

			for (int vIndex = 0; vIndex < vertices.size(); vIndex++)
				weightedCurvature[vIndex] += smoothedCurvature[lIndex][vIndex] * scale;

		}

		if (weights[4] > 0)
		{
			need_curvatures();
			float scale;

			if (scaleCur < 0.f)
				scale = -1.0/ scaleCur * weights[4] / weightSum;
			else
				scale = scaleCur * weights[4] / weightSum;

			for (int vIndex = 0; vIndex < vertices.size(); vIndex++)
				weightedCurvature[vIndex] += curv1[vIndex] * scale;
		}

	}
}; // namespace trimesh
