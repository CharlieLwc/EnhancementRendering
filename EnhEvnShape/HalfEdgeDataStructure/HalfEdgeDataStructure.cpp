
#include "stdafx.h"

#include "HalfEdgeDataStructure.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>

using namespace std;


void HE_vert::calculateFaceNum()
{
	faceNum = 1;

	HE_edge *curEdge = edge->pair->next;

	while (curEdge->face != NULL && edge->index != curEdge->index)
	{
		faceNum++;
		curEdge = curEdge->pair->next;
	}

	if (curEdge->face == NULL)
	{
		curEdge = edge->next->next->pair;
		while (curEdge->face != NULL)
		{
			faceNum++;
			edge = curEdge;
			curEdge = curEdge->next->next->pair;
		}
	}
	return;
}


void HE_edge::getCenter(float *center)
{
	center[0] = 0.5f*(vert->x + pair->vert->x);
	center[1] = 0.5f*(vert->y + pair->vert->y);
	center[2] = 0.5f*(vert->z + pair->vert->z);
}

void HE_face::getVertices(HE_vert *faceVertices[])
{
	faceVertices[0] = edge->vert;
	faceVertices[1] = edge->next->vert;
	faceVertices[2] = edge->next->next->vert;
}

void HE_face::getVertexIndexes(int *faceVertexIndexes)
{
	faceVertexIndexes[0] = edge->vert->index;
	faceVertexIndexes[1] = edge->next->vert->index;
	faceVertexIndexes[2] = edge->next->next->vert->index;
}

void HE_face::getBoundEdges(HE_edge *faceEdges[])
{
	faceEdges[0] = edge;
	faceEdges[1] = edge->next;
	faceEdges[2] = edge->next->next;
}

void HE_face::getBoundEdgeIndexes(int *faceEdgeindexes)
{
	faceEdgeindexes[0] = edge->index;
	faceEdgeindexes[1] = edge->next->index;
	faceEdgeindexes[2] = edge->next->next->index;
}

void HE_face::getAdFaces(HE_face *adFaces[])
{
	adFaces[0] = edge->pair->face;
	adFaces[1] = edge->next->pair->face;
	adFaces[2] = edge->next->next->pair->face;
}

int HE_face::findAdFace(HE_face *adFace)
{
	
	if (adFace == edge->pair->face)
		return 0;
	if (adFace == edge->next->pair->face)
		return 1;
	if (adFace == edge->next->next->pair->face)
		return 2;

	cout<<endl<<"error in findAdFace!"<<endl;
	getchar();

	return -1;
}

int HE_face::findAdFace(int faceIndex)
{

	if (faceIndex == edge->pair->face->index)
		return 0;
	if (faceIndex == edge->next->pair->face->index)
		return 1;
	if (faceIndex == edge->next->next->pair->face->index)
		return 2;

	cout<<endl<<"error in findAdFace!"<<endl;
	getchar();

	return -1;
}


int HE_face::findVertex(int vertexIndex)
{
	if (vertexIndex == edge->vert->index)
		return 0;
	if (vertexIndex == edge->next->vert->index)
		return 1;
	if (vertexIndex == edge->next->next->vert->index)
		return 2;

	cout<<endl<<"error in findVertex!"<<endl;
	getchar();

	return -1;

}

void HE_face::findVertexAdEdgeIndex(int vertexAdEdgeIndexes[][2])
{
	vertexAdEdgeIndexes[0][0] = edge->index;
	vertexAdEdgeIndexes[1][0] = edge->next->index;
	vertexAdEdgeIndexes[2][0] = edge->next->next->index;
	vertexAdEdgeIndexes[0][1] = vertexAdEdgeIndexes[2][0];
	vertexAdEdgeIndexes[1][1] = vertexAdEdgeIndexes[0][0];
	vertexAdEdgeIndexes[2][1] = vertexAdEdgeIndexes[1][0];
}

HE_mesh::HE_mesh(void)
{
	vertices = NULL;
	edges = NULL;
	faces = NULL;

	everageEdgeLength = -1;
	edgeLength = NULL;
	faceEdgeCenterDistances = NULL;
	weightsForVertexNormal = NULL;
	faceAreas = NULL;
	faceNormals = NULL;
	vertexNormals = NULL;
	distanceBetweenAdFaceCenters = NULL;

	hasCalEL = false;
	hasCalFECD = false;
	hasCalWFVN = false;
	hasCalFA = false;
	hasCalFN = false;
	hasCalVN = false;
	hasCalDBAFC = false;
}

HE_mesh::HE_mesh(char *fileName)
{
	HE_mesh();
	init(fileName);
}

HE_mesh::~HE_mesh(void)
{
	deleteEdgeLength();
	deleteFaceAreas();
	deleteFaceEdgeCenterDistances();
	deleteWeightsForVertexNormal();
	deleteFaceNormals();
	deleteVertexNormals();
	deleteDistanceBetweenAdFaceCenters();

	if (vertices)
		delete[] vertices;
	if (faces)
		delete[] faces;
	if (edges)
		delete[] edges;

}

void HE_mesh::connectingEdge(char* fileName)
{



	ifstream model_f(fileName);


	char str[80]; 

	while(strcmp(str, "Vertices:"))//取得顶点、面片个数。;
		model_f>>str;
	model_f>>vertexNum;
	while(strcmp(str, "Faces:"))
		model_f>>str;
	model_f>>faceNum;
	//顶点、法向;
	vertices = new HE_vert[vertexNum];
	faces = new HE_face[faceNum];

	int (*faceVertexIndex)[3] = new int[faceNum][3];
	int *edgeAdEdge = new int[faceNum*3];
	vector<int> *vertexFacesIndex = new vector<int>[vertexNum];


	cout<<"reading "<<vertexNum<<" vertexes.... No normal!"<<endl;

	while(strcmp(str, "v"))//读取顶点;
		model_f>>str;


	for(int vertexIndex = 0; vertexIndex < vertexNum; vertexIndex++)
	{
		model_f>>vertices[vertexIndex].x;
		model_f>>vertices[vertexIndex].y;
		model_f>>vertices[vertexIndex].z;
		vertices[vertexIndex].index = vertexIndex;
		vertices[vertexIndex].edge = NULL;
		model_f>>str;

	}
	cout<<"reading "<<faceNum<<" faces....No normal!"<<endl;
	while(strcmp(str, "f"))//读取面片;
		model_f>>str;
	for(int faceIndex = 0; faceIndex < faceNum; faceIndex++)
	{
		faces[faceIndex].index = faceIndex;
		for (int dimention = 0; dimention<3; dimention++)
		{
			int vertexIndex;
			model_f>>vertexIndex;
			vertexIndex--;
			faceVertexIndex[faceIndex][dimention] = vertexIndex;
			vertexFacesIndex[vertexIndex].push_back(faceIndex);
		}
		model_f>>str;
	}
	model_f.close();

	cout<<"connecting faces...."<<endl;

	edgeNum = faceNum * 3;
	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
		edgeAdEdge[edgeIndex] = -1;

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
	{
		for (unsigned int curFaceIndex = 0; curFaceIndex<vertexFacesIndex[vertexIndex].size(); curFaceIndex++)//遍历面片;
		{
			int curFace = vertexFacesIndex[vertexIndex][curFaceIndex];//当前面片序号;
			int dimention = 0;
			for (; dimention<3; dimention++)//找到顶点/边在面片中序号;
				if (faceVertexIndex[curFace][dimention] == vertexIndex)
					break;
			int curEdge = curFace*3 + dimention;
			if (edgeAdEdge[curEdge] >= 0)//已经找到对应边则跳过;
				continue;
			dimention++;
			int nextVertex = faceVertexIndex[curFace][dimention%3];//边的另一个顶点;

			for (unsigned int nextFaceIndex = 0; nextFaceIndex<vertexFacesIndex[vertexIndex].size(); nextFaceIndex++)//遍历面片找相邻;
			{
				int nextFace = vertexFacesIndex[vertexIndex][nextFaceIndex];//下一个面片序号;
				if (nextFace == curFace)//相同的跳过;
					continue;
				bool notHasVertex = true;
				for (dimention = 0; dimention<3; dimention++)//找到顶点/边在面片中序号;
					if (faceVertexIndex[nextFace][dimention] == nextVertex)
					{
						notHasVertex = false;
						break;
					}
				if (notHasVertex)
					continue;

				int nextEdge = nextFace*3 + dimention;
				if (edgeAdEdge[nextEdge] != -1)
					continue;

				if (edgeAdEdge[curEdge] != -1 || edgeAdEdge[nextEdge] != -1)
				{
					int a=0;
					faceVertexIndex[curEdge/3];
					faceVertexIndex[nextEdge/3];
					faceVertexIndex[edgeAdEdge[curEdge] /3];
					faceVertexIndex[edgeAdEdge[nextEdge] /3];

				}


				edgeAdEdge[curEdge] = nextEdge;
				edgeAdEdge[nextEdge] = curEdge;
				break;
			}
		}
	}

	int nullEdge = 0;
	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
		if (edgeAdEdge[edgeIndex]<0)
			nullEdge++;

	edgeNum+=nullEdge;
	nullEdge = faceNum*3;

	edges = new HE_edge[edgeNum];
	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
	{
		edges[edgeIndex].index = edgeIndex;
		edges[edgeIndex].pair = NULL;
	}
	if (nullEdge < edgeNum)
	{
		for (int edgeIndex = nullEdge; edgeIndex<edgeNum; edgeIndex++)
		{
			edges[edgeIndex].face = NULL;
			edges[edgeIndex].next = NULL;
		}
	}

	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		int edgeIndex = faceIndex * 3;

		faces[faceIndex].edge = &edges[edgeIndex];

		edges[edgeIndex].next = &edges[edgeIndex+1];
		edges[edgeIndex+1].next = &edges[edgeIndex+2];
		edges[edgeIndex+2].next = &edges[edgeIndex];

		HE_face *curFace = &faces[faceIndex];

		for (int dimention = 0; dimention<3; dimention++)
		{
			int nxtDimention = (dimention+1)%3;
			int curVertex = faceVertexIndex[faceIndex][dimention];
			edges[edgeIndex].face = curFace;
			edges[edgeIndex].vert = &vertices[curVertex];
			if (vertices[curVertex].edge == NULL)
				vertices[curVertex].edge = &edges[edgeIndex];

			int nextEdgeIndex = edgeAdEdge[edgeIndex];

			if (nextEdgeIndex >= 0)
			{
				edges[edgeIndex].pair = &edges[nextEdgeIndex];
			}
			else
			{
				int nextVertexIndex = faceVertexIndex[faceIndex][(dimention+1)%3];
				nextEdgeIndex = nullEdge++;
				edges[edgeIndex].pair = &edges[nextEdgeIndex];
				edges[nextEdgeIndex].pair = &edges[edgeIndex];
				edges[nextEdgeIndex].vert = &vertices[nextVertexIndex];
				vertices[nextVertexIndex].edge = &edges[nextEdgeIndex];
			}
			edgeIndex++;
		}
	}

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		vertices[vertexIndex].calculateFaceNum();

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
	{
		if (vertices[vertexIndex].index != vertexIndex)
			cout<<" ";
		if (vertices[vertexIndex].edge->vert != &vertices[vertexIndex])
			cout<<" ";
	}

	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
	{
		if (edges[edgeIndex].index != edgeIndex)
			cout<<" ";
		if (edges[edgeIndex].pair == NULL || edges[edgeIndex].vert == NULL )
			cout<<" ";
		if (edges[edgeIndex].pair->pair != &edges[edgeIndex])
			cout<<" ";
		if (edges[edgeIndex].face == NULL)
		{
			if (edges[edgeIndex].next != NULL)
				cout<<" ";
		}
		else
		{
			if (edges[edgeIndex].next == NULL)
				cout<<" ";
			if (edges[edgeIndex].next->next->next!= &edges[edgeIndex])
				cout<<" ";
		}
	}


	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		if (faces[faceIndex].edge->face != &faces[faceIndex])
			cout<<" ";
		if (faces[faceIndex].index != faceIndex)
			cout<<" ";
		if (faces[faceIndex].edge->vert->index != faceVertexIndex[faceIndex][0])
			cout<<" ";
		if (faces[faceIndex].edge->next->vert->index != faceVertexIndex[faceIndex][1])
			cout<<" ";
		if (faces[faceIndex].edge->next->next->vert->index != faceVertexIndex[faceIndex][2])
			cout<<" ";
	}


	delete[] faceVertexIndex;
	delete[] edgeAdEdge;
	delete[] vertexFacesIndex;




	ofstream edge("d:\\d\\HE_edge.obj");

	edge<<fileName<<endl;

	edge<<vertexNum<<" "<<edgeNum<<" "<<faceNum<<endl;

	for (int vIndex = 0; vIndex<vertexNum; vIndex++)
	{
		edge<<vertices[vIndex].x<<" "<<
			vertices[vIndex].y<<" "<<
			vertices[vIndex].z<<" "<<
			vertices[vIndex].index<<" "<<
			vertices[vIndex].faceNum<<" "<<
			vertices[vIndex].edge->index<<endl;
	}
	for (int eIndex = 0; eIndex<edgeNum; eIndex++)
	{
		edge<<edges[eIndex].vert->index<<" ";
		if (edges[eIndex].pair != NULL)
			edge<<edges[eIndex].pair->index<<" ";
		else
			edge<<-1<<" ";

		if (edges[eIndex].face != NULL)
			edge<<edges[eIndex].face->index<<" ";
		else
			edge<<-1<<" ";

		if (edges[eIndex].next != NULL)
			edge<<edges[eIndex].next->index<<" ";
		else
			edge<<-1<<" ";

		edge<<edges[eIndex].index<<endl;
	}
	for (int fIndex = 0; fIndex<faceNum; fIndex++)
	{
		edge<<faces[fIndex].edge->index<<" "<<
			faces[fIndex].index<<endl;
	}

	edge<<"end";
	edge.close();


}

void HE_mesh::init(char *fileName)
{

	ifstream edge("d:\\d\\HE_edge.obj");

	if (!edge)
	{
		edge.close();
		connectingEdge(fileName);
		return;
	}

	char str[80]; 
	edge>>str;
	if (strcmp(str, fileName)!=0)
	{
		edge.close();
		connectingEdge(fileName);
		return;
	}

	edge>>vertexNum;
	edge>>edgeNum;
	edge>>faceNum;

	vertices = new HE_vert[vertexNum];
	edges = new HE_edge[edgeNum];
	faces = new HE_face[faceNum];

	for (int vIndex = 0; vIndex<vertexNum; vIndex++)
	{
		edge>>vertices[vIndex].x;
		edge>>vertices[vIndex].y;
		edge>>vertices[vIndex].z;
		edge>>vertices[vIndex].index;
		edge>>vertices[vIndex].faceNum;
		int tempIndex;
		edge>>tempIndex;
		vertices[vIndex].edge = &edges[tempIndex];
	}



	for (int eIndex = 0; eIndex<edgeNum; eIndex++)
	{
		int tempIndex;
		edge>>tempIndex;
		edges[eIndex].vert = &vertices[tempIndex];

		edge>>tempIndex;
		if (tempIndex >= 0)
			edges[eIndex].pair = &edges[tempIndex];
		else
			edges[eIndex].pair = NULL;

		edge>>tempIndex;
		if (tempIndex >= 0)
			edges[eIndex].face = &faces[tempIndex];
		else
			edges[eIndex].face = NULL;	

		edge>>tempIndex;
		if (tempIndex >= 0)
			edges[eIndex].next = &edges[tempIndex];
		else
			edges[eIndex].next = NULL;

		edge>>edges[eIndex].index;
	}
	for (int fIndex = 0; fIndex<faceNum; fIndex++)
	{
		int tempIndex;
		edge>>tempIndex;
		faces[fIndex].edge = &edges[tempIndex];
		edge>>faces[fIndex].index;
	}
	edge.close();
}

float *HE_mesh::getEdgeLength()
{
	if (edgeLength)
		return edgeLength;

	if (hasCalEL)
	{
		cout<<endl<<"EdgeLength has been calculated!"<<endl;
		getchar();
	}
	hasCalEL = true;

	edgeLength = new float[edgeNum];

	for (int edgeIndex=0; edgeIndex<edgeNum; edgeIndex++)
		if (edges[edgeIndex].isClockwise())
		{
			edgeLength[edgeIndex] = edges[edgeIndex].getLength();
			edgeLength[edges[edgeIndex].getPairIndex()] = edgeLength[edgeIndex];
		}

	return edgeLength;
}


float HE_mesh::getEdgeLength(int edgeIndex)
{
	getEdgeLength();
	return edgeLength[edgeIndex];
}


void HE_mesh::deleteEdgeLength()
{
	if (edgeLength)
	{
		delete[] edgeLength;
		edgeLength = NULL;
	}
}

float HE_mesh::getEverageEdgeLength()
{
	if (everageEdgeLength >= 0)
		return everageEdgeLength;

	getEdgeLength();

	float totalLength = 0.f;

	for (int edgeIndex=0; edgeIndex<edgeNum; edgeIndex++)
		totalLength += edgeLength[edgeIndex];

	everageEdgeLength = totalLength/edgeNum;

	return everageEdgeLength;
}

float (*HE_mesh::getFaceEdgeCenterDistances())[3]
{
	if (faceEdgeCenterDistances)
		return faceEdgeCenterDistances;

	if (hasCalFECD)
	{
		cout<<endl<<"FaceEdgeCenterDistances has been calculated!"<<endl;
		getchar();
	}
	hasCalFECD = true;

	faceEdgeCenterDistances = new float[faceNum][3];

	float faceCenters[3];
	HE_vert *faceVertices[3];

	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].getVertices(faceVertices);
		faceCenters[0] = (faceVertices[0]->x + faceVertices[1]->x + faceVertices[2]->x)/3.f;
		faceCenters[1] = (faceVertices[0]->y + faceVertices[1]->y + faceVertices[2]->y)/3.f;
		faceCenters[2] = (faceVertices[0]->z + faceVertices[1]->z + faceVertices[2]->z)/3.f;

		HE_edge *boundEdges[3];
		float edgeCenter[3];

		faces[faceIndex].getBoundEdges(boundEdges);
		for(int adEdgeIndex=0; adEdgeIndex<3; adEdgeIndex++)
		{
			boundEdges[adEdgeIndex]->getCenter(edgeCenter);
			faceEdgeCenterDistances[faceIndex][adEdgeIndex] = DISTANCE(edgeCenter, faceCenters);
		}
	}

	return faceEdgeCenterDistances;
};

void HE_mesh::deleteFaceEdgeCenterDistances()
{
	if (faceEdgeCenterDistances)
	{
		delete[] faceEdgeCenterDistances;
		faceEdgeCenterDistances = NULL;
	}
}

double (*HE_mesh::getWeightsForVertexNormal())[3]
{

	if (weightsForVertexNormal)
		return weightsForVertexNormal;

	if (hasCalWFVN)
	{
		cout<<endl<<"weightsForVertexNormal has been calculated!"<<endl;
		getchar();
	}
	hasCalWFVN = true;

	weightsForVertexNormal = new double[faceNum][3];

	getEdgeLength();
	getFaceAreas();

	double *reciprocalSquareOfEdgeLength = new double[edgeNum];

	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
		reciprocalSquareOfEdgeLength[edgeIndex] = 1.0/((double)edgeLength[edgeIndex] * (double)edgeLength[edgeIndex]);

	int vertexAdEdgeIndexes[3][2];

	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].findVertexAdEdgeIndex(vertexAdEdgeIndexes);
		for(int vertexIndex=0; vertexIndex<3; vertexIndex++)
			weightsForVertexNormal[faceIndex][vertexIndex] = faceAreas[faceIndex]
				* reciprocalSquareOfEdgeLength[vertexAdEdgeIndexes[vertexIndex][0]]
				* reciprocalSquareOfEdgeLength[vertexAdEdgeIndexes[vertexIndex][1]];
	}

	delete[] reciprocalSquareOfEdgeLength;
	return weightsForVertexNormal;
};

void HE_mesh::deleteWeightsForVertexNormal()
{
	if (weightsForVertexNormal)
	{
		delete[] weightsForVertexNormal;
		weightsForVertexNormal = NULL;
	}
}

double *HE_mesh::getFaceAreas()
{
	if (faceAreas)
		return faceAreas;


	if (hasCalFA)
	{
		cout<<endl<<"faceAreas has been calculated!"<<endl;
		getchar();
	}

	if (!faceNormals)
	{
		getFaceNormalsCalAreas();
		return faceAreas;
	}

	hasCalFA = true;
	faceAreas = new double[faceNum];

	getEdgeLength();
	int boundEdgeIndexes[3];
	float boundEdgeLength[3];
	double perimeter;

	for(int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].getBoundEdgeIndexes(boundEdgeIndexes);
		for (int dimention = 0; dimention<3; dimention++)
			boundEdgeLength[dimention] = edgeLength[boundEdgeIndexes[dimention]];

		perimeter = 0.5 * (boundEdgeLength[0]+boundEdgeLength[1]+boundEdgeLength[2]);
		faceAreas[faceIndex] = sqrt(perimeter*(perimeter-boundEdgeLength[0])*(perimeter-boundEdgeLength[1])*(perimeter-boundEdgeLength[2]));
	}
	return faceAreas;
};

void HE_mesh::deleteFaceAreas()
{
	if (faceAreas)
	{
		delete[] faceAreas;
		faceAreas = NULL;
	}
}

float (*HE_mesh::getFaceNormalsCalAreas())[3]
{
	if (faceNormals)
		return faceNormals;


	if (hasCalFN)
	{
		cout<<endl<<"faceNormals has been calculated!"<<endl;
		getchar();
	}
	if (faceAreas)
	{
		cout<<endl<<"faceAreas has been calculated!"<<endl;
		getchar();
	}

	hasCalFN = true;
	hasCalFA = true;

	faceNormals = new float[faceNum][3];
	faceAreas = new double[faceNum];

	double v1[3], v2[3], n[3];
	double length;
	HE_vert *faceVertices[3];

	for(int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].getVertices(faceVertices);
		VECTOR(v1, faceVertices[0], faceVertices[1]);
		VECTOR(v2, faceVertices[0], faceVertices[2]);
		CROSS(n, v1, v2);
		length = LENGTH(n);
		if (length == 0.0)
			length = 0.0;
		faceAreas[faceIndex] = 0.5*length;
		length = 1.0/length;
		for (int vertexIndex = 0; vertexIndex<3; vertexIndex++)
			faceNormals[faceIndex][vertexIndex] = (float)(n[vertexIndex]*length);

	}
	return faceNormals;
};

void HE_mesh::deleteFaceNormals()
{
	if (faceNormals)
	{
		delete[] faceNormals;
		faceNormals = NULL;
	}
}

float (*HE_mesh::getVertexNormals())[3]
{
	if (vertexNormals)
		return vertexNormals;


	if (hasCalVN)
	{
		cout<<endl<<"vertexNormals has been calculated!"<<endl;
		getchar();
	}

	hasCalVN = true;
	vertexNormals = new float[vertexNum][3];
	getFaceNormalsCalAreas();
	getWeightsForVertexNormal();
	double (*tempVertexNormals)[3] = new double[vertexNum][3];

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		for (int dimention = 0; dimention<3; dimention++)
			tempVertexNormals[vertexIndex][dimention] = 0.f;

	int faceVertexIndex[3];
	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].getVertexIndexes(faceVertexIndex);
		for (int vertexIndex = 0; vertexIndex<3; vertexIndex++)
			for (int dimention = 0; dimention<3; dimention++)
				tempVertexNormals[faceVertexIndex[vertexIndex]][dimention] 
					+= (weightsForVertexNormal[faceIndex][faces[faceIndex].findVertex(faceVertexIndex[vertexIndex])] 
					* (double)faceNormals[faceIndex][dimention]);
	}

	double length;
	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
	{
		length = 1.0/LENGTH(tempVertexNormals[vertexIndex]);
		for (int dimention = 0; dimention<3; dimention++)
			vertexNormals[vertexIndex][dimention] = (float)(length*tempVertexNormals[vertexIndex][dimention]);
	}
	delete[] tempVertexNormals;
	return vertexNormals;
};

void HE_mesh::deleteVertexNormals()
{
	if (vertexNormals)
	{
		delete[] vertexNormals;
		vertexNormals = NULL;
	}
}

float (*HE_mesh::getDistanceBetweenAdFaceCenters())[3]
{
	
	if (distanceBetweenAdFaceCenters)
		return distanceBetweenAdFaceCenters;


	if (hasCalDBAFC)
	{
		cout<<endl<<"distanceBetweenAdFaceCenters has been calculated!"<<endl;
		getchar();
	}
	getFaceEdgeCenterDistances();

	HE_face *adFaces[3];
	HE_face *nextFace;
	int nextFaceIndex;

	for(int faceIndex=0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].getAdFaces(adFaces);

		for(int adFaceIndex=0; adFaceIndex<3; adFaceIndex++)
		{
			nextFace = adFaces[adFaceIndex];
			if(nextFace == NULL)
				continue;
			nextFaceIndex = nextFace->getIndex();

			if(nextFaceIndex < faceIndex)
				continue;

			int nextFaceAdFaceIndex = nextFace->findAdFace(faceIndex);

			if (nextFaceAdFaceIndex < 0)
			{
				cout<<endl<<"error in findAdFace!"<<endl;
				getchar();
			}

			distanceBetweenAdFaceCenters[faceIndex][adFaceIndex] = faceEdgeCenterDistances[faceIndex][adFaceIndex] + faceEdgeCenterDistances[nextFaceIndex][nextFaceAdFaceIndex];
			distanceBetweenAdFaceCenters[nextFaceIndex][nextFaceAdFaceIndex] = distanceBetweenAdFaceCenters[faceIndex][adFaceIndex];

		}
	}
	deleteFaceEdgeCenterDistances();
	return distanceBetweenAdFaceCenters;
}

void HE_mesh::deleteDistanceBetweenAdFaceCenters()
{
	if (distanceBetweenAdFaceCenters)
	{
		delete[] distanceBetweenAdFaceCenters;
		distanceBetweenAdFaceCenters = NULL;
	}
}


void HE_mesh::clculateSmoothedVertexNormals(float (*smoothedFaceNormal)[3], float (*smoothedVertexNormal)[3])
{
	getFaceNormalsCalAreas();
	getWeightsForVertexNormal();
	double (*tempVertexNormals)[3] = new double[vertexNum][3];

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		for (int dimention = 0; dimention<3; dimention++)
			tempVertexNormals[vertexIndex][dimention] = 0.f;

	int faceVertexIndex[3];
	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		faces[faceIndex].getVertexIndexes(faceVertexIndex);
		for (int vertexIndex = 0; vertexIndex<3; vertexIndex++)
			for (int dimention = 0; dimention<3; dimention++)
				tempVertexNormals[faceVertexIndex[vertexIndex]][dimention] 
			+= (weightsForVertexNormal[faceIndex][faces[faceIndex].findVertex(faceVertexIndex[vertexIndex])] 
			* (double)smoothedFaceNormal[faceIndex][dimention]);
	}

	double length;
	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
	{
		length = 1.0/LENGTH(tempVertexNormals[vertexIndex]);
		for (int dimention = 0; dimention<3; dimention++)
			smoothedVertexNormal[vertexIndex][dimention] = (float)(length*tempVertexNormals[vertexIndex][dimention]);
	}
	delete[] tempVertexNormals;
	return ;



}