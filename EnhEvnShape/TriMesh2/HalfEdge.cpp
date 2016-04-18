


#include "HalfEdge.h"
#include "TriMesh.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <vector>

#include "../Eigen/Sparse"
#include "../Eigen/Dense"
#include "../Eigen/SparseCore"
#include "../Eigen/OrderingMethods"
#include "../Eigen/SparseQR"
#include "../Eigen/IterativeLinearSolvers"

#define E 2.7182818284590
using namespace Eigen;
using namespace std;

namespace trimesh {

void HE_vert::calculateFaceNum()
{

	if (edge == NULL)
		cout << "no faces isolated vertex";

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

int HE_vert::findEdge(int endVertex)
{

	HE_edge *curEdge = edge;

	do
	{
		if (curEdge->pair->vert->index == endVertex)
			return curEdge->index;
		curEdge = curEdge->pair->next;

	} while (curEdge->face != NULL && curEdge->index != edge->index);

	return -1;
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

void getDirections(vec3 directions[], vec planeNormal[], vec normal)
{
	directions[0] = vec3(1.f, 2.f, 0.f);
	/*
	if (abs(normal[0]) >= 0.5f)
	{
	directions[0][0] = 0.f;
	if (abs(normal[1]) >=0.5f)
	if (normal[2] < 0.f)
	directions[0][2] = -1.f;
	else
	directions[0][2] = 1.f;
	else
	if (normal[1] < 0.f)
	directions[0][1] = -1.f;
	else
	directions[0][1] = 1.f;
	}
	else
	if (normal[0] < 0.f)
	directions[0][1] = -1;
	*/


	normalize(directions[0]);

	float scale = directions[0].dot(normal);
	directions[0] -= normal * scale;
	normalize(directions[0]);

	vec3 direction8 = -directions[0];

	directions[4] = directions[0].cross(normal);
	normalize(directions[4]);
	directions[2] = directions[0] + directions[4];
	directions[6] = direction8 + directions[4];
	normalize(directions[2]);
	normalize(directions[3]);


	directions[1] = directions[0] + directions[2];
	directions[3] = directions[2] + directions[4];
	directions[5] = directions[4] + directions[6];
	directions[7] = directions[6] + direction8;

	normalize(directions[1]);
	normalize(directions[3]);
	normalize(directions[5]);
	normalize(directions[6]);

	planeNormal[0] = directions[4];
	planeNormal[1] = directions[5];
	planeNormal[2] = directions[6];
	planeNormal[3] = directions[7];
	planeNormal[4] = -directions[0];
	planeNormal[5] = -directions[1];
	planeNormal[6] = -directions[2];
	planeNormal[7] = -directions[3];

}

void TriMesh::connectingEdge()
{


	::std::vector<HE_vert> &he_vert = halfEdge.he_vert;
	::std::vector<HE_edge> &he_edge = halfEdge.he_edge;
	::std::vector<HE_face> &he_face = halfEdge.he_face;


	int vertexNum = vertices.size();
	int faceNum = faces.size();


	//顶点、法向;
	he_vert.resize(vertexNum);
	he_face.resize(faceNum);


	int *edgeAdEdge = new int[faceNum*3];
	vector<int> *vertexFacesIndex = new vector<int>[vertexNum];



	for(int vertexIndex = 0; vertexIndex < vertexNum; vertexIndex++)
	{
		he_vert[vertexIndex].index = vertexIndex;
		he_vert[vertexIndex].edge = NULL;
	}

	for(int faceIndex = 0; faceIndex < faceNum; faceIndex++)
	{
		he_face[faceIndex].index = faceIndex;
		for (int dimention = 0; dimention<3; dimention++)
		{
			vertexFacesIndex[faces[faceIndex].v[dimention]].push_back(faceIndex);
		}
	}

	cout<<"connecting faces...."<<endl;

	int edgeNum = faceNum * 3;
	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
		edgeAdEdge[edgeIndex] = -1;

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
	{
		for (unsigned int curFaceIndex = 0; curFaceIndex<vertexFacesIndex[vertexIndex].size(); curFaceIndex++)//遍历面片;
		{
			int curFace = vertexFacesIndex[vertexIndex][curFaceIndex];//当前面片序号;
			int dimention = 0;
			for (; dimention<3; dimention++)//找到顶点/边在面片中序号;
				if (faces[curFace].v[dimention] == vertexIndex)
					break;
			int curEdge = curFace*3 + dimention;
			if (edgeAdEdge[curEdge] >= 0)//已经找到对应边则跳过;
				continue;
			dimention++;
			int nextVertex = faces[curFace].v[dimention % 3];//边的另一个顶点;

			for (unsigned int nextFaceIndex = 0; nextFaceIndex<vertexFacesIndex[vertexIndex].size(); nextFaceIndex++)//遍历面片找相邻;
			{
				int nextFace = vertexFacesIndex[vertexIndex][nextFaceIndex];//下一个面片序号;
				if (nextFace == curFace)//相同的跳过;
					continue;
				bool notHasVertex = true;
				for (dimention = 0; dimention<3; dimention++)//找到顶点/边在面片中序号;
					if (faces[nextFace].v[dimention] == nextVertex)
					{
						notHasVertex = false;
						break;
					}
				if (notHasVertex)
					continue;

				int nextEdge = nextFace*3 + dimention;
				if (edgeAdEdge[nextEdge] != -1)
					continue;
				
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

	he_edge.resize(edgeNum);
	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
	{
		he_edge[edgeIndex].index = edgeIndex;
		he_edge[edgeIndex].pair = NULL;
	}
	if (nullEdge < edgeNum)
	{
		for (int edgeIndex = nullEdge; edgeIndex<edgeNum; edgeIndex++)
		{
			he_edge[edgeIndex].face = NULL;
			he_edge[edgeIndex].next = NULL;
		}
	}

	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		int edgeIndex = faceIndex * 3;

		he_face[faceIndex].edge = &he_edge[edgeIndex];

		he_edge[edgeIndex].next = &he_edge[edgeIndex + 1];
		he_edge[edgeIndex + 1].next = &he_edge[edgeIndex + 2];
		he_edge[edgeIndex + 2].next = &he_edge[edgeIndex];

		HE_face *curFace = &he_face[faceIndex];



		for (int dimention = 0; dimention<3; dimention++)
		{


			int nxtDimention = (dimention+1)%3;
			int curVertex = faces[faceIndex].v[dimention];

			he_edge[edgeIndex].face = curFace;
			he_edge[edgeIndex].vert = &he_vert[curVertex];
			if (he_vert[curVertex].edge == NULL)
				he_vert[curVertex].edge = &he_edge[edgeIndex];





			int nextEdgeIndex = edgeAdEdge[edgeIndex];

			if (nextEdgeIndex >= 0)
			{
				he_edge[edgeIndex].pair = &he_edge[nextEdgeIndex];
			}
			else
			{

				int nextVertexIndex = faces[faceIndex].v[(dimention+1)%3];
				nextEdgeIndex = nullEdge++;
				he_edge[edgeIndex].pair = &he_edge[nextEdgeIndex];
				he_edge[nextEdgeIndex].pair = &he_edge[edgeIndex];
				he_edge[nextEdgeIndex].next = &he_edge[nextEdgeIndex];
				he_edge[nextEdgeIndex].vert = &he_vert[nextVertexIndex];

			}
			edgeIndex++;
		}

		
	}

	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
		he_vert[vertexIndex].calculateFaceNum();


	for (int vertexIndex = 0; vertexIndex<vertexNum; vertexIndex++)
	{

		if (he_vert[vertexIndex].index != vertexIndex)
			cout<<" ";
		if (he_vert[vertexIndex].edge->vert != &he_vert[vertexIndex])
			cout << "one vertex attached two same faces";
	}

	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
	{
		if (he_edge[edgeIndex].index != edgeIndex)
			cout<<" ";
		if (he_edge[edgeIndex].pair == NULL || he_edge[edgeIndex].vert == NULL)
			cout<<" ";
		if (he_edge[edgeIndex].pair->pair != &he_edge[edgeIndex])
			cout<<" ";
		if (he_edge[edgeIndex].face == NULL)
		{
			if (he_edge[edgeIndex].next != &he_edge[edgeIndex])
				cout<<" ";
		}
		else
		{
			if (he_edge[edgeIndex].next == NULL)
				cout<<" ";
			if (he_edge[edgeIndex].next->next->next != &he_edge[edgeIndex])
				cout<<" ";
		}
	}


	for (int faceIndex = 0; faceIndex<faceNum; faceIndex++)
	{
		if (he_face[faceIndex].edge->face != &he_face[faceIndex])
			cout<<" ";
		if (he_face[faceIndex].index != faceIndex)
			cout<<" ";
		if (he_face[faceIndex].edge->vert->index != faces[faceIndex].v[0])
			cout<<" ";
		if (he_face[faceIndex].edge->next->vert->index != faces[faceIndex].v[1])
			cout<<" ";
		if (he_face[faceIndex].edge->next->next->vert->index != faces[faceIndex].v[2])
			cout<<" ";
	}


	delete[] edgeAdEdge;
	delete[] vertexFacesIndex;



}

void TriMesh::need_halfEdge()
{


	::std::vector<HE_vert> &he_vert = halfEdge.he_vert;
	::std::vector<HE_edge> &he_edge = halfEdge.he_edge;
	::std::vector<HE_face> &he_face = halfEdge.he_face;

	if (!he_vert.empty())
		return;

	int vertexNum = vertices.size();
	int faceNum = faces.size();
	ifstream edge_if(dFName.getHalfEdgeDN());

	if (edge_if.is_open())
	{
		int edgeNum;

		edge_if >> edgeNum;
		
		he_vert.resize(vertexNum);
		he_edge.resize(edgeNum);
		he_face.resize(faceNum);

		for (int vIndex = 0; vIndex<vertexNum; vIndex++)
		{
			edge_if >> he_vert[vIndex].index;
			edge_if >> he_vert[vIndex].faceNum;
			int tempIndex;
			edge_if >> tempIndex;
			he_vert[vIndex].edge = &he_edge[tempIndex];
		}


		for (int eIndex = 0; eIndex<edgeNum; eIndex++)
		{
			int tempIndex;
			edge_if >> tempIndex;
			he_edge[eIndex].vert = &he_vert[tempIndex];

			edge_if >> tempIndex;
			if (tempIndex >= 0)
				he_edge[eIndex].pair = &he_edge[tempIndex];
			else
				he_edge[eIndex].pair = NULL;

			edge_if >> tempIndex;
			if (tempIndex >= 0)
				he_edge[eIndex].face = &he_face[tempIndex];
			else
				he_edge[eIndex].face = NULL;

			edge_if >> tempIndex;
			if (tempIndex >= 0)
				he_edge[eIndex].next = &he_edge[tempIndex];
			else
				he_edge[eIndex].next = NULL;

			edge_if >> he_edge[eIndex].index;
		}
		for (int fIndex = 0; fIndex<faceNum; fIndex++)
		{
			int tempIndex;
			edge_if >> tempIndex;
			he_face[fIndex].edge = &he_edge[tempIndex];
			edge_if >> he_face[fIndex].index;
		}

		edge_if.close();
		return;
	}

	edge_if.close();


	connectingEdge();



	int edgeNum = he_edge.size();
	ofstream edge_of(dFName.getHalfEdgeDN());

	edge_of << edgeNum << endl;

	for (int vIndex = 0; vIndex<vertexNum; vIndex++)
	{
		edge_of << he_vert[vIndex].index << " " <<
			he_vert[vIndex].faceNum << " " <<
			he_vert[vIndex].edge->index << endl;
	}


	for (int eIndex = 0; eIndex<edgeNum; eIndex++)
	{
		edge_of << he_edge[eIndex].vert->index << " ";
		if (he_edge[eIndex].pair != NULL)
			edge_of << he_edge[eIndex].pair->index << " ";
		else
			edge_of << -1 << " ";

		if (he_edge[eIndex].face != NULL)
			edge_of << he_edge[eIndex].face->index << " ";
		else
			edge_of << -1 << " ";

		if (he_edge[eIndex].next != NULL)
			edge_of << he_edge[eIndex].next->index << " ";
		else
			edge_of << -1 << " ";

		edge_of << he_edge[eIndex].index << endl;
	}
	for (int fIndex = 0; fIndex<faceNum; fIndex++)
	{
		edge_of << he_face[fIndex].edge->index << " " <<
			he_face[fIndex].index << endl;
	}

	edge_of << "end";
	edge_of.close();

}

void TriMesh::need_edgeLength()
{


	if (!edgeLength.empty())
		return;

	::std::vector<HE_edge> &he_edge = halfEdge.he_edge;
	need_halfEdge();
	int edgeNum = he_edge.size();
	edgeLength.resize(edgeNum);

	FILE *fin;
	fopen_s(&fin, dFName.getEdgeLengthDN(), "rb");
	if (fin)
	{
		fread(&edgeLength[0], sizeof(float), edgeNum, fin);
		fclose(fin);
		return;
	}

	cout << "Calculating EdgeLength.....";


	for (int edgeIndex = 0; edgeIndex<edgeNum; edgeIndex++)
	if (he_edge[edgeIndex].isClockwise())
	{
		edgeLength[edgeIndex] = dist(vertices[he_edge[edgeIndex].getStartVertIndex()], vertices[he_edge[edgeIndex].getEndVertIndex()]);
		edgeLength[he_edge[edgeIndex].getPairIndex()] = edgeLength[edgeIndex];
	}


	FILE *fout;
	fopen_s(&fout, dFName.getEdgeLengthDN(), "wb");

	fwrite(&edgeLength[0], sizeof(float), edgeNum, fout);

	fclose(fout);



}

bool TriMesh::need_height()
{

	//smooth normal first

	::std::vector<HE_edge> &he_edge = halfEdge.he_edge;
	::std::vector<float> &vertexHeight = halfEdge.vertexHeight;
	float &maxHeight = halfEdge.maxHeight;
	float &minHeight = halfEdge.minHeight;


	if (!vertexHeight.empty())
		return true;
	int vertexNum = vertices.size();


	FILE *fin;
	fopen_s(&fin, dFName.getVertexHeightDN(), "rb");
	if (fin)
	{
		vertexHeight.resize(vertexNum);

		fread(&minHeight, sizeof(float), 1, fin);
		fread(&maxHeight, sizeof(float), 1, fin);
		fread(&vertexHeight[0], sizeof(float), vertexNum, fin);
		fclose(fin);

		return true;
	}

	if (smoothedNormalOrigin.empty())
	{
		cout << "smooth vertexNormal first!" << endl;
		return false;
	}

	need_halfEdge();

	vertexHeight.resize(vertexNum);

	cout << "Calculating EdgeHeight.....";

	std::vector<Eigen::Triplet<double>> coefficients;

	int rowIndex = 0;
	double *bVal = new double[he_edge.size()];
	
	vector<HE_edge> *curE = &he_edge;

	for (unsigned int eIndex = 0; eIndex < he_edge.size(); eIndex++)
	{
		if (he_edge[eIndex].index > he_edge[eIndex].pair->index)
			continue;

		int startIndex = he_edge[eIndex].getStartVertIndex();
		int endIndex = he_edge[eIndex].getEndVertIndex();

		coefficients.push_back(Eigen::Triplet<double>(rowIndex, startIndex, 1.0));
		coefficients.push_back(Eigen::Triplet<double>(rowIndex, endIndex, -1.0));



		vec3 vect1 = vertices[he_edge[eIndex].getStartVertIndex()] - vertices[he_edge[eIndex].getEndVertIndex()];
		vec3 vect2 = smoothedNormalOrigin[startIndex] + smoothedNormalOrigin[endIndex];
		normalize(vect2);
		double tempHigh = vect1.dot(vect2);

		bVal[rowIndex++] = tempHigh;
	}

	coefficients.push_back(Eigen::Triplet<double>(rowIndex, 0, 1.0));

	bVal[rowIndex++] = 0.0;

	Eigen::SparseMatrix<double, Eigen::ColMajor> A(rowIndex, vertexNum);
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	A.makeCompressed();


	Eigen::SparseMatrix<double> AT = SparseMatrix<double>(A.transpose());
	SparseMatrix<double> ATA = AT*A;

	Eigen::VectorXd b(rowIndex);
	for (int eIndex = 0; eIndex < rowIndex; eIndex++)
		b(eIndex) = bVal[eIndex];
	Eigen::VectorXd ATb = AT*b;


	ConjugateGradient<SparseMatrix<double> > matCG;
	matCG.compute(ATA);
	matCG.setMaxIterations(10);
	Eigen::VectorXd x = matCG.solve(ATb);

	minHeight = FLT_MAX;
	maxHeight = FLT_MIN;


	for (int i = 0; i < vertexNum; i++)
	{
		vertexHeight[i] = (float)x(i);
		if (minHeight > vertexHeight[i])
			minHeight = vertexHeight[i];
		if (maxHeight < vertexHeight[i])
			maxHeight = vertexHeight[i];
	}

	delete[] bVal;


	FILE *fout;
	fopen_s(&fout, dFName.getVertexHeightDN(), "wb");

	fwrite(&minHeight, sizeof(float), 1, fout);
	fwrite(&maxHeight, sizeof(float), 1, fout);
	fwrite(&vertexHeight[0], sizeof(float), vertexNum, fout);

	fclose(fout);
	return true;
}

void calScaleMatrix(Eigen::Matrix4f &scaleMatrix, int scale, float sampleStep)
{
	int halfSampleNum = (int)(((float)scale - 1.f)*0.5f);
	float t = -1.0f*halfSampleNum*sampleStep;

	Eigen::Matrix4f xVTtimeXV;
	Eigen::RowVector4f xVec;
	Eigen::Vector4f xVecTrans;

	scaleMatrix.setZero();

	for (int saIndex = 0; saIndex < scale; saIndex++)
	{

		float xDimention = 1;
		for (int i = 0; i < 4; i++)
		{
			xVec(0, i) = xDimention;
			xVecTrans(i, 0) = xDimention;
			xDimention *= t;
		}
		t += sampleStep;

		xVTtimeXV = xVecTrans * xVec;
		scaleMatrix += xVTtimeXV;
	}
}

bool calHitPoint(double &t, vec3 &planeVertex, vec3 &planeNormal, vec3 &edgeStart, vec3 &edgeEnd, vec3 &hitPoint, float edgeLength)
{
	vec3 vector = edgeEnd - edgeStart;

	t = planeNormal.dot(planeVertex - edgeStart) / vector.dot(planeNormal);

	if (t<0.f || t> 1.0 || _isnan(t))
		return false;

	vector *= (float)t;
	hitPoint = edgeStart + vector;

	if (t == 0.f || t == edgeLength)
		return false;

	return true;
}

void polyFitSamples(Eigen::Matrix4f &scaleMatrix, int scale, float *parmeter, float *sampleHeight, float sampleStep)
{


	Eigen::Vector4f parmeterVec;
	Eigen::Vector4f heighthVec;

	int halfSampleNum = (int)(((float)scale - 1.f)*0.5f);
	float steps = -1.0f*halfSampleNum*sampleStep;

	heighthVec.setZero();

	for (int saIndex = 0; saIndex < scale; saIndex++)
	{

		float xDimention = sampleHeight[saIndex];


		for (int i = 0; i < 4; i++)
		{
			heighthVec(i) += xDimention;
			xDimention *= steps;
		}
		steps += sampleStep;
	}

	parmeterVec = scaleMatrix.fullPivLu().solve(heighthVec);
	

	for (int i = 0; i < 4; i++)
		parmeter[i] = parmeterVec(i);

	return;
}

void TriMesh::findSamplesHeight(vec3 planeNormal, int sampleNum, int vertexIndex, float *sampleErrorEstimate)
{

	need_averageEdgeLength();
	need_edgeLength();
	need_height();


	::std::vector<HE_vert> &he_vert = halfEdge.he_vert;
	::std::vector<float> &sampleHeight = halfEdge.sampleHeight;
	int &blackFaceNum = halfEdge.blackFaceNum;
	::std::vector<float> &vertexHeight = halfEdge.vertexHeight;
	::std::vector<int> &blackFaceIndex = halfEdge.blackFaceIndex;

	HE_vert *curHEVertex = &he_vert[vertexIndex];
	vec3 &curVertex = vertices[vertexIndex];


	float sampleStep = averageEdgeLength[0] * 0.5f;

	double t = 0.f;
	int halfSampleNum = (int)(((float)sampleNum - 1)*0.5f);
	float *halfSampleHeight[2];
	vec3 *halfNormalField[2];

	bool isVertex[2] = { false, false };

	vec3 *sampleNormal = new vec3[sampleNum];

	halfSampleHeight[0] = &sampleHeight[0];
	halfSampleHeight[1] = &sampleHeight[halfSampleNum + 1];
	halfNormalField[0] = sampleNormal;
	halfNormalField[1] = &sampleNormal[halfSampleNum + 1];


	blackFaceNum = 0;
//	blackFaceIndex[blackFaceNum++] = vertexIndex;

	for (int sIndex = 0; sIndex<halfSampleNum; sIndex++)
	{
		halfSampleHeight[0][sIndex] = 0.f;
		halfSampleHeight[1][sIndex] = 0.f;
		halfNormalField[0][sIndex] = vec3(0.f);
		halfNormalField[1][sIndex] = vec3(0.f);
	}

	sampleHeight[halfSampleNum] = vertexHeight[vertexIndex];
	sampleNormal[halfSampleNum] = normals[vertexIndex];

	HE_edge *startEdge[3];
	vec3 hitPoint[3];

	HE_edge *faceEdge = curHEVertex->edge->next;
	HE_edge *leftEdge = faceEdge->next;



	int hitPointNum = 0;
	for (int fIndex = 0; fIndex<curHEVertex->faceNum; fIndex++)
	{
		if (faceEdge->pair->next->next->vert->index == vertexIndex)
		{
			leftEdge = leftEdge->next->pair;
			faceEdge = leftEdge->next->next;
			continue;

		}
		if (calHitPoint(t, curVertex, planeNormal, vertices[faceEdge->vert->index], vertices[leftEdge->vert->index], hitPoint[hitPointNum], edgeLength[faceEdge->index]))
		{
			startEdge[hitPointNum] = faceEdge;

			hitPointNum++;
			if (hitPointNum >= 3)
			{
				cout << "too many start point" << vertexIndex << endl;
				hitPointNum = 2;
				break;
			}
		}
		leftEdge = leftEdge->next->pair;
		faceEdge = leftEdge->next->next;
	}
	for (int dIndex = 0; dIndex<hitPointNum; dIndex++)
	{
		vec3 startPoint = curVertex;
		vec3 &endPoint = hitPoint[dIndex];
		float remainLength = 0.f;
		HE_edge *curEdge = startEdge[dIndex];

		for (int sIndex = 0; sIndex<halfSampleNum;)
		{
			
			float length = dist(startPoint, endPoint);
			float totalLength = length + remainLength;

			if (totalLength >= sampleStep)
			{
				float triLength = sampleStep - remainLength;

				vec3 vector = endPoint - startPoint;
				vector *= triLength;
				vec3 samplePoint = startPoint + vector;

				float weight[3];
				double weightSum = 0.0;


//				blackFaceIndex[blackFaceNum++] = curEdge->face->index;
				for (int vIndex = 0; vIndex<3; vIndex++)
				{

					weight[vIndex] = dist(samplePoint, vertices[curEdge->vert->index]);

					halfSampleHeight[dIndex][sIndex] += vertexHeight[curEdge->getStartVertIndex()] * weight[vIndex];
					halfNormalField[dIndex][sIndex] += normals[curEdge->getStartVertIndex()] * weight[vIndex];

					curEdge = curEdge->next;
					weightSum += weight[vIndex];
				}
				halfSampleHeight[dIndex][sIndex] /= (float)weightSum;
				normalize(halfNormalField[dIndex][sIndex]);
				remainLength = totalLength - sampleStep;
				sIndex++;
			}
			else
			{
				remainLength += length;
			}


			startPoint = hitPoint[dIndex];

			if (isVertex[dIndex])
			{
				if (curEdge->pair->next == curEdge->pair)
					break;

				HE_edge *tempEdge = curEdge->pair->next;

				for (int edgeIndex = 0; edgeIndex<curEdge->vert->faceNum - 1; edgeIndex++)
				{
					faceEdge = tempEdge->next;
					if (calHitPoint(t, curVertex, planeNormal, vertices[faceEdge->vert->index], vertices[faceEdge->next->vert->index], hitPoint[dIndex], edgeLength[faceEdge->index]))
					{
						tempEdge = faceEdge;
						break;
					}
					else
						tempEdge = tempEdge->pair->next;
				}
				if (tempEdge == curEdge)
					cout << "no hit point for one ling" << vertexIndex << endl;
				else
					curEdge = tempEdge;
				isVertex[dIndex] = false;

			}
			else
			{
				if (curEdge->pair->next == curEdge->pair)
					break;
				curEdge = curEdge->pair->next;

				if (curEdge->vert->index == vertexIndex || curEdge->next->vert->index == vertexIndex)
				{
					break;
				}

				if (!calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]))
				{
					curEdge = curEdge->next;
					if (!calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]))
					{

						if (abs(t)<0.000001)
						{
							isVertex[dIndex] = true;
							hitPoint[dIndex] = vertices[curEdge->vert->index];
						}
						else
						{
//							cout << "no hit point for two lines" << vertexIndex << endl;
							while (true)
							{
								break;
								curEdge = curEdge->next;

								calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]);

							}
							break;
						}


					}
				}

			}

		}
	}

	int sIndex1 = halfSampleNum - 1;
	vec3 tempN;
	float temp;
	for (int sIndex = 0; sIndex<halfSampleNum; sIndex++)
	{
		if (sIndex >= sIndex1)
			break;
		temp = sampleHeight[sIndex1];
		sampleHeight[sIndex1] = sampleHeight[sIndex];
		sampleHeight[sIndex] = temp;


		tempN = sampleNormal[sIndex1];
		sampleNormal[sIndex1] = sampleNormal[sIndex];
		sampleNormal[sIndex] = tempN;

		sIndex1--;
	}


	for (int i = 0; i<sampleNum; i++)
		sampleErrorEstimate[i] = planeNormal.dot(sampleNormal[i]);
	delete[] sampleNormal;
}

void TriMesh::findSamplesNormal(vec3 planeNormal, int sampleNum, int vertexIndex)
{

	::std::vector<::std::vector<int>> &leftFaceIndex = halfEdge.leftFaceIndex;
	::std::vector<::std::vector<int>> &rightFaceIndex = halfEdge.rightFaceIndex;
	::std::vector<::std::vector<vec3>> &leftSampleNormal = halfEdge.leftSampleNormal;
	::std::vector<::std::vector<vec3>> &rightSampleNormal = halfEdge.rightSampleNormal;
	::std::vector<::std::vector<vec3>> &leftSamplePosition = halfEdge.leftSamplePosition;
	::std::vector<::std::vector<vec3>> &rightSamplePosition = halfEdge.rightSamplePosition;


	::std::vector<HE_vert> &he_vert = halfEdge.he_vert;


	HE_vert *curHEVertex = &he_vert[vertexIndex];
	vec3 &curVertex = vertices[vertexIndex];

	need_averageEdgeLength();
	need_edgeLength();

	float sampleStep = averageEdgeLength[0] * 0.5f;

	double t = 0.f;

	vector<vec3> *halfNormalField[2];
	vector<int> *halfFaceIndex[2];
	vector<vec3> *halfPosition[2];


	halfNormalField[0] = &leftSampleNormal[vertexIndex];
	halfNormalField[1] = &rightSampleNormal[vertexIndex];
	halfNormalField[0]->clear();
	halfNormalField[1]->clear();


	halfFaceIndex[0] = &leftFaceIndex[vertexIndex];
	halfFaceIndex[1] = &rightFaceIndex[vertexIndex];
	halfFaceIndex[0]->clear();
	halfFaceIndex[1]->clear();


	halfPosition[0] = &leftSamplePosition[vertexIndex];
	halfPosition[1] = &rightSamplePosition[vertexIndex];
	halfPosition[0]->clear();
	halfPosition[1]->clear();


	bool isVertex[2] = { false, false };


	HE_edge *startEdge[3];
	vec3 hitPoint[3];

	HE_edge *faceEdge = curHEVertex->edge->next;
	HE_edge *leftEdge = faceEdge->next;



	int hitPointNum = 0;
	for (int fIndex = 0; fIndex<curHEVertex->faceNum; fIndex++)
	{
		if (faceEdge->pair->next->next->vert->index == vertexIndex)
		{	
			leftEdge = leftEdge->next->pair;
			faceEdge = leftEdge->next->next;
			continue;

		}
		if (calHitPoint(t, curVertex, planeNormal, vertices[faceEdge->vert->index], vertices[leftEdge->vert->index], hitPoint[hitPointNum], edgeLength[faceEdge->index]))
		{
			startEdge[hitPointNum] = faceEdge;

			hitPointNum++;
			if (hitPointNum >= 3)
			{
				cout << "too many start point" << endl;
				hitPointNum = 2;
				break;
			}
		}
		leftEdge = leftEdge->next->pair;
		faceEdge = leftEdge->next->next;
	}
	for (int dIndex = 0; dIndex<hitPointNum; dIndex++)
	{
		vec3 startPoint = curVertex;
		vec3 &endPoint = hitPoint[dIndex];
		float remainLength = 0.f;
		HE_edge *curEdge = startEdge[dIndex];


		for (int sIndex = 0; sIndex<sampleNum;)
		{
			halfPosition[dIndex]->push_back(endPoint);
			float length = dist(startPoint, endPoint);
			float totalLength = length + remainLength;

			if (totalLength >= sampleStep)
			{
				float triLength = sampleStep - remainLength;

				vec3 vector = endPoint - startPoint;
				vector *= triLength;
				vec3 samplePoint = startPoint + vector;

				float weight[3];
				double weightSum = 0.0;


				halfFaceIndex[dIndex]->push_back(curEdge->face->index);
				vec3 tempNormal(0.f);
				for (int vIndex = 0; vIndex<3; vIndex++)
				{

					weight[vIndex] = dist(samplePoint, vertices[curEdge->vert->index]);

					tempNormal += normals[curEdge->getStartVertIndex()] * weight[vIndex];

					curEdge = curEdge->next;
					weightSum += weight[vIndex];
				}
				normalize(tempNormal);
				halfNormalField[dIndex]->push_back(tempNormal);
				remainLength = totalLength - sampleStep;
				sIndex++;
			}
			else
			{
				remainLength += length;
			}


			startPoint = hitPoint[dIndex];

			if (isVertex[dIndex])
			{
				if (curEdge->pair->next == curEdge->pair)
					break;

				HE_edge *tempEdge = curEdge->pair->next;

				for (int edgeIndex = 0; edgeIndex<curEdge->vert->faceNum - 1; edgeIndex++)
				{
					faceEdge = tempEdge->next;
					if (calHitPoint(t, curVertex, planeNormal, vertices[faceEdge->vert->index], vertices[faceEdge->next->vert->index], hitPoint[dIndex], edgeLength[faceEdge->index]))
					{
						tempEdge = faceEdge;
						break;
					}
					else
						tempEdge = tempEdge->pair->next;
				}
				if (tempEdge == curEdge)
					cout << "no hit point for one ling" << endl;
				else
					curEdge = tempEdge;
				isVertex[dIndex] = false;

			}
			else
			{
				if (curEdge->pair->next == curEdge->pair)
					break;
				curEdge = curEdge->pair->next;

				if (curEdge->vert->index == vertexIndex || curEdge->next->vert->index == vertexIndex)
				{
					break;
				}

				if (!calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]))
				{
					curEdge = curEdge->next;
					if (!calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]))
					{

						if (abs(t)<0.000001)
						{
							isVertex[dIndex] = true;
							hitPoint[dIndex] = vertices[curEdge->vert->index];
						}
						else
						{
							cout << "no hit point for two lines" << endl;
							while (true)
							{
								break;
								curEdge = curEdge->next;

								calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]);

							}
						}


					}
				}

			}

		}
	}
}


void analysisCurvPlan(vec3 &position, vec3 &normal, vec3 &directionDir, float curveture, vector<vec3> &positions, vector<vec3> &normals, vector<float> &curvs, vector<bool> &isCurv, float curvThreshold, float planeThreshold)
{
	isCurv.clear();
	unsigned int i = 0;

	curvThreshold = abs(curveture*curvThreshold);
	for (; i < positions.size(); i++)
	{
		if (abs(curvs[i] - curveture) <= curvThreshold)
			isCurv.push_back(true);
		else
			break;
	}

	if (i >= positions.size())
		return;

	vec3 planeNormal = normals[i];

	for (; i < positions.size(); i++)
	{
		float difference = planeNormal.dot(normals[i]);

		if (difference > planeThreshold)
			isCurv.push_back(false);
		else
			break;
	}
	return;


	bool negative = false;
	vec3 ed = positions[1] - position;
	if (ed.dot(directionDir) < 0)
		negative = true;

	for (; i < positions.size(); i++)
	{
		ed = positions[i] - position;
		float dis = len(ed) * curveture;
		float angle = asin(dis*0.5f)*2.f;

		float nFactor = cos(angle);
		float dDactor = sin(angle);
		if (negative)
			dDactor = -dDactor;


		vec3 newNormal = normal*nFactor + directionDir * dDactor;
//		normals[i] = newNormal;
//		continue;


		
		float difference = newNormal.dot(normals[i]);

		if (difference > curvThreshold)
			isCurv.push_back(true);
		else
			break;


		//normals[i] *= difference;
	}

	if (i >= positions.size())
		return;

	planeNormal = normals[i];


	for (; i < positions.size(); i++)
	{
		float difference = planeNormal.dot(normals[i]);

		if (difference > planeThreshold)
			isCurv.push_back(false);
		else
			break;
	}
}


void TriMesh::findSamplesNormalFromEdge(int sampleNum, int sampleIndex, vec3& lineDir)
{


	::std::vector<::std::vector<int>> &leftFaceIndex = halfEdge.leftFaceIndex;
	::std::vector<::std::vector<int>> &rightFaceIndex = halfEdge.rightFaceIndex;
	::std::vector<::std::vector<vec3>> &leftSampleNormal = halfEdge.leftSampleNormal;
	::std::vector<::std::vector<vec3>> &rightSampleNormal = halfEdge.rightSampleNormal;

	::std::vector<::std::vector<vec3>> &leftSamplePosition = halfEdge.leftSamplePosition;
	::std::vector<::std::vector<vec3>> &rightSamplePosition = halfEdge.rightSamplePosition;


	::std::vector<HE_vert> &he_vert = halfEdge.he_vert;




	vec3 samplePosition, sampleNormal, planeNormal, directionDir;
	int edgeIndex;
	float curveture;
	bool tempIsVertex = false;

	ridgeLines.getPointInfo(sampleIndex, samplePosition, sampleNormal, curveture, directionDir, edgeIndex, tempIsVertex);

	planeNormal = directionDir.cross(sampleNormal);

	vec3 leftNormal = lineDir.cross(sampleNormal);



	need_averageEdgeLength();
	need_edgeLength();

	float sampleStep = averageEdgeLength[0] * 0.5f;

	double t = 0.f;

	vec3 &curVertex = samplePosition;

	vector<vec3> *halfNormalField[2];
	vector<int> *halfFaceIndex[2];
	vector<vec3> *halfPosition[2];
	vector<float> *halfCurv[2];


	halfNormalField[0] = &leftSampleNormal[sampleIndex];
	halfNormalField[1] = &rightSampleNormal[sampleIndex];
	halfNormalField[0]->clear();
	halfNormalField[1]->clear();


	halfFaceIndex[0] = &leftFaceIndex[sampleIndex];
	halfFaceIndex[1] = &rightFaceIndex[sampleIndex];
	halfFaceIndex[0]->clear();
	halfFaceIndex[1]->clear();





	halfPosition[0] = &leftSamplePosition[sampleIndex];
	halfPosition[1] = &rightSamplePosition[sampleIndex];
	halfPosition[0]->clear();
	halfPosition[1]->clear();



	halfCurv[0] = &halfEdge.leftSampleCurv[sampleIndex];
	halfCurv[1] = &halfEdge.rightSampleCurv[sampleIndex];
	halfCurv[0]->clear();
	halfCurv[1]->clear();



	halfNormalField[0]->reserve(sampleNum);
	halfNormalField[1]->reserve(sampleNum);
	halfFaceIndex[0]->reserve(sampleNum);
	halfFaceIndex[1]->reserve(sampleNum);
	halfPosition[0]->reserve(sampleNum);
	halfPosition[1]->reserve(sampleNum);
	halfCurv[0]->reserve(sampleNum);
	halfCurv[1]->reserve(sampleNum);




	bool isVertex[2] = { tempIsVertex, tempIsVertex };
	

	HE_edge *startEdge[2];

	vec3 edgeDir = vertices[halfEdge.he_edge[edgeIndex].getEndVertIndex()] - vertices[halfEdge.he_edge[edgeIndex].getStartVertIndex()];
	
	startEdge[0] = &halfEdge.he_edge[edgeIndex];
	startEdge[1] = halfEdge.he_edge[edgeIndex].pair;
	


	vec3 hitPoint[2] = { samplePosition, samplePosition };


	HE_edge *faceEdge;
	int hitPointNum = 2;

	for (int dIndex = 0; dIndex<hitPointNum; dIndex++)
	{
		vec3 startPoint = curVertex;
		vec3 &endPoint = hitPoint[dIndex];
		float remainLength = 0.f;
		HE_edge *curEdge = startEdge[dIndex];


		for (int sIndex = 0; sIndex<sampleNum;)
		{
			float length = dist(startPoint, endPoint);
			float totalLength = length + remainLength;

			if (totalLength >= sampleStep)
			{
				halfPosition[dIndex]->push_back(endPoint);
				float triLength = sampleStep - remainLength;

				vec3 vector = endPoint - startPoint;
				vector *= triLength;
				vec3 samplePoint = startPoint + vector;
				float tempCurv = 0.f;

				float weight[3];
				float weightSum = 0.f;


				halfFaceIndex[dIndex]->push_back(curEdge->face->index);
				vec3 tempNormal(0.f);
				for (int vIndex = 0; vIndex<3; vIndex++)
				{

					weight[vIndex] = dist(samplePoint, vertices[curEdge->vert->index]);

					tempNormal += normals[curEdge->getStartVertIndex()] * weight[vIndex];
					tempCurv += curv1[curEdge->getStartVertIndex()] * weight[vIndex];

					curEdge = curEdge->next;
					weightSum += weight[vIndex];
				}
				tempCurv /= weightSum;
				halfCurv[dIndex]->push_back(tempCurv);
				normalize(tempNormal);
				halfNormalField[dIndex]->push_back(tempNormal);
				remainLength = totalLength - sampleStep;
				sIndex++;
			}
			else
			{
				remainLength += length;
			}


			startPoint = hitPoint[dIndex];

			if (isVertex[dIndex])
			{
				if (curEdge->pair->next == curEdge->pair)
					break;

				HE_edge *tempEdge = curEdge->pair->next;

				for (int edgeIndex = 0; edgeIndex<curEdge->vert->faceNum - 1; edgeIndex++)
				{
					faceEdge = tempEdge->next;
					if (calHitPoint(t, curVertex, planeNormal, vertices[faceEdge->vert->index], vertices[faceEdge->next->vert->index], hitPoint[dIndex], edgeLength[faceEdge->index]))
					{
						tempEdge = faceEdge;
						break;
					}
					else
						tempEdge = tempEdge->pair->next;
				}
				if (tempEdge == curEdge)
					cout << "ERROR" << endl;
				else
					curEdge = tempEdge;
				isVertex[dIndex] = false;

			}
			else
			{
				if (curEdge->pair->next == curEdge->pair)
					break;
				curEdge = curEdge->pair->next;

				if (curEdge->index == edgeIndex)
				{
					break;
				}

				if (!calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]))
				{
					curEdge = curEdge->next;
					if (!calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]))
					{

						if (abs(t)<0.000001)
						{
							isVertex[dIndex] = true;
							hitPoint[dIndex] = vertices[curEdge->vert->index];
						}
						else
						{
							cout << "ERROR" << endl;
							while (true)
							{
								break;
								curEdge = curEdge->next;

								calHitPoint(t, curVertex, planeNormal, vertices[curEdge->vert->index], vertices[curEdge->next->vert->index], hitPoint[dIndex], edgeLength[curEdge->index]);

							}
						}


					}
				}

			}

		}
	}


	//all after *********************************************************//
	if(halfPosition[0]->empty() || halfPosition[1]->empty())
		return;
	//all after *********************************************************//



	vec3 curLeft = (*halfPosition[0])[0] - samplePosition;

	if (curLeft.dot(leftNormal) > 0.f)
	{

		vector<vec3> halfNormalField1;
		vector<int> halfFaceIndex1;
		vector<vec3> halfPosition1;
		vector<float> halfCurv1;


		halfNormalField1 = leftSampleNormal[sampleIndex];
		leftSampleNormal[sampleIndex] = rightSampleNormal[sampleIndex];
		rightSampleNormal[sampleIndex] = halfNormalField1;

		halfFaceIndex1 = leftFaceIndex[sampleIndex];
		leftFaceIndex[sampleIndex] = rightFaceIndex[sampleIndex];
		rightFaceIndex[sampleIndex] = halfFaceIndex1;

		halfPosition1 = leftSamplePosition[sampleIndex];
		leftSamplePosition[sampleIndex] = rightSamplePosition[sampleIndex];
		rightSamplePosition[sampleIndex] = halfPosition1;

		halfCurv1 = halfEdge.leftSampleCurv[sampleIndex];
		halfEdge.leftSampleCurv[sampleIndex] = halfEdge.rightSampleCurv[sampleIndex];
		halfEdge.rightSampleCurv[sampleIndex] = halfCurv1;



	}








}

void analysisPoly(int scale, float* parmeter, float &weight, float &feature, float enhanceScale, float &curvature, float &featureSlop)
{


	//all after *********************************************************//
	enhanceScale = 100.0;

	//all after *********************************************************//


	float t;
	float upT, downT, inverT, Tzero;
	float center = ((float)scale - 1.f)*0.5f;

	float *u = parmeter;

	if (abs(u[3]) > 0.0001)
	{
		double u3u1 = u[1] * u[3];
		double u2u2 = u[2] * u[2];
		double temp = sqrt(3.0*sqrt(9.0*u3u1*u3u1 - 6.0*u3u1*u2u2 + u2u2*u2u2 + 5.0*u[3] * u[3]) - 6.0*u3u1 + 2.0*u2u2) / (3.0*sqrt(5.0)*u[3]);
		inverT = -1.f*u[2] / (3.f*u[3]);
		upT = inverT - temp;
		downT = inverT + temp;


		Tzero = abs(upT);
		temp = abs(downT);
		t = upT;
		if (Tzero>temp)
		{
			Tzero = temp;
			t = downT;
		}

		weight = Tzero / center;
		weight = 1.0 - clamp(weight, 0.f, 1.f);

		

		temp = 2.0*Tzero / abs(upT - downT);
		temp = clamp(temp, 0.0, 1.0);
		temp = 1.0 - temp*temp;
		feature = pow(temp, (double)enhanceScale);

	}
	else if (abs(u[2])  > 0.0001)
	{
		upT = -1.0*u[1] / (2.0*u[2]);
		t = upT;

		Tzero = abs(upT);

		weight = Tzero / center;
		weight = 1.0 - clamp(weight, 0.f, 1.f);
		double temp = 2.0*Tzero / (Tzero + center);
		temp = clamp(temp, 0.0, 1.0);
		temp = 1.0 - temp*temp;
		feature = pow(temp, (double)enhanceScale);

	}
	else
	{
		feature = 0.f;
		t = 0.f;
		weight = 0.f;
	}

	if (weight == 0.f)
	{
		feature = 0.f;
		curvature = 0.f;
	}
	else
	{
		double curve = u[1];
		curve = 1.0 + curve*curve;
		curve = (2.0*u[2]) / sqrt(curve*curve*curve);



		curve = u[1] + 2.0*u[2] * t + 3.0*u[3] * t*t;
		curve = 1.0 + curve*curve;
		curvature = -(2.0*u[2] + 6.0*u[3] * t) / sqrt(curve*curve*curve);
	}

	featureSlop = -u[1];

	


}

void TriMesh::analysisCurvture(int scaleNum, int* scale, float scaleCur)
{

	/*
	::std::vector<float> &reliefFeature1 = halfEdge.reliefFeature;
	if (!reliefFeature1.empty())
		return;
	reliefFeature1.resize(vertices.size());
	ifstream normal_if("E:\\Copy\\EnhEvnShape\\Model\\circular 35 55 75.obj");
	if (normal_if.is_open())
	{
		cout << "Reading Normal.....";
		int totalNum;
		normal_if >> totalNum;

		for (int i = 0; i < totalNum; i++)
		{
			normal_if >> reliefFeature1[i];
			reliefFeature1[i] = - reliefFeature1[i];

		}
		cout << "Done" << endl;
		normal_if.close();
		return;
		
	}
	

	*/




	int tempScale[3] = { 35, 55, 75 };
	scaleNum = 3;
	scale = tempScale;

	need_halfEdge();

	


	
	bool needCalculate = true;
	for (int scaleIndex = 0; scaleIndex < scaleNum; scaleIndex++)
	{
		if (halfEdge.scaleValue[scaleIndex] != scale[scaleIndex] || halfEdge.singleReliefFeature[0].empty())
		{

			halfEdge.scaleValue[scaleIndex] = scale[scaleIndex];
			needCalculate = true;

			halfEdge.singleReliefFeature[scaleIndex].clear();

			FILE *fin;
			fopen_s(&fin, dFName.getReliefDN(scale[scaleIndex]), "rb");
			if (fin)
			{
				cout << "Reading SurfaceRelif.....";
				int totalNum;

				fread(&totalNum, sizeof(int), 1, fin);
				if (totalNum == vertices.size())
				{
					halfEdge.singleReliefFeature[scaleIndex].resize(vertices.size());
					halfEdge.singleReliefFeatureWeight[scaleIndex].resize(vertices.size());
					fread(&halfEdge.singleReliefFeature[scaleIndex][0], sizeof(float), totalNum, fin);
					fread(&halfEdge.singleReliefFeatureWeight[scaleIndex][0], sizeof(float), totalNum, fin);
				}
				else

				cout << "Done" << endl;
				fclose(fin);
			}
		}
	}


	if (!needCalculate)
		return;


	float needFindScale = false;
	for (int scaleIndex = 0; scaleIndex < scaleNum; scaleIndex++)
	{
		if (halfEdge.singleReliefFeature[scaleIndex].empty())
		{
			halfEdge.singleReliefFeature[scaleIndex].resize(vertices.size());
			halfEdge.singleReliefFeatureWeight[scaleIndex].resize(vertices.size());

			
			needFindScale = true;
		}
	}

	::std::vector<float> &reliefFeature = halfEdge.reliefFeature;

	reliefFeature.resize(vertices.size());

	if (needFindScale)
	{
		needFindScale = false;


		::std::vector<int> &blackFaceIndex = halfEdge.blackFaceIndex;

		blackFaceIndex.resize(scale[scaleNum - 1]);


		::std::vector<float> &sampleHeight = halfEdge.sampleHeight;
		sampleHeight.resize(scale[scaleNum - 1]);


		int vertexNum = vertices.size();

		vec3 directions[8];
		vec3 planeNormal[8];

		Eigen::Matrix4f scaleMatrix[3];


		need_averageEdgeLength();
		float sampleStep = averageEdgeLength[0] * 0.5;


		float feature[8][3];
		float featureSlop[8][3];
		float weight[8][3];
		float errorEstimate[8][3];
		float curvature[8][3];



		float enhanceScale = 1.0;

		float *sampleErrorEstimate = new float[scale[scaleNum - 1]];

		for (int sIndex = 0; sIndex < scaleNum; sIndex++)
			calScaleMatrix(scaleMatrix[sIndex], scale[sIndex], sampleStep);


		for (int vIndex = 0; vIndex < vertexNum; vIndex++)
		{

			// change normals to baseNormals
			getDirections(directions, planeNormal, normals[vIndex]);



			for (int dIndex = 0; dIndex<8; dIndex++)
			{
				int sampleNum = scale[scaleNum - 1];

				findSamplesHeight(planeNormal[dIndex], sampleNum, vIndex, sampleErrorEstimate);

				for (int sIndex = 0; sIndex < scaleNum; sIndex++)
				{
					int start = (scale[scaleNum - 1] - scale[sIndex])*0.5;
					float parmeter[4];
					polyFitSamples(scaleMatrix[sIndex], scale[sIndex], parmeter, &sampleHeight[start], sampleStep);

					analysisPoly(scale[sIndex], parmeter, weight[dIndex][sIndex], feature[dIndex][sIndex], enhanceScale, curvature[dIndex][sIndex], featureSlop[dIndex][sIndex]);

					errorEstimate[dIndex][sIndex] = 0.f;

					int end = start + scale[sIndex];
					for (int i = start; i<end; i++)
						errorEstimate[dIndex][sIndex] += sampleErrorEstimate[i];

					errorEstimate[dIndex][sIndex] /= (float)scale[sIndex];
					weight[dIndex][sIndex] = (1.0 - errorEstimate[dIndex][sIndex])*weight[dIndex][sIndex];
				}

			}



			float sumWeightScale[3] = { 0.f, 0.f, 0.f };
			float alphScale[3] = { 0.f, 0.f, 0.f };
			float avgErrorScale[3] = { 0.f, 0.f, 0.f };
			float stdErrorScale[3] = { 0.f, 0.f, 0.f };

			reliefFeature[vIndex] = 0.f;

			for (int sIndex = 0; sIndex<scaleNum; sIndex++)
			{
				for (int dIndex = 0; dIndex<8; dIndex++)
				{
					sumWeightScale[sIndex] += weight[dIndex][sIndex];
					avgErrorScale[sIndex] += errorEstimate[dIndex][sIndex];
				}

				float temp;
				if (sumWeightScale[sIndex] != 0.f)
				{
					temp = 1.0 / sumWeightScale[sIndex];
					for (int dIndex = 0; dIndex<8; dIndex++)
						weight[dIndex][sIndex] *= temp;
				}



				//std error;
				avgErrorScale[sIndex] /= 8.0;
				for (int dIndex = 0; dIndex<8; dIndex++)
					stdErrorScale[sIndex] = stdErrorScale[sIndex] + (errorEstimate[dIndex][sIndex] - avgErrorScale[sIndex])
					*(errorEstimate[dIndex][sIndex] - avgErrorScale[sIndex]);
				stdErrorScale[sIndex] /= 8.0;
				stdErrorScale[sIndex] = sqrt(stdErrorScale[sIndex]);

				//alph;
				for (int dIndex = 0; dIndex<8; dIndex++)
					alphScale[sIndex] = alphScale[sIndex] + weight[dIndex][sIndex] * curvature[dIndex][sIndex];
				alphScale[sIndex] *= stdErrorScale[sIndex];
				alphScale[sIndex] = 1.0 / (1.0 + pow(E, (double)-alphScale[sIndex]));

				float curScale = 0.f;
				for (int dIndex = 0; dIndex < 8; dIndex++)
					curScale = curScale + curvature[dIndex][sIndex] * feature[dIndex][sIndex] * weight[dIndex][sIndex];// ;
				//alphScale[0] = 1;
				//alphScale[1] = 5;


				halfEdge.singleReliefFeature[sIndex][vIndex] = curScale;

				halfEdge.singleReliefFeatureWeight[sIndex][vIndex] = alphScale[sIndex] * sumWeightScale[sIndex];


			}




		}

		delete[] sampleErrorEstimate;




		for (int scaleIndex = 0; scaleIndex < scaleNum; scaleIndex++)
		{

			FILE *fout;
			fopen_s(&fout, dFName.getReliefDN(scale[scaleIndex]), "wb");

			cout << "Reading SurfaceRelif.....";
			int totalNum = vertices.size();

			fwrite(&totalNum, sizeof(int), 1, fout);

			fwrite(&halfEdge.singleReliefFeature[scaleIndex][0], sizeof(float), totalNum, fout);
			fwrite(&halfEdge.singleReliefFeatureWeight[scaleIndex][0], sizeof(float), totalNum, fout);
			cout << "Done" << endl;
			fclose(fout);
			
			
		}




	}

	if (needCalculate)
	{

		int vertexNum = vertices.size();
		for (int vIndex = 0; vIndex < vertexNum; vIndex++)
		{
			reliefFeature[vIndex] = 0.f;
			float sumWeight = 0.f;

			for (int sIndex = 0; sIndex < scaleNum; sIndex++)
			{
				sumWeight += halfEdge.singleReliefFeatureWeight[sIndex][vIndex];
				reliefFeature[vIndex] += halfEdge.singleReliefFeature[sIndex][vIndex] * halfEdge.singleReliefFeatureWeight[sIndex][vIndex];
			}

			if (sumWeight != 0.f)
				reliefFeature[vIndex] = reliefFeature[vIndex] * scaleCur / sumWeight;


			if (reliefFeature[vIndex] >= 0.9)
				reliefFeature[vIndex] = 0.9;

			if (reliefFeature[vIndex] <= -0.9)
				reliefFeature[vIndex] = -0.9;


		}



	}






	/*

	for (int vIndex = 0; vIndex<vertexNum; vIndex++)
	{
		// change normals to baseNormals
		getDirections(directions, planeNormal, normals[vIndex]);



		for (int dIndex = 0; dIndex<8; dIndex++)
		{

			int sampleNum = scale[scaleNum - 1];
			float *sampleErrorEstimate = new float[sampleNum];

			findSamplesHeight(planeNormal[dIndex], sampleNum, vIndex, sampleErrorEstimate);



			for (int sIndex = 0; sIndex < scaleNum; sIndex++)
			{
				int start = (scale[scaleNum - 1] - scale[sIndex])*0.5;
				float parmeter[4];
				polyFitSamples(scaleMatrix[sIndex], scale[sIndex], parmeter, &sampleHeight[start], sampleStep);

				analysisPoly(scale[sIndex], parmeter, weight[dIndex][sIndex], feature[dIndex][sIndex], enhanceScale, curvature[dIndex][sIndex], featureSlop[dIndex][sIndex]);

				errorEstimate[dIndex][sIndex] = 0.f;

				int end = start + scale[sIndex];
				for (int i = start; i<end; i++)
					errorEstimate[dIndex][sIndex] += sampleErrorEstimate[i];

				errorEstimate[dIndex][sIndex] /= (float)scale[sIndex];
				weight[dIndex][sIndex] = (1.0 - errorEstimate[dIndex][sIndex])*weight[dIndex][sIndex];
			}

		}

		//weight error scale;


		float sumWeight = 0.f;
		float sumWeightScale[3] = { 0.f, 0.f, 0.f };
		float alphScale[3] = { 0.f, 0.f, 0.f };
		float avgErrorScale[3] = { 0.f, 0.f, 0.f };
		float stdErrorScale[3] = { 0.f, 0.f, 0.f };

		reliefFeature[vIndex] = 0.f;

		for (int sIndex = 0; sIndex<scaleNum; sIndex++)
		{
			for (int dIndex = 0; dIndex<8; dIndex++)
			{
				sumWeightScale[sIndex] += weight[dIndex][sIndex];
				avgErrorScale[sIndex] += errorEstimate[dIndex][sIndex];
			}
			sumWeight += sumWeightScale[sIndex];

			float temp;
			if (sumWeightScale[sIndex] != 0.f)
			{
				temp = 1.0 / sumWeightScale[sIndex];
				for (int dIndex = 0; dIndex<8; dIndex++)
					weight[dIndex][sIndex] *= temp;
			}

			//std error;
			avgErrorScale[sIndex] /= 8.0;
			for (int dIndex = 0; dIndex<8; dIndex++)
				stdErrorScale[sIndex] = stdErrorScale[sIndex] + (errorEstimate[dIndex][sIndex] - avgErrorScale[sIndex])
				*(errorEstimate[dIndex][sIndex] - avgErrorScale[sIndex]);
			stdErrorScale[sIndex] /= 8.0;
			stdErrorScale[sIndex] = sqrt(stdErrorScale[sIndex]);

			//alph;
			for (int dIndex = 0; dIndex<8; dIndex++)
				alphScale[sIndex] = alphScale[sIndex] + weight[dIndex][sIndex] * curvature[dIndex][sIndex];
			alphScale[sIndex] *= stdErrorScale[sIndex];
			alphScale[sIndex] = 1.0 / (1.0 + pow(E, (double)-alphScale[sIndex]));

			float curScale = 0.f;
			for (int dIndex = 0; dIndex < 8; dIndex++)
				curScale = curScale + curvature[dIndex][sIndex] *feature[dIndex][sIndex] *weight[dIndex][sIndex];// ;
			alphScale[0] = 1;
			alphScale[1] = 5;
			reliefFeature[vIndex] = reliefFeature[vIndex] + alphScale[sIndex] * sumWeightScale[sIndex] * curScale;

		}

		if (sumWeight != 0.f)
			reliefFeature[vIndex] = reliefFeature[vIndex] / sumWeight;

		*/

		/*
		for (int sIndex = 0; sIndex<scaleNum;sIndex++)
		{

		for (int dIndex = 0; dIndex<8; dIndex++)
		{
		float tempNormal[3];
		float direc[3];

		float a = sqrt(featureSlop[dIndex][sIndex] * featureSlop[dIndex][sIndex] + 1);
		a = 1.0/a;

		for (int demention = 0; demention<3; demention++)
		{
		tempNormal[demention] = featureSlop[dIndex][sIndex] * directions[dIndex][demention] + baseNormal[vIndex][demention];
		tempNormal[demention] *=a;
		}



		for (int demention = 0; demention<3; demention++)
		normalEnhance[vIndex][demention] = normalEnhance[vIndex][demention] + alphScale[sIndex] * sumWeightScale[sIndex] * weight[dIndex][sIndex] * tempNormal[demention];

		}

		}

		NORMALIZE(normalEnhance[vIndex]);


		if (sumWeight == 0)
		for (int demention = 0; demention<3; demention++)
		normalEnhance[vIndex][demention] = vertexNormal[vIndex][demention];
		*/

//	}


	/*
	float min, max;

	min = 9999;
	max = -9999;
	for (int vIndex = 0; vIndex<vertexNum; vIndex++)
	{
	if (min > cur[vIndex])
	min = cur[vIndex];
	if (max < cur[vIndex])
	max = cur[vIndex];
	}
	max = 1.0f/(max-min);
	for (int vIndex = 0; vIndex<vertexNum; vIndex++)
	{
	cur[vIndex] -= min;
	cur[vIndex] *= max;
	}
	*/

	//	delete[] baseNormal;
}

void TriMesh::getVertexSamplePolyHeight(int vertexIndex, int scale, float angle)
{
	need_curvatures();
	need_halfEdge();


	::std::vector<int> &blackFaceIndex = halfEdge.blackFaceIndex;
	::std::vector<float> &sampleHeight = halfEdge.sampleHeight;


	blackFaceIndex.resize(scale);
	sampleHeight.resize(scale * 2);
	
	vec3 oriNormal = pdir1[vertexIndex].cross(normals[vertexIndex]);
	normalize(oriNormal);
	vec3 planeNormal = oriNormal * cos(angle) + pdir1[vertexIndex] * sin(angle);
	normalize(planeNormal);
	
	Eigen::Matrix4f scaleMatrix;
	need_averageEdgeLength();
	float sampleStep = averageEdgeLength[0] * 0.5;
	int halfSampleNum = (int)(((float)scale - 1.f)*0.5f);
	calScaleMatrix(scaleMatrix, scale, sampleStep);

	float *sampleErrorEstimate = new float[scale];
	findSamplesHeight(planeNormal, scale, vertexIndex, sampleErrorEstimate);
	delete[] sampleErrorEstimate;



	float parmeter[4];
	polyFitSamples(scaleMatrix, scale, parmeter, &sampleHeight[0], sampleStep);

	int all = scale * 2;
	float steps = -1.0*halfSampleNum*sampleStep;


	float maxSampleHeight;
	float minSampleHeight;

	float gamma;

	maxSampleHeight = FLT_MIN;
	minSampleHeight = FLT_MAX;


	for (int i = scale; i<all; i++)
	{

		sampleHeight[i] = parmeter[0] + parmeter[1] * steps + parmeter[2] * steps*steps + parmeter[3] * steps*steps*steps;
		steps += sampleStep;

		if (maxSampleHeight < sampleHeight[i])
			maxSampleHeight = sampleHeight[i];
		if (minSampleHeight > sampleHeight[i])
			minSampleHeight = sampleHeight[i];

	}

	for (int i = 0; i<scale; i++)
	{
		if (maxSampleHeight < sampleHeight[i])
			maxSampleHeight = sampleHeight[i];
		if (minSampleHeight > sampleHeight[i])
			minSampleHeight = sampleHeight[i];
	}


	gamma = 1.0 / (maxSampleHeight - minSampleHeight);
	for (int i = 0; i<all; i++)
	{
		sampleHeight[i] -= minSampleHeight;
		sampleHeight[i] *= gamma;
	}




}

void TriMesh::getVertexSampleNormal(int vertexIndex, int scale)
{


	need_curvatures();
	need_halfEdge();

	vec3 planeNormal = pdir1[vertexIndex].cross(normals[vertexIndex]);

	std::vector<vec2> &leftScrPosition = halfEdge.leftScrPosition;
	std::vector<vec2> &rightScrPosition = halfEdge.rightScrPosition;
	std::vector<vec2> &leftScrNormal = halfEdge.leftScrNormal;
	std::vector<vec2> &rightScrNormal = halfEdge.rightScrNormal;

	::std::vector<::std::vector<vec3>> &leftSampleNormal = halfEdge.leftSampleNormal;
	::std::vector<::std::vector<vec3>> &rightSampleNormal = halfEdge.rightSampleNormal;
	::std::vector<::std::vector<vec3>> &leftSamplePosition = halfEdge.leftSamplePosition;
	::std::vector<::std::vector<vec3>> &rightSamplePosition = halfEdge.rightSamplePosition;

	leftScrPosition.clear();
	leftScrNormal.clear();
	rightScrPosition.clear();
	rightScrNormal.clear();

	findSamplesNormal(planeNormal, scale, vertexIndex);


	vec2 max(FLT_MIN, FLT_MIN );
	vec2 min( FLT_MAX, FLT_MAX );
	vec2 gamma;

	leftScrPosition.push_back(vec2(0.f, 0.f));
	leftScrNormal.push_back(vec2(0.f, 1.f));

	for (unsigned int i = 0; i < leftSamplePosition[vertexIndex].size(); i++)
	{
		vec3 coord3 = leftSamplePosition[vertexIndex][i] - vertices[vertexIndex];
		leftScrPosition.push_back(vec2(coord3.dot(pdir1[vertexIndex]), coord3.dot(normals[vertexIndex])));
		vec3 &normal3 = leftSampleNormal[vertexIndex][i];
		leftScrNormal.push_back(vec2(normal3.dot(pdir1[vertexIndex]), normal3.dot(normals[vertexIndex])));
	}
	for (unsigned int i = 0; i < leftScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			if (max[d] < leftScrPosition[i][d])
				max[d] = leftScrPosition[i][d];
			if (min[d] > leftScrPosition[i][d])
				min[d] = leftScrPosition[i][d];
		}
	}



	rightScrPosition.push_back(vec2(0.f, 0.f));
	rightScrNormal.push_back(vec2(0.f, 1.f));

	for (unsigned int i = 0; i < rightSamplePosition[vertexIndex].size(); i++)
	{
		vec3 coord3 = rightSamplePosition[vertexIndex][i] - vertices[vertexIndex];
		rightScrPosition.push_back(vec2(coord3.dot(pdir1[vertexIndex]), coord3.dot(normals[vertexIndex])));
		vec3 &normal3 = rightSampleNormal[vertexIndex][i];
		rightScrNormal.push_back(vec2(normal3.dot(pdir1[vertexIndex]), normal3.dot(normals[vertexIndex])));
	}

	for (unsigned int i = 0; i < rightScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			if (max[d] < rightScrPosition[i][d])
				max[d] = rightScrPosition[i][d];
			if (min[d] > rightScrPosition[i][d])
				min[d] = rightScrPosition[i][d];
		}
	}


	gamma[0] = 1.0 / (max[0] - min[0]);
	gamma[1] = 0.8 / (max[1] - min[1]);


	for (unsigned int i = 0; i < rightScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			rightScrPosition[i][d] -= min[d];
			rightScrPosition[i][d] *= gamma[d];
		}
	}

	for (unsigned int i = 0; i < rightScrNormal.size(); i++)
	{
		rightScrNormal[i] *= 0.2;
		rightScrNormal[i] += rightScrPosition[i];
	}



	for (int i = 0; i < leftScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			leftScrPosition[i][d] -= min[d];
			leftScrPosition[i][d] *= gamma[d];
		}
	}

	for (int i = 0; i < leftScrNormal.size(); i++)
	{
		leftScrNormal[i] *= 0.2;
		leftScrNormal[i] += leftScrPosition[i];
	}



}

void TriMesh::getRidgeSampleNormal(int ridgeIndex, int scale, float curvThreshold, float planeThreshold)
{


	need_curvatures();
	need_halfEdge();

	
	std::vector<vec2> &leftScrPosition = halfEdge.leftScrPosition;
	std::vector<vec2> &rightScrPosition = halfEdge.rightScrPosition;
	std::vector<vec2> &leftScrNormal = halfEdge.leftScrNormal;
	std::vector<vec2> &rightScrNormal = halfEdge.rightScrNormal;

	::std::vector<::std::vector<vec3>> &leftSampleNormal = halfEdge.leftSampleNormal;
	::std::vector<::std::vector<vec3>> &rightSampleNormal = halfEdge.rightSampleNormal;
	::std::vector<::std::vector<vec3>> &leftSamplePosition = halfEdge.leftSamplePosition;
	::std::vector<::std::vector<vec3>> &rightSamplePosition = halfEdge.rightSamplePosition;

	leftScrPosition.clear();
	leftScrNormal.clear();
	rightScrPosition.clear();
	rightScrNormal.clear();

	findSamplesNormalFromEdge(scale, ridgeIndex, vec3(1.f));


	vec3 samplePosition, sampleNormal, pCurDir, directionDir;
	float curveture;
	int temp;
	bool isVertex;
	ridgeLines.getPointInfo(ridgeIndex, samplePosition, sampleNormal, curveture, pCurDir, temp, isVertex);
	

	directionDir = pCurDir;

	

	::std::vector<::std::vector<bool>> &leftIsSamplePlane = halfEdge.leftIsSampleCurv;
	::std::vector<::std::vector<bool>> &rightIsSamplePlane = halfEdge.rightIsSampleCurv;

	analysisCurvPlan(samplePosition, sampleNormal, directionDir, curveture, leftSamplePosition[ridgeIndex], leftSampleNormal[ridgeIndex], halfEdge.leftSampleCurv[ridgeIndex], leftIsSamplePlane[ridgeIndex], curvThreshold, planeThreshold);
	analysisCurvPlan(samplePosition, sampleNormal, directionDir, curveture, rightSamplePosition[ridgeIndex], rightSampleNormal[ridgeIndex], halfEdge.rightSampleCurv[ridgeIndex], rightIsSamplePlane[ridgeIndex], curvThreshold, planeThreshold);



	

	vec2 max =vec2(FLT_MIN, FLT_MIN);
	vec2 min =vec2(FLT_MAX, FLT_MAX);
	vec2 gamma;


	for (unsigned int i = 0; i < leftSamplePosition[ridgeIndex].size(); i++)
	{
		vec3 coord3 = leftSamplePosition[ridgeIndex][i] - samplePosition;
		leftScrPosition.push_back(vec2(coord3.dot(pCurDir), coord3.dot(sampleNormal)));
		vec3 &normal3 = halfEdge.leftSampleNormal[ridgeIndex][i];
		leftScrNormal.push_back(vec2(normal3.dot(pCurDir), normal3.dot(sampleNormal)));
	}
	for (unsigned int i = 0; i < leftScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			if (max[d] < leftScrPosition[i][d])
				max[d] = leftScrPosition[i][d];
			if (min[d] > leftScrPosition[i][d])
				min[d] = leftScrPosition[i][d];
		}
	}


	for (unsigned int i = 0; i < rightSamplePosition[ridgeIndex].size(); i++)
	{
		vec3 coord3 = rightSamplePosition[ridgeIndex][i] - samplePosition;
		rightScrPosition.push_back(vec2(coord3.dot(pCurDir), coord3.dot(sampleNormal)));
		vec3 &normal3 = rightSampleNormal[ridgeIndex][i];
		rightScrNormal.push_back(vec2(normal3.dot(pCurDir), normal3.dot(sampleNormal)));

	}
	for (unsigned int i = 0; i < rightScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			if (max[d] < rightScrPosition[i][d])
				max[d] = rightScrPosition[i][d];
			if (min[d] > rightScrPosition[i][d])
				min[d] = rightScrPosition[i][d];
		}
	}


	gamma[0] = 1.0 / (max[0] - min[0]);
	gamma[1] = 0.8 / (max[1] - min[1]);


	for (int i = 0; i < (int)rightScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			rightScrPosition[i][d] -= min[d];
			rightScrPosition[i][d] *= gamma[d];
		}
	}

	for (int i = 0; i < (int)rightScrNormal.size(); i++)
	{
		rightScrNormal[i] *= 0.2f;
		rightScrNormal[i] += rightScrPosition[i];
	}



	for (int i = 0; i < (int)leftScrPosition.size(); i++)
	{
		for (int d = 0; d < 2; d++)
		{
			leftScrPosition[i][d] -= min[d];
			leftScrPosition[i][d] *= gamma[d];
		}
	}

	for (int i = 0; i < (int)leftScrNormal.size(); i++)
	{
		leftScrNormal[i] *= 0.2f;
		leftScrNormal[i] += leftScrPosition[i];
	}



}

void TriMesh::setColorRidgeCurvPlane(int scale, float scaleCur, float ridgeThreshold, float curvThreshold, float planeThreshold)
{

	need_RidgeLines(true);
	need_RidgeLines(false);
	need_curvatures();


	int sampleNum = ridgeLines.ridgeLinesPointNormal.size() + ridgeLines.valleyLinesPointNormal.size();

	if (halfEdge.leftFaceIndex.empty())
		halfEdge.leftFaceIndex.resize(sampleNum);
	if (halfEdge.rightFaceIndex.empty())
		halfEdge.rightFaceIndex.resize(sampleNum);
	if (halfEdge.leftSampleNormal.empty())
		halfEdge.leftSampleNormal.resize(sampleNum);
	if (halfEdge.rightSampleNormal.empty())
		halfEdge.rightSampleNormal.resize(sampleNum);

	if (halfEdge.leftIsSampleCurv.empty())
		halfEdge.leftIsSampleCurv.resize(sampleNum);
	if (halfEdge.rightIsSampleCurv.empty())
		halfEdge.rightIsSampleCurv.resize(sampleNum);

	if (halfEdge.leftSamplePosition.empty())
		halfEdge.leftSamplePosition.resize(sampleNum);
	if (halfEdge.rightSamplePosition.empty())
		halfEdge.rightSamplePosition.resize(sampleNum);


	if (halfEdge.leftSampleCurv.empty())
		halfEdge.leftSampleCurv.resize(sampleNum);
	if (halfEdge.rightSampleCurv.empty())
		halfEdge.rightSampleCurv.resize(sampleNum);
	


	if (halfEdge.hasRedfile.empty())
	{
		halfEdge.oldScale = 0;
		halfEdge.oldRidgeThreshold = FLT_MAX;
		halfEdge.hasRedfile.push_back(true);
		FILE *fin;
		fopen_s(&fin, dFName.getRidgeSampleDN(), "rb");
		if (fin)
		{
			cout << "Reading ridgeSamples.....";


			fread(&halfEdge.oldScale, sizeof(int), 1, fin);
			fread(&halfEdge.oldRidgeThreshold, sizeof(float), 1, fin);

			if (halfEdge.oldScale >= scale && ridgeThreshold >= halfEdge.oldRidgeThreshold)
			{

				vector<int> leftSampleSize;
				vector<int> rightSampleSize;
				leftSampleSize.resize(sampleNum);
				rightSampleSize.resize(sampleNum);
				fread(&leftSampleSize[0], sizeof(int), sampleNum, fin);
				fread(&rightSampleSize[0], sizeof(int), sampleNum, fin);

				for (int i = 0; i < sampleNum; i++)
				{
					int tempSize = leftSampleSize[i];

					if (tempSize != 0)
					{
						halfEdge.leftFaceIndex[i].resize(tempSize);
						halfEdge.leftSampleNormal[i].resize(tempSize);
						halfEdge.leftSamplePosition[i].resize(tempSize);
						halfEdge.leftSampleCurv[i].resize(tempSize);
						fread(&halfEdge.leftFaceIndex[i][0], sizeof(int), tempSize, fin);
						fread(halfEdge.leftSampleNormal[i][0], sizeof(float), tempSize * 3, fin);
						fread(halfEdge.leftSamplePosition[i][0], sizeof(float), tempSize * 3, fin);
						fread(&halfEdge.leftSampleCurv[i][0], sizeof(float), tempSize, fin);
					}


					tempSize = rightSampleSize[i];
					if (tempSize != 0)
					{
						halfEdge.rightFaceIndex[i].resize(tempSize);
						halfEdge.rightSampleNormal[i].resize(tempSize);
						halfEdge.rightSamplePosition[i].resize(tempSize);
						halfEdge.rightSampleCurv[i].resize(tempSize);
						fread(&halfEdge.rightFaceIndex[i][0], sizeof(int), tempSize, fin);
						fread(halfEdge.rightSampleNormal[i][0], sizeof(float), tempSize * 3, fin);
						fread(halfEdge.rightSamplePosition[i][0], sizeof(float), tempSize * 3, fin);
						fread(&halfEdge.rightSampleCurv[i][0], sizeof(float), tempSize, fin);
					}


				}

				cout << "Done" << endl;
				fclose(fin);
				return;
			}
			else
			{


			}
			cout << "Done" << endl;
			fclose(fin);
		}
	}

	int totalSampleNum = ridgeLines.ridgeLinesPointNormal.size() + ridgeLines.valleyLinesPointNormal.size();

	//calculate all
	if (halfEdge.oldScale < scale)
	{
		vec3 lineDir;
		float tempTh = ridgeThreshold * averageCurv1;

		for (int i = 0; i < totalSampleNum; i++)
		{
			lineDir = ridgeLines.getLineDirection(i);


			if (abs(ridgeLines.getCurv(i)) >= tempTh)
				findSamplesNormalFromEdge(scale, i, lineDir);
		}
	}
	else if (halfEdge.oldRidgeThreshold > ridgeThreshold)
	{
		vec3 lineDir;
		float tempTh = ridgeThreshold * averageCurv1;

		for (int i = 0; i < totalSampleNum; i++)
		{
			lineDir = ridgeLines.getLineDirection(i);

			if (abs(ridgeLines.getCurv(i)) >= tempTh)
				findSamplesNormalFromEdge(halfEdge.oldScale, i, lineDir);
		}
	}

	



	vec3 samplePosition, sampleNormal, planeNormal, directionDir;
	int edgeIndex;
	float curveture;
	bool isVertex;

	halfEdge.samplePosition.clear();
	float tempTh = ridgeThreshold * averageCurv1;
	for (int i = 0; i < totalSampleNum; i++)
	{
		ridgeLines.getPointInfo(i, samplePosition, sampleNormal, curveture, directionDir, edgeIndex, isVertex);

		if (abs(curveture) >= tempTh)
		{
			halfEdge.samplePosition.push_back(ridgeLines.getPointPosition(i));


			analysisCurvPlan(samplePosition, sampleNormal, directionDir, curveture,
				halfEdge.leftSamplePosition[i], halfEdge.leftSampleNormal[i], halfEdge.leftSampleCurv[i], halfEdge.leftIsSampleCurv[i], curvThreshold, planeThreshold);
			analysisCurvPlan(samplePosition, sampleNormal, directionDir, curveture,
				halfEdge.rightSamplePosition[i], halfEdge.rightSampleNormal[i], halfEdge.rightSampleCurv[i], halfEdge.rightIsSampleCurv[i], curvThreshold, planeThreshold);

		}
	}

	
	

	if (halfEdge.oldRidgeThreshold > ridgeThreshold || halfEdge.oldScale < scale)
	{
		if (halfEdge.oldRidgeThreshold > ridgeThreshold)
			halfEdge.oldRidgeThreshold = ridgeThreshold;
		if (halfEdge.oldScale < scale)
			halfEdge.oldScale = scale;



		FILE *fout;
		fopen_s(&fout, dFName.getRidgeSampleDN(), "wb");

		fwrite(&halfEdge.oldScale, sizeof(int), 1, fout);
		fwrite(&halfEdge.oldRidgeThreshold, sizeof(float), 1, fout);


		vector<int> leftSampleSize;
		vector<int> rightSampleSize;
		leftSampleSize.resize(sampleNum);
		rightSampleSize.resize(sampleNum);

		for (int i = 0; i < sampleNum; i++)
		{
			leftSampleSize[i] = halfEdge.leftFaceIndex[i].size();
			rightSampleSize[i] = halfEdge.rightFaceIndex[i].size();

		}

		fwrite(&leftSampleSize[0], sizeof(int), sampleNum, fout);
		fwrite(&rightSampleSize[0], sizeof(int), sampleNum, fout);

		for (int i = 0; i < sampleNum; i++)
		{

			int tempSize = leftSampleSize[i];

			if (tempSize != 0)
			{
				fwrite(&halfEdge.leftFaceIndex[i][0], sizeof(int), tempSize, fout);
				fwrite(halfEdge.leftSampleNormal[i][0], sizeof(float), tempSize * 3, fout);
				fwrite(halfEdge.leftSamplePosition[i][0], sizeof(float), tempSize * 3, fout);
				fwrite(&halfEdge.leftSampleCurv[i][0], sizeof(float), tempSize, fout);


			}

			tempSize = rightSampleSize[i];

			if (tempSize != 0)
			{
				fwrite(&halfEdge.rightFaceIndex[i][0], sizeof(int), tempSize, fout);
				fwrite(halfEdge.rightSampleNormal[i][0], sizeof(float), tempSize * 3, fout);
				fwrite(halfEdge.rightSamplePosition[i][0], sizeof(float), tempSize * 3, fout);
				fwrite(&halfEdge.rightSampleCurv[i][0], sizeof(float), tempSize, fout);

			}

		}

		cout << "Done" << endl;

		fclose(fout);



	}


	



	if (colors.size() != vertices.size())
		colors.resize(vertices.size());


	setColorCurvatureGray(scaleCur);




	Color green(0.f, 1.f, 0.f);
	Color red(1.f, 0.f, 0.f);
	Color blue(0.f, 0.f, 1.f);



	for (unsigned int i = 0; i < halfEdge.rightFaceIndex.size(); i++)
	for (unsigned int fIndex = 0; fIndex < halfEdge.rightIsSampleCurv[i].size() && fIndex<scale; fIndex++)
	{
		if (halfEdge.rightIsSampleCurv[i][fIndex])
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.rightFaceIndex[i][fIndex]].v[vIndex]] = red;
		else
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.rightFaceIndex[i][fIndex]].v[vIndex]] += blue;
	}


	for (unsigned int i = 0; i < halfEdge.leftFaceIndex.size(); i++)
	for (int fIndex = 0; fIndex < (int)halfEdge.leftIsSampleCurv[i].size() && fIndex < scale; fIndex++)
	{
		if (halfEdge.leftIsSampleCurv[i][fIndex])
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.leftFaceIndex[i][fIndex]].v[vIndex]] = red;
		else
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.leftFaceIndex[i][fIndex]].v[vIndex]] += green;


	}
	return;

}

int TriMesh::setColorRidgeVertexCurvPlane(int vertexIndex, int scale, float scaleCur, float ridgeThreshold, float curvThreshold, float planeThreshold)
{

	need_RidgeLines(true);
	need_RidgeLines(false);
	need_curvatures();
	clamp(vertexIndex, 0, 100000);

	float tempTh = ridgeThreshold * averageCurv1;

	int sampleIndex = ridgeLines.changeToVaildIndex(vertexIndex, tempTh);


	int sampleNum = ridgeLines.ridgeLinesPointNormal.size() + ridgeLines.ridgeLinesPointNormal.size();


	if (halfEdge.leftFaceIndex.empty())
		halfEdge.leftFaceIndex.resize(sampleNum);
	if (halfEdge.rightFaceIndex.empty())
		halfEdge.rightFaceIndex.resize(sampleNum);
	if (halfEdge.leftSampleNormal.empty())
		halfEdge.leftSampleNormal.resize(sampleNum);
	if (halfEdge.rightSampleNormal.empty())
		halfEdge.rightSampleNormal.resize(sampleNum);
	if (halfEdge.leftIsSampleCurv.empty())
		halfEdge.leftIsSampleCurv.resize(sampleNum);
	if (halfEdge.rightIsSampleCurv.empty())
		halfEdge.rightIsSampleCurv.resize(sampleNum);

	if (halfEdge.leftSamplePosition.empty())
		halfEdge.leftSamplePosition.resize(sampleNum);
	if (halfEdge.rightSamplePosition.empty())
		halfEdge.rightSamplePosition.resize(sampleNum);

	if (halfEdge.leftSampleCurv.empty())
		halfEdge.leftSampleCurv.resize(sampleNum);
	if (halfEdge.rightSampleCurv.empty())
		halfEdge.rightSampleCurv.resize(sampleNum);



	if (colors.size() != vertices.size())
		colors.resize(vertices.size());



	setColorCurvatureGray(scaleCur);



	getRidgeSampleNormal(sampleIndex, scale, curvThreshold, planeThreshold);

	Color black(0.f);
	Color red(1.f, 0.f, 0.f);
	Color blue(0.f, 0.f, 1.f);


	for (unsigned int fIndex = 0; fIndex < halfEdge.leftFaceIndex[sampleIndex].size(); fIndex++)
	{
		if (fIndex >= halfEdge.leftIsSampleCurv[sampleIndex].size() || halfEdge.leftIsSampleCurv[sampleIndex][fIndex])
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.leftFaceIndex[sampleIndex][fIndex]].v[vIndex]] = red;
		else
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.leftFaceIndex[sampleIndex][fIndex]].v[vIndex]] = black;
	}


	for (unsigned int fIndex = 0; fIndex < halfEdge.rightFaceIndex[sampleIndex].size(); fIndex++)
	{
		if (fIndex >= halfEdge.rightIsSampleCurv[sampleIndex].size() || halfEdge.rightIsSampleCurv[sampleIndex][fIndex])
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.rightFaceIndex[sampleIndex][fIndex]].v[vIndex]] = blue;
		else
		for (int vIndex = 0; vIndex < 3; vIndex++)
			colors[faces[halfEdge.rightFaceIndex[sampleIndex][fIndex]].v[vIndex]] = black;
	}


	return sampleIndex;
}

bool TriMesh::getPlaneNormal()
{

	if (!halfEdge.leftPlaneNormal.empty())
		return true;

	int scale;
	float ridgeThreshold;
	float curvThreshold;
	float planeThreshold;


	FILE *fin;
	fopen_s(&fin, dFName.getPlaneNormalsDN(), "rb");
	if (fin)
	{

		need_RidgeLines(true);
		need_RidgeLines(false);
		int sampleNum;


		fread(&scale, sizeof(int), 1, fin);
		fread(&ridgeThreshold, sizeof(int), 1, fin);
		fread(&curvThreshold, sizeof(int), 1, fin);
		fread(&planeThreshold, sizeof(int), 1, fin);
		fread(&sampleNum, sizeof(int), 1, fin);


		halfEdge.leftPlaneNormal.resize(sampleNum);
		halfEdge.rightPlaneNormal.resize(sampleNum);
		halfEdge.samplePosition.resize(sampleNum);

		fread(&halfEdge.leftPlaneNormal[0], sizeof(int)*3, sampleNum, fin);
		fread(&halfEdge.rightPlaneNormal[0], sizeof(int)* 3, sampleNum, fin);
		fread(&halfEdge.samplePosition[0], sizeof(int)* 3, sampleNum, fin);

		fclose(fin);

		return true;

	}
	else
	{

		cout << "No file" << endl;
		fclose(fin);
		return false;
	}


}

void TriMesh::writePlaneNormal(int scale, float ridgeThreshold, float curvThreshold, float planeThreshold)
{
	if (halfEdge.leftSampleNormal.empty())
		return;

	need_RidgeLines(true);
	need_RidgeLines(false);
	need_curvatures();

	int sampleNum = ridgeLines.ridgeLinesPointNormal.size() + ridgeLines.valleyLinesPointNormal.size();

	halfEdge.leftPlaneNormal.clear();
	halfEdge.rightPlaneNormal.clear();
	halfEdge.samplePosition.clear();


	halfEdge.leftPlaneNormal.reserve(sampleNum);
	halfEdge.rightPlaneNormal.reserve(sampleNum);
	halfEdge.samplePosition.reserve(sampleNum);


	float tempTh = ridgeThreshold * averageCurv1;

	//all after *********************************************************//
	bool *isFlat = new bool[halfEdge.leftSampleNormal.size()];

	//all after *********************************************************//

	for (unsigned int rIndex = 0; rIndex < halfEdge.leftSampleNormal.size(); rIndex++)
	{
		isFlat[rIndex] = false;

		if (!halfEdge.leftSampleNormal[rIndex].empty() && !halfEdge.rightSampleNormal[rIndex].empty() && abs(ridgeLines.getCurv(rIndex)) >= tempTh)
		{
			vec3 tempLeftNormal(0.f);
			for (unsigned int sIndex = 0; sIndex < halfEdge.leftIsSampleCurv[rIndex].size(); sIndex++)
			{
				if (halfEdge.leftIsSampleCurv[rIndex][sIndex] == false)
					tempLeftNormal += halfEdge.leftSampleNormal[rIndex][sIndex];
			}
			normalize(tempLeftNormal);

			vec3 tempRightNormal = vec3(0.f);
			for (unsigned int sIndex = 0; sIndex < halfEdge.rightIsSampleCurv[rIndex].size(); sIndex++)
			{
				if (halfEdge.rightIsSampleCurv[rIndex][sIndex] == false)
					tempRightNormal += halfEdge.rightSampleNormal[rIndex][sIndex];
			}
			normalize(tempRightNormal);

			if(tempLeftNormal.dot(tempRightNormal) > 0.707)
			{
				isFlat[rIndex]  = true;
				continue;
			}

			halfEdge.leftPlaneNormal.push_back(tempLeftNormal);
			halfEdge.rightPlaneNormal.push_back(tempRightNormal);
			halfEdge.samplePosition.push_back(ridgeLines.getPointPosition(rIndex));
		}

	}





	FILE *fout;
	fopen_s(&fout, dFName.getPlaneNormalsDN(), "wb");


	sampleNum = halfEdge.leftPlaneNormal.size();

	fwrite(&scale, sizeof(int), 1, fout);
	fwrite(&ridgeThreshold, sizeof(float), 1, fout);
	fwrite(&curvThreshold, sizeof(float), 1, fout);
	fwrite(&planeThreshold, sizeof(float), 1, fout);
	fwrite(&sampleNum, sizeof(int), 1, fout);


	fwrite(&halfEdge.leftPlaneNormal[0], sizeof(float)*3, sampleNum, fout);
	fwrite(&halfEdge.rightPlaneNormal[0], sizeof(float)* 3, sampleNum, fout);
	fwrite(&halfEdge.samplePosition[0], sizeof(float)* 3, sampleNum, fout);

	fclose(fout);







	int *lineIndex = new int[halfEdge.leftPlaneNormal.size()];
	
	int* starEnd = new int[halfEdge.leftPlaneNormal.size()];

	int allIndex = 0;
	int localIndex = 0;
	int totalIndex = 0;
	int lineIndexStart = 1;
	bool isAvilable = false;
	int startIndex = 0;

	for (int i = 0; i < (int)ridgeLines.ridgeLinesPoint.size(); i++)
	{
		int start = localIndex + 1;
		isAvilable = false;
		for (int sampleIndex = 0; sampleIndex < (int)ridgeLines.ridgeLinesPoint[i].size(); sampleIndex++)
		{

			if (abs(ridgeLines.ridgeLinesPointPcurv[totalIndex]) > tempTh && isFlat[allIndex] == false)
			{
				lineIndex[localIndex++] = lineIndexStart;
				isAvilable = true;
			}
			totalIndex++;
			allIndex++;
		}
		if (isAvilable)
		{
			lineIndexStart++;
			starEnd[startIndex++] = start;
			starEnd[startIndex++] = localIndex;
		}

	}

	isAvilable = false;
	totalIndex = 0;
	for (int i = 0; i < (int)ridgeLines.valleyLinesPoint.size(); i++)
	{
		int start = localIndex + 1;
		isAvilable = false;
		for (int sampleIndex = 0; sampleIndex < (int)ridgeLines.valleyLinesPoint[i].size(); sampleIndex++)
		{
			if (abs(ridgeLines.valleyLinesPointPcurv[totalIndex]) > tempTh && isFlat[allIndex] == false)
			{
				lineIndex[localIndex++] = lineIndexStart;
				isAvilable = true;
			}
			
			allIndex++;
			totalIndex++;
		}
		if (isAvilable)
		{
			lineIndexStart++;
			starEnd[startIndex++] = start;
			starEnd[startIndex++] = localIndex;
		}
	}

	if (localIndex != halfEdge.leftPlaneNormal.size())
		cout << localIndex << endl;


	ofstream WaveCoeff_of(dFName.getPlaNorLinesDN());



	for (int i = 0; i < localIndex; i++)
	{
		WaveCoeff_of << lineIndex[i] << " ";
	}

	for (int i = 0; i < startIndex; i++)
	{
		WaveCoeff_of << starEnd[i++] << " " << starEnd[i] << endl;
	}

	WaveCoeff_of.close();

//	delete[] lineIndex;
}






/*

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
*/

}; // namespace trimesh
