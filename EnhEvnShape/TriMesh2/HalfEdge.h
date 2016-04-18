//半边数据结构

#ifndef HALFEDGE_H
#define HALFEDGE_H


#include "math.h"
#include "vec.h"

#include <vector>

namespace trimesh {


struct HE_vert;
struct HE_edge;
struct HE_face;


struct HE_vert
{

	HE_edge* edge;  // one of the half-edges emantating from the vertex;
	int index;
	int faceNum;

	void calculateFaceNum();
	int findEdge(int endVertex);
};

struct HE_edge
{

	HE_vert* vert;   // vertex at the end of the half-edge
	HE_edge* pair;   // oppositely oriented adjacent half-edge 
	HE_face* face;   // face the half-edge borders
	HE_edge* next;   // next half-edge around the face
	int index;

	int getStartVertIndex(){ return vert->index; }
	int getEndVertIndex(){ return pair->vert->index; }
	int getPairIndex(){ return pair->index; }

	HE_face *getLeftFace(){return face;}
	HE_face *getRightFace(){return pair->face;}
	bool isClockwise(){return vert>pair->vert;}
};

struct HE_face
{
	HE_edge* edge;  // one of the half-edges bordering the face
	int index;

	int getIndex(){return index;}
	void getVertices(HE_vert *[]);
	void getVertexIndexes(int *);
	void getBoundEdges(HE_edge *[]);
	void getBoundEdgeIndexes(int *);
	void getAdFaces(HE_face *[]);
	int findAdFace(HE_face *);
	int findAdFace(int);
	void findVertexAdEdgeIndex(int [][2]);
	int findVertex(int);
};



struct HalfEdge
{

	::std::vector<HE_vert> he_vert;
	::std::vector<HE_edge> he_edge;
	::std::vector<HE_face> he_face;

	::std::vector<float> vertexHeight;
	::std::vector<float> sampleHeight;
	::std::vector<float> reliefFeature;
	::std::vector<int> blackFaceIndex;



	::std::vector<float> singleReliefFeature[3];
	::std::vector<float> singleReliefFeatureWeight[3];
	int scaleValue[3];
	


	int blackFaceNum;


	float maxHeight;
	float minHeight;

	float maxReliefFeature;
	float minReliefFeature;

	int oldScale;
	float oldRidgeThreshold;

	::std::vector<::std::vector<int>> leftFaceIndex;
	::std::vector<::std::vector<int>> rightFaceIndex;
	::std::vector<::std::vector<vec3>> leftSampleNormal;
	::std::vector<::std::vector<vec3>> rightSampleNormal;
	::std::vector<::std::vector<bool>> leftIsSampleCurv;
	::std::vector<::std::vector<bool>> rightIsSampleCurv;
	::std::vector<::std::vector<vec3>> leftSamplePosition;
	::std::vector<::std::vector<vec3>> rightSamplePosition;
	::std::vector<::std::vector<float>> leftSampleCurv;
	::std::vector<::std::vector<float>> rightSampleCurv;
	std::vector<vec2> leftScrPosition;
	std::vector<vec2> rightScrPosition;
	std::vector<vec2> leftScrNormal;
	std::vector<vec2> rightScrNormal;

	std::vector<bool> hasRedfile;




	std::vector<vec3> leftPlaneNormal;
	std::vector<vec3> rightPlaneNormal;
	std::vector<vec3> samplePosition;



};


}; // namespace trimesh

#endif