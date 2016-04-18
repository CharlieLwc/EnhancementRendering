//半边数据结构



#include "stdafx.h"

#include "math.h"


struct HE_vert;
struct HE_edge;
struct HE_face;


struct HE_vert
{
	float x;
	float y;
	float z;

	HE_edge* edge;  // one of the half-edges emantating from the vertex;
	int index;
	int faceNum;

	int getIndex(){return index;}
	void calculateFaceNum();
};


static inline void VECTOR(double *v, HE_vert *v1, HE_vert *v2){
	v[0] = (double)v2->x-(double)v1->x;
	v[1] = (double)v2->y-(double)v1->y;
	v[2] = (double)v2->z-(double)v1->z;
}

static inline void VECTOR(double *v, float v1[3], HE_vert *v2){
	v[0] = (double)v2->x-(double)v1[0];
	v[1] = (double)v2->y-(double)v1[1];
	v[2] = (double)v2->z-(double)v1[2];
}

static inline void VECTOR(float *v, float v1[3], HE_vert *v2){
	v[0] = v2->x-v1[0];
	v[1] = v2->y-v1[1];
	v[2] = v2->z-v1[2];
}

static inline void VECTOR(double *v, float *v1, float *v2){
	v[0] = (double)v2[0]-(double)v1[0];
	v[1] = (double)v2[1]-(double)v1[1];
	v[2] = (double)v2[2]-(double)v1[2];
}

static inline void VECTOR(float *v, float *v1, float *v2){
	v[0] = v2[0]-v1[0];
	v[1] = v2[1]-v1[1];
	v[2] = v2[2]-v1[2];
}

static inline void VECTOR(float v[3], HE_vert *s, HE_vert *e){
	v[0] = e->x-s->x;
	v[1] = e->y-s->y;
	v[2] = e->z-s->z;
}


static inline double LENGTH(double *v){
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline float LENGTH(float *v){
	return (float)sqrt((double)v[0]*(double)v[0] + (double)v[1]*(double)v[1] + (double)v[2]*(double)v[2]);
}

static inline double DLENGTH(float *v){
	return sqrt((double)v[0]*(double)v[0] + (double)v[1]*(double)v[1] + (double)v[2]*(double)v[2]);
}

static inline float DOT(float *v1, float *v2){
	return (float)((double)v1[0]*(double)v2[0] + (double)v1[1]*(double)v2[1] + (double)v1[2]*(double)v2[2]);
}

static inline double DDOT(float v1[3], float v2[3]){
	return (double)v1[0]*(double)v2[0] + (double)v1[1]*(double)v2[1] + (double)v1[2]*(double)v2[2];
}

static inline void VECSCALVEC(float v1[3], double scal, float n[3]){
	n[0] = (float)((double)v1[0] *scal);
	n[1] = (float)((double)v1[1] *scal);
	n[2] = (float)((double)v1[2] *scal);
	return;
}

static inline void VECSCAL(float v1[3], double scal){
	v1[0] = (float)((double)v1[0] *scal);
	v1[1] = (float)((double)v1[1] *scal);
	v1[2] = (float)((double)v1[2] *scal);
	return;
}

static inline void MINUS(float v1[3], float v2[3]){
	v1[0] -= v2[0];
	v1[1] -= v2[1];
	v1[2] -= v2[2];
	return;
}

static inline void ADD(float v1[3], float v2[3], float n[3]){
	n[0] = v1[0] + v2[0];
	n[1] = v1[1] + v2[1];
	n[2] = v1[2] + v2[2];
	return;
}

static inline void ADD(HE_vert *v1, float v2[3], float n[3]){
	n[0] = v1->x + v2[0];
	n[1] = v1->y + v2[1];
	n[2] = v1->z + v2[2];
	return;
}

static inline void ADD(HE_vert *s1, HE_vert *e1, HE_vert *s2, HE_vert *e2, float n[3]){
	n[0] = e1->x + e2->x - s1->x - s2->x;
	n[1] = e1->y + e2->y - s1->y - s2->y;
	n[2] = e1->z + e2->z - s1->z - s2->z;
	return;
}

static inline void NORMALIZE(float n[3]){
	double length = DLENGTH(n);
	length = 1.0/length;
	VECSCAL(n, length);
	return;
}


static inline float DISTANCE(HE_vert *v1, HE_vert *v2){
	double vet[3];
	VECTOR(vet, v1, v2);
	return (float)LENGTH(vet);
}

static inline float DISTANCE(float v1[3], float v2[3]){
	double vet[3];
	VECTOR(vet, v1, v2);
	return (float)LENGTH(vet);
}

static inline float DISTANCE(float v1[3], HE_vert *v2){
	double vet[3];
	VECTOR(vet, v1, v2);
	return (float)LENGTH(vet);
}

static inline void CROSS(double *n, double v1[3], double v2[3]){
	n[0] = v1[1]*v2[2] - v1[2]*v2[1];
	n[1] = v1[2]*v2[0] - v1[0]*v2[2];
	n[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

static inline void CROSS(float v1[3], float v2[3], float n[3]){
	n[0] = (float)((double)v1[1]*(double)v2[2] - (double)v1[2]*(double)v2[1]);
	n[1] = (float)((double)v1[2]*(double)v2[0] - (double)v1[0]*(double)v2[2]);
	n[2] = (float)((double)v1[0]*(double)v2[1] - (double)v1[1]*(double)v2[0]);
}

struct HE_edge
{

	HE_vert* vert;   // vertex at the end of the half-edge
	HE_edge* pair;   // oppositely oriented adjacent half-edge 
	HE_face* face;   // face the half-edge borders
	HE_edge* next;   // next half-edge around the face
	int index;

	int getIndex(){return index;}
	HE_vert *getStartVert(){return vert;}
	int getStartVertIndex(){return vert->index;}
	HE_vert *getEndVert(){return pair->vert;}
	int getEndVertIndex(){return pair->vert->index;}
	int getPairIndex(){return pair->index;}

	HE_face *getLeftFace(){return face;}
	HE_face *getRightFace(){return pair->face;}
	float getLength(){return DISTANCE(vert, pair->vert);}
	bool isClockwise(){return vert>pair->vert;}
	void getCenter(float [3]);
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

class HE_mesh
{
public:

	HE_mesh(void);
	HE_mesh(char *);
	~HE_mesh();

	void init(char *);

	HE_vert *getVertices(void){return vertices;}
	HE_edge *getEdges(void){return edges;}
	HE_face *getFaces(void){return faces;}

	HE_vert *getVertex(int vertexIndex){return &vertices[vertexIndex];}
	HE_edge *getEdge(int edgeIndex){return &edges[edgeIndex];}
	HE_face *getFace(int faceIndex){return &faces[faceIndex];}
	
	int getVertexNum(){return vertexNum;}
	int getEdgeNum(){return edgeNum;}
	int getFaceNum(){return faceNum;}

	float *getEdgeLength();
	float getEdgeLength(int edgeIndex);
	void deleteEdgeLength();

	float getEverageEdgeLength();

	float (*getFaceEdgeCenterDistances())[3];
	void deleteFaceEdgeCenterDistances();

	double (*getWeightsForVertexNormal())[3];
	void deleteWeightsForVertexNormal();

	double *getFaceAreas();
	void deleteFaceAreas();

	float (*getFaceNormalsCalAreas())[3];
	void deleteFaceNormals();

	float (*getVertexNormals())[3];
	void deleteVertexNormals();

	float (*getDistanceBetweenAdFaceCenters())[3];
	void deleteDistanceBetweenAdFaceCenters();

	void clculateSmoothedVertexNormals(float (*smoothedFaceNormal)[3], float (*smoothedVertexNormal)[3]);

private:
	
	void connectingEdge(char*);

	float (*faceNormals)[3];
	bool hasCalFN;
	float (*vertexNormals)[3];
	bool hasCalVN;
	float (*faceEdgeCenterDistances)[3];
	bool hasCalFECD;
	double (*weightsForVertexNormal)[3];
	bool hasCalWFVN;
	float *edgeLength;
	bool hasCalEL;
	double *faceAreas;
	bool hasCalFA;
	float (*distanceBetweenAdFaceCenters)[3];
	bool hasCalDBAFC;


	HE_vert *vertices;
	HE_face *faces;
	HE_edge *edges;
	int vertexNum;
	int faceNum;
	int edgeNum;

	float everageEdgeLength;


};