#ifndef TRIMESH_H
#define TRIMESH_H
/*
Szymon Rusinkiewicz
Princeton University

TriMesh.h
Class for triangle meshes.
*/

#include "Vec.h"
#include "Box.h"
#include "Color.h"
#include "DFName.h"
#include <vector>
#include <string>

#include "HalfEdge.h"
#include "ridgeLines.h"


#ifndef M_PIf
# define M_PIf 3.1415927f
#endif


namespace trimesh {






class TriMesh {
public:
	//
	// Types
	//
	struct Face {
		int v[3];

		Face() {}
		Face(const int &v0, const int &v1, const int &v2)
			{ v[0] = v0; v[1] = v1; v[2] = v2; }
		Face(const int *v_)
			{ v[0] = v_[0]; v[1] = v_[1]; v[2] = v_[2]; }
		template <class S> explicit Face(const S &x)
			{ v[0] = x[0];  v[1] = x[1];  v[2] = x[2]; }
		int &operator[] (int i) { return v[i]; }
		const int &operator[] (int i) const { return v[i]; }
		operator const int * () const { return &(v[0]); }
		operator const int * () { return &(v[0]); }
		operator int * () { return &(v[0]); }
		int indexof(int v_) const
		{
			return (v[0] == v_) ? 0 :
			       (v[1] == v_) ? 1 :
			       (v[2] == v_) ? 2 : -1;
		}
	};

	struct BSphere {
		point center;
		float r;
		bool valid;
		BSphere() : valid(false)
			{}
	};

	//
	// Enums
	//
	enum TstripRep { TSTRIP_LENGTH, TSTRIP_TERM };
	enum { GRID_INVALID = -1 };
	enum StatOp { STAT_MIN, STAT_MAX, STAT_MEAN, STAT_MEANABS,
		STAT_RMS, STAT_MEDIAN, STAT_STDEV, STAT_TOTAL };
	enum StatVal { STAT_VALENCE, STAT_FACEAREA, STAT_ANGLE,
		STAT_DIHEDRAL, STAT_EDGELEN, STAT_X, STAT_Y, STAT_Z };

	//
	// Constructor
	//

	TriMesh() : grid_width(-1), grid_height(-1), flag_curr(0), tune_pdir1(false)
	{}

	TriMesh(const char* d, const char *n) : grid_width(-1), grid_height(-1), flag_curr(0), dFName(d, n), tune_pdir1(false)
	{}

	//
	// Members
	//
	DFName dFName;

	// The basics: vertices and faces
	::std::vector<point> vertices;
	::std::vector<Face> faces;

	// Triangle strips
	::std::vector<int> tstrips;

	// Grid, if present
	::std::vector<int> grid;
	int grid_width, grid_height;

	// Other per-vertex properties
	::std::vector<Color> colors;
	::std::vector<float> confidences;
	::std::vector<unsigned> flags;
	unsigned flag_curr;

	// Computed per-vertex properties
	::std::vector<vec> normals;
	::std::vector<vec> faceNormals;
	::std::vector<vec> pdir1, pdir2;
	::std::vector<float> curv1, curv2;
	float averageCurv1;
	::std::vector< Vec<4,float> > dcurv;
	::std::vector<vec> cornerareas;
	::std::vector<float> pointareas;



	// Bounding structures
	box bbox;
	BSphere bsphere;

	// Connectivity structures:
	//  For each vertex, all neighboring vertices
	::std::vector< ::std::vector<int> > neighbors;
	//  For each vertex, all neighboring faces
	::std::vector< ::std::vector<int> > adjacentfaces;
	//  For each face, the three faces attached to its edges
	//  (for example, across_edge[3][2] is the number of the face
	//   that's touching the edge opposite vertex 2 of face 3)
	::std::vector<Face> across_edge;




	//
	// Compute all this stuff...
	//
	void need_tstrips();
	void convert_strips(TstripRep rep);
	void unpack_tstrips();
	void triangulate_grid(bool remove_slivers = true);
	void need_faces()
	{
		if (!faces.empty())
			return;
		if (!tstrips.empty())
			unpack_tstrips();
		else if (!grid.empty())
			triangulate_grid();
	}
	void need_normals();
	void need_faceNormals();
	void need_pointareas();
	void need_curvatures();
	void reCalCurvatures();
	void need_dcurv();
	void need_bbox();
	void need_bsphere();
	void need_neighbors();
	void need_adjacentfaces();
	void need_across_edge();



	RidgeLines ridgeLines;
	void sortRidgeLines(::std::vector<vec> &RPI, ::std::vector<vec> &RPN, ::std::vector<float> &RPC, ::std::vector<vec> &RPD, ::std::vector<int> &RPE, ::std::vector<bool> &RPIV, bool isRidge);
	void need_RidgeLines(bool isRidge);

	void draw_face_ridges(int v0, int v1, int v2, bool do_ridge,
		::std::vector<vec> &RPI, ::std::vector<vec> &RPN, ::std::vector<float> &RPC, ::std::vector<vec> &RPD, ::std::vector<int> &RPE, ::std::vector<bool> &RPIV);
	void draw_segment_ridge(int v0, int v1, int v2, float pCurv0, float pCurv1, float pCurv2, float emax0, float emax1, float emax2, float kmax0, float kmax1, float kmax2,
		::std::vector<vec> &RPI, ::std::vector<vec> &RPN, ::std::vector<float> &RPC, ::std::vector<vec> &RPD, ::std::vector<int> &RPE, ::std::vector<bool> &RPIV, bool to_center);



	/***************************************************/
	/**************** cluster  ************************/
	/***************************************************/


	void need_cluster(int layerIndex);
	void save_cluster(int layerIndex);
	//how many vertex in each cluster
	void need_clusterSize(int layerIndex);
	//for adaptive smoothing
	void need_calHalfArc(int layerIndex);
	//the neibour of each cluster
	void need_topology(int layerIndex);
	//as it said.
	void need_averageEdgeLength();
	//as it said.
	void need_edgeLength();
	//average normal of each cluster
	void need_clusterNormal(int layerIndex);

	void need_upAndLowScale(int layer, int startIndex, int endIndex);

	void need_clusterScale(int layerIndex);

	//set random color for each cluster
	void setColorCluster(int layerIndex);

	void setNormalClusterNormal(int layerIndex);

	void setColorClusterNormal(int layerIndex);
	void setColorSmoothedNormal(int layerIndex, int radius, bool useCluster, float smoothness);


	void needSmoothNormal(int radius);
	void smoothNormal(int layerIndex, int radius, bool useCluster, float smoothness);
	void smoothNormalWhithOutCluster(int layerIndex, int radius);
	void smoothNormalWhithCluster(int layerIndex, int radius);
	void smoothNormalWithEdge(int layerIndex, int radius);
	void getColorSmoothedCurvature(int layerIndex, int radius, float curv_thre_pos, float curv_thre_neg, float scaleCur);
	void need_smoothedCurvature(int layerIndex, int radius);
	void needWeightedCurvature(float* weighted, float scaleCur);
	
	//as it said
	void creatLayer(int layerIndex, float angle, float minCLusterSize, float maxClusterSize);
	void clusterVertex(float angle);
	void reclusterVertex(float angle);
	void clusterCluster(int layerIndex, float angle);
	void calScale(float* tempClusterScale, int* tempVertexClusterIndex, int tempClusterNum);
	void calScaleVec(vector<float> &tempClusterScale, vector<int> &tempVertexClusterIndex, int tempClusterNum);


	::std::vector<int> vertexCluster[4];
	::std::vector<int> clusterSize[4];
	::std::vector<::std::vector<int>> clusterTopology[4];
	::std::vector<::std::vector<float>> clusterHalfArc[4];

	::std::vector<float> averageEdgeLength;
	::std::vector<float> edgeLength;

	::std::vector<vec> clusterNormal[4];
	::std::vector<vec> clusterNormalForVertex;

	::std::vector<float> clusterScale[4];
	::std::vector<vec> smoothedNormal[4];
	::std::vector<vec> smoothedNormalOrigin;
	::std::vector<float> smoothedCurvature[4];
	::std::vector<float> weightedCurvature;


	float lowClusterScale[4];
	float upClusterScale[4];
	::std::vector<float> averageClusterScale[4];

	int smoothRadius[4];
	int smoothRadiusForCurv[4];
	int clusterNum[4];
	float clusterAngle[4];

	vec3 eyePosition;
	/***************************************************/
	/**************** cluster  ************************/
	/***************************************************/


	void getSmoothNormal(vector<vec3>& resultNormal, float radius);
	void getSmoothNormalCluster(vector<vec3>& resultNormal, int layerIndex);



	/***************************************************/
	/**************** halfEdge  ************************/
	/***************************************************/

	HalfEdge halfEdge;

	void need_halfEdge();
	void connectingEdge();

	bool need_height();

	void setColorVertexSamplePoly(int, int, float, bool);
	void setColorSurfaceRelif(int, float, float, float);
	

	void findSamplesHeight(vec3 planeNormal, int sampleNum, int vertexIndex, float *sampleErrorEstimate);
	void findSamplesNormal(vec3 planeNormal, int sampleNum, int vertexIndex);
	void findSamplesNormalFromEdge(int sampleNum, int sampleIndex, vec3& lineDir);
	void analysisCurvture(int scaleNum, int* scale, float scaleCur);
	void getVertexSamplePolyHeight(int vertexIndex, int scale, float angle);
	
	void getVertexSampleNormal(int vertexIndex, int scale);
	void getRidgeSampleNormal(int vertexIndex, int scale, float curvThreshold, float planeThreshold);
	int setColorRidgeVertexCurvPlane(int, int, float scaleCur, float ridgeThreshold, float curvThreshold, float planeThreshold);
	void setColorRidgeCurvPlane(int scale, float scaleCur, float ridgeThreshold, float curvThreshold, float planeThreshold);


	bool getPlaneNormal();
	void writePlaneNormal(int scale, float ridgeThreshold, float curvThreshold, float planeThreshold);

	/***************************************************/
	/**************** halfEdge  ************************/
	/***************************************************/


	/***************************************************/
	/**************** color  ************************/
	/***************************************************/

	void setColorCurvature(float curv_thre_pos, float curv_thre_neg, float scaleCur);
	void setColorCurvatureGray(float scaleCur);
	void setColorDcurvature();
	void setColorBaseNormal(int radius);
	void setColorHeight();
	

	/***************************************************/
	/**************** color  ************************/
	/***************************************************/



	//user defined
	void tune_curve_direction();


	//
	// Delete everything
	//
	void clear()
	{
		vertices.clear(); faces.clear(); tstrips.clear();
		grid.clear(); grid_width = grid_height = -1;
		colors.clear(); confidences.clear();
		flags.clear(); flag_curr = 0;
		normals.clear(); pdir1.clear(); pdir2.clear();
		curv1.clear(); curv2.clear(); dcurv.clear();
		cornerareas.clear(); pointareas.clear();
		bbox.valid = bsphere.valid = false;
		neighbors.clear(); adjacentfaces.clear(); across_edge.clear();
	}

	//
	// Input and output
	//
protected:
	static bool read_helper(const char *filename, TriMesh *mesh);
	bool tune_pdir1;
public:
	static TriMesh *read(const char *direction, const char *filename);
	static TriMesh *read(const ::std::string &filename);
	bool write(const char *filename);
	bool write(const ::std::string &filename);


	//
	// Useful queries
	//

	// Is vertex v on the mesh boundary?
	bool is_bdy(int v)
	{
		if (neighbors.empty()) need_neighbors();
		if (adjacentfaces.empty()) need_adjacentfaces();
		return neighbors[v].size() != adjacentfaces[v].size();
	}

	// Centroid of face f
	vec centroid(int f)
	{
		if (faces.empty()) need_faces();
		return (1.0f / 3.0f) *
			(vertices[faces[f][0]] +
			 vertices[faces[f][1]] +
			 vertices[faces[f][2]]);
	}

	// Normal of face f
	vec trinorm(int f)
	{
		if (faces.empty()) need_faces();
		return trimesh::trinorm(vertices[faces[f][0]], vertices[faces[f][1]],
			vertices[faces[f][2]]);
	}

	// Angle of corner j in triangle i
	float cornerangle(int i, int j)
	{
		using namespace ::std;

		if (faces.empty()) need_faces();
		const point &p0 = vertices[faces[i][j]];
		const point &p1 = vertices[faces[i][(j+1)%3]];
		const point &p2 = vertices[faces[i][(j+2)%3]];
		return acos((p1 - p0) DOT (p2 - p0));
	}

	// Dihedral angle between face i and face across_edge[i][j]
	float dihedral(int i, int j)
	{
		if (across_edge.empty()) need_across_edge();
		if (across_edge[i][j] < 0) return 0.0f;
		vec mynorm = trinorm(i);
		vec othernorm = trinorm(across_edge[i][j]);
		float ang = angle(mynorm, othernorm);
		vec towards = 0.5f * (vertices[faces[i][(j+1)%3]] +
		                      vertices[faces[i][(j+2)%3]]) -
		              vertices[faces[i][j]];
		if ((towards DOT othernorm) < 0.0f)
			return M_PIf + ang;
		else
			return M_PIf - ang;
	}

	// Statistics
	float stat(StatOp op, StatVal val);
	float feature_size();

	//
	// Debugging
	//
	// Debugging printout, controllable by a "verbose"ness parameter
	static int verbose;
	static void set_verbose(int);
	static void (*dprintf_hook)(const char *);
	static void set_dprintf_hook(void (*hook)(const char *));
	static void dprintf(const char *format, ...);

	// Same as above, but fatal-error printout
	static void (*eprintf_hook)(const char *);
	static void set_eprintf_hook(void (*hook)(const char *));
	static void eprintf(const char *format, ...);

	::std::vector<Color>& getColor()
	{
		if (colors.size() != vertices.size())
			colors.clear();
		if (colors.empty())
			colors.resize(vertices.size());
		return colors;
	}

};

}; // namespace trimesh

#endif
