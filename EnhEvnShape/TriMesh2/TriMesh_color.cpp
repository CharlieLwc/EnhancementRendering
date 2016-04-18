

#include "TriMesh.h"


using namespace std;



namespace trimesh {



	void TriMesh::setColorCurvature(float curv_thre_pos, float curv_thre_neg, float scaleCur)
	{

		float scale;
		if (scaleCur < 0.f)
			scale = -scaleCur;
		else
			scale = 1.f / scaleCur;


		if (colors.size() != vertices.size())
			colors.clear();
		if (colors.empty())
			colors.resize(vertices.size());

		need_curvatures();



		for (unsigned int i = 0; i < vertices.size(); i++)
		{
			colors[i][0] = curv1[i] * scale;
			colors[i][1] = -curv1[i] * scale;
			colors[i][2] = 0.0;

			float cur = colors[i][0] + colors[i][1];
			if (colors[i][0] >= curv_thre_pos)
				colors[i][2] = 1.0;
			if (colors[i][1] >= curv_thre_neg)
				colors[i][2] = 1.0;
		}
		/*
		if (curv1[i] >= curv_thre)
		for (int d=0; d<3; d++)
		colors[i][d] = 1.0f;
		else
		for (int d=0; d<3; d++)
		colors[i][d] = 0.5f;
		*/
	}


	void TriMesh::setColorCurvatureGray(float scaleCur)
	{

		float scale;
		if (scaleCur < 0.f)
			scale = -scaleCur;
		else
			scale = 1.f / scaleCur;


		if (colors.size() != vertices.size())
			colors.clear();
		if (colors.empty())
			colors.resize(vertices.size());

		need_curvatures();

		scale *= 0.5;

		for (unsigned int i = 0; i < vertices.size(); i++)
		{
			colors[i][0] = curv1[i] * scale + 0.5f;
			colors[i][1] = colors[i][0];
			colors[i][2] = colors[i][0];
		}
		/*
		if (curv1[i] >= curv_thre)
		for (int d=0; d<3; d++)
		colors[i][d] = 1.0f;
		else
		for (int d=0; d<3; d++)
		colors[i][d] = 0.5f;
		*/
	}


	void TriMesh::setColorDcurvature()
	{
		if (colors.size() != vertices.size())
			colors.clear();
		if (colors.empty())
			colors.resize(vertices.size());

		need_dcurv();
		for (unsigned int i = 0; i < vertices.size(); i++)
		for (int d = 0; d < 3; d++)
			colors[i][d] = dcurv[i][0] * 50.0f;
	}

	void TriMesh::setColorBaseNormal(int radius)
	{


		if (colors.size() != vertices.size())
			colors.resize(vertices.size());

		needSmoothNormal(radius);


		for (unsigned int vIndex = 0; vIndex<vertices.size(); vIndex++)
		for (int d = 0; d<3; d++)
			colors[vIndex][d] = abs(smoothedNormalOrigin[vIndex][d] + 1.f) *0.5f;

	}

	void TriMesh::setColorHeight()
	{
		if (!need_height())
			return;
		
		if (colors.size() != vertices.size())
			colors.resize(vertices.size());

		float gamma = 1.f / (halfEdge.maxHeight - halfEdge.minHeight);

		for (unsigned int vIndex = 0; vIndex<vertices.size(); vIndex++)
		for (int d = 0; d<3; d++)
			colors[vIndex][d] = (halfEdge.vertexHeight[vIndex] - halfEdge.minHeight)*gamma;

	}

	void TriMesh::setColorVertexSamplePoly(int vertexIndex, int scale, float angle, bool firstTime)
	{

		if (!need_height())
			return;
		
		if (firstTime)
			setColorHeight();
		{
			float gamma = 1.f / (halfEdge.maxHeight - halfEdge.minHeight);
			for (int i = 0; i < halfEdge.blackFaceNum; i++)
			for (int vIndex = 0; vIndex < 3; vIndex++)
			for (int d = 0; d<3; d++)
				colors[faces[halfEdge.blackFaceIndex[i]].v[vIndex]][d] = (halfEdge.vertexHeight[faces[halfEdge.blackFaceIndex[i]].v[vIndex]] - halfEdge.minHeight)*gamma;
		}

		getVertexSamplePolyHeight(vertexIndex, scale, angle);


		for (int i = 0; i < halfEdge.blackFaceNum; i++)
		for (int vIndex = 0; vIndex < 3; vIndex++)
		for (int d = 0; d<3; d++)
			colors[faces[halfEdge.blackFaceIndex[i]].v[vIndex]][d] = 0.f;
	}
	void TriMesh::setColorSurfaceRelif(int scale, float scaleCur, float curv_thre_pos, float curv_thre_neg)
	{


		float scaleTemp;
		if (scaleCur < 0.f)
			scaleTemp = -scaleCur;
		else
			scaleTemp = 1.f / scaleCur;


		analysisCurvture(1, &scale, scaleTemp);


		if (colors.size() != vertices.size())
			colors.resize(vertices.size());

		float gamma = 1.f / (halfEdge.maxReliefFeature - halfEdge.minReliefFeature);

		for (unsigned int i = 0; i < vertices.size(); i++)
		{

			colors[i][0] = halfEdge.reliefFeature[i];
			colors[i][1] = -halfEdge.reliefFeature[i];
			colors[i][2] = 0.0;

			if (colors[i][0] >= curv_thre_pos)
				colors[i][2] = 1.0;
			if (colors[i][1] >= curv_thre_neg)
				colors[i][2] = 1.0;


		}

	}
	
}