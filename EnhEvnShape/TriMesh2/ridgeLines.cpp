
#include "ridgeLines.h"
#include "TriMesh.h"


using namespace std;

namespace trimesh {




	// Draw part of a ridge/valley curve on one triangle face.  v0,v1,v2
	// are the indices of the 3 vertices; this function assumes that the
	// curve connects points on the edges v0-v1 and v1-v2
	// (or connects point on v0-v1 to center if to_center is true)
	void TriMesh::draw_segment_ridge(int v0, int v1, int v2,
		float pCurv0, float pCurv1, float pCurv2,
		float emax0, float emax1, float emax2,
		float kmax0, float kmax1, float kmax2,
		::std::vector<vec> &RPI, ::std::vector<vec> &RPN, ::std::vector<float> &RPC, ::std::vector<vec> &RPD, 
		::std::vector<int> &RPE, ::std::vector<bool> &RPIV, bool to_center)
	{
		// Interpolate to find ridge/valley line segment endpoints
		// in this triangle and the curvatures there
		bool isVertex01 = false, isVertex12 = false, isVertex02 = false;
		bool isVertex011 = false, isVertex122 = false, isVertex022 = false;


		float w10 = fabs(emax0) / (fabs(emax0) + fabs(emax1));
		float w01 = 1.0f - w10;

		
		point p01 = w01 * vertices[v0] + w10 * vertices[v1];
		vec n01 = w01 * normals[v0] + w10 * normals[v1];
		normalize(n01);
		float c01 = w01 * curv1[v0] + w10 * curv1[v1];
		float k01 = fabs(w01 * kmax0 + w10 * kmax1);
		vec d01 = w01 * pdir1[v0] + w10 * pdir1[v1];
		normalize(d01);
		if (w01 < 0.00001 || w10  < 0.00001)
		{
			isVertex01 = true;
			if (w01 < 0.00001)
				isVertex011 = true;
		}
		int e01;
		if(isVertex011)
			e01 = halfEdge.he_vert[v1].findEdge(v0);
		else
			e01 = halfEdge.he_vert[v0].findEdge(v1);

		if (e01 < 0)
			return;

		if (!e01)
			cout << "error on ridge" << endl;



		float w21 = fabs(emax1) / (fabs(emax1) + fabs(emax2));
		float w12 = 1.0f - w21;
		point p12 = w12 * vertices[v1] + w21 * vertices[v2];
		vec n12 = w12 * normals[v1] + w21 * normals[v2];
		normalize(n12);
		float c12 = w12 * curv1[v1] + w21 * curv1[v2];
		float k12 = fabs(w12 * kmax1 + w21 * kmax2);
		vec d12 = w12 * pdir1[v1] + w21 * pdir1[v2];
		normalize(d12);
		if (w12 < 0.00001 || w21 < 0.00001)
		{
			isVertex12 = true;
			if (w12 < 0.00001)
				isVertex122 = true;
		}
		int e12;
		if (isVertex122)
			e12 = halfEdge.he_vert[v2].findEdge(v1);
		else
			e12 = halfEdge.he_vert[v1].findEdge(v2);

		if (!e12)
			cout << "error on ridge" << endl;


		if (e12 < 0)
			return;

		int firstTwo[3] = {0,0,0};

		if (to_center) {


			float w20 = fabs(emax2) / (fabs(emax0) + fabs(emax2));
			float w02 = 1.0f - w20;
			point p02 = w02 * vertices[v0] + w20 * vertices[v2];
			vec n02 = w02 * normals[v0] + w20 * normals[v2];
			normalize(n02);
			float c02 = w02 * curv1[v0] + w20 * curv1[v2];
			float k02 = fabs(w02 * kmax1 + w20 * kmax2);
			vec d02 = w02 * pdir1[v0] + w20 * pdir1[v2];
			normalize(d02);
			if (w02 < 0.00001 || w20 < 0.00001)
			{
				isVertex02 = true;
				if (w02 < 0.00001)
					isVertex022 = true;
			}
			int e02;
			if (isVertex022)
				e02 = halfEdge.he_vert[v2].findEdge(v0);
			else
				e02 = halfEdge.he_vert[v0].findEdge(v2);


			if (e02 < 0)
				return;

			if (!e02)
				cout << "error on ridge" << endl;

			if (c01 > c12){ firstTwo[0]++; }
			else{ firstTwo[1]++; }

			if (c02 > c12){ firstTwo[2]++; }
			else{ firstTwo[1]++; }

			if (c01 > c02){ firstTwo[0]++; }
			else{ firstTwo[2]++; }


			if (firstTwo[0] > 0)
			{
				RPI.push_back(p01);
				RPN.push_back(n01);
				RPC.push_back(c01);
				RPD.push_back(d01);
				RPE.push_back(e01);
				RPIV.push_back(isVertex01);
			}


			if (firstTwo[1] > 0)
			{
				RPI.push_back(p12);
				RPN.push_back(n12);
				RPC.push_back(c12);
				RPD.push_back(d12);
				RPE.push_back(e12);
				RPIV.push_back(isVertex12);
			}


			if (firstTwo[2] > 0)
			{
				RPI.push_back(p02);
				RPN.push_back(n02);
				RPC.push_back(c02);
				RPD.push_back(d02);
				RPE.push_back(e02);
				RPIV.push_back(isVertex02);
			}


			/*
			// Connect first point to center of triangle
			p12 = (vertices[v0] + vertices[v1] + vertices[v2]) / 3.0f;
			n12 = normals[v0] + normals[v1] + normals[v2];
			normalize(n12);
			c12 = (curv1[v0] + curv1[v1] + curv1[v2]) / 3.0f;
			k12 = fabs(kmax0 + kmax1 + kmax2) / 3.0f;
			d12 = pdir1[v0] + pdir1[v1] + pdir1[v2];
			normalize(d12);

			*/
		}
		else
		{
			RPI.push_back(p01);
			RPI.push_back(p12);
			RPN.push_back(n01);
			RPN.push_back(n12);
			RPC.push_back(c01);
			RPC.push_back(c12);
			RPD.push_back(d01);
			RPD.push_back(d12);
			RPE.push_back(e01);
			RPE.push_back(e12);
			RPIV.push_back(isVertex01);
			RPIV.push_back(isVertex12);
		}




	}



	// Draw ridges or valleys (depending on do_ridge) in a triangle v0,v1,v2
	// - uses ndotv for backface culling (enabled with do_bfcull)
	// - do_test checks for curvature maxima/minina for ridges/valleys
	//   (when off, it draws positive minima and negative maxima)
	// Note: this computes ridges/valleys every time, instead of once at the
	//   start (given they aren't view dependent, this is wasteful)
	// Algorithm based on formulas of Ohtake et al., 2004.
	void TriMesh::draw_face_ridges(int v0, int v1, int v2, bool do_ridge,
		::std::vector<vec> &RPI, ::std::vector<vec> &RPN, ::std::vector<float> &RPC, ::std::vector<vec> &RPD, 
		::std::vector<int> &RPE, ::std::vector<bool> &RPIV)
	{
		// Check if ridge possible at vertices just based on curvatures
		if (do_ridge) {
			if ((curv1[v0] <= 0.0f) ||
				(curv1[v1] <= 0.0f) ||
				(curv1[v2] <= 0.0f))
				return;
		}
		else {
			if ((curv1[v0] >= 0.0f) ||
				(curv1[v1] >= 0.0f) ||
				(curv1[v2] >= 0.0f))
				return;
		}

		// Sign of curvature on ridge/valley
		float rv_sign = do_ridge ? 1.0f : -1.0f;

		// The "tmax" are the principal directions of maximal curvature,
		// flipped to point in the direction in which the curvature
		// is increasing (decreasing for valleys).  Note that this
		// is a bit different from the notation in Ohtake et al.,
		// but the tests below are equivalent.
		const float &emax0 = dcurv[v0][0];
		const float &emax1 = dcurv[v1][0];
		const float &emax2 = dcurv[v2][0];
		vec tmax0 = rv_sign * dcurv[v0][0] * pdir1[v0];
		vec tmax1 = rv_sign * dcurv[v1][0] * pdir1[v1];
		vec tmax2 = rv_sign * dcurv[v2][0] * pdir1[v2];

		// We have a "zero crossing" if the tmaxes along an edge
		// point in opposite directions
		bool z01 = ((tmax0 DOT tmax1) <= 0.0f);
		bool z12 = ((tmax1 DOT tmax2) <= 0.0f);
		bool z20 = ((tmax2 DOT tmax0) <= 0.0f);

		if (z01 + z12 + z20 < 2)
			return;
		float do_test = true;

		if (do_test) {
			const point &p0 = vertices[v0],
				&p1 = vertices[v1],
				&p2 = vertices[v2];

			// Check whether we have the correct flavor of extremum:
			// Is the curvature increasing along the edge?
			z01 = z01 && ((tmax0 DOT(p1 - p0)) >= 0.0f ||
				(tmax1 DOT(p1 - p0)) <= 0.0f);
			z12 = z12 && ((tmax1 DOT(p2 - p1)) >= 0.0f ||
				(tmax2 DOT(p2 - p1)) <= 0.0f);
			z20 = z20 && ((tmax2 DOT(p0 - p2)) >= 0.0f ||
				(tmax0 DOT(p0 - p2)) <= 0.0f);

			if (z01 + z12 + z20 < 2)
				return;
		}

		// Draw line segment
		const float &kmax0 = curv1[v0];
		const float &kmax1 = curv1[v1];
		const float &kmax2 = curv1[v2];
		if (!z01) {
			draw_segment_ridge(v1, v2, v0,
				curv1[v0], curv1[v1], curv1[v2],
				emax1, emax2, emax0,
				kmax1, kmax2, kmax0,
				RPI, RPN, RPC, RPD, RPE, RPIV, false);
		}
		else if (!z12) {
			draw_segment_ridge(v2, v0, v1,
				curv1[v0], curv1[v1], curv1[v2],
				emax2, emax0, emax1,
				kmax2, kmax0, kmax1,
				RPI, RPN, RPC, RPD, RPE, RPIV, false);
		}
		else if (!z20) {
			draw_segment_ridge(v0, v1, v2,
				curv1[v0], curv1[v1], curv1[v2],
				emax0, emax1, emax2,
				kmax0, kmax1, kmax2,
				RPI, RPN, RPC, RPD, RPE, RPIV, false);
		}
		else {
			// All three edges have crossings -- connect all to center
			draw_segment_ridge(v1, v2, v0,
				curv1[v0], curv1[v1], curv1[v2],
				emax1, emax2, emax0,
				kmax1, kmax2, kmax0,
				RPI, RPN, RPC, RPD, RPE, RPIV, true);

			/*
			draw_segment_ridge(v2, v0, v1,
				curv1[v0], curv1[v1], curv1[v2],
				emax2, emax0, emax1,
				kmax2, kmax0, kmax1,
				RPI, RPN, RPC, RPD, RPE, true);
			draw_segment_ridge(v0, v1, v2,
				curv1[v0], curv1[v1], curv1[v2],
				emax0, emax1, emax2,
				kmax0, kmax1, kmax2,
				RPI, RPN, RPC, RPD, RPE, true);
				*/
		}
	}


	void TriMesh::sortRidgeLines(::std::vector<vec> &RPI, ::std::vector<vec> &RPN, ::std::vector<float> &RPC, ::std::vector<vec> &RPD, 
		::std::vector<int> &RPE, ::std::vector<bool> &RPIV, bool isRidge)
	{
		int *linesPerEdge = new int[halfEdge.he_edge.size()];
		int (*pointIndices)[2] = new int[halfEdge.he_edge.size()][2];



		for (int eIndex = 0; eIndex < (int)halfEdge.he_edge.size(); eIndex++)
			linesPerEdge[eIndex] = 0;



		for (int eIndex = 0; eIndex < (int)RPE.size(); eIndex++)
		{
			int currentEdge = RPE[eIndex];
			int pairEdge = halfEdge.he_edge[currentEdge].pair->index;

			if (currentEdge > pairEdge)
				currentEdge = pairEdge;
			RPE[eIndex] = currentEdge;

			if (linesPerEdge[currentEdge] <2)
				pointIndices[currentEdge][linesPerEdge[currentEdge]++] = eIndex;

		}

		if (isRidge)
		{
			for (int eIndex = 0; eIndex < (int)halfEdge.he_edge.size(); eIndex++)
			{
				if (linesPerEdge[eIndex] != 1)
					continue;

				vector<vec> ridgeLinePoint;
				ridgeLinePoint.clear();

				int currentEdge = eIndex;

				int currentPointIndex = pointIndices[currentEdge][0];
				ridgeLinePoint.push_back(RPI[currentPointIndex]);
				ridgeLines.ridgeLinesPointNormal.push_back(RPN[currentPointIndex]);
				ridgeLines.ridgeLinesPointPcurv.push_back(RPC[currentPointIndex]);
				ridgeLines.ridgeLinesPointPDir.push_back(RPD[currentPointIndex]);
				ridgeLines.ridgeLinesPointEdge.push_back(RPE[currentPointIndex]);
				ridgeLines.ridgeLinesPointIsVertex.push_back(RPIV[currentPointIndex]);

				



				while (linesPerEdge[currentEdge] > 0)
				{
					linesPerEdge[currentEdge]--;

					int nextPointIndex = currentPointIndex + 1;
					if (currentPointIndex % 2 == 1)
						nextPointIndex -= 2;

					currentEdge = RPE[nextPointIndex];


					if (pointIndices[currentEdge][0] == nextPointIndex)
						pointIndices[currentEdge][0] = pointIndices[currentEdge][1];

					currentPointIndex = pointIndices[currentEdge][0];
					linesPerEdge[currentEdge]--;

					ridgeLinePoint.push_back(RPI[nextPointIndex]);
					ridgeLines.ridgeLinesPointNormal.push_back(RPN[nextPointIndex]);
					ridgeLines.ridgeLinesPointPcurv.push_back(RPC[nextPointIndex]);
					ridgeLines.ridgeLinesPointPDir.push_back(RPD[nextPointIndex]);
					ridgeLines.ridgeLinesPointEdge.push_back(RPE[nextPointIndex]);
					ridgeLines.ridgeLinesPointIsVertex.push_back(RPIV[nextPointIndex]);
				}
				ridgeLines.ridgeLinesPoint.push_back(ridgeLinePoint);

			}

			for (int eIndex = 0; eIndex < (int)halfEdge.he_edge.size(); eIndex++)
			{

				if (linesPerEdge[eIndex] < 2 )
					continue;

				vector<vec> ridgeLinePoint;
				ridgeLinePoint.clear();

				int currentEdge = eIndex;

				int currentPointIndex = pointIndices[currentEdge][0];
				ridgeLinePoint.push_back(RPI[currentPointIndex]);
				ridgeLines.ridgeLinesPointNormal.push_back(RPN[currentPointIndex]);
				ridgeLines.ridgeLinesPointPcurv.push_back(RPC[currentPointIndex]);
				ridgeLines.ridgeLinesPointPDir.push_back(RPD[currentPointIndex]);
				ridgeLines.ridgeLinesPointEdge.push_back(RPE[currentPointIndex]);
				ridgeLines.ridgeLinesPointIsVertex.push_back(RPIV[currentPointIndex]);

				while (linesPerEdge[currentEdge] > 0)
				{
					linesPerEdge[currentEdge]--;

					int nextPointIndex = currentPointIndex + 1;
					if (currentPointIndex % 2 == 1)
						nextPointIndex -= 2;

					currentEdge = RPE[nextPointIndex];


					if (pointIndices[currentEdge][0] == nextPointIndex)
						pointIndices[currentEdge][0] = pointIndices[currentEdge][1];

					currentPointIndex = pointIndices[currentEdge][0];
					linesPerEdge[currentEdge]--;

					ridgeLinePoint.push_back(RPI[nextPointIndex]);
					ridgeLines.ridgeLinesPointNormal.push_back(RPN[nextPointIndex]);
					ridgeLines.ridgeLinesPointPcurv.push_back(RPC[nextPointIndex]);
					ridgeLines.ridgeLinesPointPDir.push_back(RPD[nextPointIndex]);
					ridgeLines.ridgeLinesPointEdge.push_back(RPE[nextPointIndex]);
					ridgeLines.ridgeLinesPointIsVertex.push_back(RPIV[nextPointIndex]);

				}
				ridgeLines.ridgeLinesPoint.push_back(ridgeLinePoint);



			}
		}
		else
		{
			for (int eIndex = 0; eIndex < (int)halfEdge.he_edge.size(); eIndex++)
			{
				if (linesPerEdge[eIndex] != 1)
					continue;

				vector<vec> ridgeLinePoint;
				ridgeLinePoint.clear();

				int currentEdge = eIndex;

				int currentPointIndex = pointIndices[currentEdge][0];
				ridgeLinePoint.push_back(RPI[currentPointIndex]);
				ridgeLines.valleyLinesPointNormal.push_back(RPN[currentPointIndex]);
				ridgeLines.valleyLinesPointPcurv.push_back(RPC[currentPointIndex]);
				ridgeLines.valleyLinesPointPDir.push_back(RPD[currentPointIndex]);
				ridgeLines.valleyLinesPointEdge.push_back(RPE[currentPointIndex]);
				ridgeLines.valleyLinesPointIsVertex.push_back(RPIV[currentPointIndex]);

				while (linesPerEdge[currentEdge] > 0)
				{
					linesPerEdge[currentEdge]--;

					int nextPointIndex = currentPointIndex + 1;
					if (currentPointIndex % 2 == 1)
						nextPointIndex -= 2;

					currentEdge = RPE[nextPointIndex];

					if (pointIndices[currentEdge][0] == nextPointIndex)
						pointIndices[currentEdge][0] = pointIndices[currentEdge][1];

					currentPointIndex = pointIndices[currentEdge][0];
					linesPerEdge[currentEdge]--;

					ridgeLinePoint.push_back(RPI[nextPointIndex]);
					ridgeLines.valleyLinesPointNormal.push_back(RPN[nextPointIndex]);
					ridgeLines.valleyLinesPointPcurv.push_back(RPC[nextPointIndex]);
					ridgeLines.valleyLinesPointPDir.push_back(RPD[nextPointIndex]);
					ridgeLines.valleyLinesPointEdge.push_back(RPE[nextPointIndex]);
					ridgeLines.valleyLinesPointIsVertex.push_back(RPIV[nextPointIndex]);

				}
				ridgeLines.valleyLinesPoint.push_back(ridgeLinePoint);
			}

			for (int eIndex = 0; eIndex < (int)halfEdge.he_edge.size(); eIndex++)
			{

				if (linesPerEdge[eIndex] < 2)
					continue;

				vector<vec> ridgeLinePoint;
				ridgeLinePoint.clear();

				int currentEdge = eIndex;

				int currentPointIndex = pointIndices[currentEdge][0];
				ridgeLinePoint.push_back(RPI[currentPointIndex]);
				ridgeLines.valleyLinesPointNormal.push_back(RPN[currentPointIndex]);
				ridgeLines.valleyLinesPointPcurv.push_back(RPC[currentPointIndex]);
				ridgeLines.valleyLinesPointPDir.push_back(RPD[currentPointIndex]);
				ridgeLines.valleyLinesPointEdge.push_back(RPE[currentPointIndex]);
				ridgeLines.valleyLinesPointIsVertex.push_back(RPIV[currentPointIndex]);

				while (linesPerEdge[currentEdge] > 0)
				{
					linesPerEdge[currentEdge]--;

					int nextPointIndex = currentPointIndex + 1;
					if (currentPointIndex % 2 == 1)
						nextPointIndex -= 2;

					currentEdge = RPE[nextPointIndex];


					if (pointIndices[currentEdge][0] == nextPointIndex)
						pointIndices[currentEdge][0] = pointIndices[currentEdge][1];

					currentPointIndex = pointIndices[currentEdge][0];
					linesPerEdge[currentEdge]--;

					ridgeLinePoint.push_back(RPI[nextPointIndex]);
					ridgeLines.valleyLinesPointNormal.push_back(RPN[nextPointIndex]);
					ridgeLines.valleyLinesPointPcurv.push_back(RPC[nextPointIndex]);
					ridgeLines.valleyLinesPointPDir.push_back(RPD[nextPointIndex]);
					ridgeLines.valleyLinesPointEdge.push_back(RPE[nextPointIndex]);
					ridgeLines.valleyLinesPointIsVertex.push_back(RPIV[nextPointIndex]);

				}
				ridgeLines.valleyLinesPoint.push_back(ridgeLinePoint);



			}
		}

		delete[] linesPerEdge;
		delete[] pointIndices;
	}

	void TriMesh::need_RidgeLines(bool isRidge)
	{


		if ((isRidge && !ridgeLines.ridgeLinesPointNormal.empty()) || (!isRidge && !ridgeLines.valleyLinesPointNormal.empty()))
			return;

		FILE *fin;
		fopen_s(&fin, dFName.getRidgeLinesDN(isRidge), "rb");
		if (fin)
		{

			int totalNum, lineNum;
			fread(&totalNum, sizeof(int), 1, fin);
			fread(&lineNum, sizeof(int), 1, fin);

			if (lineNum == 0)
			{

				fclose(fin);
				return;
			}

			int* lineLenth = new int[lineNum];
			fread(&lineLenth[0], sizeof(int), lineNum, fin);

			

			if (isRidge)
			{
				cout << "Reading ridge lines.....";

				ridgeLines.ridgeLinesPoint.resize(lineNum);

				for (int i = 0; i < lineNum; i++)
				{
					ridgeLines.ridgeLinesPoint[i].resize(lineLenth[i]);
					fread(ridgeLines.ridgeLinesPoint[i][0], sizeof(float), lineLenth[i] * 3, fin);
				}

				ridgeLines.ridgeLinesPointNormal.resize(totalNum);
				ridgeLines.ridgeLinesPointPcurv.resize(totalNum);
				ridgeLines.ridgeLinesPointPDir.resize(totalNum);
				ridgeLines.ridgeLinesPointEdge.resize(totalNum);
				ridgeLines.ridgeLinesPointIsVertex.resize(totalNum);

				fread(ridgeLines.ridgeLinesPointNormal[0], sizeof(float), totalNum * 3, fin);
				fread(&ridgeLines.ridgeLinesPointPcurv[0], sizeof(float), totalNum, fin);
				fread(ridgeLines.ridgeLinesPointPDir[0], sizeof(float), totalNum * 3, fin);
				fread(&ridgeLines.ridgeLinesPointEdge[0], sizeof(int), totalNum, fin);


				vector<int> temp;
				temp.resize(totalNum);

				fread(&temp[0], sizeof(int), totalNum, fin);


				for (unsigned int i = 0; i < ridgeLines.ridgeLinesPointIsVertex.size(); i++)
					if (temp[i] == 1)
						ridgeLines.ridgeLinesPointIsVertex[i] = true;
					else
						ridgeLines.ridgeLinesPointIsVertex[i] = false;



				

			}
			else
			{
				cout << "Reading valley lines.....";



				ridgeLines.valleyLinesPoint.resize(lineNum);

				for (int i = 0; i < lineNum; i++)
				{
					ridgeLines.valleyLinesPoint[i].resize(lineLenth[i]);
					fread(ridgeLines.valleyLinesPoint[i][0], sizeof(float), lineLenth[i] * 3, fin);
				}


				ridgeLines.valleyLinesPointNormal.resize(totalNum);
				ridgeLines.valleyLinesPointPcurv.resize(totalNum);
				ridgeLines.valleyLinesPointPDir.resize(totalNum);
				ridgeLines.valleyLinesPointEdge.resize(totalNum);
				ridgeLines.valleyLinesPointIsVertex.resize(totalNum);


				fread(ridgeLines.valleyLinesPointNormal[0], sizeof(float), totalNum * 3, fin);
				fread(&ridgeLines.valleyLinesPointPcurv[0], sizeof(float), totalNum, fin);
				fread(ridgeLines.valleyLinesPointPDir[0], sizeof(float), totalNum * 3, fin);
				fread(&ridgeLines.valleyLinesPointEdge[0], sizeof(int), totalNum, fin);
				vector<int> temp;
				temp.resize(totalNum);

				fread(&temp[0], sizeof(int), totalNum, fin);


				for (unsigned int i = 0; i < ridgeLines.valleyLinesPointIsVertex.size(); i++)
					if (temp[i] == 1)
						ridgeLines.valleyLinesPointIsVertex[i] = true;
					else
						ridgeLines.valleyLinesPointIsVertex[i] = false;

			}
			cout << "Done" << endl;
			fclose(fin);
			return;

		}

		need_halfEdge();

		::std::vector<point> ridgePoint;
		::std::vector<vec> ridgePointNormal;
		::std::vector<float> ridgePointPcurv;
		::std::vector<vec> ridgePointPDir;
		::std::vector<int> ridgePointEdge;
		::std::vector<bool> ridgePointIsVertex;

		need_neighbors();
		need_curvatures();
		need_dcurv();

		//	tune_curve_direction();
		need_normals();
		float threhould = 0.00;

		if (isRidge)
			dprintf("Computing RidgeLines... ");
		else
			dprintf("Computing ValleyLines... ");

		const int *t = &tstrips[0];
		const int *stripend = t;
		const int *end = t + tstrips.size();

		// Walk through triangle strips

		while (1) {
			if (unlikely(t >= stripend)) {
				if (unlikely(t >= end))
					break;
				// New strip: each strip is stored as
				// length followed by indices
				stripend = t + 1 + *t;
				// Skip over length plus first two indices of
				// first face
				t += 3;
			}
			draw_face_ridges(*(t - 2), *(t - 1), *t,
				isRidge, ridgePoint, ridgePointNormal, ridgePointPcurv, ridgePointPDir, ridgePointEdge, ridgePointIsVertex);
			t++;
		}

		dprintf("Done.\n");


		sortRidgeLines(ridgePoint, ridgePointNormal, ridgePointPcurv, ridgePointPDir, ridgePointEdge, ridgePointIsVertex, isRidge);



		FILE *fout;
		fopen_s(&fout, dFName.getRidgeLinesDN(isRidge), "wb");

		if (isRidge)
		{


			int lineNum = ridgeLines.ridgeLinesPoint.size();
			int *lineLength = new int[lineNum];

			for (int i = 0; i < lineNum; i++)
				lineLength[i] = ridgeLines.ridgeLinesPoint[i].size();

			cout << "Writing ridge lines.....";
			int totalNum = ridgeLines.ridgeLinesPointNormal.size();
			fwrite(&totalNum, sizeof(int), 1, fout);
			fwrite(&lineNum, sizeof(int), 1, fout);
			fwrite(&lineLength[0], sizeof(int), lineNum, fout);

			for (int i = 0; i < lineNum; i++)
				fwrite(ridgeLines.ridgeLinesPoint[i][0], sizeof(float), lineLength[i] * 3, fout);

			fwrite(ridgeLines.ridgeLinesPointNormal[0], sizeof(float), totalNum * 3, fout);
			fwrite(&ridgeLines.ridgeLinesPointPcurv[0], sizeof(float), totalNum, fout);
			fwrite(ridgeLines.ridgeLinesPointPDir[0], sizeof(float), totalNum * 3, fout);
			fwrite(&ridgeLines.ridgeLinesPointEdge[0], sizeof(int), totalNum, fout);

			vector<int> temp;
			temp.resize(ridgeLines.ridgeLinesPointIsVertex.size());
			for (unsigned int i = 0; i < ridgeLines.ridgeLinesPointIsVertex.size(); i++)
				if (ridgeLines.ridgeLinesPointIsVertex[i])
					temp[i] = 1;
				else
					temp[i] = 0;

			fwrite(&temp[0], sizeof(int), totalNum, fout);




			delete[] lineLength;
		}
		else
		{



			int lineNum = ridgeLines.valleyLinesPoint.size();

			cout << "Writing valley lines.....";
			int totalNum = ridgeLines.valleyLinesPointNormal.size();
			fwrite(&totalNum, sizeof(int), 1, fout);
			fwrite(&lineNum, sizeof(int), 1, fout);

			if (lineNum != 0)
			{


				int *lineLength = new int[lineNum];

				for (int i = 0; i < lineNum; i++)
					lineLength[i] = ridgeLines.valleyLinesPoint[i].size();


				fwrite(&lineLength[0], sizeof(int), lineNum, fout);


				for (int i = 0; i < lineNum; i++)
					fwrite(ridgeLines.valleyLinesPoint[i][0], sizeof(float), lineLength[i] * 3, fout);

				fwrite(ridgeLines.valleyLinesPointNormal[0], sizeof(float), totalNum * 3, fout);
				fwrite(&ridgeLines.valleyLinesPointPcurv[0], sizeof(float), totalNum, fout);
				fwrite(ridgeLines.valleyLinesPointPDir[0], sizeof(float), totalNum * 3, fout);
				fwrite(&ridgeLines.valleyLinesPointEdge[0], sizeof(int), totalNum, fout);



				vector<int> temp;
				temp.resize(ridgeLines.valleyLinesPointIsVertex.size());
				for (unsigned int i = 0; i < ridgeLines.valleyLinesPointIsVertex.size(); i++)
					if (ridgeLines.valleyLinesPointIsVertex[i])
						temp[i] = 1;
					else
						temp[i] = 0;

				fwrite(&temp[0], sizeof(int), totalNum, fout);


				delete[] lineLength;

			}


		}
		fclose(fout);


		/*
		ofstream ridgeLines_of(dFName.getRidgeLinesDN(isRidge));

		if (isRidge)
		{
		ridgeLines_of << ridgePoint.size() << endl;
		for (int lineIndex = 0; lineIndex<ridgePoint.size(); lineIndex++)
		{
		ridgeLines_of << ridgePoint[lineIndex] << " " << ridgePointNormal[lineIndex] << " " << ridgePointPcurv[lineIndex] << endl;
		}
		}
		else
		{
		ridgeLines_of << valleyPoint.size() << endl;
		for (int lineIndex = 0; lineIndex<valleyPoint.size(); lineIndex++)
		{
		ridgeLines_of << valleyPoint[lineIndex] << " " << valleyPointNormal[lineIndex] << " " << valleyPointPcurv[lineIndex] << endl;
		}
		}

		ridgeLines_of.close();
		*/








		cout << "Done!" << endl;

	}






};