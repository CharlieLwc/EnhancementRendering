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


namespace trimesh {

// Compute per-vertex normals
void TriMesh::need_normals()
{
	// Nothing to do if we already have normals
	int nv = vertices.size();
	if (int(normals.size()) == nv)
		return;

	normals.clear();
	normals.resize(nv);




	FILE *fin;	
	fopen_s(&fin, dFName.getNormalDN(), "rb");
	if (fin)
	{
		cout << "Reading Normal.....";
		int totalNum;

		fread(&totalNum, sizeof(int), 1, fin);
		if (totalNum == nv)
		{

			fread(normals[0], sizeof(float), totalNum*3, fin);
			cout << "Done" << endl;
			fclose(fin);
			return;
		}
		cout << "Done" << endl;
		fclose(fin);
	}

	/*

	ifstream normal_if(dFName.getNormalDN());
	if (normal_if.is_open())
	{
		cout<<"Reading Normal.....";
		int totalNum;
		normal_if>>totalNum;
		if (totalNum == nv)
		{
			for(int i = 0; i<totalNum; i++)
				normal_if>>normals[i];
			cout<<"Done"<<endl;
			normal_if.close();
			return;
		}
		cout<<"Done"<<endl;
		normal_if.close();
	}
	*/

	dprintf("Computing normals... ");

	// TODO: direct handling of grids
	if (!tstrips.empty()) {
		// Compute from tstrips
		const int *t = &tstrips[0], *end = t + tstrips.size();
		while (likely(t < end)) {
			int striplen = *t - 2;
			t += 3;
			bool flip = false;
			for (int i = 0; i < striplen; i++, t++, flip = !flip) {
				const point &p0 = vertices[*(t-2)];
				const point &p1 = vertices[*(t-1)];
				const point &p2 = vertices[* t   ];
				vec a = p0-p1, b = p1-p2, c = p2-p0;
				float l2a = len2(a), l2b = len2(b), l2c = len2(c);
				if (!l2a || !l2b || !l2c)
					continue;
				vec facenormal = flip ? (b CROSS a) : (a CROSS b);
				normals[*(t-2)] += facenormal * (1.0f / (l2a * l2c));
				normals[*(t-1)] += facenormal * (1.0f / (l2b * l2a));
				normals[* t   ] += facenormal * (1.0f / (l2c * l2b));
			}
		}
	} else if (need_faces(), !faces.empty()) {
		// Compute from faces
		int nf = faces.size();
#pragma omp parallel for
		for (int i = 0; i < nf; i++) {
			const point &p0 = vertices[faces[i][0]];
			const point &p1 = vertices[faces[i][1]];
			const point &p2 = vertices[faces[i][2]];
			vec a = p0-p1, b = p1-p2, c = p2-p0;
			float l2a = len2(a), l2b = len2(b), l2c = len2(c);
			if (!l2a || !l2b || !l2c)
				continue;
			vec facenormal = a CROSS b;
			normals[faces[i][0]] += facenormal * (1.0f / (l2a * l2c));
			normals[faces[i][1]] += facenormal * (1.0f / (l2b * l2a));
			normals[faces[i][2]] += facenormal * (1.0f / (l2c * l2b));
		}
	} else {
		// Find normals of a point cloud
		const int k = 6;
		const vec ref(0, 0, 1);
		KDtree kd(vertices);
#pragma omp parallel for
		for (int i = 0; i < nv; i++) {
			vector<const float *> knn;
			kd.find_k_closest_to_pt(knn, k, vertices[i]);
			int actual_k = knn.size();
			if (actual_k < 3) {
				dprintf("Warning: not enough points for vertex %d\n", i);
				normals[i] = ref;
				continue;
			}
			// Compute covariance
			float C[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
			// The below loop starts at 1, since element 0     
			// is just vertices[i] itself 
			for (int j = 1; j < actual_k; j++) {
				vec d = point(knn[j]) - vertices[i];
				for (int l = 0; l < 3; l++)
					for (int m = 0; m < 3; m++)
						C[l][m] += d[l] * d[m];
			}
			float e[3];
			eigdc<float,3>(C, e);
			normals[i] = vec(C[0][0], C[1][0], C[2][0]);
			if ((normals[i] DOT ref) < 0.0f)
				normals[i] = -normals[i];
		}
	}

	// Make them all unit-length
#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		normalize(normals[i]);

	dprintf("Done.\n");


	cout<<"Writing Normal.....";


	FILE *fout;
	fopen_s(&fout, dFName.getNormalDN(), "wb");

	fwrite(&nv, sizeof(int), 1, fout);
	fwrite(normals[0], sizeof(float), nv * 3, fout);

	fclose(fout);
	/*


	ofstream normal_of(dFName.getNormalDN());
	normal_of<<nv<<endl;
	for (int j=0; j<nv; j++)
		normal_of<<normals[j]<<endl;
	normal_of.close();
	cout<<"Done!"<<endl;
	*/
}

void TriMesh::need_faceNormals()
{
	// Nothing to do if we already have normals
	int fv = faces.size();
	if (int(faceNormals.size()) == fv)
		return;

	faceNormals.clear();
	faceNormals.resize(fv);




	FILE *fin;
	fopen_s(&fin, dFName.getFaceNormalDN(), "rb");
	if (fin)
	{
		cout << "Reading FacecNormal.....";
		int totalNum;

		fread(&totalNum, sizeof(int), 1, fin);
		if (totalNum == fv)
		{

			fread(faceNormals[0], sizeof(float), totalNum * 3, fin);
			cout << "Done" << endl;
			fclose(fin);
			return;
		}
		cout << "Done" << endl;
		fclose(fin);
	}

	need_tstrips();



	dprintf("Computing faceNormals... ");

	int faceIndex = 0;

	// TODO: direct handling of grids
	if (!tstrips.empty()) {
		// Compute from tstrips
		const int *t = &tstrips[0], *end = t + tstrips.size();
		while (likely(t < end)) {
			int striplen = *t - 2;
			t += 3;
			bool flip = false;
			for (int i = 0; i < striplen; i++, t++, flip = !flip) {
				const point &p0 = vertices[*(t - 2)];
				const point &p1 = vertices[*(t - 1)];
				const point &p2 = vertices[*t];
				vec a = p0 - p1, b = p1 - p2, c = p2 - p0;
				float l2a = len2(a), l2b = len2(b), l2c = len2(c);
				if (!l2a || !l2b || !l2c)
					continue;

				faceNormals[faceIndex++] = flip ? (b CROSS a) : (a CROSS b);
			}
		}
	}

	// Make them all unit-length
#pragma omp parallel for
	for (int i = 0; i < fv; i++)
		normalize(faceNormals[i]);

	dprintf("Done.\n");


	cout << "Writing Normal.....";


	FILE *fout;
	fopen_s(&fout, dFName.getFaceNormalDN(), "wb");

	fwrite(&fv, sizeof(int), 1, fout);
	fwrite(faceNormals[0], sizeof(float), fv * 3, fout);

	fclose(fout);
	/*


	ofstream normal_of(dFName.getNormalDN());
	normal_of<<nv<<endl;
	for (int j=0; j<nv; j++)
	normal_of<<normals[j]<<endl;
	normal_of.close();
	cout<<"Done!"<<endl;
	*/
}


}; // namespace trimesh
