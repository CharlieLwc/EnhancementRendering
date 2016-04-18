/*
Szymon Rusinkiewicz
Princeton University

TriMesh_pointareas.cc
Compute the area "belonging" to each vertex or each corner
of a triangle (defined as Voronoi area restricted to the 1-ring of
a vertex, or to the triangle).
*/

#include "TriMesh.h"

using namespace std;

namespace trimesh {

// Compute per-vertex point areas
void TriMesh::need_pointareas()
{
	if (pointareas.size() == vertices.size())
		return;
	need_faces();


	int nf = faces.size(), nv = vertices.size();

	FILE *fin;
	fopen_s(&fin, dFName.getPointareasDN(), "rb");
	if (fin)
	{
		cout << "Reading pointareas.....";
		int verNum, facNum;

		fread(&verNum, sizeof(int), 1, fin);
		fread(&facNum, sizeof(int), 1, fin);
		if (verNum == nv && facNum == nf)
		{
			pointareas.resize(nv); cornerareas.resize(nf);

			fread(&pointareas[0], sizeof(float), nv, fin);
			fread(&cornerareas[0], sizeof(float), nf*3, fin);


			cout << "Done" << endl;
			fclose(fin);
			return;
		}
		cout << "Done" << endl;
		fclose(fin);
	}



		/*
	ifstream pointareas_if(dFName.getPointareasDN());
	if (pointareas_if.is_open() && pointareas.empty())
	{
		cout<<"Reading pointareas.....";
		int verNum, facNum;
		pointareas_if>>verNum;
		pointareas_if>>facNum;
		if (verNum == nv && facNum == nf)
		{
			pointareas.resize(nv); cornerareas.resize(nf);
			for(int i = 0; i<nv; i++)
				pointareas_if>>pointareas[i];
			for(int i = 0; i<nf; i++)
				pointareas_if>>cornerareas[i];
			cout<<"Done"<<endl;
			pointareas_if.close();
			return;
		}
		cout<<"Done"<<endl;
		pointareas_if.close();
	}

	*/

	dprintf("Computing point areas... ");

	pointareas.clear();
	pointareas.resize(nv);
	cornerareas.clear();
	cornerareas.resize(nf);

#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		// Edges
		vec e[3] = { vertices[faces[i][2]] - vertices[faces[i][1]],
			     vertices[faces[i][0]] - vertices[faces[i][2]],
			     vertices[faces[i][1]] - vertices[faces[i][0]] };

		// Compute corner weights
		float area = 0.5f * len(e[0] CROSS e[1]);
		float l2[3] = { len2(e[0]), len2(e[1]), len2(e[2]) };
		float ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
				l2[1] * (l2[2] + l2[0] - l2[1]),
				l2[2] * (l2[0] + l2[1] - l2[2]) };
		if (ew[0] <= 0.0f) {
			cornerareas[i][1] = -0.25f * l2[2] * area /
					    (e[0] DOT e[2]);
			cornerareas[i][2] = -0.25f * l2[1] * area /
					    (e[0] DOT e[1]);
			cornerareas[i][0] = area - cornerareas[i][1] -
					    cornerareas[i][2];
		} else if (ew[1] <= 0.0f) {
			cornerareas[i][2] = -0.25f * l2[0] * area /
					    (e[1] DOT e[0]);
			cornerareas[i][0] = -0.25f * l2[2] * area /
					    (e[1] DOT e[2]);
			cornerareas[i][1] = area - cornerareas[i][2] -
					    cornerareas[i][0];
		} else if (ew[2] <= 0.0f) {
			cornerareas[i][0] = -0.25f * l2[1] * area /
					    (e[2] DOT e[1]);
			cornerareas[i][1] = -0.25f * l2[0] * area /
					    (e[2] DOT e[0]);
			cornerareas[i][2] = area - cornerareas[i][0] -
					    cornerareas[i][1];
		} else {
			float ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
				cornerareas[i][j] = ewscale * (ew[(j+1)%3] +
							       ew[(j+2)%3]);
		}
#pragma omp atomic
		pointareas[faces[i][0]] += cornerareas[i][0];
#pragma omp atomic
		pointareas[faces[i][1]] += cornerareas[i][1];
#pragma omp atomic
		pointareas[faces[i][2]] += cornerareas[i][2];
	}

	dprintf("Done.\n");

	cout<<"Writing pointareas.....";


	FILE *fout;
	fopen_s(&fout, dFName.getPointareasDN(), "wb");

	fwrite(&nv, sizeof(int), 1, fout);
	fwrite(&nf, sizeof(int), 1, fout);
	fwrite(&pointareas[0], sizeof(float), nv, fout);
	fwrite(&cornerareas[0], sizeof(float), nf*3, fout);
	fclose(fout);



	/*
	ofstream pointareas_of(dFName.getPointareasDN());
	pointareas_of<<nv<<" "<<nf<<endl;
	for (int j=0; j<nv; j++)
		pointareas_of<<pointareas[j]<<endl;
	for (int j=0; j<nf; j++)
		pointareas_of<<cornerareas[j]<<endl;
	pointareas_of.close();
	
	*/



	cout<<"Done!"<<endl;

}

}; // namespace trimesh
