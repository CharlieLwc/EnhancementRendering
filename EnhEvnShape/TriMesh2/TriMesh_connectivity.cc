/*
Szymon Rusinkiewicz
Princeton University

TriMesh_connectivity.cc
Manipulate data structures that describe connectivity between faces and verts.
*/


#include "TriMesh.h"
#include <algorithm>
using namespace std;


namespace trimesh {

// Find the direct neighbors of each vertex
void TriMesh::need_neighbors()
{
	if (!neighbors.empty())
		return;

	need_faces();
	if (faces.empty())
		return;


	int nv = vertices.size(), nf = faces.size();




	FILE *fin;
	fopen_s(&fin, dFName.getNeighborsDN(), "rb");
	if (fin)
	{
		cout << "vertex neighbors.....";
		int totalNum;

		fread(&totalNum, sizeof(int), 1, fin);
		if (totalNum == nv)
		{
			unsigned int *s = new unsigned int[nv];
			neighbors.resize(nv);
			fread(s, sizeof(int), totalNum, fin);

			for (int vertexIndex = 0; vertexIndex < nv; vertexIndex++)
			if (s[vertexIndex])
			{
				neighbors[vertexIndex].resize(s[vertexIndex]);
				fread(&neighbors[vertexIndex][0], sizeof(int), s[vertexIndex], fin);
			}

			delete[] s;
			cout << "Done" << endl;
			fclose(fin);
			return;
		}
		cout << "Done" << endl;
		fclose(fin);
	}


	/*
	ifstream neighbors_if(dFName.getNeighborsDN());
	if (neighbors_if.is_open())
	{
		cout<<"Reading vertex neighbors.....";
		int totalNum;
		neighbors_if>>totalNum;
		if (totalNum == nv)
		{
			neighbors.resize(nv);
			for (int vertexIndex=0; vertexIndex<nv; vertexIndex++)
			{
				int nn;
				neighbors_if>>nn;
				neighbors[vertexIndex].resize(nn);
				for (int neighborIndex=0; neighborIndex<nn; neighborIndex++)
					neighbors_if>>neighbors[vertexIndex][neighborIndex];
			}

			cout<<"Done"<<endl;
			neighbors_if.close();
			return;
		}
		cout<<"Done"<<endl;
		neighbors_if.close();
	}
	*/

	dprintf("Finding vertex neighbors... ");

	vector<int> numneighbors(nv);
	for (int i = 0; i < nf; i++) {
		numneighbors[faces[i][0]]++;
		numneighbors[faces[i][1]]++;
		numneighbors[faces[i][2]]++;
	}

	neighbors.resize(nv);
	for (int i = 0; i < nv; i++)
		neighbors[i].reserve(numneighbors[i]+2); // Slop for boundaries

	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < 3; j++) {
			vector<int> &me = neighbors[faces[i][j]];
			int n1 = faces[i][(j+1)%3];
			int n2 = faces[i][(j+2)%3];
			if (find(me.begin(), me.end(), n1) == me.end())
				me.push_back(n1);
			if (find(me.begin(), me.end(), n2) == me.end())
				me.push_back(n2);
		}
	}

	dprintf("Done.\n");

	cout<<"Writing vertex neighbors.....";


	FILE *fout;
	fopen_s(&fout, dFName.getNeighborsDN(), "wb");

	unsigned int *s = new unsigned int[nv];
	fwrite(&nv, sizeof(int), 1, fout);
	for (int vertexIndex = 0; vertexIndex < nv; vertexIndex++)
		s[vertexIndex] = neighbors[vertexIndex].size();
	fwrite(s, sizeof(int), nv, fout);
	for (int vertexIndex = 0; vertexIndex < nv; vertexIndex++)
	if (s[vertexIndex])
		fwrite(&neighbors[vertexIndex][0], sizeof(int), s[vertexIndex], fout);
	delete[] s;

	fclose(fout);

	/*
	ofstream neighbors_of(dFName.getNeighborsDN());
	neighbors_of<<nv<<endl;
	for (int vertexIndex=0; vertexIndex<nv; vertexIndex++)
	{
		neighbors_of<<neighbors[vertexIndex].size()<<" ";
		
		for (int neighborIndex=0; neighborIndex<neighbors[vertexIndex].size(); neighborIndex++)
			neighbors_of<<neighbors[vertexIndex][neighborIndex]<<" ";
		neighbors_of<<endl;
	}
	neighbors_of.close();
	*/



	cout<<"Done!"<<endl;



}


// Find the faces touching each vertex
void TriMesh::need_adjacentfaces()
{
	if (!adjacentfaces.empty())
		return;

	need_faces();
	if (faces.empty())
		return;

	dprintf("Finding vertex to triangle maps... ");
	int nv = vertices.size(), nf = faces.size();

	vector<int> numadjacentfaces(nv);
	for (int i = 0; i < nf; i++) {
		numadjacentfaces[faces[i][0]]++;
		numadjacentfaces[faces[i][1]]++;
		numadjacentfaces[faces[i][2]]++;
	}

	adjacentfaces.resize(vertices.size());
	for (int i = 0; i < nv; i++)
		adjacentfaces[i].reserve(numadjacentfaces[i]);

	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < 3; j++)
			adjacentfaces[faces[i][j]].push_back(i);
	}

	dprintf("Done.\n");


}


// Find the face across each edge from each other face (-1 on boundary)
// If topology is bad, not necessarily what one would expect...
void TriMesh::need_across_edge()
{
	if (!across_edge.empty())
		return;

	need_adjacentfaces();
	if (adjacentfaces.empty())
		return;

	dprintf("Finding across-edge maps... ");

	int nf = faces.size();
	across_edge.resize(nf, Face(-1,-1,-1));

#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		for (int j = 0; j < 3; j++) {
			if (across_edge[i][j] != -1)
				continue;
			int v1 = faces[i][(j+1)%3];
			int v2 = faces[i][(j+2)%3];
			const vector<int> &a1 = adjacentfaces[v1];
			const vector<int> &a2 = adjacentfaces[v2];
			for (size_t k1 = 0; k1 < a1.size(); k1++) {
				int other = a1[k1];
				if (other == i)
					continue;
				vector<int>::const_iterator it =
					find(a2.begin(), a2.end(), other);
				if (it == a2.end())
					continue;
				int ind = (faces[other].indexof(v1)+1)%3;
				if (faces[other][(ind+1)%3] != v2)
					continue;
				across_edge[i][j] = other;
				across_edge[other][ind] = i;
				break;
			}
		}
	}

	dprintf("Done.\n");
}

}; // namespace trimesh
