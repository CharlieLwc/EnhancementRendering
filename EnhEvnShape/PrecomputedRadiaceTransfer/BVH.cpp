
#include "stdafx.h"
#include "BVH.h"

		 
BVH::BVH(TriMesh* pmesh, int iMaxNodePrim): mesh(pmesh), miMaxNodePrim(iMaxNodePrim)
{


	/*

	FILE *fin;
	fopen_s(&fin, mesh->dFName.getBVHDN(), "rb");



	if (fin)
	{
		fread(&miPrimCount, sizeof(unsigned int), 1, fin);

		if (miPrimCount == mesh->faces.size())
		{
			cout << "Reading BVH.....";
			int tri4Count;



			fread(&miMaxNodePrim, sizeof(unsigned int), 1, fin);
			fread(&miPrimCount, sizeof(unsigned int), 1, fin);
			fread(&miTotalNodesCount, sizeof(unsigned int), 1, fin);
			fread(&tri4Count, sizeof(unsigned int), 1, fin);




			Vec3f_SSE* mVertices0 = new Vec3f_SSE[tri4Count];
			Vec3f_SSE* mEdges1 = new Vec3f_SSE[tri4Count];
			Vec3f_SSE* mEdges2 = new Vec3f_SSE[tri4Count];
			Vec3f_SSE* mGeomNormals = new Vec3f_SSE[tri4Count];

			fread(mVertices0, sizeof(Vec3f_SSE), tri4Count, fin);
			fread(mEdges1, sizeof(Vec3f_SSE), tri4Count, fin);
			fread(mEdges2, sizeof(Vec3f_SSE), tri4Count, fin);
			fread(mGeomNormals, sizeof(Vec3f_SSE), tri4Count, fin);

			mvPackedTriangles.resize(tri4Count);
			for (unsigned int i = 0; i<tri4Count; i++)
			{
				Triangle4& currentT4 = mvPackedTriangles[i];

				currentT4.mVertices0 = mVertices0[i];
				currentT4.mEdges1 = mEdges1[i];
				currentT4.mEdges2 = mEdges2[i];
				currentT4.mGeomNormals = mGeomNormals[i];
			}


			delete[] mVertices0;
			delete[] mEdges1;
			delete[] mEdges2;
			delete[] mGeomNormals;


			FloatSSE* fMinMaxBoundsX = new FloatSSE[miTotalNodesCount];
			FloatSSE* fMinMaxBoundsY = new FloatSSE[miTotalNodesCount];
			FloatSSE* fMinMaxBoundsZ = new FloatSSE[miTotalNodesCount];
			uint* iSecondChildOffset = new uint[miTotalNodesCount];
			uint8* iPrimCount = new uint8[miTotalNodesCount];
			uint8* iAxis = new uint8[miTotalNodesCount];
			uint16* pad = new uint16[miTotalNodesCount];



			mpNodes = AllocAligned<BVHNode>(miTotalNodesCount);
			for (unsigned int i = 0; i<miTotalNodesCount; i++)
			{
				BVHNode* currentNode = mpNodes + i;

				currentNode->fMinMaxBoundsX[0] = 0.222f;

				currentNode->fMinMaxBoundsX = fMinMaxBoundsX[i];
				currentNode->fMinMaxBoundsY = fMinMaxBoundsY[i];
				currentNode->fMinMaxBoundsZ = fMinMaxBoundsZ[i];
				currentNode->iSecondChildOffset = iSecondChildOffset[i];
				currentNode->iPrimCount = iPrimCount[i];
				currentNode->iAxis = iAxis[i];
				currentNode->pad = pad[i];
			}

			delete[] fMinMaxBoundsX;
			delete[] fMinMaxBoundsY;
			delete[] fMinMaxBoundsZ;
			delete[] iSecondChildOffset;
			delete[] iPrimCount;
			delete[] iAxis;
			delete[] pad;


			cout << "Done" << endl;
			fclose(fin);

			return;
		}
		fclose(fin);
	}
	*/

	
	ifstream BVH_if(mesh->dFName.getBVHDN());
	if (BVH_if.is_open())
	{
		BVH_if>>miPrimCount;
		if (miPrimCount == mesh->faces.size())
		{
			cout<<"Reading BVH.....";
			int tri4Count;
			BVH_if>>miMaxNodePrim;
			BVH_if>>miPrimCount;
			BVH_if>>miTotalNodesCount;
			BVH_if>>tri4Count;
			mvPackedTriangles.resize(tri4Count);


			mpNodes = AllocAligned<BVHNode>(miTotalNodesCount);

			for(unsigned int i = 0; i<tri4Count; i++)
			{
				Triangle4& currentT4 = mvPackedTriangles[i];
				BVH_if >> currentT4.mVertices0;
				BVH_if >> currentT4.mEdges1;
				BVH_if >> currentT4.mEdges2;
				BVH_if >> currentT4.mGeomNormals;
			}


			for(unsigned int i = 0; i<miTotalNodesCount; i++)
			{
				BVHNode* currentNode = mpNodes + i;

				currentNode->fMinMaxBoundsX[0] = 0.222f;

				BVH_if >> currentNode->fMinMaxBoundsX;
				BVH_if >> currentNode->fMinMaxBoundsY;
				BVH_if >> currentNode->fMinMaxBoundsZ;
				BVH_if >> currentNode->iSecondChildOffset;
				BVH_if >> currentNode->iPrimCount;
				BVH_if >> currentNode->iAxis;
				BVH_if >> currentNode->pad;
			}		   


			cout<<"Done"<<endl;
			BVH_if.close();
			
			return;
		}
		BVH_if.close();
	}


	cout<<"Computing BVH... ";
	
	
	vector<TriMesh::Face>& faces = pmesh->faces;
	vector<vec3>& vertiese = pmesh->vertices;

	mvPrimitives.resize(faces.size());

	for (int i=0; i<mvPrimitives.size(); i++)
	{
		int* vIndex = faces[i].v;
		mvPrimitives[i].init(vertiese[vIndex[0]], vertiese[vIndex[1]], vertiese[vIndex[2]]);
	}



	miPrimCount = mvPrimitives.size();

	vector<BVHPrimInfo> vBuildData;
	vBuildData.reserve(mvPrimitives.size());
	for(int i = 0; i < mvPrimitives.size(); i++)
	{
		BoundingBox bbox = mvPrimitives[i].WorldBound();
		vBuildData.push_back(BVHPrimInfo(i, bbox));
		mWorldBound = Math::Union(mWorldBound, bbox);
	}

	MemoryArena memory;
	miTotalNodesCount = 0;
	mvPackedTriangles.clear();
	mvPackedTriangles.reserve((mvPrimitives.size() + 3) / 4);
	BVHBuildNode* pRoot = Construction(vBuildData, miTotalNodesCount, 0, miPrimCount, memory);


	mpNodes = AllocAligned<BVHNode>(miTotalNodesCount);
	uint iOffset = 0;
	FlattenTree(pRoot, &iOffset);
	assert(iOffset == miTotalNodesCount);


	cout<<"Done"<<endl;

	/*
	FILE *fout;
	fopen_s(&fout, mesh->dFName.getBVHDN(), "wb");

	unsigned int tri4Count = mvPackedTriangles.size();
	cout << "Writing BVH.....";

	fwrite(&miPrimCount, sizeof(unsigned int), 1, fout);
	fwrite(&miMaxNodePrim, sizeof(unsigned int), 1, fout);
	fwrite(&miPrimCount, sizeof(unsigned int), 1, fout);
	fwrite(&miTotalNodesCount, sizeof(unsigned int), 1, fout);
	fwrite(&tri4Count, sizeof(unsigned int), 1, fout);


	Vec3f_SSE* mVertices0 = new Vec3f_SSE[tri4Count];
	Vec3f_SSE* mEdges1 = new Vec3f_SSE[tri4Count];
	Vec3f_SSE* mEdges2 = new Vec3f_SSE[tri4Count];
	Vec3f_SSE* mGeomNormals = new Vec3f_SSE[tri4Count];


	for (unsigned int i = 0; i<tri4Count; i++)
	{
		Triangle4& currentT4 = mvPackedTriangles[i];

		mVertices0[i] = currentT4.mVertices0;
		mEdges1[i] = currentT4.mEdges1;
		mEdges2[i] = currentT4.mEdges2;
		mGeomNormals[i] = currentT4.mGeomNormals;
	}

	fwrite(mVertices0, sizeof(Vec3f_SSE), tri4Count, fout);
	fwrite(mEdges1, sizeof(Vec3f_SSE), tri4Count, fout);
	fwrite(mEdges2, sizeof(Vec3f_SSE), tri4Count, fout);
	fwrite(mGeomNormals, sizeof(Vec3f_SSE), tri4Count, fout);

	delete[] mVertices0;
	delete[] mEdges1;
	delete[] mEdges2;
	delete[] mGeomNormals;


	FloatSSE* fMinMaxBoundsX = new FloatSSE[miTotalNodesCount];
	FloatSSE* fMinMaxBoundsY = new FloatSSE[miTotalNodesCount];
	FloatSSE* fMinMaxBoundsZ = new FloatSSE[miTotalNodesCount];
	uint* iSecondChildOffset = new uint[miTotalNodesCount];
	uint8* iPrimCount = new uint8[miTotalNodesCount];
	uint8* iAxis = new uint8[miTotalNodesCount];
	uint16* pad = new uint16[miTotalNodesCount];


	for (unsigned int i = 0; i<miTotalNodesCount; i++)
	{
		BVHNode& currentNode = mpNodes[i];

		fMinMaxBoundsX[i] = currentNode.fMinMaxBoundsX;
		fMinMaxBoundsY[i] = currentNode.fMinMaxBoundsY;
		fMinMaxBoundsZ[i] = currentNode.fMinMaxBoundsZ;
		iSecondChildOffset[i] = currentNode.iSecondChildOffset;
		iPrimCount[i] = currentNode.iPrimCount;
		iAxis[i] = currentNode.iAxis;
		pad[i] = currentNode.pad;
	}

	fwrite(fMinMaxBoundsX, sizeof(FloatSSE), miTotalNodesCount, fout);
	fwrite(fMinMaxBoundsY, sizeof(FloatSSE), miTotalNodesCount, fout);
	fwrite(fMinMaxBoundsZ, sizeof(FloatSSE), miTotalNodesCount, fout);
	fwrite(iSecondChildOffset, sizeof(uint), miTotalNodesCount, fout);
	fwrite(iPrimCount, sizeof(uint8), miTotalNodesCount, fout);
	fwrite(iAxis, sizeof(uint8), miTotalNodesCount, fout);
	fwrite(pad, sizeof(uint16), miTotalNodesCount, fout);

	delete[] fMinMaxBoundsX;
	delete[] fMinMaxBoundsY;
	delete[] fMinMaxBoundsZ;
	delete[] iSecondChildOffset;
	delete[] iPrimCount;
	delete[] iAxis;
	delete[] pad;

	fclose(fout);
	*/


	ofstream BVH_of(mesh->dFName.getBVHDN());

	int tri4Count = mvPackedTriangles.size();
	cout<<"Writing BVH.....";
	BVH_of<<miPrimCount<<" "<<
		miMaxNodePrim<<" "<<
		miPrimCount<<" "<<
		miTotalNodesCount<<" "<<
		tri4Count<<endl;
	
	for(unsigned int i = 0; i<tri4Count; i++)
	{
		Triangle4& currentT4 = mvPackedTriangles[i];
		BVH_of << currentT4.mVertices0<<endl;
		BVH_of << currentT4.mEdges1<<endl;
		BVH_of << currentT4.mEdges2<<endl;
		BVH_of << currentT4.mGeomNormals<<endl;
	}


	for(unsigned int i = 0; i<miTotalNodesCount; i++)
	{
		BVHNode& currentNode = mpNodes[i];

		BVH_of << currentNode.fMinMaxBoundsX<<endl;
		BVH_of << currentNode.fMinMaxBoundsY<<endl;
		BVH_of << currentNode.fMinMaxBoundsZ<<endl;
		BVH_of << currentNode.iSecondChildOffset<<endl;
		BVH_of << currentNode.iPrimCount<<endl;
		BVH_of << currentNode.iAxis<<endl;
		BVH_of << currentNode.pad<<endl;
	}


	cout<<"Done"<<endl;
	

	BVH_of.close();

	vertexOccludedOutside = NULL;
	vertexOccludedInside = NULL;



}

BVHBuildNode* BVH::Construction(vector<BVHPrimInfo>& vBuildInfo,
								uint& iTotalNodes,
								uint iStart, uint iEnd,
								MemoryArena& memory)
{
	assert(iStart < iEnd);
	iTotalNodes++; // Increment total node count

	// Allocate memory from arena
	BVHBuildNode* pNode = memory.Alloc<BVHBuildNode>();
	*pNode = BVHBuildNode();

	// Compute bounding box for the current node
	BoundingBox nodeBBox;
	for(auto i = iStart; i < iEnd; i++)
		nodeBBox = Math::Union(nodeBBox, vBuildInfo[i].bbox);

	uint iNodePrimCount = iEnd - iStart;
	auto InitLeaf = [&] // Lambda for leaf creation
	{
		uint iFirstPrimOffset = mvPackedTriangles.size();
		uint iPackedCount = (iNodePrimCount + 3) / 4;
		uint iInfoIdx = iStart;
		for(auto i = 0; i < iPackedCount; i++)
		{
			Triangle aTris[4];
			uint iCount = 0;
			for(; iCount < 4; iCount++)
			{
				if(iInfoIdx >= iEnd)
					break;

				aTris[iCount] = mvPrimitives[vBuildInfo[iInfoIdx++].iPrimIdx];
			}

			mvPackedTriangles.push_back(Triangle4((Triangle*)aTris, iCount));
		}
		pNode->InitLeaf(iFirstPrimOffset, iPackedCount, nodeBBox);
	};

	if(iNodePrimCount <= miMaxNodePrim) // If current node has only 1 primitive, create leaf node
	{
		InitLeaf();
	}
	else // Split the node
	{
		// Get the centroid bounds of current node's primitive
		BoundingBox centroidBounds;
		for(auto i = iStart; i < iEnd; i++)
			centroidBounds = Math::Union(centroidBounds, vBuildInfo[i].vCentroid);
		// Get split dimension
		int iDim = centroidBounds.MaximumExtent();
		// Create leaf node if maximum extent is zero
		if(centroidBounds.mptMin[iDim] == centroidBounds.mptMax[iDim])
		{
			InitLeaf();
			return pNode;
		}

		// Partition primitives
		uint iMid = 0;
		auto EqualCountSplit = [&]
		{
			iMid = (iStart + iEnd) / 2;
			std::nth_element(&vBuildInfo[iStart], &vBuildInfo[iMid], &vBuildInfo[iEnd - 1] + 1,
				[iDim](const BVHPrimInfo& lhs, const BVHPrimInfo& rhs)
			{
				return lhs.vCentroid[iDim] < rhs.vCentroid[iDim];
			});
		};

		// 			if(iNodePrimCount <= 4)
		// 			{
		// 				EqualCountSplit();
		// 			}
		// 			else // Apply SAH
		{
			struct Bucket
			{
				int iCount;
				BoundingBox bounds;
				Bucket() : iCount(0) {}
			};
			const int iBucketCount = 128;
			Bucket aBuckets[iBucketCount];

			// Get primitive count and bounds for each bucket
			for(auto i = iStart; i < iEnd; i++)
			{
				int iBukIdx = iBucketCount * ((vBuildInfo[i].vCentroid[iDim] - nodeBBox.mptMin[iDim]) /
					(nodeBBox.mptMax[iDim] - nodeBBox.mptMin[iDim]));

				iBukIdx = Math::Min(iBukIdx, iBucketCount - 1);
				aBuckets[iBukIdx].iCount++;
				aBuckets[iBukIdx].bounds = Math::Union(aBuckets[iBukIdx].bounds, vBuildInfo[i].bbox);
			}

			float afCost[iBucketCount - 1] = {0};
			float fMinCost = Math::EDX_INFINITY;
			int iBestSplit = -1;
			for(auto i = 0; i < iBucketCount - 1; i++)
			{
				BoundingBox boundsLow, boundsHigh;
				uint iCountLow = 0, iCountHigh = 0;
				for(auto j = 0; j <= i; j++)
				{
					iCountLow += aBuckets[j].iCount;
					boundsLow = Math::Union(boundsLow, aBuckets[j].bounds);
				}
				for(auto j = i + 1; j < iBucketCount; j++)
				{
					iCountHigh += aBuckets[j].iCount;
					boundsHigh = Math::Union(boundsHigh, aBuckets[j].bounds);
				}

				afCost[i] = 0.5f + (iCountLow * boundsLow.Area() + iCountHigh * boundsHigh.Area()) / nodeBBox.Area();
				if(afCost[i] < fMinCost)
				{
					fMinCost = afCost[i];
					iBestSplit = i;
				}
			}

			if(fMinCost < iNodePrimCount)
			{
				BVHPrimInfo* pMid = std::partition(&vBuildInfo[iStart], &vBuildInfo[iEnd - 1] + 1,
					[&](const BVHPrimInfo& par) -> bool
				{
					int iBukIdx = iBucketCount * ((par.vCentroid[iDim] - nodeBBox.mptMin[iDim]) /
						(nodeBBox.mptMax[iDim] - nodeBBox.mptMin[iDim]));

					iBukIdx = Math::Min(iBukIdx, iBucketCount - 1);
					return iBukIdx <= iBestSplit;
				});

				iMid = pMid - &vBuildInfo[0];
			}
			else
			{
				InitLeaf();
				return pNode;
			}

		}

		pNode->InitInterior(iDim,
			Construction(vBuildInfo, iTotalNodes, iStart, iMid, memory),
			Construction(vBuildInfo, iTotalNodes, iMid, iEnd, memory));
	}

	return pNode;
}

uint BVH::FlattenTree(BVHBuildNode* pBuildNode, uint* piOffset)
{
	BVHNode* pNode = mpNodes + *piOffset;

	pNode->fMinMaxBoundsX[0] = pBuildNode->boundsLeft.mptMin.x;
	pNode->fMinMaxBoundsX[1] = pBuildNode->boundsRight.mptMin.x;
	pNode->fMinMaxBoundsX[2] = pBuildNode->boundsLeft.mptMax.x;
	pNode->fMinMaxBoundsX[3] = pBuildNode->boundsRight.mptMax.x;

	pNode->fMinMaxBoundsY[0] = pBuildNode->boundsLeft.mptMin.y;
	pNode->fMinMaxBoundsY[1] = pBuildNode->boundsRight.mptMin.y;
	pNode->fMinMaxBoundsY[2] = pBuildNode->boundsLeft.mptMax.y;
	pNode->fMinMaxBoundsY[3] = pBuildNode->boundsRight.mptMax.y;

	pNode->fMinMaxBoundsZ[0] = pBuildNode->boundsLeft.mptMin.z;
	pNode->fMinMaxBoundsZ[1] = pBuildNode->boundsRight.mptMin.z;
	pNode->fMinMaxBoundsZ[2] = pBuildNode->boundsLeft.mptMax.z;
	pNode->fMinMaxBoundsZ[3] = pBuildNode->boundsRight.mptMax.z;

	uint iCurrOffset = (*piOffset)++;

	if(pBuildNode->iPrimCount > 0) // Create leaf node
	{
		pNode->iPrimOffset = pBuildNode->iFirstPrimOffset;
		pNode->iPrimCount = pBuildNode->iPrimCount;
	}
	else // Interior
	{
		pNode->iPrimCount = 0;
		pNode->iAxis = pBuildNode->iSplitAxis;
		FlattenTree(pBuildNode->pChildren[0], piOffset);
		pNode->iSecondChildOffset = FlattenTree(pBuildNode->pChildren[1], piOffset);
	}

	return iCurrOffset;
}

bool BVH::Intersect(const Ray& ray, Intersection* pIsect) const
{
	const IntSSE identity = _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2, 1, 0);
	const IntSSE swap     = _mm_set_epi8( 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10, 9, 8);
	const IntSSE shuffleX = ray.mvDir.x >= 0 ? identity : swap;
	const IntSSE shuffleY = ray.mvDir.y >= 0 ? identity : swap;
	const IntSSE shuffleZ = ray.mvDir.z >= 0 ? identity : swap;

	const IntSSE pn = IntSSE(0x00000000,0x00000000,0x80000000,0x80000000);
	const Vec3f_SSE norg(-ray.mptOrg.x,-ray.mptOrg.y,-ray.mptOrg.z);
	const Vec3f_SSE rdir = Vec3f_SSE(FloatSSE(1.0f / ray.mvDir.x) ^ pn, FloatSSE(1.0f / ray.mvDir.y) ^ pn, FloatSSE(1.0f / ray.mvDir.z) ^ pn);
	FloatSSE nearFar(ray.mfMin, ray.mfMin, -ray.mfMax, -ray.mfMax);

	struct TravStackItem
	{
		float fDist;
		int iIndex;
	};
	TravStackItem TodoStack[64];
	uint iStackTop = 0, iNodeIdx = 0;

	bool bHit = false;
	while(true)
	{
		const BVHNode* pNode = &mpNodes[iNodeIdx];

		// Interior node
		if(pNode->iPrimCount == 0)
		{
			const FloatSSE tNearFarX = (SSE::Shuffle8(pNode->fMinMaxBoundsX, shuffleX) + norg.x) * rdir.x;
			const FloatSSE tNearFarY = (SSE::Shuffle8(pNode->fMinMaxBoundsY, shuffleY) + norg.y) * rdir.y;
			const FloatSSE tNearFarZ = (SSE::Shuffle8(pNode->fMinMaxBoundsZ, shuffleZ) + norg.z) * rdir.z;
			const FloatSSE tNearFar = SSE::Max(SSE::Max(tNearFarX, tNearFarY), SSE::Max(tNearFarZ, nearFar)) ^ pn;
			const BoolSSE lrhit = tNearFar <= SSE::Shuffle8(tNearFar, swap);

			if(lrhit[0] != 0 && lrhit[1] != 0)
			{
				if(tNearFar[0] < tNearFar[1]) // First child first
				{
					TodoStack[iStackTop].iIndex = pNode->iSecondChildOffset;
					TodoStack[iStackTop].fDist = tNearFar[1];
					iStackTop++;
					iNodeIdx = iNodeIdx + 1;
				}
				else
				{
					TodoStack[iStackTop].iIndex = iNodeIdx + 1;
					TodoStack[iStackTop].fDist = tNearFar[0];
					iStackTop++;
					iNodeIdx = pNode->iSecondChildOffset;
				}
			}
			else if(lrhit[0] != 0)
			{
				iNodeIdx = iNodeIdx + 1;
			}
			else if(lrhit[1] != 0)
			{
				iNodeIdx = pNode->iSecondChildOffset;
			}
			else // If miss the node's bounds
			{
				do 
				{
					if(iStackTop == 0)
						return bHit;
					iStackTop--;
					iNodeIdx = TodoStack[iStackTop].iIndex;
				} while (TodoStack[iStackTop].fDist > pIsect->mfDist);
			}
		}
		else // Leaf node
		{
			for(auto i = 0; i < pNode->iPrimCount; i++)
			{
				if(mvPackedTriangles[pNode->iPrimOffset + i].Intersect(ray, pIsect))
					bHit = true;
			}
			nearFar = SSE::Shuffle<0,1,2,3>(nearFar, -pIsect->mfDist);

			do 
			{
				if(iStackTop == 0)
					return bHit;
				iStackTop--;
				iNodeIdx = TodoStack[iStackTop].iIndex;
			} while (TodoStack[iStackTop].fDist > pIsect->mfDist);
		}
	}

	return bHit;
}

bool BVH::Occluded(const Ray& ray) const
{
	const IntSSE identity = _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2, 1, 0);
	const IntSSE swap     = _mm_set_epi8( 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10, 9, 8);
	const IntSSE shuffleX = ray.mvDir.x >= 0 ? identity : swap;
	const IntSSE shuffleY = ray.mvDir.y >= 0 ? identity : swap;
	const IntSSE shuffleZ = ray.mvDir.z >= 0 ? identity : swap;

	const IntSSE pn = IntSSE(0x00000000,0x00000000,0x80000000,0x80000000);
	const Vec3f_SSE norg(-ray.mptOrg.x,-ray.mptOrg.y,-ray.mptOrg.z);
	const Vec3f_SSE rdir = Vec3f_SSE(FloatSSE(1.0f / ray.mvDir.x) ^ pn, FloatSSE(1.0f / ray.mvDir.y) ^ pn, FloatSSE(1.0f / ray.mvDir.z) ^ pn);
	FloatSSE nearFar(ray.mfMin, ray.mfMin, -ray.mfMax, -ray.mfMax);

	uint aiTodoStack[64];
	uint iStackTop = 0, iNodeIdx = 0;

	while(true)
	{
		const BVHNode* pNode = &mpNodes[iNodeIdx];

		// If is a leaf node
		if(pNode->iPrimCount == 0)
		{
			const FloatSSE tNearFarX = (SSE::Shuffle8(pNode->fMinMaxBoundsX, shuffleX) + norg.x) * rdir.x;
			const FloatSSE tNearFarY = (SSE::Shuffle8(pNode->fMinMaxBoundsY, shuffleY) + norg.y) * rdir.y;
			const FloatSSE tNearFarZ = (SSE::Shuffle8(pNode->fMinMaxBoundsZ, shuffleZ) + norg.z) * rdir.z;
			const FloatSSE tNearFar = SSE::Max(SSE::Max(tNearFarX, tNearFarY), SSE::Max(tNearFarZ, nearFar)) ^ pn;
			const BoolSSE lrhit = tNearFar <= SSE::Shuffle8(tNearFar, swap);

			if(lrhit[0] != 0 && lrhit[1] != 0)
			{
				if(tNearFar[0] < tNearFar[1]) // First child first
				{
					aiTodoStack[iStackTop++] = pNode->iSecondChildOffset;
					iNodeIdx = iNodeIdx + 1;
				}
				else
				{
					aiTodoStack[iStackTop++] = iNodeIdx + 1;
					iNodeIdx = pNode->iSecondChildOffset;
				}
			}
			else if(lrhit[0] != 0)
			{
				iNodeIdx = iNodeIdx + 1;
			}
			else if(lrhit[1] != 0)
			{
				iNodeIdx = pNode->iSecondChildOffset;
			}
			else // If miss the node's bounds
			{
				if(iStackTop == 0)
					break;
				iNodeIdx = aiTodoStack[--iStackTop];
			}
		}
		else // Interior node
		{
			for(auto i = 0; i < pNode->iPrimCount; i++)
			{
				if(mvPackedTriangles[pNode->iPrimOffset + i].Occluded(ray))
					return true;
			}

			if(iStackTop == 0)
				break;
			iNodeIdx = aiTodoStack[--iStackTop];
		}
	}

	return false;
}

vector<uint> BVH::IntersectDebug(const Ray& ray, Intersection* pIsect) const
{
	const IntSSE identity = _mm_set_epi8(15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2, 1, 0);
	const IntSSE swap     = _mm_set_epi8( 7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10, 9, 8);
	const IntSSE shuffleX = ray.mvDir.x >= 0 ? identity : swap;
	const IntSSE shuffleY = ray.mvDir.y >= 0 ? identity : swap;
	const IntSSE shuffleZ = ray.mvDir.z >= 0 ? identity : swap;

	const IntSSE pn = IntSSE(0x00000000,0x00000000,0x80000000,0x80000000);
	const Vec3f_SSE norg(-ray.mptOrg.x,-ray.mptOrg.y,-ray.mptOrg.z);
	const Vec3f_SSE rdir = Vec3f_SSE(FloatSSE(1.0f / ray.mvDir.x) ^ pn, FloatSSE(1.0f / ray.mvDir.y) ^ pn, FloatSSE(1.0f / ray.mvDir.z) ^ pn);
	FloatSSE nearFar(ray.mfMin, ray.mfMin, -ray.mfMax, -ray.mfMax);

	struct TravStackItem
	{
		float fDist;
		int iIndex;
	};
	TravStackItem TodoStack[64];
	uint iStackTop = 0, iNodeIdx = 0;

	vector<uint> vIndices;
	while(true)
	{
		const BVHNode* pNode = &mpNodes[iNodeIdx];
		vIndices.push_back(iNodeIdx);

		// Interior node
		if(pNode->iPrimCount == 0)
		{
			const FloatSSE tNearFarX = (SSE::Shuffle8(pNode->fMinMaxBoundsX, shuffleX) + norg.x) * rdir.x;
			const FloatSSE tNearFarY = (SSE::Shuffle8(pNode->fMinMaxBoundsY, shuffleY) + norg.y) * rdir.y;
			const FloatSSE tNearFarZ = (SSE::Shuffle8(pNode->fMinMaxBoundsZ, shuffleZ) + norg.z) * rdir.z;
			const FloatSSE tNearFar = SSE::Max(SSE::Max(tNearFarX, tNearFarY), SSE::Max(tNearFarZ, nearFar)) ^ pn;
			const BoolSSE lrhit = tNearFar <= SSE::Shuffle8(tNearFar, swap);

			if(lrhit[0] != 0 && lrhit[1] != 0)
			{
				if(tNearFar[0] < tNearFar[1]) // First child first
				{
					TodoStack[iStackTop].iIndex = pNode->iSecondChildOffset;
					TodoStack[iStackTop].fDist = tNearFar[1];
					iStackTop++;
					iNodeIdx = iNodeIdx + 1;
				}
				else
				{
					TodoStack[iStackTop].iIndex = iNodeIdx + 1;
					TodoStack[iStackTop].fDist = tNearFar[0];
					iStackTop++;
					iNodeIdx = pNode->iSecondChildOffset;
				}
			}
			else if(lrhit[0] != 0)
			{
				iNodeIdx = iNodeIdx + 1;
			}
			else if(lrhit[1] != 0)
			{
				iNodeIdx = pNode->iSecondChildOffset;
			}
			else // If miss the node's bounds
			{
				do 
				{
					if(iStackTop == 0)
						return vIndices;
					iStackTop--;
					iNodeIdx = TodoStack[iStackTop].iIndex;
				} while (TodoStack[iStackTop].fDist > pIsect->mfDist);
			}
		}
		else // Leaf node
		{
			do 
			{
				if(iStackTop == 0)
					return vIndices;
				iStackTop--;
				iNodeIdx = TodoStack[iStackTop].iIndex;
			} while (TodoStack[iStackTop].fDist > pIsect->mfDist);
		}
	}

	return vIndices;
}

bool BVH::IsVertVertOutsideOccluded(int startVertexIndex, int endtVertexIndex)
{

	if (vertexOccludedOutside == NULL)
	{
		int vertexNum;
		ifstream VertexOcl_if(mesh->dFName.getVertexVertexOutsideOclDN());
		if (VertexOcl_if.is_open())
		{
			VertexOcl_if>>vertexNum;
			if (vertexNum == mesh->vertices.size())
			{
				cout<<"Reading vertices visibility...."<<endl;
				int totalNum = vertexNum*(vertexNum-1)/2;
				vertexOccludedOutside = new bool[totalNum];
				for (int i=0; i<totalNum; i++)
					VertexOcl_if >> vertexOccludedOutside[i];
			}
			VertexOcl_if.close();
		}

		if (vertexOccludedOutside == NULL)
		{

			cout<<"Calculating vertices visibility...."<<endl;
			vertexNum = mesh->vertices.size();
			int totalNum = vertexNum*(vertexNum-1)/2;
			vertexOccludedOutside = new bool[totalNum];

			vector<vec3> &Normals = mesh->normals;
			vector<vec3> &Positions = mesh->vertices;
			int currentVInd = 0;
			for (int sVInd=1; sVInd<vertexNum; sVInd++)
			{
				vec3 &sVNor = Normals[sVInd];
				vec3 &sVPos = Positions[sVInd];
				for (int eVInd=0; eVInd<sVInd; eVInd++)
				{
					vec3 &eVNor = Normals[eVInd];
					vec3 &eVPos = Positions[eVInd];
					vec3 lightDir = eVPos - sVPos;
					if (sVNor.dot(lightDir)<=0.f || eVNor.dot(-lightDir) <= 0.f)
						vertexOccludedOutside[currentVInd++] = true;
					else
					{
						float max = len(lightDir) - 4*EPSILON;
						normalize(lightDir);
						Ray ray(sVPos+2*EPSILON*sVNor,lightDir, max);
						vertexOccludedOutside[currentVInd++] = Occluded(ray);
					}

				}
			}

			cout<<"Writing vertices visibility...."<<endl;
			ofstream VertexOcl_of(mesh->dFName.getVertexVertexOutsideOclDN());
			VertexOcl_of<<vertexNum<<endl;
			for (int i=0; i<totalNum; i++)
				VertexOcl_of << vertexOccludedOutside[i] << " ";
			VertexOcl_of.close();
		}
	}

	
	if (startVertexIndex == endtVertexIndex)
		return true;

	if (startVertexIndex<endtVertexIndex)
	{
		int temp = startVertexIndex;
		startVertexIndex = endtVertexIndex;
		endtVertexIndex = temp;
	}

	return vertexOccludedOutside[startVertexIndex*(startVertexIndex - 1) / 2 + endtVertexIndex];


}

bool BVH::IsVertVertInsideOccluded(int startVertexIndex, int endtVertexIndex)
{

	if (vertexOccludedInside == NULL)
	{
		int vertexNum;
		ifstream VertexOcl_if(mesh->dFName.getVertexVertexInsideOclDN());
		if (VertexOcl_if.is_open())
		{
			VertexOcl_if >> vertexNum;
			if (vertexNum == mesh->vertices.size())
			{
				cout << "Reading vertices visibility...." << endl;
				int totalNum = vertexNum*(vertexNum - 1) / 2;
				vertexOccludedInside = new bool[totalNum];
				for (int i = 0; i<totalNum; i++)
					VertexOcl_if >> vertexOccludedInside[i];
			}
			VertexOcl_if.close();
		}

		if (vertexOccludedInside == NULL)
		{

			cout << "Calculating vertices visibility...." << endl;
			vertexNum = mesh->vertices.size();
			int totalNum = vertexNum*(vertexNum - 1) / 2;
			vertexOccludedInside = new bool[totalNum];

			vector<vec3> &Normals = mesh->normals;
			vector<vec3> &Positions = mesh->vertices;
			int currentVInd = 0;
			for (int sVInd = 1; sVInd<vertexNum; sVInd++)
			{
				vec3 &sVNor = Normals[sVInd];
				vec3 &sVPos = Positions[sVInd];
				for (int eVInd = 0; eVInd<sVInd; eVInd++)
				{
					vec3 &eVNor = Normals[eVInd];
					vec3 &eVPos = Positions[eVInd];
					vec3 lightDir = eVPos - sVPos;
					if (sVNor.dot(lightDir) >= 0.f || eVNor.dot(-lightDir) >= 0.f)
						vertexOccludedInside[currentVInd++] = true;
					else
					{
						float max = len(lightDir) - 4 * EPSILON;
						normalize(lightDir);
						Ray ray(sVPos + 2 * EPSILON*lightDir, lightDir, max);
						vertexOccludedInside[currentVInd++] = Occluded(ray);
					}

				}
			}

			cout << "Writing vertices visibility...." << endl;
			ofstream VertexOcl_of(mesh->dFName.getVertexVertexInsideOclDN());
			VertexOcl_of << vertexNum << endl;
			for (int i = 0; i<totalNum; i++)
				VertexOcl_of << vertexOccludedInside[i] << " ";
			VertexOcl_of.close();
		}
	}


	if (startVertexIndex == endtVertexIndex)
		return true;

	if (startVertexIndex<endtVertexIndex)
	{
		int temp = startVertexIndex;
		startVertexIndex = endtVertexIndex;
		endtVertexIndex = temp;
	}

	return vertexOccludedInside[startVertexIndex*(startVertexIndex - 1) / 2 + endtVertexIndex];


}

float BVH::ShortesDistanceInside(int startVertexIndex, int endtVertexIndex)
{

	if (vertexDistanceInside == NULL)
	{
		int vertexNum;
		ifstream VertexOcl_if(mesh->dFName.getVertexVertexInsideDisDN());
		if (VertexOcl_if.is_open())
		{
			VertexOcl_if >> vertexNum;
			if (vertexNum == mesh->vertices.size())
			{
				cout << "Reading vertices visibility...." << endl;
				int totalNum = vertexNum*(vertexNum - 1) / 2;
				vertexDistanceInside = new float[totalNum];
				for (int i = 0; i<totalNum; i++)
					VertexOcl_if >> vertexDistanceInside[i];
			}
			VertexOcl_if.close();
		}

		if (vertexDistanceInside == NULL)
		{

			/*
			cout << "Calculating vertices distance...." << endl;


			CRichModel CRmesh(mesh->dFName.getObjDN());

			cout << "----------------------model info begin----------------\n";
			cout << "File name:\t" << CRmesh.GetFileName() << endl;
			try
			{
				CRmesh.LoadModel();
			}
			catch (const char* msg)
			{
				cout << "ERRORS happen!\n" << msg << endl;
				getchar();
			}

			CRmesh.Preprocess();

			if (!CRmesh.HasBeenLoad() || !CRmesh.HasBeenProcessed())
			{
				cout << "The model fails to be handled." << endl;
				getchar();
			}


			cout << "Face number:\t" << CRmesh.GetNumOfFaces() << endl;
			cout << "Vertex number:\t" << CRmesh.GetNumOfVerts() << endl;
			cout << "Edge number:\t" << CRmesh.GetNumOfEdges() << endl;
			if (CRmesh.IsClosedModel())
				cout << "It is a closed model." << endl;
			else
				cout << "It is an open model with " << CRmesh.GetNumOfHoles() << " holes." << endl;
			if (CRmesh.GetNumOfIsolated() > 1)
			{
				cout << "Perhaps it is composed of several components." << endl;
			}
			cout << "----------------------model info end------------------\n";



			CExactMethodForDGP * algorithm = new CXinWangImprovedCH(CRmesh, 0);

			algorithm->Execute();


			*/



			vertexNum = mesh->vertices.size();
			int totalNum = vertexNum*(vertexNum - 1) / 2;
			vertexDistanceInside = new float[totalNum];

			vector<vec3> &Normals = mesh->normals;
			vector<vec3> &Positions = mesh->vertices;
			int currentVInd = 0;
			for (int sVInd = 1; sVInd<vertexNum; sVInd++)
			{
				vec3 &sVNor = Normals[sVInd];
				vec3 &sVPos = Positions[sVInd];
				for (int eVInd = 0; eVInd<sVInd; eVInd++)
				{
					vec3 &eVNor = Normals[eVInd];
					vec3 &eVPos = Positions[eVInd];
					vec3 lightDir = eVPos - sVPos;
					if (sVNor.dot(lightDir) >= 0.f || eVNor.dot(-lightDir) >= 0.f)
						vertexDistanceInside[currentVInd++] = true;
					else
					{
						float max = len(lightDir) - 4 * EPSILON;
						normalize(lightDir);
						Ray ray(sVPos + 2 * EPSILON*lightDir, lightDir, max);
						vertexDistanceInside[currentVInd++] = Occluded(ray);
					}

				}
			}

			cout << "Writing vertices visibility...." << endl;
			ofstream VertexOcl_of(mesh->dFName.getVertexVertexInsideOclDN());
			VertexOcl_of << vertexNum << endl;
			for (int i = 0; i<totalNum; i++)
				VertexOcl_of << vertexDistanceInside[i] << " ";
			VertexOcl_of.close();
		}
	}


	if (startVertexIndex == endtVertexIndex)
		return true;

	if (startVertexIndex<endtVertexIndex)
	{
		int temp = startVertexIndex;
		startVertexIndex = endtVertexIndex;
		endtVertexIndex = temp;
	}

	return vertexDistanceInside[startVertexIndex*(startVertexIndex - 1) / 2 + endtVertexIndex];


}

void BVHVisualizer::Initialize(const BVH* pBVH)
{
	mpBVH = pBVH;

	const BVHNode* pRoot = pBVH->GetRootPtr();
	for(auto i = 0; i < pBVH->miTotalNodesCount; i++)
	{
		BoundingBox boundsLeft, boundsRight;

		boundsLeft.mptMin.x = pRoot[i].fMinMaxBoundsX[0];
		boundsRight.mptMin.x = pRoot[i].fMinMaxBoundsX[1];
		boundsLeft.mptMax.x = pRoot[i].fMinMaxBoundsX[2];
		boundsRight.mptMax.x = pRoot[i].fMinMaxBoundsX[3];

		boundsLeft.mptMin.y = pRoot[i].fMinMaxBoundsY[0];
		boundsRight.mptMin.y = pRoot[i].fMinMaxBoundsY[1];
		boundsLeft.mptMax.y = pRoot[i].fMinMaxBoundsY[2];
		boundsRight.mptMax.y = pRoot[i].fMinMaxBoundsY[3];

		boundsLeft.mptMin.z = pRoot[i].fMinMaxBoundsZ[0];
		boundsRight.mptMin.z = pRoot[i].fMinMaxBoundsZ[1];
		boundsLeft.mptMax.z = pRoot[i].fMinMaxBoundsZ[2];
		boundsRight.mptMax.z = pRoot[i].fMinMaxBoundsZ[3];

		mvBBoxList.push_back(Math::Union(boundsLeft, boundsRight));
	}
}

void BVHVisualizer::GenRayTraverseIndices(const Ray& ray)
{
	Intersection isect;
	mvRayTravIndices = mpBVH->IntersectDebug(ray, &isect);
}



