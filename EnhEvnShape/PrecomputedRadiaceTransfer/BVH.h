#pragma once

#include "stdafx.h"
#include "assert.h"
#include "SSE.h"


using namespace EDX;

#define RAY_EPSI 1e-3f


///////////////////////////////////////////////////////////////////////
//                               Ray                                 //
///////////////////////////////////////////////////////////////////////

namespace EDX
{
#define RAY_EPSI 1e-3f
	class Ray
	{
	public:
		Vector3 mptOrg;
		Vector3 mvDir;

		mutable float mfMin, mfMax;
		int miDepth;

	public:
		Ray()
			: mptOrg(Vector3::ZERO), mvDir(Vector3::UNIT_Z), 
			miDepth(0), mfMin(RAY_EPSI), mfMax(float(Math::EDX_INFINITY)) {}
		Ray(const Vector3& ptOrig, const Vector3& vDir, float fMax = float(Math::EDX_INFINITY), float fMin = 0.0f, int iDepth = 0)
			: mptOrg(ptOrig), mvDir(vDir), miDepth(iDepth), mfMin(fMin + RAY_EPSI), mfMax(fMax - RAY_EPSI) {}
		Ray(const vec3& ptOrig, const vec3& vDir, float fMax = float(Math::EDX_INFINITY), float fMin = 0.0f, int iDepth = 0)
			: mptOrg(ptOrig[0], ptOrig[1], ptOrig[2]), mvDir(vDir[0], vDir[1], vDir[2]), miDepth(iDepth), mfMin(fMin + RAY_EPSI), mfMax(fMax - RAY_EPSI) {}
		// 		Ray(float fOrigX, float fOrigY, float fOrigZ, 
		// 			float fDirX, float fDirY, float fDirZ, 
		// 			float fMin = 0.0f, float fMax = float(Math::EDX_INFINITY))
		// 			: mptOrg(fOrigX, fOrigY, fOrigZ), mvDir(fDirX, fDirY, fDirZ), 
		// 			mfMin(fMin + RAY_EPSI), mfMax(fMax - RAY_EPSI) {}

		~Ray() {}

		inline Vector3 CalcPoint(float fDist) const { return mptOrg + mvDir * fDist; }
	};

	class RayDifferential : public Ray
	{
	public:
		Vector3 mvDxOrg, mvDyOrg;
		Vector3 mvDxDir, mvDyDir;

		bool mbHasDifferential;

	public:
		RayDifferential()
			: Ray(), mbHasDifferential(false) {}
		RayDifferential(const Vector3& ptOrig, const Vector3& vDir, float fMax = float(Math::EDX_INFINITY), float fMin = 0.0f, int iDepth = 0)
			: Ray(ptOrig, vDir, fMax, fMin, iDepth), mbHasDifferential(false) {}
		RayDifferential(const Ray& ray)
			: Ray(ray), mbHasDifferential(false) {}
	};

	namespace Math
	{
		inline float Distance(const Vector3& p1, const Vector3& p2);
	}
	// Axis Aligned Bounding Box definition


	///////////////////////////////////////////////////////////////////////
	//                           BoundingBox                             //
	///////////////////////////////////////////////////////////////////////


	class BoundingBox
	{
	public:
		Vector3 mptMin, mptMax;

	public:
		BoundingBox()
			: mptMin(float(Math::EDX_INFINITY), float(Math::EDX_INFINITY), float(Math::EDX_INFINITY))
			, mptMax(float(Math::EDX_NEG_INFINITY), float(Math::EDX_NEG_INFINITY), float(Math::EDX_NEG_INFINITY))
		{}

		BoundingBox(const Vector3& pt)
			: mptMin(pt)
			, mptMax(pt)
		{}

		BoundingBox(const vec3& pt)
			: mptMin(pt[0], pt[1], pt[2])
			, mptMax(pt[0], pt[1], pt[2])
		{}

		BoundingBox(
			const Vector3& pt1,
			const Vector3& pt2)
		{
			mptMin = Vector3(Math::Min(pt1.x, pt2.x), Math::Min(pt1.y, pt2.y), Math::Min(pt1.z, pt2.z));
			mptMax = Vector3(Math::Max(pt1.x, pt2.x), Math::Max(pt1.y, pt2.y), Math::Max(pt1.z, pt2.z));
		}

		BoundingBox(
			const vec3& pt1,
			const vec3& pt2)
		{
			mptMin = Vector3(Math::Min(pt1[0], pt2[0]), Math::Min(pt1[1], pt2[1]), Math::Min(pt1[2], pt2[2]));
			mptMax = Vector3(Math::Max(pt1[0], pt2[0]), Math::Max(pt1[1], pt2[1]), Math::Max(pt1[2], pt2[2]));
		}

		bool Occluded(const Ray& ray, float* pfHitNear = NULL, float* pfHitFar = NULL) const;

		inline bool Inside(const Vector3& point) const
		{
			return point.x >= mptMin.x && point.x <= mptMax.x &&
				point.y >= mptMin.y && point.y <= mptMax.y &&
				point.z >= mptMin.z && point.z <= mptMax.z;
		}

		inline bool Overlaps(const BoundingBox& bbox) const
		{
			return mptMin.x <= bbox.mptMax.x && mptMax.x >= bbox.mptMin.x &&
				mptMin.y <= bbox.mptMax.y && mptMax.y >= bbox.mptMin.y &&
				mptMin.z <= bbox.mptMax.z && mptMax.z >= bbox.mptMin.x;
		}

		inline int MaximumExtent() const
		{
			Vector3 vDiag = mptMax - mptMin;
			if (vDiag.x > vDiag.y && vDiag.x > vDiag.z)
				return 0;
			else if (vDiag.y > vDiag.z)
				return 1;
			else
				return 2;
		}

		inline Vector3 Centroid() const
		{
			return 0.5f * (mptMin + mptMax);
		}

		inline Vector3 Offset(const Vector3& pt) const
		{
			return Vector3((pt.x - mptMin.x) / (mptMax.x - mptMin.x),
				(pt.y - mptMin.y) / (mptMax.y - mptMin.y),
				(pt.z - mptMin.z) / (mptMax.z - mptMin.z));
		}

		inline float Area() const
		{
			Vector3 vDiag = mptMax - mptMin;
			return 2.0f * (vDiag.x * vDiag.y + vDiag.y * vDiag.z + vDiag.x * vDiag.z);
		}

		inline float Volume() const
		{
			Vector3 vDiag = mptMax - mptMin;
			return vDiag.x * vDiag.y * vDiag.z;
		}

		inline void BoundingSphere(Vector3* pvCenter, float* pfRadius)
		{
			*pvCenter = 0.5f * (mptMin + mptMax);
			*pfRadius = Inside(*pvCenter) ? Math::Distance(*pvCenter, mptMax) : 0.0f;
		}

		inline void Expand(const float fDelta)
		{
			mptMin -= Vector3(fDelta);
			mptMax += Vector3(fDelta);
		}

		inline bool operator == (const BoundingBox& bbox)
		{
			return mptMin == bbox.mptMin && mptMax == bbox.mptMax;
		}

		inline bool operator != (const BoundingBox& bbox)
		{
			return mptMin != bbox.mptMin || mptMax != bbox.mptMax;
		}

		inline const Vector3& operator [] (int i) const
		{
			assert(i == 0 || i == 1);
			return (&mptMin)[i];
		}

		inline Vector3 operator [] (int i)
		{
			assert(i == 0 || i == 1);
			return (&mptMin)[i];
		}

		void RenderInGL() const;
	};

	namespace Math
	{
		template<class T>
		inline T Dot(const Vec<2, T>& vVec1, const Vec<2, T>& vVec2) { return vVec1.x * vVec2.x + vVec1.y * vVec2.y; }
		template<class T>
		inline T Dot(const Vec<3, T>& vVec1, const Vec<3, T>& vVec2) { return vVec1.x * vVec2.x + vVec1.y * vVec2.y + vVec1.z * vVec2.z; }

		template<class T>
		inline T AbsDot(const Vec<2, T>& vVec1, const Vec<2, T>& vVec2) { T ret = Dot(vVec1, vVec2); return ret >= 0 ? ret : -ret; }
		template<class T>
		inline T AbsDot(const Vec<3, T>& vVec1, const Vec<3, T>& vVec2) { T ret = Dot(vVec1, vVec2); return ret >= 0 ? ret : -ret; }

		template<class T>
		inline T Cross(const Vec<2, T>& vVec1, const Vec<2, T>& vVec2)
		{
			return vVec1.x * vVec2.y - vVec1.y * vVec2.x;
		}
		template<class T>
		inline Vec<3, T> Cross(const Vec<3, T>& vVec1, const Vec<3, T>& vVec2)
		{
			return Vec<3, T>(vVec1.y * vVec2.z - vVec1.z * vVec2.y,
				vVec1.z * vVec2.x - vVec1.x * vVec2.z,
				vVec1.x * vVec2.y - vVec1.y * vVec2.x);
		}

		template<class T>
		inline T Curl(const Vec<2, T>& vDvdx, const Vec<2, T>& vDvdy)
		{
			return vDvdx.y - vDvdy.x;
		}
		template<class T>
		inline Vec<3, T> Curl(const Vec<3, T>& vDvdx, const Vec<3, T>& vDvdy, const Vec<3, T>& vDvdz)
		{
			return Vec<3, T>(vDvdy.z - vDvdz.y, vDvdz.x - vDvdx.z, vDvdx.y - vDvdy.x);
		}

		template<class T>
		inline T LengthSquared(const Vec<2, T>& vVec) { return Dot(vVec, vVec); }
		template<class T>
		inline T LengthSquared(const Vec<3, T>& vVec) { return Dot(vVec, vVec); }

		template<uint Dimension>
		inline Vec<Dimension, int> FloorToInt(const Vec<Dimension, float>& vec)
		{
			Vec<Dimension, int> vRet;
			for(auto d = 0; d < Dimension; d++)
				vRet[d] = FloorToInt(vec[d]);
			return vRet;
		}

		template<uint Dimension>
		inline Vec<Dimension, int> RoundToInt(const Vec<Dimension, float>& vec)
		{
			Vec<Dimension, int> vRet;
			for(auto d = 0; d < Dimension; d++)
				vRet[d] = RoundToInt(vec[d]);
			return vRet;
		}

		inline float Length(const Vector2& vVec) { return Math::Sqrt(LengthSquared(vVec)); }
		inline float Length(const Vector3& vVec) { return Math::Sqrt(LengthSquared(vVec)); }

		inline Vector2 Normalize(const Vector2& v) { return v / Length(v); }
		inline Vector3 Normalize(const Vector3& v) { return v / Length(v); }

		inline void CoordinateSystem(const Vector3& vVec1, Vector3* vVec2, Vector3* vVec3)
		{
			if(Math::Abs(vVec1.x) > Math::Abs(vVec1.y))
			{
				float fInvLen = 1.0f / Math::Sqrt(vVec1.x * vVec1.x + vVec1.z * vVec1.z);
				*vVec2 = Vector3(-vVec1.z * fInvLen, 0.0f, vVec1.x * fInvLen);
			}
			else
			{
				float fInvLen = 1.0f / Math::Sqrt(vVec1.y * vVec1.y + vVec1.z * vVec1.z);
				*vVec2 = Vector3(0.0f, vVec1.z * fInvLen, -vVec1.y * fInvLen);
			}
			*vVec3 = Cross(vVec1, *vVec2);
		}


		inline float Distance(const Vector2& p1, const Vector2& p2) { return Length(p1 - p2); }
		inline float Distance(const Vector3& p1, const Vector3& p2) { return Length(p1 - p2); }
		inline float DistanceSquared(const Vector2& p1, const Vector2& p2) { return LengthSquared(p1 - p2); }
		inline float DistanceSquared(const Vector3& p1, const Vector3& p2) { return LengthSquared(p1 - p2); }

		inline Vector3 Reflect(const Vector3& vInci, const Vector3& nNorm) { return vInci + Vector3(2 * Dot(-vInci, nNorm) * nNorm); }
		inline Vector3 Refract(const Vector3& vInci, const Vector3& nNorm, float eta)
		{
			float NDotI = Dot(nNorm, vInci);
			float k = 1.0f - eta * eta * (1.0f - NDotI * NDotI); 
			if(k < 0.0f)
				return Vector3::ZERO;
			else
				return eta * vInci - (eta * NDotI + Math::Sqrt(k)) * Vector3(nNorm);
		}

		inline Vector3 FaceForward(const Vector3& n, const Vector3& v) { return (Dot(n, v) < 0.0f) ? -n : n; }

		inline BoundingBox Union(const BoundingBox& bbox1, const BoundingBox& bbox2)
		{
			BoundingBox ret;

			ret.mptMin.x = Math::Min(bbox1.mptMin.x, bbox2.mptMin.x);
			ret.mptMin.y = Math::Min(bbox1.mptMin.y, bbox2.mptMin.y);
			ret.mptMin.z = Math::Min(bbox1.mptMin.z, bbox2.mptMin.z);
			ret.mptMax.x = Math::Max(bbox1.mptMax.x, bbox2.mptMax.x);
			ret.mptMax.y = Math::Max(bbox1.mptMax.y, bbox2.mptMax.y);
			ret.mptMax.z = Math::Max(bbox1.mptMax.z, bbox2.mptMax.z);

			return ret;
		}

		inline BoundingBox Union(const BoundingBox& bbox, const Vector3& point)
		{
			BoundingBox ret;

			ret.mptMin.x = Math::Min(bbox.mptMin.x, point.x);
			ret.mptMin.y = Math::Min(bbox.mptMin.y, point.y);
			ret.mptMin.z = Math::Min(bbox.mptMin.z, point.z);
			ret.mptMax.x = Math::Max(bbox.mptMax.x, point.x);
			ret.mptMax.y = Math::Max(bbox.mptMax.y, point.y);
			ret.mptMax.z = Math::Max(bbox.mptMax.z, point.z);

			return ret;
		}

		inline BoundingBox Union(const BoundingBox& bbox, const vec3& point)
		{
			BoundingBox ret;

			ret.mptMin.x = Math::Min(bbox.mptMin.x, point[0]);
			ret.mptMin.y = Math::Min(bbox.mptMin.y, point[1]);
			ret.mptMin.z = Math::Min(bbox.mptMin.z, point[2]);
			ret.mptMax.x = Math::Max(bbox.mptMax.x, point[0]);
			ret.mptMax.y = Math::Max(bbox.mptMax.y, point[1]);
			ret.mptMax.z = Math::Max(bbox.mptMax.z, point[2]);

			return ret;
		}
		inline Vector3 SphericalDirection(float fSinTheta, float fCosTheta, float fPhi)
		{
			return Vector3(fSinTheta * Math::Cos(fPhi),
				fCosTheta,
				fSinTheta * Math::Sin(fPhi));
		}


		inline Vector3 SphericalDirection(float fSinTheta, float fCosTheta,
			float fPhi, const Vector3& vX,
			const Vector3& vY, const Vector3& vZ)
		{
			return fSinTheta * Math::Cos(fPhi) * vX +
				fCosTheta * vY + fSinTheta * Math::Sin(fPhi) * vZ;
		}


		inline float SphericalTheta(const Vector3& vVec)
		{
			return Math::Sin(Math::Clamp(vVec.y, -1.0f, 1.0f));
		}


		inline float SphericalPhi(const Vector3& vVec)
		{
			float p = Math::Atan2(vVec.z, vVec.x);
			return (p < 0.0f) ? p + 2.0f * float(float(Math::EDX_PI)) : p;
		}
	};

}



///////////////////////////////////////////////////////////////////////
//                          Intersection                             //
///////////////////////////////////////////////////////////////////////
class Intersection
{
public:
	Vector3 mptPosition;
	Vector3 mvDpdu, mvDpdv;
	Vector3 mnDndu, mnDndv;

	// Differentials
	mutable Vector3 mvDpdx, mvDpdy;
	mutable float mfDudx, mfDudy, mfDvdx, mfDvdy;
	
	mutable float mfDist;
	
public:
	Intersection() 
		: mptPosition(Vector3::ZERO), 
		mfDist(float(Math::EDX_INFINITY)),
		mvDpdu(Vector3::ZERO), mvDpdv(Vector3::ZERO), 
		mnDndu(Vector3::ZERO), mnDndv(Vector3::ZERO),
		mvDpdx(Vector3::ZERO), mvDpdy(Vector3::ZERO),
		mfDudx(0.0f), mfDudy(0.0f), mfDvdx(0.0f), mfDvdy(0.0f)
	{
	}

	Intersection(const Vector3& vPos, float fDist, 
		const Vector3& vDpdu = Vector3::ZERO, const Vector3& vDpdv = Vector3::ZERO,
		const Vector3& nDndu = Vector3::ZERO, const Vector3& nDndv = Vector3::ZERO) 
		: mptPosition(vPos), 
		mfDist(fDist), 
		mvDpdu(vDpdu), mvDpdv(vDpdv), 
		mnDndu(nDndu), mnDndv(nDndv),
		mvDpdx(Vector3::ZERO), mvDpdy(Vector3::ZERO),
		mfDudx(0.0f), mfDudy(0.0f), mfDvdx(0.0f), mfDvdy(0.0f)
	{
	}

	~Intersection() {}
};

///////////////////////////////////////////////////////////////////////
//                            Triangle                               //
///////////////////////////////////////////////////////////////////////
class Triangle
{
private:
	vec3* pt[3];

public:
	Triangle(){}
	void init(vec3& pt1, vec3& pt2, vec3& pt3)
	{
		pt[0] = &pt1;
		pt[1] = &pt2;
		pt[2] = &pt3;
	}
	~Triangle()
	{
	}

	/*
	inline BoundingBox ObjectBound() const
	{
//		const Vector3& pt1 = mpMesh->vertices[mpIndex[0]];
//		const Vector3& pt2 = mpMesh->vertices[mpIndex[1]];
//		const Vector3& pt3 = mpMesh->vertices[mpIndex[2]];
//		return BoundingBox(Math::Union(BoundingBox(pt1, pt2), pt3));
	}
	*/
	inline BoundingBox WorldBound() const
	{
		return BoundingBox(Math::Union(BoundingBox(*pt[0], *pt[1]), *pt[2]));
	}
	

	inline const vec3& GetVertex(uint iIdx) const { return *pt[iIdx]; }
};

class Triangle4
{
public:
	Vec3f_SSE mVertices0;
	Vec3f_SSE mEdges1;
	Vec3f_SSE mEdges2;
	Vec3f_SSE mGeomNormals;

	
public:
	Triangle4()
	{
	}

	Triangle4(const Triangle* pTris, const size_t iCount)
	{
		assert(iCount <= 4);

		for(auto i = 0; i < iCount; i++)
		{
			// Get triangle vertices in pt1, pt2 and pt3
			const vec3& vertex0 = pTris[i].GetVertex(0);
			const vec3& vertex1 = pTris[i].GetVertex(1);
			const vec3& vertex2 = pTris[i].GetVertex(2);

			// Pack 4 triangles into SOA form
			mVertices0.x[i] = vertex0[0];
			mVertices0.y[i] = vertex0[1];
			mVertices0.z[i] = vertex0[2];

			vec3 vEdge1 = vertex0 - vertex1;
			vec3 vEdge2 = vertex2 - vertex0;

			mEdges1.x[i] = vEdge1[0];
			mEdges1.y[i] = vEdge1[1];
			mEdges1.z[i] = vEdge1[2];

			mEdges2.x[i] = vEdge2[0];
			mEdges2.y[i] = vEdge2[1];
			mEdges2.z[i] = vEdge2[2];

			// Pack normals
			vec3 vGeomN = vEdge1.cross(vEdge2);
			mGeomNormals.x[i] = vGeomN[0];
			mGeomNormals.y[i] = vGeomN[1];
			mGeomNormals.z[i] = vGeomN[2];
		}
	}

	__forceinline bool Intersect(const Ray& ray, Intersection* pIsect) const
	{
		// Calculate determinant
		const Vec3f_SSE vOrgs = Vec3f_SSE(ray.mptOrg);
		const Vec3f_SSE vDirs = Vec3f_SSE(ray.mvDir);
		const Vec3f_SSE vC = mVertices0 - vOrgs;
		const Vec3f_SSE vR = Math::Cross(vDirs, vC);
		const FloatSSE det = Math::Dot(mGeomNormals, vDirs);
		const FloatSSE absDet = SSE::Abs(det);
		const FloatSSE signedDet = SSE::SignMask(det);

		// Edge tests
		const FloatSSE U = Math::Dot(vR, mEdges2) ^ signedDet;
		const FloatSSE V = Math::Dot(vR, mEdges1) ^ signedDet;
		BoolSSE valid = (det != FloatSSE(Math::EDX_ZERO)) & (U >= 0.0f) & (V >= 0.0f) & (U + V <= absDet);
		if(SSE::None(valid))
			return false;

		const FloatSSE T = Math::Dot(mGeomNormals, vC) ^ signedDet;
		valid &= (T > absDet * FloatSSE(ray.mfMin)) & (T < absDet * FloatSSE(pIsect->mfDist));
		if(SSE::None(valid))
			return false;

		const FloatSSE invAbsDet = SSE::Rcp(absDet);
		const FloatSSE u = U * invAbsDet;
		const FloatSSE v = V * invAbsDet;
		const FloatSSE t = T * invAbsDet;

		const size_t idx = SSE::SelectMin(valid, t);
		float fT = t[idx];

		Vector3 vPos = ray.CalcPoint(fT);

		Vector3 vDpdu, vDpdv;
		const Vector3 pE1 = Vector3(mEdges1.x[idx], mEdges1.y[idx], mEdges1.z[idx]);
		const Vector3 pE2 = Vector3(mEdges2.x[idx], mEdges2.y[idx], mEdges2.z[idx]);

		ray.mfMax = fT;
		*pIsect = Intersection(vPos, fT, vDpdu, vDpdv);

		return true;
	}

	__forceinline bool Occluded(const Ray& ray) const
	{
		// Calculate determinant
		const Vec3f_SSE vOrgs = Vec3f_SSE(ray.mptOrg);
		const Vec3f_SSE vDirs = Vec3f_SSE(ray.mvDir);
		const Vec3f_SSE vC = mVertices0 - vOrgs;
		const Vec3f_SSE vR = Math::Cross(vDirs, vC);
		const FloatSSE det = Math::Dot(mGeomNormals, vDirs);
		const FloatSSE absDet = SSE::Abs(det);
		const FloatSSE signedDet = SSE::SignMask(det);

		// Edge tests
		Vec3f_SSE me = mEdges2;
		Vec3f_SSE me1 = mEdges1;
		FloatSSE temp = Math::Dot(vR, me);
		const FloatSSE U = temp ^ signedDet;
		const FloatSSE V = Math::Dot(vR, me1) ^ signedDet;
		BoolSSE valid = (det != FloatSSE(Math::EDX_ZERO)) & (U >= 0.0f) & (V >= 0.0f) & (U + V <= absDet);
		if(SSE::None(valid))
			return false;

		const FloatSSE T = Math::Dot(mGeomNormals, vC) ^ signedDet;
		valid &= (T > absDet * FloatSSE(ray.mfMin)) & (T < absDet * FloatSSE(ray.mfMax));
		if(SSE::None(valid))
			return false;

		return true;
	}
};



///////////////////////////////////////////////////////////////////////
//                              BVH                                  //
///////////////////////////////////////////////////////////////////////

struct BVHPrimInfo
{
	int iPrimIdx;
	Vector3 vCentroid;
	BoundingBox bbox;

	BVHPrimInfo(int idx, const BoundingBox& box)
		: iPrimIdx(idx)
		, bbox(box)
	{
		vCentroid = bbox.Centroid();
	}
};

struct BVHBuildNode
{
	BoundingBox boundsLeft;
	BoundingBox boundsRight;
	BVHBuildNode* pChildren[2];
	int iSplitAxis;
	int iFirstPrimOffset;
	int iPrimCount;

	BVHBuildNode()
	{
		pChildren[0] = pChildren[1] = NULL;
	}
	void InitLeaf(int iFirst, int iCount, const BoundingBox& box)
	{
		iFirstPrimOffset = iFirst;
		iPrimCount = iCount;
		boundsLeft = box;
		boundsRight = box;
	}
	void InitInterior(int iAx, BVHBuildNode* pChild0, BVHBuildNode* pChild1)
	{
		pChildren[0] = pChild0;
		pChildren[1] = pChild1;
		boundsLeft = Math::Union(pChild0->boundsLeft, pChild0->boundsRight);
		boundsRight = Math::Union(pChild1->boundsLeft, pChild1->boundsRight);
		iSplitAxis = iAx;
		iPrimCount = 0;
	}
};

struct BVHNode
{
	FloatSSE fMinMaxBoundsX;
	FloatSSE fMinMaxBoundsY;
	FloatSSE fMinMaxBoundsZ;
	union
	{
		uint iPrimOffset;
		uint iSecondChildOffset;
	};
	uint8 iPrimCount;
	uint8 iAxis;
	uint16 pad;
};

class BVH
{
private:
	friend class BVHVisualizer;

	// Node tree top ptr
	BVHNode* mpNodes;

	// Primitives storage
	vector<Triangle> mvPrimitives;
	vector<Triangle4> mvPackedTriangles;

	// Bounding box
	BoundingBox mWorldBound;

	// Total primitives count
	uint miPrimCount;
	uint miMaxNodePrim;

	uint miTotalNodesCount;
public:
	BVH(TriMesh* pmesh, int iMaxNodePrim = 16);
	~BVH()
	{
		FreeAligned(mpNodes);
		if (vertexOccludedInside)
			delete[] vertexOccludedInside;
		if (vertexOccludedOutside)
			delete[] vertexOccludedOutside;
	}

public:
	// Accelerated ray-primitive intersection routine
	bool Intersect(const Ray& ray, Intersection* pIsect) const;
	bool Occluded(const Ray& ray) const;
	bool IsVertVertOutsideOccluded(int startVertexIndex, int endtVertexIndex);
	bool IsVertVertInsideOccluded(int startVertexIndex, int endtVertexIndex);

	float ShortesDistanceInside(int startVertexIndex, int endtVertexIndex);

	BVHNode* GetRootPtr() const { return mpNodes; }

	// For debugging the intersection routine, returns the index of bounding boxes traversed
	vector<uint> IntersectDebug(const Ray& ray, Intersection* pIsect) const;

	BoundingBox WorldBound() const { return mWorldBound; }

private:
	BVHBuildNode* Construction(vector<BVHPrimInfo>& vBuildData,
		uint& iTotalNodes,
		uint iStart, uint iEnd,
		MemoryArena& memory);
	TriMesh *mesh;
	uint FlattenTree(BVHBuildNode* pNode, uint* piOffset);
	bool *vertexOccludedInside;
	bool *vertexOccludedOutside;
	float *vertexDistanceInside;
};

class BVHVisualizer
{
private:
	const BVH* mpBVH;
	vector<BoundingBox> mvBBoxList;
	vector<uint> mvRayTravIndices;

public:
	BVHVisualizer()
	{
	}
	~BVHVisualizer()
	{
	}

	void Initialize(const BVH* pBVH);
	void GenRayTraverseIndices(const Ray& ray);
	void RenderIntersection(const Ray& ray, int iIndex) const;
};
