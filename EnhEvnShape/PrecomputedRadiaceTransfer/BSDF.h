#include "stdafx.h"
#include <ctime>
#include "assert.h"
/*



inline void CoordinateSystem(const vec3& vVec1, vec3* vVec2, vec3* vVec3)
{
	if(abs(vVec1[0]) > abs(vVec1[1]))
	{
		float fInvLen = 1.0f / sqrt(vVec1[0] * vVec1[0] + vVec1[2] * vVec1[2]);
		*vVec2 = vec3(-vVec1[2] * fInvLen, 0.0f, vVec1[0] * fInvLen);
	}
	else
	{
		float fInvLen = 1.0f / sqrt(vVec1[1] * vVec1[1] + vVec1[2] * vVec1[2]);
		*vVec2 = vec3(0.0f, vVec1[2] * fInvLen, -vVec1[1] * fInvLen);
	}
	*vVec3 = vVec1.cross(vVec2);
}

class RandomGen
{
private:
	__m128i* mpCurSeed;

public:
	RandomGen(unsigned int seed = (unsigned int)time(NULL))
	{
		mpCurSeed = (__m128i*)_aligned_malloc(sizeof(__m128i), 16);
		*mpCurSeed =_mm_set_epi32(seed, seed+1, seed, seed+1);
	}

	~RandomGen()
	{
		_aligned_free(mpCurSeed);
	}

	inline unsigned int UnsignedInt() 
	{
		__declspec(align(16)) unsigned int result[4];
		__declspec(align(16)) __m128i cur_seed_split;
		__declspec(align(16)) __m128i multiplier;
		__declspec(align(16)) __m128i adder;
		__declspec(align(16)) __m128i mod_mask;
		__declspec(align(16)) __m128i sra_mask;
		__declspec(align(16)) static const unsigned int mult[4] =
		{ 214013, 17405, 214013, 69069 };
		__declspec(align(16)) static const unsigned int gadd[4] =
		{ 2531011, 10395331, 13737667, 1 };
		__declspec(align(16)) static const unsigned int mask[4] =
		{ 0xFFFFFFFF, 0, 0xFFFFFFFF, 0 };
		__declspec(align(16)) static const unsigned int masklo[4] =
		{ 0x00007FFF, 0x00007FFF, 0x00007FFF, 0x00007FFF };

		adder =_mm_load_si128((__m128i*) gadd);
		multiplier =_mm_load_si128((__m128i*) mult);
		mod_mask =_mm_load_si128((__m128i*) mask);
		sra_mask =_mm_load_si128((__m128i*) masklo);
		cur_seed_split =_mm_shuffle_epi32(*mpCurSeed, _MM_SHUFFLE(2, 3, 0, 1));

		*mpCurSeed =_mm_mul_epu32(*mpCurSeed, multiplier);
		multiplier =_mm_shuffle_epi32(multiplier, _MM_SHUFFLE(2, 3, 0, 1));
		cur_seed_split =_mm_mul_epu32(cur_seed_split, multiplier);
		*mpCurSeed =_mm_and_si128(*mpCurSeed, mod_mask);
		cur_seed_split =_mm_and_si128(cur_seed_split, mod_mask);
		cur_seed_split =_mm_shuffle_epi32(cur_seed_split, _MM_SHUFFLE(2, 3, 0, 1));
		*mpCurSeed =_mm_or_si128(*mpCurSeed, cur_seed_split);
		*mpCurSeed =_mm_add_epi32(*mpCurSeed, adder);

#ifdef COMPATABILITY 
		__declspec(align(16)) __m128i sseresult;
		// Add the lines below if you wish to reduce your results to 16-bit vals... 
		sseresult =_mm_srai_epi32(*mpCurSeed, 16);
		sseresult =_mm_and_si128(sseresult, sra_mask);
		_mm_storeu_si128((__m128i*) result, sseresult);
		return;
#endif 
		_mm_storeu_si128((__m128i*) result, *mpCurSeed);
		return result[0];
	}

	inline float Float()
	{
		__declspec(align(16)) unsigned int temp;
		temp = UnsignedInt();
#ifdef COMPATABILITY
		return temp / (float)(0xFFFF);
#else
		return temp / (float)(0xFFFFFFFF);
#endif
	}

	inline float GaussFloat()
	{
		float u1 = Float();
		float u2 = Float();
		if (u1 < 1e-6f)
			u1 = 1e-6f;
		return sqrt(-2.0f * logf(u1)) * cos(2 * M_PI * u2);
	}
};


struct Sample3D
{
	float fDir1, fDir2;
	float fComponent;
	Sample3D()
		: fDir1(0.0f)
		, fDir2(0.0f)
		, fComponent(0.0f)
	{
	}
	Sample3D(float f1, float f2, float f3)
		: fDir1(f1)
		, fDir2(f2)
		, fComponent(f3)
	{
	}
	Sample3D(RandomGen& random) 
	{
		fDir1 = random.Float();
		fDir2 = random.Float();
		fComponent = random.Float();
	}
};





class Frame
{
public:
	vec3 mvX, mvY, mvZ;

public:
	Frame()
	{
		mvX = vec3(1.f, 0.f, 0.f);
		mvY = vec3(0.f, 1.f, 0.f);
		mvZ = vec3(0.f, 0.f, 1.f);
	};

	Frame(const vec3& x,
		const vec3& y,
		const vec3& z
		)
		: mvX(x)
		, mvY(y)
		, mvZ(z)
	{}

	Frame(const vec3& vNormal)
	{
		mvZ = vNormal;
		normalize(mvZ);

		CoordinateSystem(mvZ, &mvX, &mvY);
	}

	vec3 LocalToWorld(const vec3& vVec) const
	{
		return mvX * vVec[0] + mvY * vVec[1] + mvZ * vVec[2];
	}

	vec3 WorldToLocal(const vec3& vVec) const
	{
		return vec3(vVec.dot(mvX), vVec.dot(mvY), vVec.dot(mvZ));
	}

	const vec3& Binormal() const { return mvX; }
	const vec3& Tangent () const { return mvY; }
	const vec3& Normal  () const { return mvZ; }

};

class Intersection
{
public:
	vec3 mptPosition;
	vec3 mnNormal;
	vec3 mnGeomNormal;
	vec2 mvTexCoord;
	vec3 mvDpdu, mvDpdv;
	vec3 mnDndu, mnDndv;

	// Differentials
	mutable vec3 mvDpdx, mvDpdy;
	mutable float mfDudx, mfDudy, mfDvdx, mfDvdy;

	Frame mShadingFrame;
	Frame mGeomFrame;

	mutable float mfDist;

	const BSDF* mpBSDF;

public:
	Intersection() 
		: mptPosition(vec3(0.f)), mnNormal(vec3(0.f, 1.f, 0.f)), mnGeomNormal(vec3(0.f, 1.f, 0.f)), 
		mfDist(FLT_MAX), mpBSDF(NULL), 
		mvDpdu(vec3(0.f)), mvDpdv(vec3(0.f)), 
		mnDndu(vec3(0.f)), mnDndv(vec3(0.f)),
		mvDpdx(vec3(0.f)), mvDpdy(vec3(0.f)),
		mfDudx(0.0f), mfDudy(0.0f), mfDvdx(0.0f), mfDvdy(0.0f),
		mShadingFrame(mnNormal), mGeomFrame(mnGeomNormal)
	{
	}

	Intersection(const vec3& vPos, const vec3& nNorm, const vec3& nGeomNorm, const vec2& vTex, 
		float fDist, const BSDF* pBSDF,
		const vec3& vDpdu = vec3(0.f), const vec3& vDpdv = vec3(0.f),
		const vec3& nDndu = vec3(0.f), const vec3& nDndv = vec3(0.f)) 
		: mptPosition(vPos), mnNormal(nNorm), mnGeomNormal(nGeomNorm), mvTexCoord(vTex),
		mfDist(fDist), mpBSDF(pBSDF), 
		mvDpdu(vDpdu), mvDpdv(vDpdv), 
		mnDndu(nDndu), mnDndv(nDndv),
		mvDpdx(vec3(0.f)), mvDpdy(vec3(0.f)),
		mfDudx(0.0f), mfDudy(0.0f), mfDvdx(0.0f), mfDvdy(0.0f),
		mShadingFrame(mnNormal), mGeomFrame(mnGeomNormal)
	{
	}

	~Intersection() {}

	inline const unsigned int GetPrimitiveID() const { return 0; }
	inline const BSDF* GetBSDF() const { return mpBSDF; }
	Color Emit(const vec3& vOut) const;

	inline vec3 WorldToLocal(const vec3& vVec) const 
	{
		return mShadingFrame.WorldToLocal(vVec);
	}
	inline vec3 LocalToWorld(const vec3& vVec) const 
	{
		return mShadingFrame.LocalToWorld(vVec);
	}

	void ComputeDifferentials(const RayDifferential& ray) const;
};

template<class T>
class Texture
{
public:
	virtual ~Texture() {}
	virtual T Sample(const Intersection& isect) const = 0;
};

template<class T>
class ConstantTexture : public Texture<T>
{
private:
	T mVal;

public:
	ConstantTexture(const T& val)
		: mVal(val) {}

	T Sample(const Intersection& isect) const
	{
		return mVal;
	}
};

template<class T>
class ImageTexture : public Texture<T>
{
private:
	int miTexWidth;
	int miTexHeight;

	T* mpTexels;

public:
	ImageTexture(const char* strFile);
	~ImageTexture()
	{
		SafeDeleteArray(mpTexels);
	}

	T Sample(const Intersection& isect) const;

	static T GammaCorrect(T tIn, float fGamma = 2.2f)
	{
		return Math::Pow(tIn, fGamma);
	}
	static T ConvertIn(byte _In)
	{
		return GammaCorrect(_In * 0.00390625f);
	}
};


enum ScatterType
{
	BSDF_REFLECTION   = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE      = 1 << 2,
	BSDF_GLOSSY       = 1 << 3,
	BSDF_SPECULAR     = 1 << 4,
	BSDF_ALL_TYPES			= BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
	BSDF_ALL_REFLECTION		= BSDF_REFLECTION | BSDF_ALL_TYPES,
	BSDF_ALL_TRANSMISSION	= BSDF_TRANSMISSION | BSDF_ALL_TYPES,
	BSDF_ALL				= BSDF_ALL_REFLECTION | BSDF_ALL_TRANSMISSION
};

enum class BSDFType
{
	Diffuse, Mirror, Glass
};

class BSDF
{
public:
	const ScatterType mType;
	const BSDFType mBSDFType;
	Texture<Color>* mpColorTex;

public:
	BSDF(ScatterType t, BSDFType t2, const Color& cColor);
	BSDF(ScatterType t, BSDFType t2, const char* strTexPath);
	virtual ~BSDF();

	bool MatchesTypes(ScatterType flags) const { return (mType & flags) == mType; }
	bool IsSpecular() const { return (ScatterType(BSDF_SPECULAR | BSDF_DIFFUSE | BSDF_GLOSSY) & mType) == ScatterType(BSDF_SPECULAR); }
	virtual Color Eval(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types = BSDF_ALL) const;
	virtual float PDF(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types = BSDF_ALL) const;
	virtual Color SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF,
		ScatterType types = BSDF_ALL, ScatterType* pSampledTypes = NULL) const = 0;
	Color GetColor(const Intersection& isect) const;
	BSDFType GetBSDFType() const { return mBSDFType; }

	static BSDF* CreateBSDF(const BSDFType type, const Color color);
	static BSDF* CreateBSDF(const BSDFType type, const char* strTexPath);

private:
	virtual Color Eval(const vec3& vOut, const vec3& vIn, ScatterType types = BSDF_ALL) const = 0;
	virtual float PDF(const vec3& vOut, const vec3& vIn, ScatterType types = BSDF_ALL) const = 0;
};

class LambertianDiffuse : public BSDF
{
public:
	LambertianDiffuse(const Color cColor = Color(1.0))
		: BSDF(ScatterType(BSDF_REFLECTION | BSDF_DIFFUSE), BSDFType::Diffuse, cColor)
	{
	}
	LambertianDiffuse(const char* strTexPath)
		: BSDF(ScatterType(BSDF_REFLECTION | BSDF_DIFFUSE), BSDFType::Diffuse, strTexPath)
	{
	}

	Color SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF,
		ScatterType types = BSDF_ALL, ScatterType* pSampledTypes = NULL) const;

private:
	float PDF(const vec3& vIn, const vec3& vOut, ScatterType types = BSDF_ALL) const;
	Color Eval(const vec3& vOut, const vec3& vIn, ScatterType types = BSDF_ALL) const;
};


class Mirror : public BSDF
{
public:
	Mirror(const Color cColor = Color(1.0))
		: BSDF(ScatterType(BSDF_REFLECTION | BSDF_SPECULAR), BSDFType::Mirror, cColor)
	{
	}
	Mirror(const char* strTexPath)
		: BSDF(ScatterType(BSDF_REFLECTION | BSDF_SPECULAR), BSDFType::Mirror, strTexPath)
	{
	}

	Color Eval(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types = BSDF_ALL) const;
	float PDF(const vec3& vIn, const vec3& vOut, const Intersection& isect, ScatterType types = BSDF_ALL) const;
	Color SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF,
		ScatterType types = BSDF_ALL, ScatterType* pSampledTypes = NULL) const;

private:
	Color Eval(const vec3& vOut, const vec3& vIn, ScatterType types = BSDF_ALL) const;
	float PDF(const vec3& vIn, const vec3& vOut, ScatterType types = BSDF_ALL) const;
};

class Glass : public BSDF
{
private:
	float mfEtai, mfEtat;
	ScatterType mTypeRef, mTypeRfr;

public:
	Glass(const Color cColor =  Color(1.0), float fetai = 1.0f, float fetat = 1.5f) 
		: BSDF(ScatterType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR), BSDFType::Glass, cColor)
		, mTypeRef(ScatterType(BSDF_REFLECTION | BSDF_SPECULAR))
		, mTypeRfr(ScatterType(BSDF_TRANSMISSION | BSDF_SPECULAR))
		, mfEtai(fetai)
		, mfEtat(fetat)
	{
	}
	Glass(const char* strTexPath, float fetai = 1.0f, float fetat = 1.5f) 
		: BSDF(ScatterType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR), BSDFType::Glass, strTexPath)
		, mTypeRef(ScatterType(BSDF_REFLECTION | BSDF_SPECULAR))
		, mTypeRfr(ScatterType(BSDF_TRANSMISSION | BSDF_SPECULAR))
		, mfEtai(fetai)
		, mfEtat(fetat)
	{
	}

	Color Eval(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types = BSDF_ALL) const;
	float PDF(const vec3& vIn, const vec3& vOut, const Intersection& isect, ScatterType types = BSDF_ALL) const;
	Color SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF,
		ScatterType types = BSDF_ALL, ScatterType* pSampledTypes = NULL) const;

private:
	Color Eval(const vec3& vOut, const vec3& vIn, ScatterType types = BSDF_ALL) const;
	float PDF(const vec3& vIn, const vec3& vOut, ScatterType types = BSDF_ALL) const;
	float Fresnel(float fCosi) const;
};


namespace BSDFCoordinate
{
	inline float CosTheta(const vec3& vVec) { return vVec[2]; }
	inline float AbsCosTheta(const vec3& vVec) { return abs(vVec[2]); }
	inline float SinTheta2(const vec3& vVec) { return max(0.0f, 1.0f - CosTheta(vVec) * CosTheta(vVec)); }
	inline float SinTheta(const vec3& vVec) { return sqrt(SinTheta2(vVec)); }
	inline float CosPhi(const vec3& vVec)
	{
		float sintheta = SinTheta(vVec);
		if (sintheta == 0.0f) return 1.0f;
		return clamp(vVec[0] / sintheta, -1.0f, 1.0f);
	}

	inline float SinPhi(const vec3& vVec)
	{
		float sintheta = SinTheta(vVec);
		if(sintheta == 0.0f)
			return 0.0f;
		return clamp(vVec[1] / sintheta, -1.0f, 1.0f);
	}

	inline bool SameHemisphere(const vec3& vVec1, const vec3& vVec2) { return vVec1[2] * vVec2[2] > 0.0f; }
};