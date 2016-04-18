#include "stdafx.h"
#include "BSDF.h"
/*
// -----------------------------------------------------------------------------------------------------------------------
// BSDF interface implementation
// -----------------------------------------------------------------------------------------------------------------------
BSDF::BSDF(ScatterType t, BSDFType t2, const Color& cColor)
	: mType(t), mBSDFType(t2)
{
	mpColorTex = new ConstantTexture<Color>(cColor);
}
BSDF::BSDF(ScatterType t, BSDFType t2, const char* strTexPath)
	: mType(t), mBSDFType(t2)
{
	mpColorTex = new ImageTexture<Color>(strTexPath);
}

BSDF* BSDF::CreateBSDF(const BSDFType type, const Color color)
{
	switch(type)
	{
	case BSDFType::Diffuse:
		return new LambertianDiffuse(color);
	case BSDFType::Mirror:
		return new Mirror(color);
	case BSDFType::Glass:
		return new Glass(color);
	}

	assert(0);
	return NULL;
}
BSDF* BSDF::CreateBSDF(const BSDFType type, const char* strTexPath)
{
	switch(type)
	{
	case BSDFType::Diffuse:
		return new LambertianDiffuse(strTexPath);
	case BSDFType::Mirror:
		return new Mirror(strTexPath);
	case BSDFType::Glass:
		return new Glass(strTexPath);
	}

	assert(0);
	return NULL;
}

BSDF::~BSDF()
{
	delete[] mpColorTex;
}

Color BSDF::Eval(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types) const
{
	if(vOut.dot(isect.mnGeomNormal) * vIn.dot(isect.mnGeomNormal) > 0.0f)
		types = ScatterType(types & ~BSDF_TRANSMISSION);
	else
		types = ScatterType(types & ~BSDF_REFLECTION);

	if(!MatchesTypes(types))
	{
		return Color(0.0);
	}

	vec3 vWo = isect.WorldToLocal(vOut);
	vec3 vWi = isect.WorldToLocal(vIn);

	return GetColor(isect) * Eval(vWo, vWi, types);
}

float BSDF::PDF(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types ) const
{
	if(!MatchesTypes(types))
	{
		return 0.0f;
	}

	vec3 vWo = isect.WorldToLocal(vOut);
	vec3 vWi = isect.WorldToLocal(vIn);

	return PDF(vWo, vWi, types);
}

Color BSDF::GetColor(const Intersection& isect) const
{
	return mpColorTex->Sample(isect);
}

// -----------------------------------------------------------------------------------------------------------------------
// Lambertian BRDF Implementation
// -----------------------------------------------------------------------------------------------------------------------
Color LambertianDiffuse::Eval(const vec3& vOut, const vec3& vIn, ScatterType types) const
{
	if(!BSDFCoordinate::SameHemisphere(vOut, vIn))
	{
		return Color(0.0);
	}

	return Color(FLT_MAX);
}

float LambertianDiffuse::PDF(const vec3& vOut, const vec3& vIn, ScatterType types ) const
{
	if(!BSDFCoordinate::SameHemisphere(vOut, vIn))
	{
		return 0.0f;
	}

	return BSDFCoordinate::AbsCosTheta(vIn) * FLT_MAX;
}

Color LambertianDiffuse::SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF,
										 ScatterType types, ScatterType* pSampledTypes) const
{
	if(!MatchesTypes(types))
	{
		*pfPDF = 0.0f;
		return Color(0.0);
	}

	vec3 vWo = isect.WorldToLocal(vOut), vWi;
	vWi = Sampling::CosineSampleHemisphere(sample.fDir1, sample.fDir2);

	if(vWo[2] < 0.0f)
		vWi[2] *= -1.0f;

	*pvIn = isect.LocalToWorld(vWi);
	*pfPDF = PDF(vWo, vWi, types);
	if(pSampledTypes != NULL)
	{
		*pSampledTypes = mType;
	}

	return GetColor(isect) * Eval(vWo, vWi, types);
}

// -----------------------------------------------------------------------------------------------------------------------
// Mirror Implementation
// -----------------------------------------------------------------------------------------------------------------------
Color Mirror::Eval(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types) const
{
	return Color(0.0);
}

Color Mirror::Eval(const vec3& vOut, const vec3& vIn, ScatterType types) const
{
	return Color(0.0);
}

float Mirror::PDF(const vec3& vIn, const vec3& vOut, const Intersection& isect, ScatterType types ) const
{
	return 0.0f;
}

float Mirror::PDF(const vec3& vIn, const vec3& vOut, ScatterType types ) const
{
	return 0.0f;
}

Color Mirror::SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF, ScatterType types , ScatterType* pSampledTypes ) const
{
	if(!MatchesTypes(types))
	{
		*pfPDF = 0.0f;
		return Color(0.0);
	}

	vec3 vWo = isect.WorldToLocal(vOut), vWi;
	vWi = vec3(-vWo[0], -vWo[1], vWo[2]);

	*pvIn = isect.LocalToWorld(vWi);
	*pfPDF = 1.0f;
	if(pSampledTypes != NULL)
	{
		*pSampledTypes = mType;
	}

	return GetColor(isect) / BSDFCoordinate::AbsCosTheta(vWi);
}

// -----------------------------------------------------------------------------------------------------------------------
// Glass Implementation
// -----------------------------------------------------------------------------------------------------------------------
Color Glass::Eval(const vec3& vOut, const vec3& vIn, const Intersection& isect, ScatterType types) const
{
	return Color(0.0);
}

Color Glass::Eval(const vec3& vOut, const vec3& vIn, ScatterType types) const
{
	return Color(0.0);
}

float Glass::PDF(const vec3& vIn, const vec3& vOut, const Intersection& isect, ScatterType types ) const
{
	return 0.0f;
}

float Glass::PDF(const vec3& vIn, const vec3& vOut, ScatterType types ) const
{
	return 0.0f;
}

Color Glass::SampleScattered(const vec3& vOut, const Sample3D& sample, const Intersection& isect, vec3* pvIn, float* pfPDF, ScatterType types , ScatterType* pSampledTypes ) const
{
	bool bSampleReflect = (types & mTypeRef) == mTypeRef;
	bool bSampleRefract = (types & mTypeRfr) == mTypeRfr;

	if(!bSampleReflect && !bSampleRefract)
	{
		*pfPDF = 0.0f;
		return Color(0.0);
	}

	bool bSampleBoth = bSampleReflect == bSampleRefract;

	vec3 vWo = isect.WorldToLocal(vOut), vWi;

	float fFresnel = Fresnel(BSDFCoordinate::CosTheta(vWo));
	float fProb = fFresnel;//0.5f * fFresnel + 0.25f;

	if(sample.fComponent <= fProb && bSampleBoth || (bSampleReflect && !bSampleBoth)) // Sample reflection
	{
		vWi = vec3(-vWo[0], -vWo[1], vWo[2]);

		*pvIn = isect.LocalToWorld(vWi);
		*pfPDF = !bSampleBoth ? 1.0f : fProb;
		if(pSampledTypes != NULL)
		{
			*pSampledTypes = mTypeRef;
		}

		return fFresnel * GetColor(isect) / BSDFCoordinate::AbsCosTheta(vWi);
	}
	else if(sample.fComponent > fProb && bSampleBoth || (bSampleRefract && !bSampleBoth)) // Sample refraction
	{
		bool bEntering = BSDFCoordinate::CosTheta(vWo) > 0.0f;
		float etai = mfEtai, etat = mfEtat;
		if(!bEntering)
			swap(etai, etat);

		float fSini2 = BSDFCoordinate::SinTheta2(vWo);
		float fEta = etai / etat;
		float fSint2 = fEta * fEta * fSini2;

		if(fSint2 > 1.0f)
			return Color(0.0);

		float fCost = sqrt(max(0.0f, 1.0f - fSint2));
		if(bEntering)
			fCost = -fCost;
		float fSintOverSini = fEta;

		vWi = vec3(fSintOverSini * -vWo[0], fSintOverSini * -vWo[1], fCost);

		*pvIn = isect.LocalToWorld(vWi);
		*pfPDF = !bSampleBoth ? 1.0f : 1.0f - fProb;
		if(pSampledTypes != NULL)
		{
			*pSampledTypes = mTypeRfr;
		}

		return (1.0f - fFresnel) * GetColor(isect) / BSDFCoordinate::AbsCosTheta(vWi);
	}

	return Color(0.0);
}

float Glass::Fresnel(float fCosi) const
{
	fCosi = clamp(fCosi, -1.0f, 1.0f);

	bool bEntering = fCosi > 0.0f;
	float ei = mfEtai, et = mfEtat;
	if (!bEntering)
		swap(ei, et);

	float fSint = ei / et * sqrt(max(0.0f, 1.0f - fCosi * fCosi));

	if(fSint >= 1.0f)
		return 1.0f;
	else
	{
		float fCost = sqrt(max(0.0f, 1.0f - fSint * fSint));
		fCosi = abs(fCosi);

		float fPara = ((mfEtat * fCosi) - (mfEtai * fCost)) /
			((mfEtat * fCosi) + (mfEtai * fCost));
		float fPerp = ((mfEtai * fCosi) - (mfEtat * fCost)) /
			((mfEtai * fCosi) + (mfEtat * fCost));

		return (fPara * fPara + fPerp * fPerp) / 2.0f;
	}
}
*/