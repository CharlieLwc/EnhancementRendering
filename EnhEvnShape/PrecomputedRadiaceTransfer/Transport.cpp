
#include "stdafx.h"

#include "Transport.h"
#include <iomanip>

#include"../GL/glew.h"
#include"../GL/glut.h"

// for demonstration
//#define usingOcluded false



#define PI 3.141592654f

float shiness = 5.0f;
float kd = 0.3f;
float ks = 0.5f;

float eta = 1.3;
vec3 sigma_a(0.0011,0.0024, 0.014);
vec3 sigma_s(2.55, 3.12, 3.77);

float alph = 0.5;


vec3 iosDir = vec3(0.0f, 1.0f, 0.f);
float alphaX = 0.5f;
float alphaY = 0.13f;
float diffuseValue = 0.5f;
float specularValue = 0.5f;




// Clamp a floating point number to [0,1]
inline float Clamp(float num)
{
	return num < 0 ? 0 : num > 1 ? 1 : num;
}

// Quantize and gamma correct a floating point
// color component
inline unsigned char QuantizeAndGC(float comp)
{
	return (unsigned char)(pow(Clamp(comp), 1.0f/2.2f)*255+.5);
}



void calc_d(float *sumDiff, int lightNum, VectorXf& data)
{



	vec3 LABmin(0.f, -85.f, -106.f);
	vec3 LABmax(98.f, 96.6f, 92.9f);
	data[0] = 0.f;

	for (int i = 1; i<lightNum; i++) {
		data[i] = sumDiff[i] - sumDiff[i - 1] + lightNum*data[i - 1];
		data[i] /= (float)lightNum;
	}


	float max = FLT_MIN;
	float min = FLT_MAX;

	for (int i = 0; i < lightNum; i++)
	{
		if (data[i] > max)
			max = data[i];
		if (data[i] < min)
			min = data[i];

	}

	float scale = 180.f / (max - min);
	float trans = -min*scale - 90.f;

	for (int i = 0; i < lightNum; i++)
	{
		data[i] = data[i] * scale + trans;
	}



}

float diffuseShader(vec3 normal, vec3 lightDir, vec3 viewDir)
{
	//Calculate cosine term for this sample
	float LdotN = max(lightDir.dot(normal), 0.f);
	return LdotN;
}

float wardShader(vec3 normal, vec3 lightDir, vec3 viewDir)
{
	// renormalize

	
	vec3 tangent = normalize(iosDir.cross(normal));
	vec3 bitangent = normalize(tangent.cross(normal));


	float NdotL = max(normal.dot(lightDir), 0.f);

	float VdotN = viewDir.dot(normal);

	

	// calculate intermediary values
	vec3 halfwayVector = lightDir + viewDir;
	halfwayVector = normalize(halfwayVector);

	float HdotN = halfwayVector.dot(normal);
	float HdotTAlphaX = halfwayVector.dot(tangent) / alphaX;
	float HdotBAlphaY = halfwayVector.dot(bitangent) / alphaY;


	// calculate N dot L
	
	float diffuse = diffuseValue * NdotL;


	// calculate the specular value
	float exponent = -2.0 * (HdotTAlphaX * HdotTAlphaX + HdotBAlphaY * HdotBAlphaY) / (1.0 + HdotN);
	// Evaluate the specular denominator
	float s_den = 4.0f * 3.14159f;
	s_den *= alphaX;
	s_den *= alphaY;
	s_den *= sqrt(NdotL * VdotN)+0.0001;

	float specular = exp(exponent) / s_den  * specularValue;

	specular = max(specular, 0.f);
	float finalValue = (diffuseValue + specular) * NdotL;
	
	if (finalValue < 0.f)
		cout << "negtive " << endl;

	return finalValue;


	
}

float translucent(point &startPoint, point &endPoint, vec &startNormal, vec &endNormal, vec &lightDirection, vec &eyeDirection, bool isGlossy, RefractionParemeter &rp, float lengthUnitization)
{
	if (isGlossy)
	{


		return 0.f;



	}
	else
	{
		vec3 InsideLightDir = startPoint - endPoint;


		float dpos, dneg, Rd1r, Rd2r, Rd1v, Rd2v, Rd;

		float d2 = len(InsideLightDir) * lengthUnitization * 100;

		dpos = sqrt(d2 + rp.zpos2[2]);
		dneg = sqrt(d2 + rp.zneg2[2]);

		float trDpos = rp.sigma_tr[2] * dpos;
		float trDneg = rp.sigma_tr[2] * dneg;

		Rd1r = rp.zpos[2] * (trDpos + 1.0);
		Rd2r = exp(-trDpos) / (dpos * dpos * dpos);
		Rd1v = rp.zneg[2] * (trDneg + 1.0);
		Rd2v = exp(-trDneg) / (dneg * dneg * dneg);

		Rd = rp.Rd0[2] * (Rd1r * Rd2r - Rd1v * Rd2v);

		return Clamp(Rd);



	}
}

float radianceDistLightChange(vec3 normal, vec3 lightDir, vec3 viewDir, vec3 crossResult, vec3 curvDir, float sharpness, bool isGlossy)
{
	if (isGlossy)
	{
		vec proLightDir = lightDir - lightDir.dot(crossResult)*crossResult;
		float cos = proLightDir.dot(normal);

		cos *= cos;
		cos *= cos;
		cos *= proLightDir.dot(curvDir);

		if (cos > 0.00001)
			cos = pow(cos, sharpness);
		else if (cos < 0.00001)
			cos = -pow(abs(cos), sharpness);
		return cos;
	}
	else
	{
		vec proLightDir = lightDir - lightDir.dot(crossResult)*crossResult;
		float sin = proLightDir.dot(normal);


		if (sin > 0.00001)
			sin = pow(sin, sharpness);
		else if (sin < 0.00001)
			sin = -pow(abs(sin), sharpness);
		return sin;
	}

}

///////////////////////////////////////////////////////////////////////
//                          Directions                               //
///////////////////////////////////////////////////////////////////////

CubeDirections::~CubeDirections()
{
	if (isInit)
	{
		for (int i=0; i<cubeSideLength; i++)
		{
			delete[] a[i];
			delete[] b[i];
			delete[] c[i];
		}
		delete[] a;
		delete[] b;
		delete[] c;
	}
}

vec3 CubeDirections::getDirection(int f, int x, int y)
{
	if (!isInit)
	{
		isInit = true;
		a = new float*[cubeSideLength];
		b = new float*[cubeSideLength];
		c = new float*[cubeSideLength];
		float *uv = new float[cubeSideLength];
		float *uuvv = new float[cubeSideLength];

		for (int x = 0; x < cubeSideLength; x++)
		{
			a[x] = new float[cubeSideLength];
			b[x] = new float[cubeSideLength];
			c[x] = new float[cubeSideLength];
			uv[x] = 2.f*((float)x+0.5f)/float(cubeSideLength)-1.f;
			uuvv[x] = uv[x]*uv[x];
		}

		for (int x = 0; x < cubeSideLength; x++)
		{
			for (int y = 0; y < cubeSideLength; y++)
			{
				float u = ((float)x+0.5)/float(cubeSideLength);
				float v = ((float)y+0.5)/float(cubeSideLength);
				a[x][y] = 1.f / sqrt(1.f + uuvv[x] + uuvv[y]);
				b[x][y] = a[x][y]*uv[x];
				c[x][y] = a[x][y]*uv[y];
			}
		}
		delete[] uuvv;
		delete[] uv;
	}
	
	switch (f)
	{
	case 0:	return vec3(b[x][y],a[x][y],c[x][y]);
	case 1:	return -vec3(a[x][y],c[x][y],-b[x][y]);
	case 2:	return vec3(b[x][y],-c[x][y],a[x][y]);
	case 3:	return vec3(a[x][y],-c[x][y],-b[x][y]);
	case 4:	return -vec3(-b[x][y],a[x][y],c[x][y]);
	case 5:	return -vec3(-b[x][y],-c[x][y],a[x][y]);
	default: return vec3(0.f, 0.f, 0.f);
	}
	

}



///////////////////////////////////////////////////////////////////////
//                       WAVELETCOMPRESSOR                           //
///////////////////////////////////////////////////////////////////////

WaveletCompressor::WaveletCompressor(int dim)
	: dimension(dim)
{
	xform.resize(dim, dim);
	inverse.resize(dim, dim);
	areaWeights.resize(dim, dim);

	// Construct normalized Haar wavelet transform (and its
	// inverse/transpose) in the specified dimension
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			xform(i, j) = (i == j ? 1.0f : 0.0f);
		}
	}

	// integer logarithm base 2.
	int iters = 0;
	int temp = dim;
	for (; temp != 1; iters++, temp >>= 1);


	for (int i = 0; i < iters; i++)
	{
		Iterate(i, xform);
	}

	for (int i = 0; i < xform.cols(); i++)
	{
		xform.col(i).normalize();
	}


	if (!xform.isUnitary())
	{
		cout << "Haar wavelet transform is not properly normalized." << endl;
		exit(1);
	}
	inverse = xform.transpose();

	// Compute wavelet area weights
	Diagonal(0, 0, dim, dim, 1);
	// Normalize them
	areaWeights *= (1.0f / (dim*dim)); 
}

void WaveletCompressor::Compress(const MatrixXf& dataIn, MatrixXf& dataOut)
{
	dataOut = inverse * dataIn * xform;
}

void WaveletCompressor::Decompress(const MatrixXf& dataIn, MatrixXf& dataOut)
{
	dataOut = xform * dataIn * inverse;
}

void WaveletCompressor::Iterate(int iter, MatrixXf& mat)
{
	MatrixXf next(dimension, dimension);
	// Set up identity matrix

	for (int i = 0; i < dimension; i++)
	for (int j = 0; j < dimension; j++)
		next(i, j) = (i == j ? 1.0f : 0.0f);
	
	
	// Set up Haar transform step
	int partdim = (iter > 0 ? dimension / (2 << (iter - 1)) : dimension);
	int pdby2 = partdim / 2;
	// First, the averaging
	for (int i = 0; i < partdim; i++)
	for (int j = 0; j < pdby2; j++)
		if (i == 2 * j || i == 2 * j + 1)
			next(i, j) = 0.5f;
		else next(i, j) = 0.0f;
	
	
	// Then, the differencing

	for (int i = 0; i < partdim; i++)
	for (int j = 0; j < pdby2; j++)
		if (i == 2 * j)
			next(i, j + pdby2) = 0.5f;
		else if (i == 2 * j + 1)
			next(i, j + pdby2) = -0.5f;
		else next(i, j + pdby2) = 0.0f;
	
	
	// Then accumulate
	mat = mat * next;
}

void WaveletCompressor::Diagonal(int imin, int jmin, int imax, int jmax, float weight)
{
	// Weight the lower-right corner
	for (int i = (imax + imin) / 2; i < imax; i++)
	{
		for (int j = (jmax + jmin) / 2; j < jmax; j++)
		{
			areaWeights(i, j) = weight;
		}
	}

	// Base case
	if (imin == imax && jmin == jmax) return;

	// Recurse on the other corners
	Horizontal(imin, (jmax + jmin) / 2, (imax + imin) / 2, jmax, weight * 2);
	Vertical((imax + imin) / 2, jmin, imax, (jmax + jmin) / 2, weight * 2);
	Diagonal(imin, jmin, (imax + imin) / 2, (jmax + jmin) / 2, weight * 4);
}

void WaveletCompressor::Horizontal(int imin, int jmin, int imax, int jmax, float weight)
{
	// Weight the bottom
	for (int i = (imax + imin) / 2; i < imax; i++)
	{
		for (int j = jmin; j < jmax; j++)
		{
			areaWeights(i, j) = weight;
		}
	}

	// Base case
	if (imin == imax) return;

	// Recurse
	Horizontal(imin, jmin, (imax + imin) / 2, jmax, weight * 2);
}

void WaveletCompressor::Vertical(int imin, int jmin, int imax, int jmax, float weight)
{
	// Weight the right side
	for (int i = imin; i < imax; i++)
	{
		for (int j = (jmin + jmax) / 2; j < jmax; j++)
		{
			areaWeights(i, j) = weight;
		}
	}

	// Base case
	if (jmin == jmax) return;

	// Recurse
	Vertical(imin, jmin, imax, (jmin + jmax) / 2, weight * 2);
}




///////////////////////////////////////////////////////////////////////
//                       WavTranEnvironment                          //
///////////////////////////////////////////////////////////////////////


void WavTranEnvironment::initRaw()
{
	if (isRawInit)
		return;

	isRawInit = true;
	initBase();


	int csl2 = cubeSideLength*cubeSideLength;
	rawLightVec[0].resize(6 * csl2);
	rawLightVec[1].resize(6 * csl2);
	rawLightVec[2].resize(6 * csl2);

	// Reset the rotation
	rot.setIdentity();
	RotateRaw(Eigen::Quaternionf::Identity());
}

void WavTranEnvironment::initWav()
{
	if (isWavInit)
		return;
	isWavInit = true;
	initBase();

	// Need to set size of lightVec
	int csl2 = cubeSideLength*cubeSideLength;
	wavLightVec[0].resize(6 * csl2);
	wavLightVec[1].resize(6 * csl2);
	wavLightVec[2].resize(6 * csl2);

	for (int f = 0; f < 6; f++)
	{
		waveletFaces[f][0].resize(cubeSideLength, cubeSideLength);
		waveletFaces[f][1].resize(cubeSideLength, cubeSideLength);
		waveletFaces[f][2].resize(cubeSideLength, cubeSideLength);
	}

	// Reset the rotation
	rot.setIdentity();
	RotateWav(Eigen::Quaternionf::Identity());
}

void WavTranEnvironment::RotateRaw(Eigen::Quaternionf& q)
{
	initRaw();
	rawHasRotted = true;
	rot = q * rot;
	Matrix3f rotMat = rot.toRotationMatrix();


	// Do rotation, write back to light vector!
	int csl2 = cubeSideLength*cubeSideLength;
	float gamma = 1.0 / (samplesPerPixel*(6 * csl2));
	for (int f = 0; f < 6; f++)
	{
#pragma omp parallel for
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				int vecI = f*csl2 + y*cubeSideLength + x;
				rawLightVec[0][vecI] = rawLightVec[1][vecI] = rawLightVec[2][vecI] = 0.0f;

				// Use each precomputed sample direction, index
				// into high-res cube map (super sampling)
				for (int s = 0; s < samplesPerPixel; s++)
				{
					Vector3f dir;
					dir = rotMat * sampleDirs[f][y][x][s];
					int newf;
					float newu, newv;
					DirToCubeCoord(newf, newu, newv, dir);
					int j = (int)(newu*(sourceSideLength - 1));
					int i = (int)(newv*(sourceSideLength - 1));
					// Normalized for numerical cubature
					rawLightVec[0][vecI] += faces[newf][0](i, j);
					rawLightVec[1][vecI] += faces[newf][1](i, j);
					rawLightVec[2][vecI] += faces[newf][2](i, j);
				}

				rawLightVec[0][vecI] *= gamma;
				rawLightVec[1][vecI] *= gamma;
				rawLightVec[2][vecI] *= gamma;
			}
		}
	}
}

void WavTranEnvironment::RotateWav(Eigen::Quaternionf& q)
{
	initWav();
	wavHasRotted = true;
	rot = q * rot;
	Matrix3f rotMat = rot.toRotationMatrix();

	// Do rotation, write to waveletFaces
	int csl2 = cubeSideLength*cubeSideLength;
	int colorIndex = 0;


	float gamma = 1.0 / samplesPerPixel;

	for (int f = 0; f < 6; f++)
	{
#pragma omp parallel for
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++, colorIndex++)
			{
				waveletFaces[f][0](y, x) = 0.0;
				waveletFaces[f][1](y, x) = 0.0;
				waveletFaces[f][2](y, x) = 0.0;

				// Use each precomputed sample direction, index
				// into high-res cube map (super sampling)
				for (int s = 0; s < samplesPerPixel; s++)
				{
					Vector3f dir;
					dir = rotMat * sampleDirs[f][y][x][s];
					int newf;
					float newu, newv;
					DirToCubeCoord(newf, newu, newv, dir);
					int j = (int)(newu*(sourceSideLength - 1));
					int i = (int)(newv*(sourceSideLength - 1));
					// Normalized for numerical cubature

					waveletFaces[f][0](y, x) += faces[newf][0](i, j);
					waveletFaces[f][1](y, x) += faces[newf][1](i, j);
					waveletFaces[f][2](y, x) += faces[newf][2](i, j);
				}
			}
		}
		// Then compress
#pragma omp parallel for 
		for (int c = 0; c < 3; c++)
		{
			waveletFaces[f][c] *= gamma;
			comp.Compress(waveletFaces[f][c], waveletFaces[f][c]);
		}
	}




	// Then write to lightVec
	//FacesToLightVec

	wavLightVec[0].setZero();
	wavLightVec[1].setZero();
	wavLightVec[2].setZero();

	gamma = 1.0 / (6 * csl2);

	if (transWeightSelect)
	{
		// Select numWaveletLights values, sorted by absolute value
		// and weighted by transport energy.
		int numSelected = 0;
		float** vals = new float*[numWaveletLights];
		int* indices = new int[numWaveletLights];
		for (int i = 0; i < numWaveletLights; i++)
		{
			vals[i] = new float[3];
		}

		for (int f = 0; f < 6; f++)
		{
			for (int y = 0; y < cubeSideLength; y++)
			{
				for (int x = 0; x < cubeSideLength; x++)
				{
					int i = f*csl2 + y*cubeSideLength + x;
					if (numSelected < numWaveletLights)
					{
						vals[numSelected][0] = waveletFaces[f][0](y, x);
						vals[numSelected][1] = waveletFaces[f][1](y, x);
						vals[numSelected][2] = waveletFaces[f][2](y, x);
						indices[numSelected] = i;
						numSelected++;
					}
					else for (int j = 0; j < numWaveletLights; j++)
					{
						if (IsBetterWaveletLight(waveletFaces[f], i, j, x, y, vals[j]))
						{
							vals[j][0] = waveletFaces[f][0](y, x);
							vals[j][1] = waveletFaces[f][1](y, x);
							vals[j][2] = waveletFaces[f][2](y, x);
							indices[j] = i;
							break;
						}
					}
				}
			}
		}

		// Finally, write stuff into the lightVec
		for (int i = 0; i < numWaveletLights; i++)
		{
			wavLightVec[0].fill(indices[i]) = vals[i][0];
			wavLightVec[1].fill(indices[i]) = vals[i][1];
			wavLightVec[2].fill(indices[i]) = vals[i][2];
		}

		// Clean up
		for (int i = 0; i < numWaveletLights; i++)
		{
			delete[] vals[i];
		}
		delete[] vals;
		delete[] indices;
	}
	else
	{
		int vecI = 0;
		// Sparsify by rejecting elements with fabs < threshold
		for (int f = 0; f < 6; f++)
		{
			for (int y = 0; y < cubeSideLength; y++)
			{
				for (int x = 0; x < cubeSideLength; x++)
				{
					float red, green, blue;
					red = waveletFaces[f][0](y, x);
					green = waveletFaces[f][1](y, x);
					blue = waveletFaces[f][2](y, x);

					// (Optionally), do area weighting
					float weight = 1.0f;
					if (useAreaWeights)
					{
						weight = comp.areaWeights(y, x);
					}

					if (weight*ColorNorm(red, green, blue) >= threshold)
					{
						wavLightVec[0].fill(vecI) = red * gamma;
						wavLightVec[1].fill(vecI) = green * gamma;
						wavLightVec[2].fill(vecI) = blue * gamma;
					}
					else
					{
						waveletFaces[f][0](y, x) = 0.f;
						waveletFaces[f][1](y, x) = 0.f;
						waveletFaces[f][2](y, x) = 0.f;
					}

					//if (fabs(red) >= threshold) lightVec[0].fill(vecI) = red;
					//if (fabs(green) >= threshold) lightVec[1].fill(vecI) = green;
					//if (fabs(blue) >= threshold) lightVec[2].fill(vecI) = blue;
					vecI++;
				}
			}

			// Then Decompress
#pragma omp parallel for
			for (int c = 0; c < 3; c++)
			{
				comp.Decompress(waveletFaces[f][c], waveletFaces[f][c]);
			}
		}
	}


	int numNonZero = wavLightVec[0].nonZeros();
	int percentNonZero = (int)(100 * numNonZero / (float)(6 * csl2));
	cout << "Num Wavelet Lights: " << numNonZero << " (" << percentNonZero << "%)" << endl;

}

void WavTranEnvironment::initTextureRaw(unsigned int texCubeMap[])
{

	initRaw();
	int csl2 = cubeSideLength*cubeSideLength * 3;
	float gamma = 2 * csl2;
	float* pixel = new float[csl2];

	int vecInd = 0;
	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int rgbInd = 0; rgbInd < 3; rgbInd++)
					pixel[pixelIndex++] = rawLightVec[rgbInd][vecInd] * gamma;
				vecInd++;
			}
		}
		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, cubeSideLength, cubeSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}
	delete[] pixel;
}

void WavTranEnvironment::initTextureWav(unsigned int texCubeMap[])
{
	initWav();

	int csl2 = cubeSideLength*cubeSideLength * 3;
	float* pixel = new float[csl2];

	int colorIndex = 0;
	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int rgbInd = 0; rgbInd < 3; rgbInd++)
					pixel[pixelIndex++] = waveletFaces[f][rgbInd](y, x);
				colorIndex++;
			}
		}

		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, cubeSideLength, cubeSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}

	delete[] pixel;
}

void WavTranEnvironment::initTextureAvg(unsigned int texCubeMap[])
{
	initAverageLightVec();
	int csl2 = cubeSideLength*cubeSideLength * 3;
	float gamma = 2 * csl2;
	float* pixel = new float[csl2];

	int vecInd = 0;
	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int rgbInd = 0; rgbInd < 3; rgbInd++)
					pixel[pixelIndex++] = rawAverageLightVec[rgbInd][vecInd] * gamma;
				vecInd++;
			}
		}
		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, cubeSideLength, cubeSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}
	delete[] pixel;
}

void WavTranEnvironment::initTextureScaledRaw(unsigned int texCubeMap[])
{

	initRaw();
	int csl2 = cubeSideLength*cubeSideLength * 3;
	float gamma = 2 * csl2;
	float* pixel = new float[csl2];

	int vecInd = 0;
	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int rgbInd = 0; rgbInd < 3; rgbInd++)
					pixel[pixelIndex++] = rawScaledLightVec[rgbInd][vecInd] *gamma;
				vecInd++;
			}
		}
		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, cubeSideLength, cubeSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}
	delete[] pixel;
}

void WavTranEnvironment::initTextureScaledWav(unsigned int texCubeMap[])
{
	initWav();

	int csl2 = cubeSideLength*cubeSideLength * 3;
	float* pixel = new float[csl2];

	int colorIndex = 0;
	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int rgbInd = 0; rgbInd < 3; rgbInd++)
					pixel[pixelIndex++] = waveletFaces[f][rgbInd](y, x);
				colorIndex++;
			}
		}

		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, cubeSideLength, cubeSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}

	delete[] pixel;
}


void WavTranEnvironment::initAverageLightVec()
{

	if (isAverageLightInit)
		return;
	isAverageLightInit = true;


	int fcsl2 = cubeSideLength*cubeSideLength * 6;

	initRaw();

	for (int i = 0; i < 3; i++)
		envScale[i] = vec3(0.f);
	float distance[3] = { 0.f };

	/*
	float *luminance = new float[fcsl2];

	for (int vecI = 0; vecI < fcsl2; vecI++)
	{
		luminance[vecI] = 0.2989*rawLightVec[0][vecI] + 0.587*rawLightVec[1][vecI] + 0.114*rawLightVec[2][vecI];

		rawAverageLight += luminance[vecI];
	}

	rawAverageLight = 0.5f;
	rawAverageLight /= (float)fcsl2;
	*/


	rawAverageLightVec[0].resize(fcsl2);
	rawAverageLightVec[1].resize(fcsl2);
	rawAverageLightVec[2].resize(fcsl2);
	


	// Do rotation, write back to light vector!
	int csl2 = cubeSideLength*cubeSideLength;
	float gamma = 1.0 / (samplesPerPixel*(6 * csl2));
	for (int f = 0; f < 6; f++)
	{
#pragma omp parallel for
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				int vecI = f*csl2 + y*cubeSideLength + x;
				rawAverageLightVec[0][vecI] = rawAverageLightVec[1][vecI] = rawAverageLightVec[2][vecI] = 0.0f;

				// Use each precomputed sample direction, index
				// into high-res cube map (super sampling)
				for (int s = 0; s < samplesPerPixel; s++)
				{

					vec3 tempColor;


					Vector3f dir;
					dir = sampleDirs[f][y][x][s];
					int newf;
					float newu, newv;
					DirToCubeCoord(newf, newu, newv, dir);
					int j = (int)(newu*(sourceSideLength - 1));
					int i = (int)(newv*(sourceSideLength - 1));
					// Normalized for numerical cubature

					tempColor[0] = faces[newf][0](i, j);
					tempColor[1] = faces[newf][1](i, j);
					tempColor[2] = faces[newf][2](i, j);




					float curDis = 0.f;
					float maxCol = 0.f;
					int curInd = 0;

					for (int cIndex = 0; cIndex < 3; cIndex++)
					{
						if (maxCol < tempColor[cIndex])
						{
							curInd = cIndex;
							maxCol = tempColor[cIndex];
						}
					}

					maxCol = 1.f / maxCol;
					tempColor *= maxCol;

					for (int cIndex = 0; cIndex < 3; cIndex++)
					{
						if (cIndex == curInd)
							continue;

						curDis += (1.f - tempColor[cIndex]) * (1.f - tempColor[cIndex]);
					}

					if (distance[curInd] < curDis)
					{
						distance[curInd] = curDis;
						envScale[curInd] = tempColor;
					}

					rawAverageLightVec[0][vecI] += faces[newf][0](i, j);
					rawAverageLightVec[1][vecI] += faces[newf][1](i, j);
					rawAverageLightVec[2][vecI] += faces[newf][2](i, j);
				}



				rawAverageLightVec[0][vecI] *= gamma;
				rawAverageLightVec[1][vecI] *= gamma;
				rawAverageLightVec[2][vecI] *= gamma;
			}
		}
	}


	/*
	for (int i = 0; i < 3; i++)
	{
		envScale[i] = vec3(0.f);
		distance[i] = 0.f;
	}


	for (int vecI = 0; vecI < csl2 * 6; vecI++)
	{

		vec3 tempColor;
		tempColor[0] = rawAverageLightVec[0][vecI];
		tempColor[1] = rawAverageLightVec[1][vecI];
		tempColor[2] = rawAverageLightVec[2][vecI];

		float curDis = 0.f;
		float maxCol = 0.f;
		int curInd = 0;

		for (int cIndex = 0; cIndex < 3; cIndex++)
		{
			if (maxCol < tempColor[cIndex])
			{
				curInd = cIndex;
				maxCol = tempColor[cIndex];
			}
		}

		maxCol = 1.f / maxCol;
		tempColor *= maxCol;


		for (int cIndex = 0; cIndex < 3; cIndex++)
		{
			if (cIndex == curInd)
				continue;

			curDis += (1.f - tempColor[cIndex]) * (1.f - tempColor[cIndex]);
		}

		if (distance[curInd] < curDis)
		{
			distance[curInd] = curDis;
			envScale[curInd] = tempColor;
		}

	}

	*/


//	delete[] luminance;

	return ;

}

void WavTranEnvironment::outputEnvScale()
{
	initAverageLightVec();

	char str[180];

	strcpy_s(str, 180, "..\\..\\matlabSoveNonlinearEquations\\");
	strcat_s(str, 180, fileName);
	strcat_s(str, 180, "_scale.obj");

	ofstream cluster_of(str);

	if (cluster_of.is_open())
	{
		for (int curInd = 0; curInd < 3; curInd++)
		{
			for (int cIndex = 0; cIndex < 3; cIndex++)
				cluster_of << envScale[curInd][cIndex] << " ";
			cluster_of << endl;
		}
	}

	cluster_of.close();
}

void WavTranEnvironment::initScaledLights()
{
	
	if (isScaledLightInit)
		return;

	isScaledLightInit = true;

	int csl2 = cubeSideLength*cubeSideLength;
	rawScaledLightVec[0].resize(6 * csl2);
	rawScaledLightVec[1].resize(6 * csl2);
	rawScaledLightVec[2].resize(6 * csl2);

}

void WavTranEnvironment::initScaledLightsRefract(int vertexNum)
{


	if (isScaledLightRefractInit)
		return;

	isScaledLightRefractInit = true;

	int csl2 = vertexNum * vertexNum;
	rawScaledRefractLightVecRefract[0].resize(vertexNum);
	rawScaledRefractLightVecRefract[1].resize(vertexNum);
	rawScaledRefractLightVecRefract[2].resize(vertexNum);

	rawMidLightVecRefract.resize(vertexNum);
	lightNearDistance.resize(vertexNum);


}

float WavTranEnvironment::getRawLightAvg()
{
	if (rawAverageLight != 0.f)
		return rawAverageLight;

	initAverageLightVec();

	return rawAverageLight;


}



void WavTranEnvironment::changeWaveEnvironment(char* newEvn)
{
	changeEnvironmentMap(newEvn);

	isRawInit = false;
	isWavInit = false;
	isAverageLightInit = false;
	isScaledLightInit = false;
	isScaledLightRefractInit = false;


}



///////////////////////////////////////////////////////////////////////
//                             TRANSPORT                             //
///////////////////////////////////////////////////////////////////////

Transport::~Transport()
{

	if (!wavNeedInitRGB)
	{
		for (int rgbInd = 1; rgbInd<RGBnum; rgbInd++)
			delete[] colEnergy[rgbInd];
	}
	if (isWavInit)
	{
		delete[] colEnergy[0];
		delete[] colEnergy;
		delete[] resultColor;
	}
}

void Transport::getVertexLightOcl(bool*** occluded, int face, BVH& meshBVH)
{

	int testVN;

	FILE *fin = fopen(meshes[0]->dFName.getVertexLightOclDN(cubeSideLength, face), "rb");
	if (fin != NULL)
	{
		fread(&testVN, sizeof(int), 1, fin);
		if (testVN == vertexNum)
		{
			std::cout << "Reading visibility..." << endl;
			for (int p = 0; p < vertexNum; p++)
			for (int x = 0; x < cubeSideLength; x++)
				fread(occluded[p][x], sizeof(bool), cubeSideLength, fin);
			fclose(fin);
			return;
		}
		fclose(fin);
	}

	/*
	ifstream VertOcl_if(meshes[0]->dFName.getVertexLightOclDN(cubeSideLength, face));
	if (VertOcl_if.is_open())
	{
		VertOcl_if>>testVN;
		if (testVN == vertexNum)
		{
			std::cout << "Reading visibility..."<<endl;
			for (int x = 0; x < cubeSideLength; x++)
				for (int y = 0; y < cubeSideLength; y++)
					for (int p = 0; p < vertexNum; p++)
						VertOcl_if>>occluded[p][x][y];
			VertOcl_if.close();
			return;
		}
		VertOcl_if.close();
	}
	*/


	std::cout << "Calculating visibility..." <<endl;

	vector<vec3> &currentObjectNormal=meshes[0]->normals;
	vector<vec3> &currentObject=meshes[0]->vertices;

	for (int x = 0; x < cubeSideLength; x++)
	{
		for (int y = 0; y < cubeSideLength; y++)
		{
			vec3& lightDir = cubeDirections.getDirection(face,x,y);
			for (int p = 0; p < vertexNum; p++)
			{
				//Calculate cosine term for this sample
				if(lightDir.dot(currentObjectNormal[p]) <= 0.f)
					occluded[p][x][y] = true;
				else
				{
					//Fill in a RAY structure for this sample
					//See if the ray is blocked by any object
					Ray ray(currentObject[p]+2*EPSILON*currentObjectNormal[p],lightDir);
					occluded[p][x][y] = meshBVH.Occluded(ray);
				}
			}
		}
		std::cout<<x+1<<"/"<<cubeSideLength<<"\r";
	}


	std::cout << "Writing visibility..." << endl;

	FILE *fout = fopen(meshes[0]->dFName.getVertexLightOclDN(cubeSideLength, face), "wb");
	float a = sizeof(bool);

	fwrite(&vertexNum, sizeof(int), 1, fout);
	for (int p = 0; p < vertexNum; p++)
	for (int x = 0; x < cubeSideLength; x++)
		fwrite(occluded[p][x], sizeof(bool), cubeSideLength, fout);

	fclose(fout);

	/*
	ofstream VertOcl_of(meshes[0]->dFName.getVertexLightOclDN(cubeSideLength, face));
	VertOcl_of<<vertexNum<<endl;
	for (int x = 0; x < cubeSideLength; x++)
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int p = 0; p < vertexNum; p++)
				VertOcl_of<<occluded[p][x][y]<<" ";
			VertOcl_of<<endl;
		}
	VertOcl_of.close();
	*/
	return;


}

void Transport::calDirectMat(MatrixXf(*matrices)[3], int face, BVH& meshBVH, bool*** vertLitOcl, bool isUnOclReflect)
{




	float(*materialShader)(vec3 normal, vec3 lightDir, vec3 viewDir);

	switch (curMaterial)
	{
	case diffuse:
		materialShader = diffuseShader;
		break;
	case ward:
		materialShader = wardShader;
		break;

	default:
		materialShader = diffuseShader;
		break;
	}



	// For each pixel, form color channel matrices for the face
	std::cout << "Constructing face matrices..." << face+1 << "/6..."<<endl;

	if (!isUnOclReflect)
		getVertexLightOcl(vertLitOcl, face, meshBVH);

	vector<vec3> &currentObjectNormal=meshes[0]->normals;
	vector<vec3> &currentObject=meshes[0]->vertices;
	vec3 eye = meshes[0]->eyePosition;

	for (int x = 0; x < cubeSideLength; x++)
	{
		for (int y = 0; y < cubeSideLength; y++)
			for (int p = 0; p < vertexNum; p++)
				if ((!isUnOclReflect) && vertLitOcl[p][x][y])
					for (int rgbInd=0; rgbInd<RGBnum; rgbInd++)
						matrices[p][rgbInd](y,x) = 0.f;
				else
				{
					point& pos = currentObject[p];
					vec3 viewDir = eye - pos;
					normalize(viewDir);
				
					for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
						matrices[p][rgbInd](y, x) = materialShader(currentObjectNormal[p], cubeDirections.getDirection(face, x, y), viewDir);
				}
		std::cout<<x+1<<"/"<<cubeSideLength<<"\r";
	}
}

void Transport::calIndirectMat(MatrixXf(*indirectMatrices)[3], int face, BVH& meshBVH)
{
	// For each pixel, form color channel matrices for the face
	std::cout << "Constructing face matrices..." << face+1 << "/6..."<<endl;

	vector<vec3> &currentObjectNormal=meshes[0]->normals;
	vector<vec3> &currentObject=meshes[0]->vertices;


	MatrixXf (*matrices)[3] = new MatrixXf[vertexNum][3];
	for (int i = 0; i < vertexNum; i++)
		for (int rgbInd=0; rgbInd<RGBnum; rgbInd++)
			matrices[i][rgbInd].resize(cubeSideLength, cubeSideLength);

		getTranMat(matrices, face, false, false, false);
	for (int i = 0; i < vertexNum; i++)
		for (int x = 0; x < cubeSideLength; x++)
			for (int y = 0; y < cubeSideLength; y++)
				for (int rgbInd=0; rgbInd<RGBnum; rgbInd++)
					indirectMatrices[i][rgbInd](y,x)=matrices[i][rgbInd](y,x);







	float(*materialShader)(vec3 normal, vec3 lightDir, vec3 viewDir);

	switch (curMaterial)
	{
	case diffuse:
		materialShader = diffuseShader;
		break;
	case ward:
		materialShader = wardShader;
		break;

	default:
		materialShader = diffuseShader;
		break;
	}



	vec3 eye = meshes[0]->eyePosition;


	for (int x = 0; x < cubeSideLength; x++)
	{
		for (int y = 0; y < cubeSideLength; y++)
		{
			vec3& lightDir = cubeDirections.getDirection(face,x,y);
			for (int sVInd = 0; sVInd < vertexNum; sVInd++)
			{
				float radiance = matrices[sVInd][0](y, x);
				point& sPos = currentObject[sVInd];
				vec3& sNor = currentObjectNormal[sVInd];

				if (matrices[sVInd][0](y, x) <= 0.f)
					continue;
				for (int eVInd = 0; eVInd < vertexNum; eVInd++)
				{
					if (eVInd == sVInd || meshBVH.IsVertVertOutsideOccluded(sVInd, eVInd))
						continue;
					{
						vec3& eNor = currentObjectNormal[eVInd];
						point& ePos = currentObject[eVInd];
						vec3 indirectLightDir = ePos - sPos;
						normalize(indirectLightDir);
						float indirectRadiace = materialShader(sNor, lightDir, indirectLightDir);


						vec3 viewDir = eye - ePos;
						normalize(viewDir);


						float viewInRa = materialShader(eNor, -indirectLightDir, viewDir)*indirectRadiace*radiance;
						for (int rgbInd=0; rgbInd<RGBnum; rgbInd++)
							indirectMatrices[eVInd][rgbInd](y,x) += viewInRa;
					}
				}
			}
			std::cout<<y+1<<"/"<<cubeSideLength<<"::"<<x+1<<"/"<<cubeSideLength<<"\r";
		}
	}
	delete[] matrices;
}

void Transport::calRefractionMat(MatrixXf(*refractionMatrices)[3], int face, BVH& meshBVH, bool isGlossy)
{
	RefractionParemeter rp(eta, sigma_a, sigma_s);

	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;

	// For each pixel, form color channel matrices for the face
	std::cout << "Constructing face matrices..." << face + 1 << "/6..." << endl;

	vector<vec3> &currentObjectNormal = meshes[0]->normals;
	vector<vec3> &currentObject = meshes[0]->vertices;

	MatrixXf(*matrices)[3] = new MatrixXf[vertexNum][3];
	for (int i = 0; i < vertexNum; i++)
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		matrices[i][rgbInd].resize(cubeSideLength, cubeSideLength);

	getTranMat(matrices, face, false, false, false);
	for (int i = 0; i < vertexNum; i++)
	for (int x = 0; x < cubeSideLength; x++)
	for (int y = 0; y < cubeSideLength; y++)
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		refractionMatrices[i][rgbInd](y, x) = 0.f;


	vec3 temp;
	float a = 0.5/meshes[0]->bsphere.r;


	for (int x = 0; x < cubeSideLength; x++)
	{
		for (int y = 0; y < cubeSideLength; y++)
		{
			vec3& lightDir = cubeDirections.getDirection(face, x, y);
			for (int sVInd = 0; sVInd < vertexNum; sVInd++)
			{
				float lightRadiance = matrices[sVInd][0](y, x);
				point& sPos = currentObject[sVInd];
				vec3& sNor = currentObjectNormal[sVInd];

//				if (matrices[sVInd][0](y, x) <= 0.f)
//					continue;

				float vertexRadiance = lightDir.dot(sNor) * lightRadiance;
								
				for (int eVInd = 0; eVInd < vertexNum; eVInd++)
				{
					if (false && meshBVH.IsVertVertInsideOccluded(sVInd, eVInd))
						continue;
					{
						vec3& eNor = currentObjectNormal[eVInd];
						point& ePos = currentObject[eVInd];

						float surfaceRadiance = translucent(sPos, ePos, temp, temp, temp, temp, false, rp, a);
						
						for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
							refractionMatrices[eVInd][rgbInd](y, x) += surfaceRadiance;
					}
				}
				
			}
			std::cout << y + 1 << "/" << cubeSideLength << "::" << x + 1 << "/" << cubeSideLength << "\r";
		}
	}
	delete[] matrices;
}

void Transport::getRefractionVertexAsLightMat(MatrixXf *refractionMatrices, bool isGlossy)
{
	bool needCal = true;
	int testVN;


	FILE *fin = fopen(meshes[0]->dFName.getTranMatRefractVertexAsLightDN(RGBnum, isGlossy), "rb");
	if (fin != NULL)
	{
		fread(&testVN, sizeof(int), 1, fin);
		if (testVN == vertexNum)
			needCal = false;
		fclose(fin);
	}

	/*
	ifstream TranMat_if(meshes[0]->dFName.getTranMatRefractVertexAsLightDN(RGBnum, isGlossy));
	if (TranMat_if.is_open())
	{
		TranMat_if >> testVN;
		if (testVN == vertexNum)
			needCal = false;
		TranMat_if.close();
	}
	*/

	vec3 temp;
	float *s = new float[vertexNum];
	//calculate all
	if (needCal)
	{
		RefractionParemeter rp(eta, sigma_a, sigma_s);

		// For each pixel, form color channel matrices for the face

		vector<vec3> &currentObject = meshes[0]->vertices;
		
		float a = 0.5 / meshes[0]->bsphere.r;


		for (int lInd = 0; lInd < vertexNum; lInd++)
		{
			point& lPos = currentObject[lInd];
			for (int vInd = 0; vInd < vertexNum; vInd++)
			{

				point& vPos = currentObject[vInd];

				float surfaceRadiance = translucent(lPos, vPos, temp, temp, temp, temp, false, rp, a);

				for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
					refractionMatrices[rgbInd](vInd, lInd) = surfaceRadiance;

			}
		}

		// For each pixel, form color channel matrices for the face
		std::cout << "Writing..." << endl;

		FILE *fout = fopen(meshes[0]->dFName.getTranMatRefractVertexAsLightDN(RGBnum, isGlossy), "wb");
		fwrite(&vertexNum, sizeof(int), 1, fout);

		for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		for (int lInd = 0; lInd < vertexNum; lInd++)
		{
			for (int vInd = 0; vInd < vertexNum; vInd++)
					s[vInd] = refractionMatrices[rgbInd](vInd, lInd);
			fwrite(s, sizeof(float), vertexNum, fout);
			std::cout << lInd << "/" << vertexNum << "\r";
		}
		fclose(fout);
		/*
		ofstream TranMat_of(meshes[0]->dFName.getTranMatRefractVertexAsLightDN(RGBnum, isGlossy));
		TranMat_of << vertexNum << endl;
		for (int lInd = 0; lInd < vertexNum; lInd++)
		{
			for (int vInd = 0; vInd < vertexNum; vInd++)
			{
				for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
					TranMat_of << setiosflags(ios::fixed) << refractionMatrices[rgbInd](vInd, lInd) << " ";
			}
			TranMat_of << endl;
			std::cout << lInd << "/" << vertexNum << "\r";
		}
		TranMat_of.close();
		*/
	}

	//read current
	// For each pixel, form color channel matrices for the face


	fin = fopen(meshes[0]->dFName.getTranMatRefractVertexAsLightDN(RGBnum, isGlossy), "rb");


	fread(&testVN, sizeof(int), 1, fin);
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
	for (int lInd = 0; lInd < vertexNum; lInd++)
	{
		fread(s, sizeof(float), vertexNum, fin);
		for (int vInd = 0; vInd < vertexNum; vInd++)
			refractionMatrices[rgbInd](vInd, lInd) = s[vInd];
	}
	fclose(fin);

	delete[] s;
	/*
	TranMat_if.open(meshes[0]->dFName.getTranMatRefractVertexAsLightDN(RGBnum, isGlossy));
	TranMat_if >> testVN;
	for (int lInd = 0; lInd < vertexNum; lInd++)
	for (int vInd = 0; vInd < vertexNum; vInd++)
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		TranMat_if >> refractionMatrices[rgbInd](vInd, lInd);
	*/
}

void Transport::getTranMat(MatrixXf(*matrices)[3], int face, bool isIndirect, bool isRefraction, bool isUnOclReflect)
{

	bool needCal = true;
	int testVN;


	FILE *fin = fopen(meshes[0]->dFName.getTranMatDN(cubeSideLength, face, isIndirect, isRefraction, RGBnum, isUnOclReflect, curMaterial, curEye), "rb");
	if (fin!=NULL)
	{
		fread(&testVN, sizeof(int), 1, fin);
		if (testVN == vertexNum)
			needCal = false;
		fclose(fin);
	}
	/*
	ifstream TranMat_if(meshes[0]->dFName.getTranMatDN(cubeSideLength, face, isIndirect, isRefraction, RGBnum, isGlossy));
	if (TranMat_if.is_open())
	{
		TranMat_if>>testVN;
		if (testVN == vertexNum)
			needCal = false;
		TranMat_if.close();
	}
	*/

	float *s = new float[vertexNum];
	//calculate all
	if (needCal)
	{

		bool*** vertLitOcl = NULL;
		if (!isIndirect && !isRefraction && !isUnOclReflect)
		{
			vertLitOcl = new bool**[vertexNum];
			for (int p = 0; p < vertexNum; p++)
			{
				vertLitOcl[p] = new bool*[cubeSideLength];
				for (int x = 0; x < cubeSideLength; x++)
					vertLitOcl[p][x] = new bool[cubeSideLength];
			}
		}

		BVH meshBVH(meshes[0]);
		for (int f=0; f<6; f++)
		{
			if (isIndirect)
				calIndirectMat(matrices, f, meshBVH);
			else if (isRefraction)
				calRefractionMat(matrices, f, meshBVH, false);
			else
				calDirectMat(matrices, f, meshBVH, vertLitOcl, isUnOclReflect);

			// For each pixel, form color channel matrices for the face
			std::cout << "Writing..." << endl;

			FILE *fout = fopen(meshes[0]->dFName.getTranMatDN(cubeSideLength, f, isIndirect, isRefraction, RGBnum, isUnOclReflect, curMaterial, curEye), "wb");
			fwrite(&vertexNum, sizeof(int), 1, fout);
			for (int rgbInd = 0; rgbInd < RGBnum; rgbInd++)
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int y = 0; y < cubeSideLength; y++)
				{
					for (int p = 0; p < vertexNum; p++)
						s[p] = matrices[p][rgbInd](y, x);
					fwrite(s, sizeof(float), vertexNum, fout);
				}
				std::cout << x << "/" << cubeSideLength << "\r";
			}
			fclose(fout);
			/*
			
			ofstream TranMat_of(meshes[0]->dFName.getTranMatDN(cubeSideLength, f+6, isIndirect, isRefraction, RGBnum, isGlossy));
			TranMat_of<<vertexNum<<endl;
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int y = 0; y < cubeSideLength; y++)
				{
					for (int p = 0; p < vertexNum; p++)
						for (int rgbInd=0; rgbInd<RGBnum; rgbInd++)
		
							TranMat_of << setiosflags(ios::fixed) << matrices[p][rgbInd](y, x) << " ";
					TranMat_of<<endl;
				}
				std::cout<<x<<"/"<<cubeSideLength<<"\r";
			}
			TranMat_of.close();
			*/
		}

		if (!isIndirect && !isRefraction && !isUnOclReflect)
		{
			for (int p = 0; p < vertexNum; p++)
			{
				for (int x = 0; x < cubeSideLength; x++)
					delete[] vertLitOcl[p][x];
				delete[] vertLitOcl[p];
			}
			delete[] vertLitOcl;
		}
	}

	//read current
	// For each pixel, form color channel matrices for the face
	std::cout << "Reading face"<<face<<" matrices..." << endl;


	
	fin = fopen(meshes[0]->dFName.getTranMatDN(cubeSideLength, face, isIndirect, isRefraction, RGBnum, isUnOclReflect, curMaterial, curEye), "rb");


	fread(&vertexNum, sizeof(int), 1, fin);
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
	for (int x = 0; x < cubeSideLength; x++)
	for (int y = 0; y < cubeSideLength; y++)
	{
		fread(s, sizeof(float), vertexNum, fin);
		for (int p = 0; p < vertexNum; p++)
			matrices[p][rgbInd](y, x) = s[p];
	}
	fclose(fin);


	delete[] s;

	/*
	ifstream TranMat_if;
	TranMat_if.open(meshes[0]->dFName.getTranMatDN(cubeSideLength, face + 6, isIndirect, isRefraction, RGBnum, isGlossy));
	TranMat_if>>testVN;
	for (int x = 0; x < cubeSideLength; x++)
		for (int y = 0; y < cubeSideLength; y++)
			for (int p = 0; p < vertexNum; p++)
			for (int rgbInd = 0; rgbInd < RGBnum; rgbInd++)
			{
				cout << x << y << p << " ";
				cout << matrices[p][rgbInd](y, x)<<" ";
				TranMat_if >> matrices[p][rgbInd](y, x);
				cout << matrices[p][rgbInd](y, x) << " ";
				getchar();
			}
	*/
}

void Transport::compressTranMat(bool isGlossy)
{

	std::cout << "COMPRESSING Wavelet BASIS" << endl;

	MatrixXf(*matrices)[3] = new MatrixXf[vertexNum][3];
	for (int i = 0; i < vertexNum; i++)
		for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
			matrices[i][rgbInd].resize(cubeSideLength, cubeSideLength);

	// Wavelet transform the matrices
	WaveletCompressor comp(cubeSideLength);

	ofstream WaveCoeff_of(meshes[0]->dFName.getWaveCofDN(cubeSideLength, wavUsingIndirect, RGBnum, isGlossy));
	WaveCoeff_of << vertexNum << endl;

	// Compress each face, one at a time
	for (int f = 0; f < 6; f++)
	{
		std::cout << "Compressing face " << f + 1 << "/6..." << endl;

		getTranMat(matrices, f, wavUsingIndirect, wavUsingRefraction, false);

		std::cout << "Compressing... " << endl;
		for (int p = 0; p < vertexNum; p++)
			for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
				comp.Compress(matrices[p][rgbInd], matrices[p][rgbInd]);

		std::cout << "Writing file..." << endl;
		// Overwrite old floatmaps;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				for (int p = 0; p < vertexNum; p++)
					for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
						WaveCoeff_of << setiosflags(ios::fixed) << matrices[p][rgbInd](y, x) << " ";
				WaveCoeff_of << endl;
			}
			std::cout << y + 1 << "/" << cubeSideLength << "\r";
		}
		WaveCoeff_of << endl;
	}

	delete[] matrices;

	std::cout << "All done!" << endl;

}

//get the transport matrix
void Transport::LoadRaw(bool isRefracetVertex, bool isRefracet, bool isReflect, bool isIndirect, bool isRGB, bool isRefractGlossy, bool isUnOclReflect, MaterialType useMtaerial, EyeType useEye)
{
	if (isIndirect == rawUsingIndirect && isRGB == rawUsingRGB && rawUsingRefVertex == isRefracetVertex && isRawRefractInit && isRawReflectInit && useMtaerial == curMaterial && curEye == useEye)
		return;


	if (isRefracetVertex != rawUsingRefVertex || rawRefractGlossy != isRefractGlossy)
	{
		isRawRefractInit = false;
	}

	if (useMtaerial != curMaterial)
	{
		needRelight = true;
		isRawReflectInit = false;
		curMaterial = useMtaerial;
	}
	
	if (curMaterial == diffuse)
	{
		useEye = current;
		curEye = current;
	}

	if (useEye != curEye)
	{
		useEye = curEye;
		needRelight = true;
		isRawReflectInit = false;
	}


	bool needInit = false;
	bool isRefract = false; 
	rawUsingIndirect = isIndirect;
	rawUsingRGB = isRGB;
	rawUsingRefVertex = isRefracetVertex;
	rawUsingRefLight = isRefracet;
	rawRefractGlossy = isRefractGlossy;
	MatrixXf *tempMat = rawReflectTransportMat;
	bool tempGlossy = false;

	if (rawUsingRGB)
		RGBnum = col;
	else
		RGBnum = rad;

	int csl2 = cubeSideLength*cubeSideLength;

	if (!isRawReflectInit && isReflect)
	{
		vertexNum = meshes[0]->vertices.size();
		rawReflectTransportMat[0].resize(vertexNum, 6 * cubeSideLength*cubeSideLength);
		isRawReflectInit = true; 
		needInit = true;
		isRefract = false;
		tempMat = rawReflectTransportMat;
	}

	if (RGBnum == col && rawNeedInitRGB)
	{
		rawNeedInitRGB = false;
		for (int rgbInd = 1; rgbInd<RGBnum; rgbInd++)
			rawReflectTransportMat[rgbInd].resize(vertexNum, 6 * cubeSideLength*cubeSideLength);
	}

	if (!isRawRefractInit && isRefracet)
	{
		vertexNum = meshes[0]->vertices.size();

		if (isRefracetVertex)
		{
			needRelight = true;
			rawRefractTransportMat[0].resize(vertexNum, vertexNum);
			isRawRefractInit = true;
			getRefractionVertexAsLightMat(rawRefractTransportMat, rawRefractGlossy);
			return;
		}
		else
		{
			rawRefractTransportMat[0].resize(vertexNum, 6 * cubeSideLength*cubeSideLength);
			isRawRefractInit = true;
			needInit = true;
			isRefract = true;
			tempMat = rawRefractTransportMat;
		}
	}


	if (needInit)
	{

		needRelight = true;
		MatrixXf(*matrices)[3] = new MatrixXf[vertexNum][3];
		for (int i = 0; i < vertexNum; i++)
		for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
			matrices[i][rgbInd].resize(cubeSideLength, cubeSideLength);



		// Load the basis

		for (int f = 0; f < 6; f++)
		{
			getTranMat(matrices, f, rawUsingIndirect, isRefract, false);
			for (int y = 0; y < cubeSideLength; y++)
			for (int x = 0; x < cubeSideLength; x++)
			{
				int matCol = f*csl2 + y*cubeSideLength + x;
				// Extract pixels, populate transportMat
				for (int p = 0; p < vertexNum; p++)
				for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
					tempMat[rgbInd](p, matCol) = matrices[p][rgbInd](y, x);
			}
		}
		delete[] matrices;
	}
	/*
	if (isUnOclInit == false && isUnOclReflect)
	{
		isUnOclInit = true;
		needRelight = true;
		MatrixXf (*matrices)[3] = new MatrixXf[vertexNum][3];
		for (int i = 0; i < vertexNum; i++)
			matrices[i][0].resize(cubeSideLength, cubeSideLength);



		// Load the basis
		rawReflectTransportMatNoOcl.resize(vertexNum, 6 * cubeSideLength*cubeSideLength);
		for (int f = 0; f < 6; f++)
		{
			getTranMat(matrices, f, rawUsingIndirect, isRefract, isUnOclReflect);
			for (int y = 0; y < cubeSideLength; y++)
			for (int x = 0; x < cubeSideLength; x++)
			{
				int matCol = f*csl2 + y*cubeSideLength + x;
				// Extract pixels, populate transportMat
				for (int p = 0; p < vertexNum; p++)
				rawReflectTransportMatNoOcl(p, matCol) = matrices[p][0](y, x);
			}
		}
		delete[] matrices;


	}
	*/
}

void Transport::LoadWav(bool isIndirect, bool isRefraction, bool isRGB, bool isGlossy)
{
	if (isIndirect == wavUsingIndirect && isRGB == wavUsingRGB && wavUsingRefraction == isRefraction && isWavInit)
		return;


	wavUsingIndirect = isIndirect;
	wavUsingRGB = isRGB;
	wavUsingRefraction = isRefraction;
	if (wavUsingRGB)
		RGBnum = col;
	else
		RGBnum = rad;

	int csl2 = cubeSideLength*cubeSideLength;
	if (!isWavInit)
	{
		vertexNum = meshes[0]->vertices.size();
		colEnergy[0] = new float[6 * csl2];
		wavTransportMat[0].resize(vertexNum, 6 * csl2);
		resultColor = new Color[6 * csl2];
		isWavInit = true;
	}
	if (RGBnum == col && wavNeedInitRGB)
	{
		wavNeedInitRGB = false;
		for (int rgbInd = 1; rgbInd<RGBnum; rgbInd++)
		{
			colEnergy[rgbInd] = new float[6 * csl2];
			wavTransportMat[rgbInd].resize(vertexNum, 6 * csl2);
		}
	}


	// Alloc space for per-column energies.
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		colEnergy[rgbInd] = new float[6 * csl2];
	int numNonZero[3] = { 0, 0, 0 };

	// Define size of transport matrix
	// Load the basis
	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		wavTransportMat[rgbInd].resize(vertexNum, 6 * csl2);

	ifstream WaveCoeffTh_if(meshes[0]->dFName.getWaveCofThreDN(cubeSideLength, threshold, wavUsingIndirect, RGBnum, isGlossy));
	if (WaveCoeffTh_if.is_open())
	{
		int testVN;
		WaveCoeffTh_if >> testVN;
		if (testVN == vertexNum)
		{
			cout << "Loading Wavelet coffees..." << endl;
			for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
			{
				wavTransportMat[rgbInd].startFill();
				WaveCoeffTh_if >> numNonZero[rgbInd];
				for (int valInd = 0, matCol = 0; valInd<numNonZero[rgbInd]; matCol++)
				{
					int colPoint;
					WaveCoeffTh_if >> colPoint;
					for (; valInd<colPoint; valInd++)
					{
						float rgbVal;
						int p;
						WaveCoeffTh_if >> rgbVal;
						WaveCoeffTh_if >> p;
						wavTransportMat[rgbInd].fill(p, matCol) = rgbVal;
					}
				}
				wavTransportMat[rgbInd].endFill();
				for (int matCol = 0; matCol<6 * csl2; matCol++)
					WaveCoeffTh_if >> colEnergy[rgbInd][matCol];
			}
			WaveCoeffTh_if.close();
			for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
				std::cout << "Transport Matrix Nonzero Elements:" << numNonZero[rgbInd] << " (" << (int)(100 * numNonZero[rgbInd] / (vertexNum * 6 * csl2)) << "%)" << endl;
			return;
		}
		WaveCoeffTh_if.close();
	}

	bool needRecal = true;
	ifstream WaveCoeff_if(meshes[0]->dFName.getWaveCofDN(cubeSideLength, wavUsingIndirect, RGBnum, isGlossy));
	if (WaveCoeff_if.is_open())
	{
		int testVN;
		WaveCoeff_if >> testVN;
		if (testVN == vertexNum)
			needRecal = false;
	}

	if (needRecal)
	{
		WaveCoeff_if.close();
		compressTranMat(isGlossy);
		WaveCoeff_if.open(meshes[0]->dFName.getWaveCofDN(cubeSideLength, wavUsingIndirect, RGBnum, isGlossy));
		WaveCoeff_if >> vertexNum;
	}

	for (int f = 0; f < 6; f++)
	{
		std::cout << "Loading Wavelet coffees face " << f + 1 << "/6..." << endl;
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				int matCol = f*csl2 + y*cubeSideLength + x;
				for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
					colEnergy[rgbInd][matCol] = 0.f;
				// Extract pixels, populate transportMat
				// Sparsify by dropping elements with fabs < threshold;
				for (int p = 0; p < vertexNum; p++)
				{
					for (int rgbInd = 0; rgbInd < RGBnum; rgbInd++)
					{
						float rgbVal;
						WaveCoeff_if >> rgbVal;
						if (fabs(rgbVal) >= threshold)
						{
							wavTransportMat[rgbInd].fill(p, matCol) = rgbVal;
							colEnergy[rgbInd][matCol] += rgbVal*rgbVal;
							numNonZero[rgbInd]++;
						}
					}
				}
				for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
					colEnergy[rgbInd][matCol] = sqrt(colEnergy[rgbInd][matCol]);
			}
		}
	}


	std::cout << "Writing  file..." << endl;
	ofstream WaveCoeffTh_of(meshes[0]->dFName.getWaveCofThreDN(cubeSideLength, threshold, wavUsingIndirect, RGBnum, isGlossy));
	WaveCoeffTh_of << vertexNum << endl;

	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
	{
		const int* indices = wavTransportMat[rgbInd]._innerIndexPtr();
		const float* values = wavTransportMat[rgbInd]._valuePtr();

//		wavTransportMat[rgbInd].endFill();
		WaveCoeffTh_of << numNonZero[rgbInd] << endl;
		int outerSize = wavTransportMat[rgbInd].outerSize();
		int valInd = 0;
		for (int outInd = 1; outInd<outerSize && valInd<numNonZero[rgbInd]; outInd++)
		{
			int colPoint = wavTransportMat[rgbInd]._outerIndexPtr()[outInd];
			WaveCoeffTh_of << colPoint << " ";
			for (; valInd<colPoint; valInd++)
			{
				WaveCoeffTh_of << values[valInd] << " " << indices[valInd] << " ";
			}
			WaveCoeffTh_of << endl;
		}
		if (valInd<numNonZero[rgbInd])
		{
			WaveCoeffTh_of << numNonZero[rgbInd] << " ";
			for (; valInd<numNonZero[rgbInd]; valInd++)
			{
				WaveCoeffTh_of << values[valInd] << " " << indices[valInd] << " ";
			}
			WaveCoeffTh_of << endl;
		}
		for (int matCol = 0; matCol<6 * csl2; matCol++)
			WaveCoeffTh_of << colEnergy[rgbInd][matCol] << " ";
	}


	WaveCoeffTh_of.close();

	for (int rgbInd = 0; rgbInd<RGBnum; rgbInd++)
		std::cout << "Transport Matrix Nonzero Elements:" << numNonZero[rgbInd] << " (" << (int)(100 * numNonZero[rgbInd] / (vertexNum * 6 * csl2)) << "%)" << endl;

}

void Transport::designLightReflectAllvertexRadiance(vector<float> &cur, bool isIndirect, bool isRGB, int evnLightIndex, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy)
{
	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;

	bool needCal = true;


	if (oldSharpness == sharpness && oldReflectGlossy == isGlossy)
		needCal = false;

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy)
	{

		int testVN;
		FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");
		if (fin != NULL)
		{
			fread(&testVN, sizeof(int), 1, fin);
			if (testVN == vertexNum)
				needCal = false;
			fclose(fin);
		}

		/*
		ifstream DesignedLights_if(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy));
		if (DesignedLights_if.is_open())
		{
		DesignedLights_if >> testVN;
		if (testVN == vertexNum)
		needCal = false;

		DesignedLights_if.close();
		}
		*/
	}

	VectorXf* maxLightRadiance = wavTranEnvironment.getRawScaledLightVec();

	if (designReflectLightMid.size() == 0)
		designReflectLightMid.resize(meshes[0]->vertices.size(), fcsl2);

	meshes[0]->need_curvatures();
	meshes[0]->tune_curve_direction();



	//calculate all
	if (needCal)
	{
		Vector3f**** sampleDirs = wavTranEnvironment.getSampleDirs();
		int samplesPerPixel = wavTranEnvironment.getSamplePP();


		std::cout << "Writing Designed Lighting..." << endl;

		float *s = new float[fcsl2];
		/*
		ofstream DesignedLights_of(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy));
		DesignedLights_of << vertexNum << endl;
		*/
		FILE *fout = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "wb");
		fwrite(&vertexNum, sizeof(int), 1, fout);

		for (int verIndex = 0; verIndex < vertexNum; verIndex++)
		{
			vec curvDir = meshes[0]->pdir1[verIndex];
			vec normal = meshes[0]->normals[verIndex];
			vec crossResult = curvDir.cross(normal);
			float weight = exp(meshes[0]->curv1[verIndex] * enhScale);


			for (int f = 0; f < 6; f++)
			{
				for (int y = 0; y < cubeSideLength; y++)
				{
					for (int x = 0; x < cubeSideLength; x++)
					{
						int vecI = f*csl2 + y*cubeSideLength + x;
						float lightRadiance = 0.f;

						// Use each precomputed sample direction, index
						// into high-res cube map (super sampling)
						for (int sIndex = 0; sIndex < samplesPerPixel; sIndex++)
						{
							vec lightDir(sampleDirs[f][y][x][sIndex][0], sampleDirs[f][y][x][sIndex][1], sampleDirs[f][y][x][sIndex][2]);
							if (lightDir.dot(normal) < 0.f)
								continue;


							// Normalized for numerical cubature
							lightRadiance += radianceDistLightChange(normal, lightDir, lightDir, crossResult, curvDir, sharpness, isGlossy);
						}

						designReflectLightMid(verIndex, vecI) = lightRadiance;

						s[vecI] = lightRadiance;
						//						DesignedLights_of << setiosflags(ios::fixed) << lightRadiance << " ";


					}

				}
			}

			fwrite(s, sizeof(float), fcsl2, fout);
			//			DesignedLights_of << endl;
		}
		fclose(fout);
		delete[] s;
		//		DesignedLights_of.close();
	}



	if ((oldSharpness != sharpness || oldReflectGlossy != isGlossy) && !needCal)
	{
		std::cout << "Reading Designed Lighting..." << endl;
		int testVN;


		float *s = new float[fcsl2];
		FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");


		fread(&testVN, sizeof(int), 1, fin);
		for (int verIndex = 0; verIndex < vertexNum; verIndex++)
		{
			fread(s, sizeof(float), fcsl2, fin);
			for (int vecI = 0; vecI < fcsl2; vecI++)
				designReflectLightMid(verIndex, vecI) = s[vecI];
		}
		fclose(fin);
		delete[] s;

		/*
		ifstream DesignedLights_if(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy));

		DesignedLights_if >> testVN;
		for (int verIndex = 0; verIndex < vertexNum; verIndex++)
		for (int vecI = 0; vecI < fcsl2; vecI++)
		DesignedLights_if >> designReflectLightMid(verIndex, vecI);
		*/
	}


	VectorXf& lightMidRawVec = wavTranEnvironment.getRawMidLightVec();

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale)
	{

		//		LoadRaw(false, false, true, isIndirect, isRGB, isGlossy);

		for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
			lightMidRawVec[lightIndex] = 0.f;

		float tempEnhScale = enhScale / curv_scale;


#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			maxLightRadiance[0][vecI] = 0.f;
			maxLightRadiance[1][vecI] = 0.f;
		}

#pragma omp parallel for
		//for (int verIndex = evnLightIndex; verIndex < evnLightIndex+1; verIndex++)
		for (int verIndex = 0; verIndex < vertexNum; verIndex++)
		{

			float weight = exp(meshes[0]->curv1[verIndex] * tempEnhScale);

			for (int vecI = 0; vecI < fcsl2; vecI++)
			{
				float lightRadiance = designReflectLightMid(verIndex, vecI) * weight;

				lightMidRawVec[vecI] += lightRadiance;

				if (lightRadiance > maxLightRadiance[0][vecI])
					maxLightRadiance[0][vecI] = lightRadiance;
				if (lightRadiance < maxLightRadiance[1][vecI])
					maxLightRadiance[1][vecI] = lightRadiance;

			}
		}


		float gamma = 0.f;

#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			if (lightMidRawVec[vecI] > maxLightRadiance[0][vecI])
				lightMidRawVec[vecI] = maxLightRadiance[0][vecI];
			else if (lightMidRawVec[vecI] < maxLightRadiance[1][vecI])
				lightMidRawVec[vecI] = maxLightRadiance[1][vecI];

			if (gamma < abs(lightMidRawVec[vecI]))
				gamma = abs(lightMidRawVec[vecI]);
		}

		gamma = 1.0f / gamma;

#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			lightMidRawVec[vecI] *= gamma;
		}

	}

	VectorXf* lightScaledRawVec = wavTranEnvironment.getRawScaledLightVec();
	VectorXf* lightRawVec = wavTranEnvironment.getRawLightVec();

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale || oldAverage != average || oldScale != scale)
	{

		float averageLight = wavTranEnvironment.getRawLightAvg();

		float gamma = 0.f;

		average *= averageLight;
		scale = scale * averageLight;
#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			Color temp(lightRawVec[0][vecI] * fcsl2, lightRawVec[1][vecI] * fcsl2, lightRawVec[2][vecI] * fcsl2);

			temp = temp.convert(temp.RGB, temp.CIELAB);

			temp[0] = average + lightMidRawVec[vecI] * scale;

			temp = temp.convert(temp.CIELAB, temp.RGB);
			lightScaledRawVec[0][vecI] = temp[0] / fcsl2;
			lightScaledRawVec[1][vecI] = temp[1] / fcsl2;
			lightScaledRawVec[2][vecI] = temp[2] / fcsl2;
		}

		needRelight = true;

		oldSharpness = sharpness;
		oldReflectGlossy = isGlossy;
		oldEnhScale = enhScale;
		oldAverage = average;
		oldScale = scale;
	}







	/*

	for (int vertexIndex = 0; vertexIndex < vertexNum; vertexIndex++)
	{
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] += rawTransportMat[0](vertexIndex, lightIndex)*cur[vertexIndex];
	}
	}
	}
	/*
	float temp = 1.0 / vertexNum;
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] *= temp;

	lightScaledRawVec[rgbInd][lightIndex] = exp(lightScaledRawVec[rgbInd][lightIndex] * 5.0)*0.01;
	lightScaledRawVec[rgbInd][lightIndex] = 0.0;

	}
	}



	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	float a = rawTransportMat[0](evnLightIndex, lightIndex);
	if (a)
	a = a;
	lightScaledRawVec[rgbInd][lightIndex] = a;
	}
	}


	float gammaS = 1.0 / fcsl2;

	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] *= gammaS;
	}
	}



	needRelight = true;

	#pragma omp parallel for
	for (int i = 0; i < fcsl2; i++)
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	rawScaledLightVec[rgbInd][i] = lightVec[rgbInd][i] * 2.0;
	*/


	/*


	Eigen::SparseVector<float>* lightWavVec = wavTranEnvironment.getWavLightVec();
	Eigen::SparseVector<float>* lightScaledWavVec = wavTranEnvironment.getWavScaledLightVec();





	lightScaledWavVec[0].setZero();
	lightScaledWavVec[1].setZero();
	lightScaledWavVec[2].setZero();



	for (int rgbInd = 0; rgbInd<3; rgbInd++)
	{

	const int* indices = lightWavVec[rgbInd]._innerIndexPtr();
	const float* values = lightWavVec[rgbInd]._valuePtr();

	int numNonZero = lightWavVec[rgbInd].nonZeros();

	for (int valInd = 0; valInd<numNonZero; valInd++)
	{
	lightScaledWavVec[rgbInd].fill(indices[valInd]) = values[valInd] * 2.0;
	}
	}
	*/
}

void Transport::designLightReflectRidgeRadiance(vector<float> &cur, bool isIndirect, bool isRGB, int evnLightIndex, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy)
{
	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;


	bool needCal = true;



	meshes[0]->need_RidgeLines(true);
	meshes[0]->need_RidgeLines(false);


	/*
	vector<point> &ridgePoint = meshes[0]->ridgeLines.ridgePoint;
	vector<vec> &ridgePointNormal = meshes[0]->ridgeLines.ridgePointNormal;
	vector<float> &ridgePointPcurv = meshes[0]->ridgeLines.ridgePointPcurv;
	vector<vec> &ridgePointPDir = meshes[0]->ridgeLines.ridgePointPDir;

	vector<point> &valleyPoint = meshes[0]->ridgeLines.valleyPoint;
	vector<vec> &valleyPointNormal = meshes[0]->ridgeLines.valleyPointNormal;
	vector<float> &valleyPointPcurv = meshes[0]->ridgeLines.valleyPointPcurv;
	vector<vec> &valleyPointPDir = meshes[0]->ridgeLines.valleyPointPDir;


	int featureNum = valleyPoint.size() + ridgePoint.size();

	if (oldSharpness == sharpness && oldReflectGlossy == isGlossy)
		needCal = false;

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy)
	{

		int testVN;
		FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");
		if (fin != NULL)
		{
			fread(&testVN, sizeof(int), 1, fin);
			if (testVN == featureNum)
				needCal = false;
			fclose(fin);
		}

	}

	VectorXf* maxLightRadiance = wavTranEnvironment.getRawScaledLightVec();

	if (designReflectLightMid.size() == 0)
		designReflectLightMid.resize(featureNum, fcsl2);




	//calculate all
	if (needCal)
	{
		Vector3f**** sampleDirs = wavTranEnvironment.getSampleDirs();
		int samplesPerPixel = wavTranEnvironment.getSamplePP();


		std::cout << "Writing Designed Lighting..." << endl;


		for (int featureIndex = 0; featureIndex < ridgePoint.size(); featureIndex++)
		{
			vec curvDir = ridgePointPDir[featureIndex];
			vec normal = ridgePointNormal[featureIndex];
			vec crossResult = curvDir.cross(normal);
			float weight = exp(ridgePointPcurv[featureIndex] * enhScale);


			for (int f = 0; f < 6; f++)
			{
				for (int y = 0; y < cubeSideLength; y++)
				{
					for (int x = 0; x < cubeSideLength; x++)
					{
						int vecI = f*csl2 + y*cubeSideLength + x;
						float lightRadiance = 0.f;

						// Use each precomputed sample direction, index
						// into high-res cube map (super sampling)
						for (int sIndex = 0; sIndex < samplesPerPixel; sIndex++)
						{
							vec lightDir(sampleDirs[f][y][x][sIndex][0], sampleDirs[f][y][x][sIndex][1], sampleDirs[f][y][x][sIndex][2]);
							if (lightDir.dot(normal) < 0.f)
								continue;


							// Normalized for numerical cubature
							lightRadiance += radianceDistLightChange(normal, lightDir, lightDir, crossResult, curvDir, sharpness, isGlossy);
						}

						designReflectLightMid(featureIndex, vecI) = lightRadiance;

//						DesignedLights_of << setiosflags(ios::fixed) << lightRadiance << " ";
					}
				}
				//			DesignedLights_of << endl;
			}
		}

		int totalFeatureIndex = ridgePoint.size() - 1;
		for (int featureIndex = 0; featureIndex < valleyPoint.size(); featureIndex++)
		{
			totalFeatureIndex++;
			vec curvDir = valleyPointPDir[featureIndex];
			vec normal = valleyPointNormal[featureIndex];
			vec crossResult = curvDir.cross(normal);
			float weight = exp(valleyPointPcurv[featureIndex] * enhScale);


			for (int f = 0; f < 6; f++)
			{
				for (int y = 0; y < cubeSideLength; y++)
				{
					for (int x = 0; x < cubeSideLength; x++)
					{
						int vecI = f*csl2 + y*cubeSideLength + x;
						float lightRadiance = 0.f;

						// Use each precomputed sample direction, index
						// into high-res cube map (super sampling)
						for (int sIndex = 0; sIndex < samplesPerPixel; sIndex++)
						{
							vec lightDir(sampleDirs[f][y][x][sIndex][0], sampleDirs[f][y][x][sIndex][1], sampleDirs[f][y][x][sIndex][2]);
							if (lightDir.dot(normal) < 0.f)
								continue;


							// Normalized for numerical cubature
							lightRadiance += radianceDistLightChange(normal, lightDir, lightDir, crossResult, curvDir, sharpness, isGlossy);
						}

						designReflectLightMid(totalFeatureIndex, vecI) = lightRadiance;

						//						DesignedLights_of << setiosflags(ios::fixed) << lightRadiance << " ";
					}
				}
				//			DesignedLights_of << endl;
			}
		}


		FILE *fout = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "wb");
		fwrite(&featureNum, sizeof(int), 1, fout);
		float *s = new float[fcsl2];


		for (int featureIndex = 0; featureIndex < featureNum; featureIndex++)
		{
			for (int vecI = 0; vecI < fcsl2; vecI++)
				s[vecI] = designReflectLightMid(featureIndex, vecI);
			fwrite(s, sizeof(float), fcsl2, fout);
		}

		fclose(fout);
		delete[] s;
//		DesignedLights_of.close();
	}



	if ((oldSharpness != sharpness || oldReflectGlossy != isGlossy) && !needCal)
	{
		std::cout << "Reading Designed Lighting..." << endl;
		int testVN;


		float *s = new float[fcsl2];
		FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");


		fread(&testVN, sizeof(int), 1, fin);
		for (int featureIndex = 0; featureIndex < featureNum; featureIndex++)
		{
			fread(s, sizeof(float), fcsl2, fin);
			for (int vecI = 0; vecI < fcsl2; vecI++)
				designReflectLightMid(featureIndex, vecI) = s[vecI];
		}
		fclose(fin);
		delete[] s;

	}


	VectorXf& lightMidRawVec = wavTranEnvironment.getRawMidLightVec();

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale)
	{

//		LoadRaw(false, false, true, isIndirect, isRGB, isGlossy);

		for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
			lightMidRawVec[lightIndex] = 0.f;

		float tempEnhScale = enhScale / curv_scale;

		
#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			maxLightRadiance[0][vecI] = 0.f;
			maxLightRadiance[1][vecI] = 0.f;
		}

#pragma omp parallel for
		//for (int verIndex = evnLightIndex; verIndex < evnLightIndex+1; verIndex++)

		for (int featureIndex = 0; featureIndex < ridgePoint.size(); featureIndex++)
		{

			float weight = exp(ridgePointPcurv[featureIndex] * tempEnhScale);

			for (int vecI = 0; vecI < fcsl2; vecI++)
			{
				float lightRadiance = designReflectLightMid(featureIndex, vecI) * weight;

				lightMidRawVec[vecI] += lightRadiance;

				if (lightRadiance > maxLightRadiance[0][vecI])
					maxLightRadiance[0][vecI] = lightRadiance;
				if (lightRadiance < maxLightRadiance[1][vecI])
					maxLightRadiance[1][vecI] = lightRadiance;

			}
		}


		int totalFeatureIndex = ridgePoint.size() - 1;
		for (int featureIndex = 0; featureIndex < valleyPoint.size(); featureIndex++)
		{
			totalFeatureIndex++;
			float weight = exp(valleyPointPcurv[featureIndex] * tempEnhScale);

			for (int vecI = 0; vecI < fcsl2; vecI++)
			{
				float lightRadiance = designReflectLightMid(totalFeatureIndex, vecI) * weight;

				lightMidRawVec[vecI] += lightRadiance;

				if (lightRadiance > maxLightRadiance[0][vecI])
					maxLightRadiance[0][vecI] = lightRadiance;
				if (lightRadiance < maxLightRadiance[1][vecI])
					maxLightRadiance[1][vecI] = lightRadiance;

			}
		}


		float gamma = 0.f;

#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			if (lightMidRawVec[vecI] > maxLightRadiance[0][vecI])
				lightMidRawVec[vecI] = maxLightRadiance[0][vecI];
			else if (lightMidRawVec[vecI] < maxLightRadiance[1][vecI])
				lightMidRawVec[vecI] = maxLightRadiance[1][vecI];

			if (gamma < abs(lightMidRawVec[vecI]))
				gamma = abs(lightMidRawVec[vecI]);
		}

		gamma = 1.0f / gamma;

#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			lightMidRawVec[vecI] *= gamma;
		}

	}

	VectorXf* lightScaledRawVec = wavTranEnvironment.getRawScaledLightVec();
	VectorXf* lightRawVec = wavTranEnvironment.getRawLightVec();

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale || oldAverage != average || oldScale != scale)
	{

		float averageLight = wavTranEnvironment.getRawLightAvg();


		average *= averageLight;
		scale = scale * averageLight;
#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			float gamma = (average + lightMidRawVec[vecI] * scale) / averageLight;

			for (int cInd = 0; cInd < 3; cInd++)
				lightScaledRawVec[cInd][vecI] = lightRawVec[cInd][vecI] * gamma;
		}

		needRelight = true;

		oldSharpness = sharpness;
		oldReflectGlossy = isGlossy;
		oldEnhScale = enhScale;
		oldAverage = average;
		oldScale = scale;
	}
	

	*/




	/*

	for (int vertexIndex = 0; vertexIndex < vertexNum; vertexIndex++)
	{
		for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
		{
			for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			{
				lightScaledRawVec[rgbInd][lightIndex] += rawTransportMat[0](vertexIndex, lightIndex)*cur[vertexIndex];
			}
		}
	}
	/*
	float temp = 1.0 / vertexNum;
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			lightScaledRawVec[rgbInd][lightIndex] *= temp;

			lightScaledRawVec[rgbInd][lightIndex] = exp(lightScaledRawVec[rgbInd][lightIndex] * 5.0)*0.01;
			lightScaledRawVec[rgbInd][lightIndex] = 0.0;

		}
	}


	
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			float a = rawTransportMat[0](evnLightIndex, lightIndex);
			if (a)
				a = a;
			lightScaledRawVec[rgbInd][lightIndex] = a;
		}
	}


		float gammaS = 1.0 / fcsl2;

		for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			lightScaledRawVec[rgbInd][lightIndex] *= gammaS;
		}
	}



	needRelight = true;

#pragma omp parallel for
	for (int i = 0; i < fcsl2; i++)
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		rawScaledLightVec[rgbInd][i] = lightVec[rgbInd][i] * 2.0;
	*/

	
	/*


	Eigen::SparseVector<float>* lightWavVec = wavTranEnvironment.getWavLightVec();
	Eigen::SparseVector<float>* lightScaledWavVec = wavTranEnvironment.getWavScaledLightVec();





	lightScaledWavVec[0].setZero();
	lightScaledWavVec[1].setZero();
	lightScaledWavVec[2].setZero();

		

	for (int rgbInd = 0; rgbInd<3; rgbInd++)
	{

		const int* indices = lightWavVec[rgbInd]._innerIndexPtr();
		const float* values = lightWavVec[rgbInd]._valuePtr();

		int numNonZero = lightWavVec[rgbInd].nonZeros();

		for (int valInd = 0; valInd<numNonZero; valInd++)
		{
			lightScaledWavVec[rgbInd].fill(indices[valInd]) = values[valInd] * 2.0;
		}
	}
	*/
}

void Transport::designLightReflectRidgeVertexColor(vec3 eyeDirection, int vertexIndex, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy)
{
	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;

	bool needCal = true;

	meshes[0]->need_RidgeLines(true);
	meshes[0]->need_RidgeLines(false);
	if (!meshes[0]->getPlaneNormal())
		return;
	int featureNum = meshes[0]->halfEdge.leftPlaneNormal.size();


	//if (oldSharpness == sharpness && oldReflectGlossy == isGlossy)
	//	needCal = false;
	/*
	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy)
	{

		int testVN;
		FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");
		if (fin != NULL)
		{
			fread(&testVN, sizeof(int), 1, fin);
			if (testVN == featureNum)
				needCal = false;
			fclose(fin);
		}

	}
	*/

	if (designReflectLightMid.size() == 0)
		designReflectLightMid.resize(fcsl2, fcsl2);


	//calculate all
	if (needCal)
	{
		vector<float> lightWeight;
		lightWeight.resize(fcsl2);
		Vector3f**** sampleDirs = wavTranEnvironment.getSampleDirs();
		int samplesPerPixel = wavTranEnvironment.getSamplePP();

		vector<vec3> &leftPlaneNormal = meshes[0]->halfEdge.leftPlaneNormal;
		vector<vec3> &rightPlaneNormal = meshes[0]->halfEdge.rightPlaneNormal;

		std::cout << "Writing Designed Lighting..." << endl;

		/*
		ofstream DesignedLights_of(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy));
		DesignedLights_of << vertexNum << endl;
		*/


		for (int posVecI = 0; posVecI < fcsl2; posVecI++)
		for (int negVecI = posVecI; negVecI < fcsl2; negVecI++)
		{
			designReflectLightMid(posVecI, negVecI) = 0.f;
			designReflectLightMid(negVecI, posVecI) = 0.f;
		}


		for (int featureIndex = vertexIndex; featureIndex < vertexIndex+1; featureIndex++)
		{
			vec3 leftNormal = leftPlaneNormal[featureIndex];
			vec3 rightNormal = rightPlaneNormal[featureIndex];

			int vecI = 0;
			for (int f = 0; f < 6; f++)
			for (int y = 0; y < cubeSideLength; y++)
			for (int x = 0; x < cubeSideLength; x++)
			{
				int vecI = f*csl2 + y*cubeSideLength + x;
				float lightRadiance = 0.f;
				// Use each precomputed sample direction, index
				// into high-res cube map (super sampling)
				for (int sIndex = 0; sIndex < samplesPerPixel; sIndex++)
				{
					vec lightDir(sampleDirs[f][y][x][sIndex][0], sampleDirs[f][y][x][sIndex][1], sampleDirs[f][y][x][sIndex][2]);
					
					float leftWeight = lightDir.dot(leftNormal);
					float rightWeight = lightDir.dot(rightNormal);

					if (leftWeight>0.1f && rightWeight <= 0.1f)
					{
						leftWeight = pow(leftWeight, sharpness * 2.f);
						if (leftWeight > 0.1)
							lightRadiance += leftWeight;
					}

					if (rightWeight>0.1f && leftWeight <= 0.1f)
					{
						rightWeight = pow(rightWeight, sharpness * 2.f);
						if (rightWeight > 0.1)
							lightRadiance -= rightWeight;
					}
				}
				lightWeight[vecI++] = lightRadiance;
			}

			for (int posVecI = 0; posVecI < fcsl2; posVecI++)
			{
				if (lightWeight[posVecI] <= 0.f)
					continue;
				for (int negVecI = 0; negVecI < fcsl2; negVecI++)
				{
					if (lightWeight[negVecI] >= 0.f)
						continue;

					float difference = lightWeight[posVecI] - lightWeight[negVecI];
					if (designReflectLightMid(posVecI, negVecI) < difference)
					{
						//difference should be positive;
						designReflectLightMid(posVecI, negVecI) = difference;
						designReflectLightMid(negVecI, posVecI) = difference;
					}
				}
			}
		}

		

		vec2* coord2D = new vec2[fcsl2];

		vec3 tempAixeZ(eyeDirection[0], 0.f, eyeDirection[2]);
		float xLength = len(tempAixeZ);
		float yLength = eyeDirection[1];
		vec3 tempAixeY = tempAixeZ * (-yLength / xLength);
		tempAixeY[1] = xLength;
		tempAixeZ[1] = eyeDirection[1];
		vec3 tempAixeX = tempAixeY.cross(tempAixeZ);

		int cenSampleInd = samplesPerPixel*0.5;

		int vecI = 0;
		for (int f = 0; f < 6; f++)
		for (int y = 0; y < cubeSideLength; y++)
		for (int x = 0; x < cubeSideLength; x++)
		{
			vec lightDir(sampleDirs[f][y][x][cenSampleInd][0], sampleDirs[f][y][x][cenSampleInd][1], sampleDirs[f][y][x][cenSampleInd][2]);
			
			float xCoord = lightDir.dot(tempAixeX);
			float yCoord = lightDir.dot(tempAixeY);
			float zCoord = lightDir.dot(tempAixeZ);

			coord2D[vecI][0] = xCoord >= 0 ? PI - acos(zCoord) : acos(zCoord) + PI;
			coord2D[vecI][1] = yCoord >= 0 ? PI - acos(zCoord) : acos(zCoord) + PI;

			coord2D[vecI][0] *= len(vec2(xCoord, zCoord));
			coord2D[vecI][1] *= len(vec2(yCoord, zCoord));

			vecI++;
		}


		float* sumDifA = new float[fcsl2];
		float* sumDifB = new float[fcsl2];


		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			sumDifA[vecI] = 0.f;
			sumDifB[vecI] = 0.f;
		}


		for (int posVecI = 0; posVecI < fcsl2; posVecI++)
		{
			vec2& coord2DPos = coord2D[posVecI];

			for (int negVecI = posVecI+1; negVecI < fcsl2; negVecI++)
			{
				if (designReflectLightMid(posVecI, negVecI) <= 0.f)
					continue;

				vec2& coord2DNeg = coord2D[negVecI];

				vec2 vecCoord2D = coord2DPos - coord2DNeg;
				normalize(vecCoord2D);

				float diffA = designReflectLightMid(posVecI, negVecI) * vecCoord2D[0];
				float diffB = designReflectLightMid(posVecI, negVecI) * vecCoord2D[1];

				sumDifA[posVecI] += diffA;
				sumDifB[posVecI] += diffB;
				sumDifA[negVecI] += -diffA;
				sumDifB[negVecI] += -diffB;
			}
		}


		VectorXf* lightLabRawVec = wavTranEnvironment.getLabLightVecColor();

		calc_d(sumDifA, fcsl2, lightLabRawVec[1]);
		calc_d(sumDifB, fcsl2, lightLabRawVec[2]);


		float gamma = 98.f / 143.f;
#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
			lightLabRawVec[0][vecI] = sqrt(lightLabRawVec[1][vecI] * lightLabRawVec[1][vecI] + lightLabRawVec[2][vecI] * lightLabRawVec[2][vecI])*gamma;


		delete[] sumDifA;
		delete[] sumDifB;


		/*

		FILE *fout = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "wb");
		fwrite(&featureNum, sizeof(int), 1, fout);
		float *s = new float[fcsl2];


		for (int vecJ = 0; vecJ < fcsl2; vecJ++)
		{
			for (int vecI = 0; vecI < fcsl2; vecI++)
				s[vecI] = designReflectLightMid(vecJ, vecI);
			fwrite(s, sizeof(float), fcsl2, fout);
		}

		fclose(fout);
		delete[] s;
		//		DesignedLights_of.close();

		*/
	}


	/*
	if ((oldSharpness != sharpness || oldReflectGlossy != isGlossy) && !needCal)
	{
		std::cout << "Reading Designed Lighting..." << endl;
		int testVN;


		float *s = new float[fcsl2];
		FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");


		fread(&testVN, sizeof(int), 1, fin);
		for (int featureIndex = 0; featureIndex < featureNum; featureIndex++)
		{
			fread(s, sizeof(float), fcsl2, fin);
			for (int vecI = 0; vecI < fcsl2; vecI++)
				designReflectLightMid(featureIndex, vecI) = s[vecI];
		}
		fclose(fin);
		delete[] s;

	}
	*/






	

	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale)
	{

	}

	VectorXf* lightScaledRawVec = wavTranEnvironment.getRawScaledLightVec();
	VectorXf* lightLabRawVec = wavTranEnvironment.getLabLightVecColor();

	if (true || oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale || oldAverage != average || oldScale != scale)
	{




		needRelight = true;

		oldSharpness = sharpness;
		oldReflectGlossy = isGlossy;
		oldEnhScale = enhScale;
		oldAverage = average;
		oldScale = scale;
	}







	/*

	for (int vertexIndex = 0; vertexIndex < vertexNum; vertexIndex++)
	{
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] += rawTransportMat[0](vertexIndex, lightIndex)*cur[vertexIndex];
	}
	}
	}
	/*
	float temp = 1.0 / vertexNum;
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] *= temp;

	lightScaledRawVec[rgbInd][lightIndex] = exp(lightScaledRawVec[rgbInd][lightIndex] * 5.0)*0.01;
	lightScaledRawVec[rgbInd][lightIndex] = 0.0;

	}
	}



	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	float a = rawTransportMat[0](evnLightIndex, lightIndex);
	if (a)
	a = a;
	lightScaledRawVec[rgbInd][lightIndex] = a;
	}
	}


	float gammaS = 1.0 / fcsl2;

	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] *= gammaS;
	}
	}



	needRelight = true;

	#pragma omp parallel for
	for (int i = 0; i < fcsl2; i++)
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	rawScaledLightVec[rgbInd][i] = lightVec[rgbInd][i] * 2.0;
	*/


	/*


	Eigen::SparseVector<float>* lightWavVec = wavTranEnvironment.getWavLightVec();
	Eigen::SparseVector<float>* lightScaledWavVec = wavTranEnvironment.getWavScaledLightVec();





	lightScaledWavVec[0].setZero();
	lightScaledWavVec[1].setZero();
	lightScaledWavVec[2].setZero();



	for (int rgbInd = 0; rgbInd<3; rgbInd++)
	{

	const int* indices = lightWavVec[rgbInd]._innerIndexPtr();
	const float* values = lightWavVec[rgbInd]._valuePtr();

	int numNonZero = lightWavVec[rgbInd].nonZeros();

	for (int valInd = 0; valInd<numNonZero; valInd++)
	{
	lightScaledWavVec[rgbInd].fill(indices[valInd]) = values[valInd] * 2.0;
	}
	}
	*/
}

void Transport::writeTransferMatrixForLRNormals()
{


	float(*materialShader)(vec3 normal, vec3 lightDir, vec3 viewDir);

	switch (curMaterial)
	{
	case diffuse:
		materialShader = diffuseShader;
		break;
	case ward:
		materialShader = wardShader;
		break;

	default:
		materialShader = diffuseShader;
		break;
	}



	MatrixXf color;
	
	meshes[0]->need_RidgeLines(true);
	meshes[0]->need_RidgeLines(false);
	if (!meshes[0]->getPlaneNormal())
		return;

	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;
	int featureNum = meshes[0]->halfEdge.leftPlaneNormal.size();
	vector<vec3> &leftPlaneNormal = meshes[0]->halfEdge.leftPlaneNormal;
	vector<vec3> &rightPlaneNormal = meshes[0]->halfEdge.rightPlaneNormal;
	Vector3f**** sampleDirs = wavTranEnvironment.getSampleDirs();

	vec3 eye = meshes[0]->eyePosition;
	vector<vec3> &position = meshes[0]->halfEdge.samplePosition;

	int samplesPerPixel = wavTranEnvironment.getSamplePP();
	float sPPrevers = 1.f / samplesPerPixel;
	
	color.resize(featureNum*2, fcsl2);
	
	int lightNum = fcsl2 * 2;

	for (int featureIndex = 0; featureIndex < featureNum; featureIndex++)
	{
		vec3 leftNormal = leftPlaneNormal[featureIndex];
		vec3 rightNormal = rightPlaneNormal[featureIndex];

		int vecI = 0;
		for (int f = 0; f < 6; f++)
		for (int y = 0; y < cubeSideLength; y++)
		for (int x = 0; x < cubeSideLength; x++)
		{
			float lightRadianceLeft = 0.f;
			float lightRadianceRight = 0.f;
			// Use each precomputed sample direction, index
			// into high-res cube map (super sampling)
			for (int sIndex = 0; sIndex < samplesPerPixel; sIndex++)
			{
				vec lightDir(sampleDirs[f][y][x][sIndex][0], sampleDirs[f][y][x][sIndex][1], sampleDirs[f][y][x][sIndex][2]);

				point& pos = position[featureIndex];
				vec3 viewDir = eye - pos;
				normalize(viewDir);

				float leftWeight = materialShader(leftNormal, lightDir, viewDir);
					
				float rightWeight = materialShader(rightNormal, lightDir, viewDir);

				if (leftWeight>=0.f)
					lightRadianceLeft += leftWeight;
				if (rightWeight >= 0.f)
					lightRadianceRight += rightWeight;
			}
			lightRadianceLeft *= sPPrevers;
			lightRadianceRight *= sPPrevers;
			
			color(featureIndex, vecI) = lightRadianceLeft;
			color(featureIndex+featureNum, vecI) = lightRadianceRight;


			vecI++;
		}

	}
	
	ofstream VertOcl_col(meshes[0]->dFName.getFeatureToLightDN(cubeSideLength, featureNum));

	for (int featureIndex = 0; featureIndex < featureNum + featureNum; featureIndex++)
	{
		for (int i = 0; i < fcsl2; i++)
		{
			VertOcl_col << color(featureIndex, i) << " ";
		}
		VertOcl_col << endl;
	}

	VertOcl_col.close();
}

void Transport::readingExistingLight()
{


	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;


	VectorXf* lightLabRawVec = wavTranEnvironment.getLabLightVecColor();


	FILE *fin;
	fopen_s(&fin, meshes[0]->dFName.getOptLightDN(cubeSideLength, curMaterial, curEye, wavTranEnvironment.getFileName()), "rb");
	if (fin)
	{
		cout << "Reading optimetic light.....";

		fread(&designReflectLightTransA, sizeof(float), 1, fin);

		float *tempLight = new float[fcsl2];

		for (int labIndex = 0; labIndex < 3; labIndex++)
		{
			fread(&tempLight[0], sizeof(float), fcsl2, fin);
			for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
				lightLabRawVec[labIndex][lightIndex] = tempLight[lightIndex];
		}
		cout << "Done" << endl;
		fclose(fin);

		delete[] tempLight;
	}
	else
	{
		char str[180];

		strcpy_s(str, 180, "..\\..\\matlabSoveNonlinearEquations\\");
		strcat_s(str, 180, wavTranEnvironment.getFileName());
		strcat_s(str, 180, "_bestLight.txt");

		ifstream VertOcl_col(str);
		VertOcl_col >> designReflectLightTransA;


		for (int labIndex = 0; labIndex < 3; labIndex++)
		for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
			VertOcl_col >> lightLabRawVec[labIndex][lightIndex];

		VertOcl_col.close();

		FILE *fout;
		fopen_s(&fout, meshes[0]->dFName.getOptLightDN(cubeSideLength, curMaterial, curEye, wavTranEnvironment.getFileName()), "wb");


		fwrite(&designReflectLightTransA, sizeof(float), 1, fout);

		float *tempLight = new float[fcsl2];
		for (int labIndex = 0; labIndex < 3; labIndex++)
		{
			for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
				tempLight[lightIndex] = lightLabRawVec[labIndex][lightIndex];
			fwrite(&tempLight[0], sizeof(float), fcsl2, fout);
		}

		delete[] tempLight;
		fclose(fout);

	}
	

	
	VectorXf* rawScaledLightVec = wavTranEnvironment.getRawScaledLightVec();

	needRelight = true;

	float gammaS = 1.0 / fcsl2;

	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
		rawScaledLightVec[0][lightIndex] = lightLabRawVec[0][lightIndex] * gammaS;
		rawScaledLightVec[1][lightIndex] = lightLabRawVec[1][lightIndex] * gammaS;
		rawScaledLightVec[2][lightIndex] = lightLabRawVec[2][lightIndex] * gammaS;
	}
	


	/*

	for (int vertexIndex = 0; vertexIndex < vertexNum; vertexIndex++)
	{
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] += rawTransportMat[0](vertexIndex, lightIndex)*cur[vertexIndex];
	}
	}
	}
	/*
	float temp = 1.0 / vertexNum;
	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] *= temp;

	lightScaledRawVec[rgbInd][lightIndex] = exp(lightScaledRawVec[rgbInd][lightIndex] * 5.0)*0.01;
	lightScaledRawVec[rgbInd][lightIndex] = 0.0;

	}
	}



	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	float a = rawTransportMat[0](evnLightIndex, lightIndex);
	if (a)
	a = a;
	lightScaledRawVec[rgbInd][lightIndex] = a;
	}
	}


	float gammaS = 1.0 / fcsl2;

	for (int lightIndex = 0; lightIndex < fcsl2; lightIndex++)
	{
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	{
	lightScaledRawVec[rgbInd][lightIndex] *= gammaS;
	}
	}



	needRelight = true;

	#pragma omp parallel for
	for (int i = 0; i < fcsl2; i++)
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	rawScaledLightVec[rgbInd][i] = lightVec[rgbInd][i] * 2.0;
	*/


	/*


	Eigen::SparseVector<float>* lightWavVec = wavTranEnvironment.getWavLightVec();
	Eigen::SparseVector<float>* lightScaledWavVec = wavTranEnvironment.getWavScaledLightVec();





	lightScaledWavVec[0].setZero();
	lightScaledWavVec[1].setZero();
	lightScaledWavVec[2].setZero();



	for (int rgbInd = 0; rgbInd<3; rgbInd++)
	{

	const int* indices = lightWavVec[rgbInd]._innerIndexPtr();
	const float* values = lightWavVec[rgbInd]._valuePtr();

	int numNonZero = lightWavVec[rgbInd].nonZeros();

	for (int valInd = 0; valInd<numNonZero; valInd++)
	{
	lightScaledWavVec[rgbInd].fill(indices[valInd]) = values[valInd] * 2.0;
	}
	}
	*/
}

void Transport::designLightReflectRidgeColor(vec3 eyeDirection, float sharpness, float enhScale, float curv_scale, float average, float scale, bool isGlossy)
{
	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;

	bool needCal = true;

	meshes[0]->need_RidgeLines(true);
	meshes[0]->need_RidgeLines(false);
	if (!meshes[0]->getPlaneNormal())
		return;
	int featureNum = meshes[0]->halfEdge.leftPlaneNormal.size();


	//if (oldSharpness == sharpness && oldReflectGlossy == isGlossy)
	//	needCal = false;
	/*
	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy)
	{

	int testVN;
	FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");
	if (fin != NULL)
	{
	fread(&testVN, sizeof(int), 1, fin);
	if (testVN == featureNum)
	needCal = false;
	fclose(fin);
	}

	}
	*/

	if (designReflectLightMid.size() == 0)
		designReflectLightMid.resize(fcsl2, fcsl2);


	//calculate all
	if (needCal)
	{
		vector<float> lightWeight;
		lightWeight.resize(fcsl2);
		Vector3f**** sampleDirs = wavTranEnvironment.getSampleDirs();
		int samplesPerPixel = wavTranEnvironment.getSamplePP();

		vector<vec3> &leftPlaneNormal = meshes[0]->halfEdge.leftPlaneNormal;
		vector<vec3> &rightPlaneNormal = meshes[0]->halfEdge.rightPlaneNormal;

		std::cout << "Writing Designed Lighting..." << endl;

		/*
		ofstream DesignedLights_of(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy));
		DesignedLights_of << vertexNum << endl;
		*/

		VectorXf* lightLabRawVec = wavTranEnvironment.getLabLightVecColor();

		for (int posVecI = 0; posVecI < fcsl2; posVecI++)
		{
			lightLabRawVec[0][posVecI] = 0.f;
			for (int negVecI = posVecI; negVecI < fcsl2; negVecI++)
			{
				designReflectLightMid(posVecI, negVecI) = 0.f;
				designReflectLightMid(negVecI, posVecI) = 0.f;
			}

		}
		float maxDif = 0.f;

		for (int featureIndex = 0; featureIndex < featureNum; featureIndex++)
		{
			vec3 leftNormal = leftPlaneNormal[featureIndex];
			vec3 rightNormal = rightPlaneNormal[featureIndex];

			int vecI = 0;
			for (int f = 0; f < 6; f++)
			for (int y = 0; y < cubeSideLength; y++)
			for (int x = 0; x < cubeSideLength; x++)
			{
				float lightRadiance = 0.f;
				// Use each precomputed sample direction, index
				// into high-res cube map (super sampling)
				for (int sIndex = 0; sIndex < samplesPerPixel; sIndex++)
				{
					vec lightDir(sampleDirs[f][y][x][sIndex][0], sampleDirs[f][y][x][sIndex][1], sampleDirs[f][y][x][sIndex][2]);

					float leftWeight = lightDir.dot(leftNormal);
					float rightWeight = lightDir.dot(rightNormal);

					if (leftWeight>0.1f && rightWeight <= 0.1f)
					{
						leftWeight = pow(leftWeight, sharpness * 2.f);
						if (leftWeight > 0.1)
							lightRadiance += leftWeight;
					}

					if (rightWeight>0.1f && leftWeight <= 0.1f)
					{
						rightWeight = pow(rightWeight, sharpness * 2.f);
						if (rightWeight > 0.1)
							lightRadiance -= rightWeight;
					}
				}
				lightWeight[vecI++] = lightRadiance;
			}

			for (int posVecI = 0; posVecI < fcsl2; posVecI++)
			{
				if (lightWeight[posVecI] <= 0.f)
					continue;
				for (int negVecI = 0; negVecI < fcsl2; negVecI++)
				{
					if (lightWeight[negVecI] >= 0.f)
						continue;

					float difference = lightWeight[posVecI] - lightWeight[negVecI];
					if (designReflectLightMid(posVecI, negVecI) < difference)
					{
						//difference should be positive;
						designReflectLightMid(posVecI, negVecI) = difference;
						designReflectLightMid(negVecI, posVecI) = difference;

						if (lightLabRawVec[0][posVecI] < difference)
						{
							lightLabRawVec[0][posVecI] = difference;
							if (maxDif < difference)
								maxDif = difference;
						}
						if (lightLabRawVec[0][negVecI] < difference)
							lightLabRawVec[0][negVecI] = difference;




					}
				}
			}
		}



		vec2* coord2D = new vec2[fcsl2];

		vec3 tempAixeZ(eyeDirection[0], 0.f, eyeDirection[2]);
		float xLength = len(tempAixeZ);
		float yLength = eyeDirection[1];
		vec3 tempAixeY = tempAixeZ * (-yLength / xLength);
		tempAixeY[1] = xLength;
		tempAixeZ[1] = eyeDirection[1];
		vec3 tempAixeX = tempAixeY.cross(tempAixeZ);

		int cenSampleInd = samplesPerPixel*0.5;

		int vecI = 0;
		for (int f = 0; f < 6; f++)
		for (int y = 0; y < cubeSideLength; y++)
		for (int x = 0; x < cubeSideLength; x++)
		{
			vec lightDir(sampleDirs[f][y][x][cenSampleInd][0], sampleDirs[f][y][x][cenSampleInd][1], sampleDirs[f][y][x][cenSampleInd][2]);

			float xCoord = lightDir.dot(tempAixeX);
			float yCoord = lightDir.dot(tempAixeY);
			float zCoord = lightDir.dot(tempAixeZ);

			coord2D[vecI][0] = xCoord >= 0 ? PI - acos(zCoord) : acos(zCoord) + PI;
			coord2D[vecI][1] = yCoord >= 0 ? PI - acos(zCoord) : acos(zCoord) + PI;

			coord2D[vecI][0] *= len(vec2(xCoord, zCoord));
			coord2D[vecI][1] *= len(vec2(yCoord, zCoord));

			vecI++;
		}


		float* sumDifA = new float[fcsl2];
		float* sumDifB = new float[fcsl2];


		for (int vecI = 0; vecI < fcsl2; vecI++)
		{
			sumDifA[vecI] = 0.f;
			sumDifB[vecI] = 0.f;
		}


		for (int posVecI = 0; posVecI < fcsl2; posVecI++)
		{
			vec2& coord2DPos = coord2D[posVecI];

			for (int negVecI = posVecI + 1; negVecI < fcsl2; negVecI++)
			{
				if (designReflectLightMid(posVecI, negVecI) <= 0.f)
					continue;

				vec2& coord2DNeg = coord2D[negVecI];

				vec2 vecCoord2D = coord2DPos - coord2DNeg;
				normalize(vecCoord2D);

				float diffA = designReflectLightMid(posVecI, negVecI) * vecCoord2D[0];
				float diffB = designReflectLightMid(posVecI, negVecI) * vecCoord2D[1];

				diffA = designReflectLightMid(posVecI, negVecI);

				sumDifA[posVecI] += diffA;
				sumDifB[posVecI] += diffB;
				sumDifA[negVecI] += -diffA;
				sumDifB[negVecI] += -diffB;
			}
		}



		calc_d(sumDifA, fcsl2, lightLabRawVec[1]);
		calc_d(sumDifB, fcsl2, lightLabRawVec[2]);


//		float gamma = 98.f / 143.f;
		float gamma = 98.f / maxDif;
#pragma omp parallel for
		for (int vecI = 0; vecI < fcsl2; vecI++)
			lightLabRawVec[0][vecI] *= gamma;// sqrt(lightLabRawVec[1][vecI] * lightLabRawVec[1][vecI] + lightLabRawVec[2][vecI] * lightLabRawVec[2][vecI])*gamma;

		for (int vecI = 0; vecI < fcsl2; vecI++)
			lightLabRawVec[2][vecI] = 0.f;
		delete[] sumDifA;
		delete[] sumDifB;


		/*

		FILE *fout = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "wb");
		fwrite(&featureNum, sizeof(int), 1, fout);
		float *s = new float[fcsl2];


		for (int vecJ = 0; vecJ < fcsl2; vecJ++)
		{
		for (int vecI = 0; vecI < fcsl2; vecI++)
		s[vecI] = designReflectLightMid(vecJ, vecI);
		fwrite(s, sizeof(float), fcsl2, fout);
		}

		fclose(fout);
		delete[] s;
		//		DesignedLights_of.close();

		*/
	}


	/*
	if ((oldSharpness != sharpness || oldReflectGlossy != isGlossy) && !needCal)
	{
	std::cout << "Reading Designed Lighting..." << endl;
	int testVN;


	float *s = new float[fcsl2];
	FILE *fin = fopen(meshes[0]->dFName.getDesignedLightDN(cubeSideLength, sharpness, isGlossy), "rb");


	fread(&testVN, sizeof(int), 1, fin);
	for (int featureIndex = 0; featureIndex < featureNum; featureIndex++)
	{
	fread(s, sizeof(float), fcsl2, fin);
	for (int vecI = 0; vecI < fcsl2; vecI++)
	designReflectLightMid(featureIndex, vecI) = s[vecI];
	}
	fclose(fin);
	delete[] s;

	}
	*/








	if (oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale)
	{

	}

	VectorXf* lightScaledRawVec = wavTranEnvironment.getRawScaledLightVec();
	VectorXf* lightLabRawVec = wavTranEnvironment.getLabLightVecColor();

	if (true || oldSharpness != sharpness || oldReflectGlossy != isGlossy || oldEnhScale != enhScale || oldAverage != average || oldScale != scale)
	{




		needRelight = true;

		oldSharpness = sharpness;
		oldReflectGlossy = isGlossy;
		oldEnhScale = enhScale;
		oldAverage = average;
		oldScale = scale;
	}

}

void Transport::designLightRefract(vector<float> &cur, bool isIndirect, bool isRGB, int evnLightIndex, float sharpness, float enhScale, float curv_scale, float scale, bool isGlossy)// , vec3 eye)
{
	
	/*
	point &light = meshes[0]->vertices[evnLightIndex];

	VectorXf* lightScaledRawVec = wavTranEnvironment.getRawScaledLightVecRefract(vertexNum);


	RefractionParemeter rp(eta, sigma_a, sigma_s);
	float a = 0.5 / meshes[0]->bsphere.r;

	float average111 = 0.f;

	float teScale = scale;

	for (int vInd = 0; vInd < vertexNum; vInd++)
	{
		if (vInd == evnLightIndex)
		{

			lightScaledRawVec[0][vInd] = 0.f;
			lightScaledRawVec[1][vInd] = 0.f;
			lightScaledRawVec[2][vInd] = 0.f;
			continue;
		}
		float length = len(meshes[0]->vertices[vInd] - light) *a *teScale;

		lightScaledRawVec[0][vInd] = radianceTranslucent(length, rp);

		average111 += lightScaledRawVec[0][vInd];
	}

	average111 = 0.5 / average111;

	for (int vInd = 0; vInd < vertexNum; vInd++)
	{
		lightScaledRawVec[0][vInd] *= average111;
		lightScaledRawVec[1][vInd] = lightScaledRawVec[0][vInd];
		lightScaledRawVec[2][vInd] = lightScaledRawVec[0][vInd];
	}


	return;

	*/








	/*



	bool needCal = true;


	if (oldRefractSharpness == sharpness && oldRefractGlossy == isGlossy)
		needCal = false;
	
	meshes[0]->need_RidgeLines(false);

	vec3 eyeDir(0.f, 0.f, 1.f);

	vector<point> &valleyPoint = meshes[0]->ridgeLines.valleyPoint;
	vector<vec> &valleyPointNormal = meshes[0]->ridgeLines.valleyPointNormal;
	vector<float> &valleyPointPcurv = meshes[0]->ridgeLines.valleyPointPcurv;
	int featureNum = valleyPoint.size();

	if (designRefractLightMid.size() == 0)
		designRefractLightMid.resize(featureNum, vertexNum);

	if (oldRefractSharpness != sharpness || oldRefractGlossy != isGlossy)
	{

		int testVN;
		FILE *fin = fopen(meshes[0]->dFName.getDesignedRefractLightDN(sharpness, isGlossy), "rb");
		if (fin != NULL)
		{
			fread(&testVN, sizeof(int), 1, fin);
			if (testVN == vertexNum)
				needCal = false;
			fclose(fin);
		}
	}


	//calculate all
	if (needCal)
	{

		float distance2 = meshes[0]->feature_size() * 2;
		distance2 *= distance2;


		int *neibourNum = new int[featureNum];

		for (int vVer = 0; vVer < vertexNum; vVer++)
		{
			for (int fVer = 0; fVer < featureNum; fVer++)
				designRefractLightMid(fVer, vVer) = 0.f;
		}
		
		for (int fVer = 0; fVer < featureNum; fVer++)
			neibourNum[fVer] = 0;


		//find near vertex
		for (int vVer = 0; vVer < vertexNum; vVer++)
		{


			point &currentVertex = meshes[0]->vertices[vVer];

			//find valley vertex
			for (int fVer = 0; fVer < featureNum; fVer++)
			{

				point &valleyVertex = valleyPoint[fVer];

				//find near feature point
				float length2 = len2(currentVertex - valleyVertex);

				if (length2 > distance2)
					continue;

				neibourNum[fVer]++;
				//cal light weights
				for (int lVer = 0; lVer < vertexNum; lVer++)
				{
					float distance = len(currentVertex - meshes[0]->vertices[lVer]);

					designRefractLightMid(fVer, lVer) += distance;

				}
			}
		}

		std::cout << "Writing Designed Lighting..." << endl;



		FILE *fout = fopen(meshes[0]->dFName.getDesignedRefractLightDN(sharpness, isGlossy), "wb");
		fwrite(&vertexNum, sizeof(int), 1, fout);
		fwrite(&featureNum, sizeof(int), 1, fout);


		vector<vec3> &currentObject = meshes[0]->vertices;



		float *s1 = new float[vertexNum];
		float *s2 = new float[vertexNum];

	
		float chageDisToWeight = 0.5f / meshes[0]->bsphere.r;

		//find valley vertex
		for (int fVer = 0; fVer < featureNum; fVer++)
		{

			point &valleyVertex = valleyPoint[fVer];

			//cal light weights
			for (int lVer = 0; lVer < vertexNum; lVer++)
			{

				float distance = len(valleyVertex - meshes[0]->vertices[lVer]);

				float weight = 1.0 - distance * chageDisToWeight;


				if (weight > 0.00001)
					weight = pow(weight, sharpness);
				else if (weight < 0.00001)
					weight = -pow(abs(weight), sharpness);

				designRefractLightMid(fVer, lVer) -= distance * neibourNum[fVer];

				designRefractLightMid(fVer, lVer) *= weight;
				
				s1[lVer] = designRefractLightMid(fVer, lVer);
				s2[lVer] = distance  * chageDisToWeight;

			}


			fwrite(s1, sizeof(float), vertexNum, fout);
			fwrite(s2, sizeof(float), vertexNum, fout);
			//			DesignedLights_of << endl;


//			DesignedLights_of << endl;
		}
		

//		DesignedLights_of.close();

		fclose(fout);
		delete[] s1;
		delete[] s2;
		delete[] neibourNum;
	}




	VectorXf& lightMidRefractRawVec = wavTranEnvironment.getRawMidRefractLightVec(vertexNum);
	VectorXf& lightNearDistance = wavTranEnvironment.getRawRefractNearDistance(vertexNum);



	if ((oldRefractSharpness != sharpness || oldRefractGlossy != isGlossy) && !needCal)
	{



		bool *isVertexVisible = new bool[vertexNum];
		bool *isFeatureVisible = new bool[featureNum];


		for (int fVer = 0; fVer < featureNum; fVer++)
		{
			if (eyeDir.dot(valleyPointNormal[fVer]) < 0.f)
				isFeatureVisible[fVer] = false;
			else
				isFeatureVisible[fVer] = true;

		}

		for (int vVer = 0; vVer < vertexNum; vVer++)
		{
			if (eyeDir.dot(meshes[0]->normals[vVer]) < 0.f)
				isVertexVisible[vVer] = false;
			else
				isVertexVisible[vVer] = true;

			lightNearDistance[vVer] = 9999999.f;
		}




		std::cout << "Reading Designed Lighting..." << endl;
		int testVN;



		float *s1 = new float[vertexNum];
		float *s2 = new float[vertexNum];
		FILE *fin = fopen(meshes[0]->dFName.getDesignedRefractLightDN(sharpness, isGlossy), "rb");

		fread(&testVN, sizeof(int), 1, fin);
		fread(&featureNum, sizeof(int), 1, fin);




		for (int fVer = 0; fVer < featureNum; fVer++)
		{

			fread(s1, sizeof(float), vertexNum, fin);
			fread(s2, sizeof(float), vertexNum, fin);


			for (int lVer = 0; lVer < vertexNum; lVer++)
			{
				designRefractLightMid(fVer, lVer) = s1[lVer];

				if (!isFeatureVisible[fVer] || isVertexVisible[lVer])
					designRefractLightMid(fVer, lVer) = 0.f;
				else
				if (s2[lVer] < lightNearDistance[lVer])
					lightNearDistance[lVer] = s2[lVer];
			}
		}

		for (int lVer = 0; lVer < vertexNum; lVer++)
		if (isVertexVisible[lVer])
			lightNearDistance[lVer] = 1.f;



		fclose(fin);
		delete[] s1;
		delete[] s2;


		
		
		delete[] isVertexVisible;
		delete[] isFeatureVisible;
	}

	if (oldRefractSharpness != sharpness || oldRefractGlossy != isGlossy || oldRefractEnhScale != enhScale || oldEvnLightIndex != evnLightIndex)
	{


		VectorXf* maxRefractLightRadiance = wavTranEnvironment.getRawScaledLightVecRefract(vertexNum);
#pragma omp parallel for
		for (int lVer = 0; lVer < vertexNum; lVer++)
		{
			maxRefractLightRadiance[0][lVer] = 0.f;
			maxRefractLightRadiance[1][lVer] = 9999999.f;
		}


		for (int lVer = 0; lVer < vertexNum; lVer++)
			lightMidRefractRawVec[lVer] = 0.f;

		float tempEnhScale = enhScale / curv_scale;

#pragma omp parallel for
		//for (int fVer = evnLightIndex; fVer < evnLightIndex+1; fVer++)
		for (int fVer = 0; fVer < featureNum; fVer++)
		{

			float weight = exp(valleyPointPcurv[fVer] * tempEnhScale);

			for (int lVer = 0; lVer < vertexNum; lVer++)
			{
				if (designRefractLightMid(fVer, lVer) == 0.f)
					continue;

				float lightRadiance = designRefractLightMid(fVer, lVer) * weight;

				lightMidRefractRawVec[lVer] += lightRadiance;

				if (lightRadiance > maxRefractLightRadiance[0][lVer])
					maxRefractLightRadiance[0][lVer] = lightRadiance;
				if (lightRadiance < maxRefractLightRadiance[1][lVer])
					maxRefractLightRadiance[1][lVer] = lightRadiance;

			}
		}


		float gamma = 0.f;



		
#pragma omp parallel for
		for (int vecI = 0; vecI < vertexNum; vecI++)
		{
			if (gamma < abs(lightMidRefractRawVec[vecI]))
			gamma = abs(lightMidRefractRawVec[vecI]);
			
		}

		gamma = 100.f / gamma /vertexNum;

		#pragma omp parallel for
		for (int vecI = 0; vecI < vertexNum; vecI++)
		{
			lightMidRefractRawVec[vecI] *= gamma;
		}
		
	}


	VectorXf* lightScaledRawVecRefract = wavTranEnvironment.getRawScaledLightVecRefract(vertexNum);






	if (!isVertexScaledLightInited)
	{
		LoadRaw(true, true, true, isIndirect, isRGB, false, isGlossy, false);

		VectorXf* lightVec = wavTranEnvironment.getRawLightVec();

		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			int tIndex = rgbInd<RGBnum ? rgbInd : 0;
			lightScaledRawVecRefract[rgbInd] = rawReflectTransportMat[tIndex] * lightVec[rgbInd];
		}
		isVertexScaledLightInited = true;
	}




	if (oldRefractSharpness != sharpness || oldRefractGlossy != isGlossy || oldRefractEnhScale != enhScale || oldRefractScale != scale || oldEvnLightIndex != evnLightIndex)
	{

		float averageLight = wavTranEnvironment.getRawLightAvg();

		float gamma = 0.f;

//		scale = scale * averageLight;
#pragma omp parallel for
		for (int vecI = 0; vecI < vertexNum; vecI++)
		{


			Color temp(lightScaledRawVecRefract[0][vecI], lightScaledRawVecRefract[1][vecI], lightScaledRawVecRefract[2][vecI]);

			temp = temp.convert(temp.RGB, temp.CIELAB);

			temp[0] = lightMidRefractRawVec[vecI] * pow(lightNearDistance[vecI], scale);

			temp = temp.convert(temp.CIELAB, temp.RGB);
			lightScaledRawVecRefract[0][vecI] = temp[0];
			lightScaledRawVecRefract[1][vecI] = temp[1];
			lightScaledRawVecRefract[2][vecI] = temp[2];
		}

		needRelight = true;

		oldRefractSharpness = sharpness;
		oldRefractGlossy = isGlossy;
		oldRefractEnhScale = enhScale;
		oldRefractScale = scale;
		oldEvnLightIndex = evnLightIndex;
	}
	*/
	
}

void Transport::setColorWav(vector<Color>& color, bool isIndirect, bool isRefraction, bool enhCur, bool isRGB, bool isGlossy)
{

	// Do the matrix multiplication


	LoadWav(isIndirect, isRefraction, isRGB, isGlossy);




	// Do the matrix multiplication
	if (enhCur)
	{
		if (wavScaledResult[0].size() == 0)
		{

			Eigen::SparseVector<float>* lightVec = wavTranEnvironment.getWavScaledLightVec();

#pragma omp parallel for
			for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			{
				int tIndex = rgbInd<RGBnum ? rgbInd : 0;
				wavScaledResult[rgbInd] = wavTransportMat[tIndex] * lightVec[rgbInd];
			}
		}


		// Quantize/Gamma correct values, write to pixBuf
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			const int* indices = wavScaledResult[rgbInd]._innerIndexPtr();
			const float* values = wavScaledResult[rgbInd]._valuePtr();
			int numCoeffs = wavScaledResult[rgbInd].nonZeros();

			for (int ind = 0, sparseInd = 0; ind < vertexNum; ind++)
			{
				// Check if this value exists in the sparse result

				if (sparseInd >= numCoeffs || ind != indices[sparseInd])
					color[ind][rgbInd] = 0.f;
				else if (ind == indices[sparseInd])
				{
					color[ind][rgbInd] = (values[sparseInd]);
					sparseInd++;
				}
			}

		}
	}
	else
	{


		if (wavResult[0].size() == 0)
		{
			Eigen::SparseVector<float>* lightVec = wavTranEnvironment.getWavLightVec();

#pragma omp parallel for
			for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			{
				int tIndex = rgbInd<RGBnum ? rgbInd : 0;
				wavResult[rgbInd] = wavTransportMat[tIndex] * lightVec[rgbInd];
			}
		}  


		// Quantize/Gamma correct values, write to pixBuf
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			const int* indices = wavResult[rgbInd]._innerIndexPtr();
			const float* values = wavResult[rgbInd]._valuePtr();
			int numCoeffs = wavResult[rgbInd].nonZeros();

			for (int ind = 0, sparseInd = 0; ind < vertexNum; ind++)
			{
				// Check if this value exists in the sparse result

				if (sparseInd >= numCoeffs || ind != indices[sparseInd])
					color[ind][rgbInd] = 0.f;
				else if (ind == indices[sparseInd])
				{
					color[ind][rgbInd] = (values[sparseInd]);
					sparseInd++;
				}
			}

		}


	}



}

void Transport::initTexture(unsigned int texCubeMap[], bool isOri, bool isAvg, bool isWavelet, bool isEnhReflect, bool isEnhRefract)
{ 
	if (isOri || isEnhRefract)
		wavTranEnvironment.initTexture(texCubeMap);
	else if (isEnhReflect)
		wavTranEnvironment.initTextureScaledRaw(texCubeMap);
	else if (isWavelet)
		wavTranEnvironment.initTextureWav(texCubeMap);
	else if (isAvg)
		wavTranEnvironment.initTextureAvg(texCubeMap);
	else
		wavTranEnvironment.initTextureRaw(texCubeMap);

	return;
}

//calculate the result
void Transport::setColor(vector<Color>& color, VectorXf *result, VectorXf* lightVec, MatrixXf *transportMat, float radianceScale)
{


	if (result[0].size() == 0 || needRelight)
	{
		needRelight = false;
#pragma omp parallel for
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			int tIndex = rgbInd<RGBnum ? rgbInd : 0;
			result[rgbInd] = transportMat[tIndex] * lightVec[rgbInd];
		}
	}
#pragma omp parallel for
	for (int i = 0; i < vertexNum; i++)
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		color[i][rgbInd] += result[rgbInd][i] * radianceScale;
}

//set ab of lab
void Transport::setAB(VectorXf *result, VectorXf* lightVec, MatrixXf *transportMat)
{
	int csl2 = cubeSideLength*cubeSideLength;
	int fcsl2 = csl2 * 6;
	float gamma = 100.f/ fcsl2;

	if (result[0].size() == 0 || needRelight)
	{
		result[0] = transportMat[0] * lightVec[0];
		result[1] = transportMat[0] * lightVec[1];
		result[2] = transportMat[0] * lightVec[2];
	}



}
//set ab of lab
void Transport::setL(VectorXf *result, VectorXf* lightVec, MatrixXf *transportMat)
{
	if (result[2].size() == 0 || needRelight)
	{
		result[2] = transportMat[0] * lightVec[2];
	}
}

void Transport::transHVSToRGB(vector<Color>& color, VectorXf *source, float transX, float transY, float rotAngle, float scale)
{

	if (scale < 0)
		scale = -1.f / scale;

	if (scale >= 1.0)
		scale = 0.9;
	scale = 0.0f;
	float gamma = designReflectLightTransA;
	float gamma2 = 1.0/ (1.0 - scale);
	needRelight = true;

#pragma omp parallel for
	for (int vecI = 0; vecI < vertexNum; vecI++)
	{
		color[vecI][0] = (source[0][vecI] * gamma - scale) * gamma2;
		color[vecI][1] = (source[1][vecI] * gamma - scale)* gamma2;
		color[vecI][2] = (source[2][vecI] * gamma - scale)* gamma2;
	}

	/*
	if (scale != 1.f)
	{
#pragma omp parallel for
		for (int vecI = 0; vecI < vertexNum; vecI++)
		{
			color[vecI][2] *= scale;
		}

	}

	if (rotAngle != 0.f)
	{
#pragma omp parallel for
		for (int vecI = 0; vecI < vertexNum; vecI++)
		{
			color[vecI][0] += rotAngle;
		}
	}

	
	if (rotAngle != 0.f)
	{
#pragma omp parallel for
		for (int vecI = 0; vecI < vertexNum; vecI++)
		{
			float length = sqrt(color[vecI][1] * color[vecI][1] + color[vecI][2] * color[vecI][2]);
			float oriAngle = atan2(color[vecI][2], color[vecI][1]);
			oriAngle += rotAngle;
			color[vecI][1] = cos(oriAngle) * length;
			color[vecI][2] = sin(oriAngle) * length;

		}

	}


	for (int vecI = 0; vecI < vertexNum; vecI++)
	{
		color[vecI] = color[vecI].convert(color[vecI].HSV, color[vecI].RGB);
	}
	*/

	return;

}


// pixBuf is assumed to be of size 3*imageSideLength*imageSideLength
//method for calling
void Transport::setColorRaw(vector<Color>& color, bool isIndirect, bool isRGB, bool isReflect, bool isRefract, 
	bool enhReflect, bool enhRefract, int evnLightIndex, float radianceScale, bool isRefractGlossy, MaterialType useMtaerial, EyeType useEye)
{

	LoadRaw(enhRefract, isRefract, isReflect, isIndirect, isRGB, isRefractGlossy, enhReflect, useMtaerial, useEye);


#pragma omp parallel for
	for (int i = 0; i < vertexNum; i++)
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		color[i][rgbInd] = 0.f;


	if (isReflect)
	{

		if (enhReflect)
		{

			setAB(rawScaledReflectResult, wavTranEnvironment.getLabLightVecColor(), rawReflectTransportMat);// &rawReflectTransportMatNoOcl);
//			setL(rawScaledReflectResult, wavTranEnvironment.getLabLightVecColor(), rawReflectTransportMat);
			transHVSToRGB(color, rawScaledReflectResult, 0, 0, 0.f, radianceScale);
		}
		else
		{
			setColor(color, rawReflectResult, wavTranEnvironment.getRawLightVec(), rawReflectTransportMat, radianceScale);
		}
	}

	if (isRefract)
	{
		if (enhRefract)
			setColor(color, rawScaledRefractResult, wavTranEnvironment.getRawScaledLightVecRefract(vertexNum), rawRefractTransportMat, radianceScale);
		else
			setColor(color, rawRefractResult, wavTranEnvironment.getRawLightVec(), rawRefractTransportMat, radianceScale);
	}

			/*

			for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			for (int v = 0; v < vertexNum; v++)
			{
			float a = lightVec[0][v];
			float b = lightVec[1][v];
			float c = lightVec[2][v];
			rawResult[rgbInd][v] = lightVec[rgbInd][v];


			}



			meshes[0]->need_RidgeLines(false);

			float distance2 = meshes[0]->feature_size() * 2;
			distance2 *= distance2;


			vector<point> &valleyPoint = meshes[0]->valleyPoint;

			//find near vertex
			for (int vVer = 0; vVer < vertexNum; vVer++)
			{


			point &currentVertex = meshes[0]->vertices[vVer];

			//find valley vertex


			int fVer = evnLightIndex;

			point &valleyVertex = valleyPoint[fVer];

			//find near feature point
			float length2 = len2(currentVertex - valleyVertex);

			if (length2 > distance2)
			continue;

			rawResult[0][vVer] = 1.0;
			}
			*/



	/*


	// Do the matrix multiplication
	if (enhCur)
	{




		vec3 eye = meshes[0]->bsphere.center;
		eye[2] += meshes[0]->bsphere.r*2.f
		;

		vec curvDir = meshes[0]->pdir1[evnLightIndex];
		vec normal = meshes[0]->normals[evnLightIndex];
		vec lightDir = curvDir*(float)sin(0.464) + normal*(float)cos(0.464);

		normalize(lightDir);

		#pragma omp parallel for
		for (int i = 0; i < vertexNum; i++)
		{
		color[i][0] = meshes[0]->normals[i].dot(lightDir);


		float lightT = color[i][0] * color[i][0];


		lightT *= lightT;

		//			color[i][0] = lightT * sqrt(1.0 - color[i][0] * color[i][0])*3.3f;


		color[i][0] = color[i][0] * lightT * 0.6f + 0.4f;
		color[i][1] = color[i][0];
		color[i][2] = color[i][0];
		}
		color[evnLightIndex][0] = 1.f;

		color[evnLightIndex][1] = 0.f;
		color[evnLightIndex][2] = 0.f;


	}
	
		*/


}

