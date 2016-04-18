
#include "stdafx.h"
#include "Environment.h"
#include "FreeImage/FreeImage.h"
#include <omp.h>
#pragma comment(lib, "FreeImage.lib")

#include"../GL/glew.h"
#include"../GL/glut.h"

///////////////////////////////////////////////////////////////////////
//                           Environment                             //
///////////////////////////////////////////////////////////////////////

Environment::~Environment()
{
	if (!isBaseInit)
		return;
	// Clean up sampleDirs
	for (int f = 0; f < 6; f++)
	{
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				delete[] sampleDirs[f][y][x];
			}
			delete[] sampleDirs[f][y];
		}
		delete[] sampleDirs[f];
	}
	
}

void Environment::LoadImage()
{


	FIBITMAP* env = FreeImage_Load(FIF_PFM, fileDirectionName);
	if (!env)
	{
		cout << "ERROR: Failed to open WavTranEnvironment \"" << fileDirectionName << "\"" << endl;
		exit(1);
	}
	FIRGBF* pixels = (FIRGBF*)FreeImage_GetBits(env);

	sourceSideLength = FreeImage_GetWidth(env) / 3;

	// Define size of faces
	for (int f = 0; f < 6; f++)
	{
		faces[f][0].resize(sourceSideLength, sourceSideLength);
		faces[f][1].resize(sourceSideLength, sourceSideLength);
		faces[f][2].resize(sourceSideLength, sourceSideLength);
	}

	// Populate faces
	int width = 3 * sourceSideLength;

	maxGamma = 0.f;

	for (int f = 0; f < 6; f++)
		for (int y = 0; y < sourceSideLength; y++)
			for (int x = 0; x < sourceSideLength; x++)
			{
				// Which pixels should we grab out of the image?
				int i;
				switch (f)
				{
				case 0:
					i = y*width + (x + sourceSideLength); break;
				case 1:
					i = (y + sourceSideLength)*width + x; break;
				case 2:
					i = (y + sourceSideLength)*width + (x + sourceSideLength); break;
				case 3:
					i = (y + sourceSideLength)*width + (x + 2 * sourceSideLength); break;
				case 4:
					i = (y + 2 * sourceSideLength)*width + (x + sourceSideLength); break;
				case 5:
					i = (y + 3 * sourceSideLength)*width + (x + sourceSideLength); break;
				default: break;
				}

				// Write pixels to appropriate place in face matrix
				faces[f][0](y, x) = pixels[i].red;
				faces[f][1](y, x) = pixels[i].green;
				faces[f][2](y, x) = pixels[i].blue;

				if (faces[f][0](y, x) > 1.0)
				{
					float a = faces[f][0](y, x);
					a = faces[f][0](y, x);
				}

				for (int rgbInd = 0; rgbInd < 3; rgbInd++)
					if (faces[f][rgbInd](y, x)>maxGamma)
						maxGamma = faces[f][rgbInd](y, x);

			}



}

void Environment::initBase()
{
	if (isBaseInit)
		return;
	isBaseInit = true;
	
	LoadImage();


	//ComputeSampleDirs;

	// Allocate space

	for (int f = 0; f < 6; f++)
	{
		sampleDirs[f] = new Vector3f** [cubeSideLength];
		for (int y = 0; y < cubeSideLength; y++)
		{
			sampleDirs[f][y] = new Vector3f* [cubeSideLength];
			for (int x = 0; x < cubeSideLength; x++)
			{
				sampleDirs[f][y][x] = new Vector3f [samplesPerPixel];
			}
		}
	}

	// Generate directions
	for (int f = 0; f < 6; f++)
	{
#pragma omp parallel for
		for (int y = 0; y < cubeSideLength; y++)
		for (int x = 0; x < cubeSideLength; x++)
		{
			float minU, minV, maxU, maxV;
			minU = (float)x/cubeSideLength;
			minV = (float)y/cubeSideLength;
			maxU = (float)(x+1)/cubeSideLength;
			maxV = (float)(y+1)/cubeSideLength;
			for (int s = 0; s < samplesPerPixel; s++)
			{
				float ru = ((float)rand())/RAND_MAX;
				float rv = ((float)rand())/RAND_MAX;
				float u = ru*minU + (1-ru)*maxU;
				float v = rv*minV + (1-rv)*maxV;
				CubeCoordToDir(f, u, v, sampleDirs[f][y][x][s]);
			}
		}
		
	}
}

void Environment::getDirectionColor(int f, int x, int y, vec3& outdir, Color& color)
{
	initBase();

	float u = ((float)x + 0.5) / cubeSideLength;
	float v = ((float)y + 0.5) / cubeSideLength;

	float a, b, c;

	a = 1.0f / sqrt(1 + (2 * u - 1)*(2 * u - 1) + (2 * v - 1)*(2 * v - 1));
	// Negative root for "negative" faces
	if (f == 1 || f == 4 || f == 5) a = -a;
	b = a*(2 * u - 1);
	c = a*(2 * v - 1);

	switch (f)
	{
	case 0: outdir = vec3(b, a, c); break;
	case 1: outdir = vec3(a, c, -b); break;
	case 2: outdir = vec3(b, -c, a); break;
	case 3: outdir = vec3(a, -c, -b); break;
	case 4: outdir = vec3(-b, a, c); break;
	case 5: outdir = vec3(-b, -c, a); break;
	default: break;
	}

	color = Color(0.0);

	// Do rotation, write back to light vector!
	int csl2 = cubeSideLength*cubeSideLength;
	float gamma = 1.0 / samplesPerPixel;
	
	// Use each precomputed sample direction, index
	// into high-res cube map (super sampling)
	for (int s = 0; s < samplesPerPixel; s++)
	{
		Vector3f dir;
		dir = sampleDirs[f][y][x][s];
		int newf;
		float newu, newv;
		DirToCubeCoord(newf, newu, newv, dir);
		int i = (int)(newu*(sourceSideLength - 1));
		int j = (int)(newv*(sourceSideLength - 1));
		// Normalized for numerical cubature

		Color a(
			faces[newf][0](i, j),
			faces[newf][1](i, j),
			faces[newf][2](i, j));

		color += a;
	}

	color *= gamma;


}

void Environment::CubeCoordToDir(int f, float u, float v, Vector3f& outdir)
{
	float a, b, c;

	a = 1.0f / sqrt(1 + (2*u-1)*(2*u-1) + (2*v-1)*(2*v-1));
	// Negative root for "negative" faces
	if (f == 1 || f == 4 || f == 5) a = -a;
	b = a*(2*u-1);
	c = a*(2*v-1);

	switch(f)
	{
	case 0: outdir = Vector3f(b,a,c); break;
	case 1: outdir = Vector3f(a,c,-b); break;
	case 2: outdir = Vector3f(b,-c,a); break;
	case 3: outdir = Vector3f(a,-c,-b); break;
	case 4: outdir = Vector3f(-b,a,c); break;
	case 5: outdir = Vector3f(-b,-c,a); break;
	default: break;
	}
}

void Environment::DirToCubeCoord(int& f, float& u, float& v, Vector3f& indir)
{
	float a, b, c;

	// Which face to look into?
	float maxcomp = max(fabs(indir[0]), max(fabs(indir[1]), fabs(indir[2])));
	if (maxcomp == fabs(indir[0])) f = (indir[0] > 0 ? 3 : 1);
	if (maxcomp == fabs(indir[1])) f = (indir[1] > 0 ? 0 : 4);
	if (maxcomp == fabs(indir[2])) f = (indir[2] > 0 ? 2 : 5);

	switch(f)
	{
	case 0: a = indir[1]; b = indir[0]; c = indir[2]; break;
	case 1: a = indir[0]; b = -indir[2]; c = indir[1]; break;
	case 2: a = indir[2]; b = indir[0]; c = -indir[1]; break;
	case 3: a = indir[0]; b = -indir[2]; c = -indir[1]; break;
	case 4: a = indir[1]; b = -indir[0]; c = indir[2]; break;
	case 5: a = indir[2]; b = -indir[0]; c = -indir[1]; break;
	}

	u = (a+b)/(2*a);
	v = (a+c)/(2*a);
}

void Environment::Export(string filename)
{
	initBase();
	int width = 3 * cubeSideLength;
	int height = 4 * cubeSideLength;
	FIBITMAP* output = FreeImage_AllocateT(FIT_RGBF, width, height);
	FIRGBF* pixels = (FIRGBF*)FreeImage_GetBits(output);

	// First fill up with zeros
	for (int p = 0; p < width*height; p++)
	{
		pixels[p].red = pixels[p].green = pixels[p].blue = 0.0f;
	}

	// Write faces
	for (int f = 0; f < 6; f++)
	{
		for (int y = 0; y < cubeSideLength; y++)
		{
			for (int x = 0; x < cubeSideLength; x++)
			{
				// Which pixels should we write into the image?
				int i;
				switch (f)
				{
				case 0:
					i = y*width + (x + cubeSideLength); break;
				case 1:
					i = (y + cubeSideLength)*width + x; break;
				case 2:
					i = (y + cubeSideLength)*width + (x + cubeSideLength); break;
				case 3:
					i = (y + cubeSideLength)*width + (x + 2 * cubeSideLength); break;
				case 4:
					i = (y + 2 * cubeSideLength)*width + (x + cubeSideLength); break;
				case 5:
					i = (y + 3 * cubeSideLength)*width + (x + cubeSideLength); break;
				default: break;
				}

				// Write pixels from face matrices
				pixels[i].red = faces[f][0](y, x);
				pixels[i].green = faces[f][1](y, x);
				pixels[i].blue = faces[f][2](y, x);
			}
		}
	}

	FreeImage_Save(FIF_PFM, output, filename.c_str());
}

bool Environment::getSourceColor(int f, int x, int y, Color &sourceColor, float weight)
{
	bool xChange = (x < 0 || x >= sourceSideLength) ? true : false;
	bool yChange = (y < 0 || y >= sourceSideLength) ? true : false;

	if (xChange && yChange)
		return false;

	if (x < 0)
	{
		x = -x - 1;
		if (f == 2 || f == 3 || f == 4) x = sourceSideLength - x - 1;
		if (f == 1 || f == 4 || f == 5) y = sourceSideLength - y - 1;

		if (f == 0 || f == 4)
		{
			int temp = x;
			x = y;
			y = temp;
		}

		switch (f)
		{
			case 0:	f = 1; break;
			case 1:	f = 5; break;
			case 2:	f = 1; break;
			case 3:	f = 2; break;
			case 4:	f = 1; break;
			case 5:	f = 1; break;
			default: break;
		}
	}

	if (x >= sourceSideLength)
	{
		x = x - sourceSideLength;
		if (f == 3 || f == 4 || f == 5) x = sourceSideLength - x - 1;
		if (f == 0 || f == 3 || f == 5) y = sourceSideLength - y - 1;

		if (f == 0 || f == 4)
		{
			int temp = x;
			x = y;
			y = temp;
		}

		switch (f)
		{
			case 0:	f = 3; break;
			case 1:	f = 2; break;
			case 2:	f = 3; break;
			case 3:	f = 5; break;
			case 4:	f = 3; break;
			case 5:	f = 3; break;
			default: break;
		}

	}

	if (y < 0)
	{
		y = -y - 1;
		if (f == 3) x = sourceSideLength - x - 1;
		if (f != 1) y = sourceSideLength - y - 1;

		if (f == 1 || f == 3)
		{
			int temp = x;
			x = y;
			y = temp;
		}

		switch (f)
		{
			case 0:	f = 5; break;
			case 1:	f = 0; break;
			case 2:	f = 0; break;
			case 3:	f = 0; break;
			case 4:	f = 2; break;
			case 5:	f = 4; break;
			default: break;
		}
	}


	if (y >= sourceSideLength)
	{
		y = y - sourceSideLength;
		if (f == 1) x = sourceSideLength - x - 1;
		if (f == 3) y = sourceSideLength - y - 1;

		if (f == 1 || f == 3)
		{
			int temp = x;
			x = y;
			y = temp;
		}

		switch (f)
		{
			case 0:	f = 2; break;
			case 1:	f = 4; break;
			case 2:	f = 4; break;
			case 3:	f = 4; break;
			case 4:	f = 5; break;
			case 5:	f = 0; break;
			default: break;
		}

	}	
	sourceColor[0] = faces[f][0](y, x);
	sourceColor[1] = faces[f][1](y, x);
	sourceColor[2] = faces[f][2](y, x);
	return true;
}

Color Environment::getAverageColor(vec3 direction)
{
	initBase();



	int newf;
	float newu, newv;


	Vector3f dir(direction[0], direction[1], direction[2]);


	DirToCubeCoord(newf, newu, newv, dir);

	int i = (int)(newu*(sourceSideLength - 1));
	int j = (int)(newv*(sourceSideLength - 1));
	// Normalized for numerical cubature

	int halfSuperPixelSize = sourceSideLength / cubeSideLength * 0.5f;


	Color sampleColor;

	Color a(0.f);
	float weightSum = 0.f;

	int startx = i - halfSuperPixelSize;
	int endx = i + halfSuperPixelSize;
	int starty = j - halfSuperPixelSize;
	int endy = j + halfSuperPixelSize;
	for (int x = startx; x < endx; x++)
	{
		for (int y = starty; y < endy; y++)
		{
			float weight = max(halfSuperPixelSize - sqrt((x - i)*(x - i) + (y - j)*(y - j)), 0.0);

			if (getSourceColor(newf, x, y, sampleColor, weight))
			{
				a += sampleColor * weight;
				weightSum += weight;
			}

		}

	}
	a /= weightSum;
	return a;

}

void Environment::initTextureHole(unsigned int texCubeMap[], vec3 direction)
{
	initBase();
	int ssl2 = sourceSideLength*sourceSideLength * 3;
	float* pixel = new float[ssl2];

	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < sourceSideLength; y++)
		for (int x = 0; x < sourceSideLength; x++)
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			faces[f][rgbInd](y, x) = 0.5;
	}

	int newf;
	float newu, newv;


	Vector3f dir(direction[0], direction[1], direction[2]);


	DirToCubeCoord(newf, newu, newv, dir);

	int i = (int)(newu*(sourceSideLength - 1));
	int j = (int)(newv*(sourceSideLength - 1));
	// Normalized for numerical cubature

	int halfSuperPixelSize = sourceSideLength / cubeSideLength * 0.5f;


	Color sampleColor;

	Color a(0.f);

	int startx = i - halfSuperPixelSize;
	int endx = i + halfSuperPixelSize;
	int starty = j - halfSuperPixelSize;
	int endy = j + halfSuperPixelSize;
	for (int x = startx; x < endx; x++)
	{
		for (int y = starty; y < endy; y++)
		{

			float weight = max(halfSuperPixelSize - sqrt((x - i)*(x - i) + (y - j)*(y - j)), 0.0);


			getSourceColor(newf, x, y, sampleColor, weight);
		}
	}



	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < sourceSideLength; y++)
		for (int x = 0; x < sourceSideLength; x++)
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			pixel[pixelIndex++] = faces[f][rgbInd](y, x);

		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, sourceSideLength, sourceSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}
	delete[] pixel;


	return;
}

void Environment::initTexture(unsigned int texCubeMap[])
{
	initBase();
	int ssl2 = sourceSideLength*sourceSideLength * 3;
	float* pixel = new float[ssl2];


	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < sourceSideLength; y++)
		for (int x = 0; x < sourceSideLength; x++)
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
			pixel[pixelIndex++] = faces[f][rgbInd](y, x);

		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, sourceSideLength, sourceSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}
	delete[] pixel;


	return;
	


	float *luminance = new float[sourceSideLength*sourceSideLength*6];


	int lightIndex = 0;
	float rawAverageLight = 0.f;
	for (int f = 0; f < 6; f++)
	for (int y = 0; y < sourceSideLength; y++)
	for (int x = 0; x < sourceSideLength; x++)
	{
		luminance[lightIndex] = 0.2989*faces[f][0](y, x) + 0.587*faces[f][1](y, x) + 0.114*faces[f][2](y, x);
		rawAverageLight += luminance[lightIndex];
		lightIndex++;
	}
	rawAverageLight /= (float)(sourceSideLength*sourceSideLength * 6);
	rawAverageLight = 0.25;

	lightIndex = 0;
	for (int f = 0; f < 6; f++)
	{
		int pixelIndex = 0;
		for (int y = 0; y < sourceSideLength; y++)
		for (int x = 0; x < sourceSideLength; x++)
		{
			luminance[lightIndex] = rawAverageLight / luminance[lightIndex];

			for (int rgbInd = 0; rgbInd < 3; rgbInd++)
				pixel[pixelIndex++] = faces[f][rgbInd](y, x) * luminance[lightIndex];

			lightIndex++;
		}
		glBindTexture(GL_TEXTURE_2D, texCubeMap[f]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, sourceSideLength, sourceSideLength, 0, GL_RGB, GL_FLOAT, pixel);
	}
	delete[] pixel;




	delete[] luminance;
}

void Environment::changeEnvironmentMap(char *newFile)
{

	char* filepath = "..\\Image\\";
	char* fileAfter = ".pfm";

	strcpy_s(fileDirectionName, 80, filepath);
	strcat_s(fileDirectionName, 80, newFile);
	strcat_s(fileDirectionName, 80, fileAfter);
	strcpy_s(fileName, 80, newFile);


	LoadImage();




}