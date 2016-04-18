
#include "stdafx.h"

#include "SHLighting.h"

///////////////////////////////////////////////////////////////////////
//                 Spherical Harmonic Equation                       //
///////////////////////////////////////////////////////////////////////

//Evaluate an Associated Legendre Polynomial P(l, m) at x
double P(int l, int m, double x)
{
	//First generate the value of P(m, m) at x
	double pmm=1.0;

	if(m>0)
	{
		double sqrtOneMinusX2=sqrt(1.0-x*x);

		double fact=1.0;

		for(int i=1; i<=m; ++i)
		{
			pmm*=(-fact)*sqrtOneMinusX2;
			fact+=2.0;
		}
	}

	//If l==m, P(l, m)==P(m, m)
	if(l==m)
		return pmm;

	//Use rule 3 to calculate P(m+1, m) from P(m, m)
	double pmp1m=x*(2.0*m+1.0)*pmm;

	//If l==m+1, P(l, m)==P(m+1, m)
	if(l==m+1)
		return pmp1m;

	//Otherwise, l>m+1.
	//Iterate rule 1 to get the result
	double plm=0.0;

	for(int i=m+2; i<=l; ++i)
	{
		plm=((2.0*i-1.0)*x*pmp1m-(i+m-1.0)*pmm)/(i-m);
		pmm=pmp1m;
		pmp1m=plm;
	}

	return plm;
}

//Calculate the normalisation constant for an SH function
//No need to use |m| since SH always passes positive m
double K(int l, int m)
{
	double temp=((2.0*l+1.0)*Factorial(l-m))/((4.0*M_PI)*Factorial(l+m));

	return sqrt(temp);
}

//Sample a spherical harmonic basis function Y(l, m) at a point on the unit sphere
double SH(int l, int m, double theta, double phi)
{
	const double sqrt2=sqrt(2.0);

	if(m==0)
		return K(l, 0)*P(l, m, cos(theta));

	if(m>0)
		return sqrt2*K(l, m)*cos(m*phi)*P(l, m, cos(theta));

	//m<0
	return sqrt2*K(l,-m)*sin(-m*phi)*P(l, -m, cos(theta));
}

//Calculate n! (n>=0)
int Factorial(int n)
{
	if(n<=1)
		return 1;

	int result=n;

	while(--n > 1)
		result*=n;

	return result;
}

///////////////////////////////////////////////////////////////////////
//                          Samples                                  //
///////////////////////////////////////////////////////////////////////

SampleInfo::SampleInfo(Environment &environment, int cSL, int numB = 4) :cubeSideLength(cSL), numBands(numB)
{
	if (cubeSideLength <= 0)
	{
		cout << "The number of sample is negative!" << cubeSideLength << endl;
	}
	numSamples = 6 * cubeSideLength*cubeSideLength;
	samples.resize(numSamples);
	GenerateSamples(environment);

}

SampleInfo::~SampleInfo()
{
}


bool SampleInfo::GenerateSamples(Environment &environment)
{
	int numFunctions=numBands*numBands;

	//Create space for the SH values in each sample
	for(int i=0; i<numSamples; ++i)
	{
		samples[i].shValues=new double[numFunctions];
		if(!samples[i].shValues)
		{
			cout<<"Unable to allocate space for SH values in samples"<<endl;
			return false;
		}
	}

	int index=0;
	
	for (int f = 0; f < 6; f++)
	for (int y = 0; y < cubeSideLength; y++)
	for (int x = 0; x < cubeSideLength; x++)
	{

		environment.getDirectionColor(f, x, y, samples[index].direction, samples[index].color);
		

		vec3 &d = samples[index].direction;

		samples[index].theta = acos(d[2]);

		if (d[2] >= 1.f)
			samples[index].phi = 0.f;
		else
		{
			vec2 dir(d[0], d[1]);
			normalize(dir);

			if (dir[1] < 0)
				samples[index].phi = M_PI*2.0 - acos(dir[0]);
			else
				samples[index].phi = acos(dir[0]);

		}
		
		//Compute SH coefficients for this sample

		for (int l = 0; l<numBands; ++l)
		for (int m = -l; m <= l; ++m)
		{
			int index2 = l*(l + 1) + m;

			samples[index].shValues[index2] = SH(l, m, samples[index].theta, samples[index].phi);
		}
		
		index++;

	}
	
	

	return true;
}


///////////////////////////////////////////////////////////////////////
//                      SHDirectCoeffs                               //
///////////////////////////////////////////////////////////////////////

SHDirectCoeffs::SHDirectCoeffs(vector<TriMesh *> &objects, int numB = 4):numBands(numB), meshes(objects)
{

}

SHDirectCoeffs::~SHDirectCoeffs()
{
}


void SHDirectCoeffs::calculateSHDirectCoeffs(SampleInfo &sampleInfo)
{


	int numObjects=meshes.size();
	vertexSHCoffs.resize(numObjects);

	vector<SAMPLE> &samples = sampleInfo.getSamples();
	const int numFunctions=numBands*numBands;
	for(int i=0; i<numObjects; ++i)
	{
		TriMesh* mesh = meshes[i];
		int numVertices = mesh->vertices.size();
		vertexSHCoffs[i].resize(numVertices);
		vector<SHVertex> &currentObejct = vertexSHCoffs[i];

		for(int j=0; j<numVertices; ++j)
		{
			SHVertex &currentVertex = currentObejct[j];
			currentVertex.unshadowedCoeffs=new double[numFunctions];
			currentVertex.shadowedCoeffs=new double[numFunctions];

			if(!currentVertex.unshadowedCoeffs || !currentVertex.shadowedCoeffs)
				cout<<"Unable to create space for vertex SH coefficients"<<endl;
		}

		calculateSHDirectCoeffsSingle(samples, mesh, i);


	}


}

void SHDirectCoeffs::calculateSHDirectCoeffsSingle(vector<SAMPLE> &samples, TriMesh* object, int objectIndex)
{

	const int numFunctions=numBands*numBands;
	//Is there a file containing the coefficients, or do they need to be regenerated?
	bool regenerateCoeffs=false;
	int numSamples = samples.size();
	vector<SHVertex> &currentObjectSH=vertexSHCoffs[objectIndex];
	const int numVertices = currentObjectSH.size();

	//Does the file exist?
	std::ifstream inFile(object->dFName.getSHDirectCoeffsDN(numSamples), std::ios::in | std::ios::binary);
	if(!inFile.is_open())
	{
		cout<<"Unable to open directcoeffs.dat, regenerating coefficients..."<<endl;
		regenerateCoeffs=true;
	}


	//Are the number of bands and samples in the file correct?
	if(!regenerateCoeffs)
	{
		int numFileBands, numFileSamples;

		inFile.read((char *)&numFileBands, sizeof(int));
		inFile.read((char *)&numFileSamples, sizeof(int));

		if(numFileBands!=numBands || numFileSamples!=numSamples)
		{
			cout<<"directcoeffs.dat has different number of bands/samples, regenerating coefficients..."<<endl;
			regenerateCoeffs=true;
			inFile.close();
		}
	}

	//If the file is good, read in the coefficients
	if(!regenerateCoeffs)
	{

		cout<<"Reading Coeffs....";

		for(int j=0; j<numVertices; ++j)
		{

			SHVertex &currentVertexSH = currentObjectSH[j];

			inFile.read((char *)currentVertexSH.unshadowedCoeffs,
				numFunctions*sizeof(double));
			inFile.read((char *)currentVertexSH.shadowedCoeffs,
				numFunctions*sizeof(double));
		}

		cout<<"Done"<<endl;

		inFile.close();
		return;
	}


	BVH meshBVH(object);

	cout<<"Computing Coeffs...."<<endl;

	//Otherwise, regenerate the coefficients
	//Loop through the vertices and project the transfer function into SH space

	vector<vec3> &currentObjectNormal=object->normals;
	vector<vec3> &currentObject=object->vertices;

	cout<<"Calculating Transfer Function SH Coefficients"<<endl;
	cout<<"Direct Illumination"<<endl;
	cout<<"Object "<<objectIndex+1<<" of "<<vertexSHCoffs.size()<<endl;
	float numHun = 0.01f*(float)numVertices;
	int curPer = 0;

	for(int j=0; j<numVertices; ++j)
	{
		SHVertex &currentVertexSH=currentObjectSH[j];
		vec3 &currentVertexNormal=currentObjectNormal[j];
		point &currentVertex=currentObject[j];

		//Clear the coefficients for this vertex
		for(int k=0; k<numFunctions; ++k)
		{
			currentVertexSH.unshadowedCoeffs[k]=0.0;
			currentVertexSH.shadowedCoeffs[k]=0.0;
		}

		//Loop through samples
		for(int k=0; k<numSamples; ++k)
		{
			//Calculate cosine term for this sample
			double dot=(double)samples[k].direction.dot(currentVertexNormal);

			//Clamp to [0, 1]
			if(dot>0.0)
			{


				Ray ray(currentVertex+2*EPSILON*currentVertexNormal,
					samples[k].direction);



				//Fill in a RAY structure for this sample
				//				RAY ray(currentVertex+2*EPSILON*currentVertexNormal,
				//					samples[k].direction);

				//See if the ray is blocked by any object
				bool rayBlocked = meshBVH.Occluded(ray);

				//Add the contribution of this sample to the coefficients
				for(int l=0; l<numFunctions; ++l)
				{
					double contribution=dot*samples[k].shValues[l];
					currentVertexSH.unshadowedCoeffs[l]+=contribution;

					if(!rayBlocked)
						currentVertexSH.shadowedCoeffs[l]+=contribution;
				}
			}
		}

		//Rescale the coefficients
		for(int k=0; k<numFunctions; ++k)
		{
			currentVertexSH.unshadowedCoeffs[k]*=4*M_PI/numSamples;
			currentVertexSH.shadowedCoeffs[k]*=4*M_PI/numSamples;
		}
		if (j>=numHun)
		{
			cout<<curPer<<"/100\r";
			numHun+=0.01f*(float)numVertices;
			curPer++;
		}
	}


	cout<<"Done"<<endl;

	cout<<"Writing Coffiens....";

	//Save the coefficients to a file
	std::ofstream outFile(object->dFName.getSHDirectCoeffsDN(numSamples), std::ios::out | std::ios::binary | std::ios::trunc);

	//First save the number of bands and samples
	outFile.write((const char *)&numBands, sizeof(int));
	outFile.write((const char *)&numSamples, sizeof(int));

	//Loop through the vertices and save the coefficients for each

	for(int j=0; j<numVertices; ++j)
	{
		SHVertex &currentVertexSH=currentObjectSH[j];

		outFile.write(	(const char *)currentVertexSH.unshadowedCoeffs,
			numFunctions*sizeof(double));
		outFile.write(	(const char *)currentVertexSH.shadowedCoeffs,
			numFunctions*sizeof(double));
	}


	cout<<"Done"<<endl;
	outFile.close();
}


///////////////////////////////////////////////////////////////////////
//                          SHLighting                               //
///////////////////////////////////////////////////////////////////////

//SHLighting::SHLighting(vector<TriMesh *> &ms, int numB = 4, int sqrtNumSamples = 50)
SHLighting::~SHLighting()
{
	delete[] lightCoeffs;
	delete[] rotatedLightCoeffs;
}

void SHLighting::init()
{
	if (isInit)
		return;


	SampleInfo sampleInfo(environment, cubeSideLength, numBands);

	int numFunctions = numBands * numBands;
	lightCoeffs = new double[numFunctions][3];
	rotatedLightCoeffs = new double[numFunctions][3];


	getSHCoefficients(sampleInfo);

	shDiredctCoeffs.calculateSHDirectCoeffs(sampleInfo);
	isInit = true;
}

void SHLighting::setColor(vector<Color>& color)
{
	init();

	if (needRote)
	{
		RotateSHCoefficients();
		needRote = false;
	}


	const int numFunctions=numBands*numBands;

	int vertexNum = color.size();
	vector<SHVertex>& MeshSH = shDiredctCoeffs.vertexSHCoffs[0];
	for (int i=0; i<vertexNum; i++)
	{
		for (int rgbInd = 0; rgbInd < 3; rgbInd++)
		{
			color[i][rgbInd] = 0.f;
			double* unShSH = MeshSH[i].shadowedCoeffs;
			for (int l = 0; l < numFunctions; ++l)
			{
				color[i][rgbInd] += rotatedLightCoeffs[l][rgbInd] * unShSH[l];
			}
		}
	}



}

double SHLighting::light(double theta, double phi)
{


	vec3 direction = vec3(float(sin(theta)*cos(phi)),
		float(sin(theta)*sin(phi)),
		float(cos(theta)));
	
	/*
	float = 


	float y = 1.0 - acos(direction[1])/M_PI;
	vec2 coordn = normalize(vec2(direction[0], direction[2]));
	coordn[0] = acos(coordn[0])/M_PI;
	coordn[1] = acos(coordn[1])/M_PI;
	float x = coordn[0]>0.5 ? 0.5*coordn[1] : 1.0-0.5*coordn[1];

	gl_FragData[0] = vec4(texture2D(texEvn,vec2(x,y)).xyz, 1.0);       

	*/
	return 1.f;

	return (theta<M_PI/6) ? 1.f : 0.f;
}

void SHLighting::getSHCoefficients(SampleInfo& sampleInfo)
{

	int numFunctions = numBands * numBands;
	vector<SAMPLE> &samples = sampleInfo.getSamples();
	int numSamples = samples.size();

	float gamma = 4 * M_PI / numSamples / 20.f;

	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	for(int i=0; i<numFunctions; ++i)
	{
		lightCoeffs[i][rgbInd] = 0.0;

		for(int j=0; j<numSamples; ++j)
			lightCoeffs[i][rgbInd] += samples[j].color[rgbInd] * samples[j].shValues[i];
//lightCoeffs[i][rgbInd] += light(samples[j].theta, samples[j].phi) * samples[j].shValues[i];

		lightCoeffs[i][rgbInd] *= gamma;
	}

}

void SHLighting::rote(int mouseXInc, int mouseYInc)
{
	theta += 0.1*mouseXInc;
	phi += 0.1*mouseYInc;
	//Clamp theta to [0, 180]
	if(theta<0.0)
		theta+=360.0;

	if(theta>360)
		theta-=360;
	needRote = true;
}

void SHLighting::GetZRotationMatrix(int band, double * entries, double angle)
{
	//Calculate the size of the matrix
	int size=2*band+1;

	//Convert angle to radians
	angle*=M_PI/180.0;

	//Entry index
	int currentEntry=0;

	//Loop through the rows and columns of the matrix
	for(int i=0; i<size; ++i)
	for(int j=0; j<size; ++j, ++currentEntry)
	{
		//Initialise this entry to zero
		entries[currentEntry]=0.0;

		//For the central row (i=(size-1)/2), entry is 1 if j==i, else zero
		if(i==(size-1)/2)
		{
			if(j==i)
				entries[currentEntry]=1.0;

			continue;
		}

		//For i<(size-1)/2, entry is cos if j==i or sin if j==size-i-1
		//The angle used is k*angle where k=(size-1)/2-i
		if(i<(size-1)/2)
		{
			int k=(size-1)/2-i;

			if(j==i)
				entries[currentEntry]=cos(k*angle);

			if(j==size-i-1)
				entries[currentEntry]=sin(k*angle);

			continue;
		}

		//For i>(size-1)/2, entry is cos if j==i or -sin if j==size-i-1
		//The angle used is k*angle where k=i-(size-1)/2
		if(i>(size-1)/2)
		{
			int k=i-(size-1)/2;

			if(j==i)
				entries[currentEntry]=cos(k*angle);

			if(j==size-i-1)
				entries[currentEntry]=-sin(k*angle);

			continue;
		}
	}
	

	return;
}

void SHLighting::GetX90DegreeRotationMatrix(int band, double *entries)
{
	//Ensure that 0<=band<=3
	if(band>3)
	{
		cout<<"X rotation matrices are only known for bands 0-3"<<endl;
		return;
	}

	if(band==0)
	{
		entries[0]= 1.0;
	}

	if(band==1)
	{
		entries[0]= 0.0;
		entries[1]= 1.0;
		entries[2]= 0.0;
		entries[3]=-1.0;
		entries[4]= 0.0;
		entries[5]= 0.0;
		entries[6]= 0.0;
		entries[7]= 0.0;
		entries[8]= 1.0;
	}

	if(band==2)
	{
		entries[ 0]= 0.0;
		entries[ 1]= 0.0;
		entries[ 2]= 0.0;
		entries[ 3]= 1.0;
		entries[ 4]= 0.0;
		entries[ 5]= 0.0;
		entries[ 6]=-1.0;
		entries[ 7]= 0.0;
		entries[ 8]= 0.0;
		entries[ 9]= 0.0;
		entries[10]= 0.0;
		entries[11]= 0.0;
		entries[12]=-0.5;
		entries[13]= 0.0;
		entries[14]=-sqrt(3.0)/2;
		entries[15]=-1.0;
		entries[16]= 0.0;
		entries[17]= 0.0;
		entries[18]= 0.0;
		entries[19]= 0.0;
		entries[20]= 0.0;
		entries[21]= 0.0;
		entries[22]=-sqrt(3.0)/2;
		entries[23]= 0.0;
		entries[24]= 0.5;
	}

	if(band==3)
	{
		//Initialise all entries to 0
		for(int i=0; i<49; ++i)
			entries[i]=0.0;

		entries[ 3]=-sqrt(0.625);
		entries[ 5]= sqrt(0.375);

		entries[ 8]=-1.0;

		entries[17]=-sqrt(0.375);
		entries[19]=-sqrt(0.625);

		entries[21]= sqrt(0.625);
		entries[23]= sqrt(0.375);

		entries[32]=-0.25;
		entries[34]=-sqrt(15.0)/4;

		entries[35]=-sqrt(0.375);
		entries[37]= sqrt(0.625);

		entries[46]=-sqrt(15.0)/4;
		entries[48]= 0.25;
	}

	return;
}

void SHLighting::ApplyMatrix(	int size, double* matrix, bool transpose,
				 double* inVector, double* outVector)
{
	//Loop through entries
	for(int i=0; i<size; ++i)
	{
		//Clear this entry of outVector
		outVector[i]=0.0;

		//Loop through matrix row/column
		for(int j=0; j<size; ++j)
		{
			if(transpose)
				outVector[i]+=matrix[j*size+i]*inVector[j];
			else
				outVector[i]+=matrix[i*size+j]*inVector[j];
		}
	}
}

void SHLighting::RotateSHCoefficients()
{
	int numFunctions=numBands*numBands;

	for(int rgbInd = 0; rgbInd < 3; rgbInd++)
	for(int i=0; i<numFunctions; ++i)
		rotatedLightCoeffs[i][rgbInd] = lightCoeffs[i][rgbInd];

	//Band 0 coefficient is unchanged

	//Rotate band 1 coefficients
	for(int rgbInd = 0; rgbInd < 3; rgbInd++)
	if(numBands>1)
	{
		//Get the rotation matrices for band 1 (want to apply Z1*Xt*Z2*X)
		double band1X[9];
		double band1Z1[9];
		double band1Z2[9];

		GetZRotationMatrix(1, band1Z1, phi);
		GetX90DegreeRotationMatrix(1, band1X);
		GetZRotationMatrix(1, band1Z2, theta);

		//Create space to hold the intermediate results
		double band1A[3], band1B[3], band1C[3];

		//Apply the matrices
		ApplyMatrix(3, band1X, false, &lightCoeffs[1][rgbInd], band1A);
		ApplyMatrix(3, band1Z2, false, band1A, band1B);
		ApplyMatrix(3, band1X, true, band1B, band1C);

		ApplyMatrix(3, band1Z1, false, band1C, &rotatedLightCoeffs[1][rgbInd]);
	}

	//Rotate band 2 coefficients
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	if(numBands>2)
	{
		double band2X[25];
		double band2Z1[25];
		double band2Z2[25];

		GetZRotationMatrix(2, band2Z1, phi);
		GetX90DegreeRotationMatrix(2, band2X);
		GetZRotationMatrix(2, band2Z2, theta);

		//Create space to hold the intermediate results
		double band2A[5], band2B[5], band2C[5];

		//Apply the matrices
		ApplyMatrix(5, band2X, false, &lightCoeffs[4][rgbInd], band2A);
		ApplyMatrix(5, band2Z2, false, band2A, band2B);
		ApplyMatrix(5, band2X, true, band2B, band2C);

		ApplyMatrix(5, band2Z1, false, band2C, &rotatedLightCoeffs[4][rgbInd]);
	}

	//Rotate band 3 coefficients
	for (int rgbInd = 0; rgbInd < 3; rgbInd++)
	if(numBands>3)
	{
		double band3X[49];
		double band3Z1[49];
		double band3Z2[49];

		GetZRotationMatrix(3, band3Z1, phi);
		GetX90DegreeRotationMatrix(3, band3X);
		GetZRotationMatrix(3, band3Z2, theta);

		//Create space to hold the intermediate results
		double band3A[7], band3B[7], band3C[7];

		//Apply the matrices
		ApplyMatrix(7, band3X, false, &lightCoeffs[9][rgbInd], band3A);
		ApplyMatrix(7, band3Z2, false, band3A, band3B);
		ApplyMatrix(7, band3X, true, band3B, band3C);

		ApplyMatrix(7, band3Z1, false, band3C, &rotatedLightCoeffs[9][rgbInd]);
	}
}

void SHLighting::initTexture(unsigned int texCubeMap[], bool isOri)
{
		environment.initTexture(texCubeMap);
	return;
}
