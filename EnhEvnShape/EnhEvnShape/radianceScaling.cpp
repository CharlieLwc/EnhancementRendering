#include "stdafx.h"

#include "radianceScaling.h"

#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#define PI 3.141592654f
//calculate time 
clock_t cstart, finish; 
using namespace std;

//control the number of environment light samples
int sampleInter = 4;

//just for debug
int mode = 0;
float enhanceScale2 = 0.2f;

//creat sample by divide a octahedral
void EvnLight::getSamplePoints(GLuint texEnvDir, ::std::vector<vec> &directions)
{


	int pointNum = 6;
	int edgeNum = 12;
	int faceNum = 8;

	//calculate the number of points, edges and faces
	for (int iIndex=1; iIndex<interNum; iIndex++)
	{
		pointNum += edgeNum;
		faceNum *=4;
		edgeNum *=4;
	}

	directions.reserve(pointNum);
	::std::vector<ivec2> edges;
	edges.reserve(edgeNum);
	::std::vector<ivec3> faceEdges;
	edges.reserve(faceNum);

	edges.reserve(faceNum);

//inite location of the points of the octahedral
	float a = 1 / sqrt(2.0f);


	directions.push_back(point(0.0, 1.0, 0.0));
	directions.push_back(point(0.0, -1.0, 0.0));
	directions.push_back(point(-a, 0.0, -a));
	directions.push_back(point( -a, 0.0, a));
	directions.push_back(point( a, 0.0,  a));
	directions.push_back(point(a, 0.0,  -a));
	
//the indexes of points on edges 

	edges.push_back(ivec2(0, 3));
	edges.push_back(ivec2(3, 4));
	edges.push_back(ivec2(0, 4));
	edges.push_back(ivec2(0, 5));
	edges.push_back(ivec2(4, 5));
	edges.push_back(ivec2(0, 2));
	edges.push_back(ivec2(2, 5));
	edges.push_back(ivec2(2, 3));
	edges.push_back(ivec2(1, 3));
	edges.push_back(ivec2(1, 4));
	edges.push_back(ivec2(1, 5));
	edges.push_back(ivec2(1, 2));
	
	//p1-e2
	//p2-e0
	//p0-e1

	//the indexes of edges on triangle faces 
	faceEdges.push_back(ivec3(0, 1, 2  ));
	faceEdges.push_back(ivec3(2, 4, 3  ));
	faceEdges.push_back(ivec3(3, 6, 5  ));
	faceEdges.push_back(ivec3(0, 5, 7  ));
	faceEdges.push_back(ivec3(1, 8, 9  ));
	faceEdges.push_back(ivec3(4, 9, 10 ));
	faceEdges.push_back(ivec3(6, 10, 11));
	faceEdges.push_back(ivec3(7, 11, 8 ));

	//the indexes of points on triangle faces 

	facePoints.push_back(ivec3(0, 3, 4));
	facePoints.push_back(ivec3(0, 4, 5));
	facePoints.push_back(ivec3(0, 5, 2));
	facePoints.push_back(ivec3(3, 0, 2));
	facePoints.push_back(ivec3(4, 3, 1));
	facePoints.push_back(ivec3(5, 4, 1));
	facePoints.push_back(ivec3(2, 5, 1));
	facePoints.push_back(ivec3(3, 2, 1));
	
	pointNum = 6;
	edgeNum = 12;
	faceNum = 8;

	//divide one face into 4 faces
	//by the middle points of the edges
	for (int iIndex=1; iIndex<interNum; iIndex++)
	{
		for (int eIndex=0; eIndex<edgeNum; eIndex++)
		{

			if (edges[eIndex][1] <= edges[eIndex][0])
				cout << endl;

			//find middle points
			//normalize point location
			directions.push_back(directions[edges[eIndex][0]] + directions[edges[eIndex][1]]);
			int middlePointDirInd = directions.size()-1;
			normalize(directions[middlePointDirInd]);

			//divide edges into two
			edges.push_back(ivec2(edges[eIndex][1],middlePointDirInd));
			edges[eIndex][1] = middlePointDirInd;


			if (edges[eIndex][1] <= edges[eIndex][0])
				cout << endl;

		}

		for (int fIndex=0; fIndex<faceNum; fIndex++)
		{

			if (fIndex == 32)
				fIndex = 32;
			//creat center face
			facePoints.push_back(ivec3(edges[faceEdges[fIndex][0]][1], edges[faceEdges[fIndex][1]][1], edges[faceEdges[fIndex][2]][1]));
			int centerFaceIndex = facePoints.size()-1;



			if (edges[faceEdges[fIndex][0]][1] < edges[faceEdges[fIndex][1]][1])
				edges.push_back(ivec2(edges[faceEdges[fIndex][0]][1], edges[faceEdges[fIndex][1]][1]));
			else
				edges.push_back(ivec2(edges[faceEdges[fIndex][1]][1], edges[faceEdges[fIndex][0]][1]));


			if (edges[faceEdges[fIndex][1]][1] < edges[faceEdges[fIndex][2]][1])
				edges.push_back(ivec2(edges[faceEdges[fIndex][1]][1], edges[faceEdges[fIndex][2]][1]));
			else
				edges.push_back(ivec2(edges[faceEdges[fIndex][2]][1], edges[faceEdges[fIndex][1]][1]));


			if (edges[faceEdges[fIndex][0]][1] < edges[faceEdges[fIndex][2]][1])

				edges.push_back(ivec2(edges[faceEdges[fIndex][0]][1], edges[faceEdges[fIndex][2]][1]));
			else
				edges.push_back(ivec2(edges[faceEdges[fIndex][2]][1], edges[faceEdges[fIndex][0]][1]));


			int lastEdgeIndex = edges.size()-1;

			faceEdges.push_back(ivec3(lastEdgeIndex-2, lastEdgeIndex-1, lastEdgeIndex));
			
			int edgeStart[3] = {faceEdges[fIndex][0], faceEdges[fIndex][1], faceEdges[fIndex][2]};
			int edgeEnd[3] = {faceEdges[fIndex][0], faceEdges[fIndex][1], faceEdges[fIndex][2]};

			facePoints[fIndex][0]<facePoints[fIndex][1] ? edgeEnd[0] += edgeNum : edgeStart[0] += edgeNum;
			facePoints[fIndex][0]<facePoints[fIndex][2] ? edgeEnd[2] += edgeNum : edgeStart[2] += edgeNum;
			facePoints[fIndex][1]<facePoints[fIndex][2] ? edgeEnd[1] += edgeNum : edgeStart[1] += edgeNum;


			//creat face1
			faceEdges.push_back(ivec3(edgeEnd[0], edgeStart[1], faceEdges[centerFaceIndex][0]));
			facePoints.push_back(ivec3(facePoints[centerFaceIndex][0], facePoints[fIndex][1], facePoints[centerFaceIndex][1]));


			//creat face2
			faceEdges.push_back(ivec3(edgeEnd[1], edgeEnd[2], faceEdges[centerFaceIndex][1]));
			facePoints.push_back(ivec3(facePoints[centerFaceIndex][1], facePoints[fIndex][2], facePoints[centerFaceIndex][2]));


			//change face to face3
			faceEdges[fIndex] = ivec3(edgeStart[2], edgeStart[0], faceEdges[centerFaceIndex][2]);
			facePoints[fIndex] = ivec3(facePoints[centerFaceIndex][2], facePoints[fIndex][0], facePoints[centerFaceIndex][0]);
			
		}


		for (int fIndex = 0; fIndex<faceNum; fIndex++)
		for (int d = 0; d < 3; d++)
		if (facePoints[fIndex][d] != edges[faceEdges[fIndex][d]][1] && facePoints[fIndex][d] != edges[faceEdges[fIndex][d]][0])
			cout << endl;


		pointNum += edgeNum;
		faceNum *=4;
		edgeNum *=4;

	}
	directionNum = pointNum;


	int sampleScale  = pow(2,sampleInter);
	glBindTexture(GL_TEXTURE_2D,texEnvDir);  
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F_ARB,sampleScale,sampleScale,0,GL_RGB,GL_FLOAT,directions[0]); 
}

void EvnLight::readBmp(GLuint texEvn, char* fileName, ::std::vector<Color> &PixelDataf)
{

	static GLint PixelLength;//像素数据长度;
	static GLubyte* PixelData; //像素数据;
	static GLuint texName;

	FILE* pFile;
	fopen_s(&pFile, fileName, "rb");//打开BMP图像文件
	if (pFile == 0)
		exit(0);//如果文件不存在则退出;
	fseek(pFile, 0x0012, SEEK_SET);//找到文件记录长宽的位置;
	fread(&ImageWidth, sizeof(ImageWidth), 1, pFile);//读取宽度;
	fread(&ImageHeight, sizeof(ImageHeight), 1, pFile);//读取高度;
	PixelLength = ImageWidth * 3; //读取像素数据长度;
	while (PixelLength % 4 != 0)
		++PixelLength;//读取像素数据长度;
	PixelLength *= ImageHeight;//读取像素数据长度;
	PixelData = (GLubyte*) malloc(PixelLength); //读取像素数据;

	int pixelNum = ImageHeight*ImageWidth;
	PixelDataf.reserve(pixelNum);

	if (PixelData == 0)
		exit(0); //读取像素数据;
	fseek(pFile, 54, SEEK_SET); //读取像素数据;
	fread(PixelData, PixelLength, 1, pFile); //读取像素数据;
	fclose(pFile);//关闭文件;
	

	int j = 0;
	for (int i=0; i<pixelNum; i++)
	{
		PixelDataf.push_back(Color((GLfloat)PixelData[j+2]/255.0f, (GLfloat)PixelData[j+1]/255.0f, (GLfloat)PixelData[j]/255.0f));
		j++;
		j++;
		j++;
	}

	glBindTexture(GL_TEXTURE_2D,texEvn);  
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F_ARB,ImageWidth,ImageHeight,0,GL_RGB,GL_FLOAT,PixelDataf[0]); 

	delete[] PixelData;

}

void EvnLight::getSampleColor(GLuint texEnvCol, ::std::vector<Color> &PixelDataf, ::std::vector<vec> &directions)
{
	//the color of sample;
	::std::vector<Color> colors;

	colors.reserve(directionNum);
	//find the color of the sample on the environment light
	for (int dIndex=0; dIndex<directionNum; dIndex++)
	{
		float y = 1.0f-acos(directions[dIndex][1])/PI;	

	//	y = abs(0.5f-y)+0.5f;

		float length = sqrt(directions[dIndex][0]*directions[dIndex][0] + directions[dIndex][2]*directions[dIndex][2]);
		if (length != 0.f)
			length = 1.0f/length;

		float x = acos(directions[dIndex][0]*length)/PI>0.5 ? 0.5f*acos(directions[dIndex][2]*length)/PI : 1.0f-0.5f*acos(directions[dIndex][2]*length)/PI;
		while (x>= 1.0)
			x -= 1.0;

		int ss = (int)(y*(float)ImageHeight)*ImageWidth + (int)(x*(float)ImageWidth);
		if (ss >= ImageHeight*ImageWidth)
			ss -= ImageHeight*ImageWidth;
		float sum = PixelDataf[ss][0] + PixelDataf[ss][1] + PixelDataf[ss][2];
		if (sum > 2.4 || sum<0.6)
		{
			if (sum <= 0.0)
				colors.push_back(Color(0.2));
			else if (sum > 2.4)
			{
				float scale = 2.4f/sum;
				colors.push_back(Color(PixelDataf[ss][0] * scale, PixelDataf[ss][1] * scale, PixelDataf[ss][2] * scale));
			}
			else if (sum < 0.6)
			{
				float scale = 0.6/sum;
				colors.push_back(Color(PixelDataf[ss][0] * scale, PixelDataf[ss][1] * scale, PixelDataf[ss][2] * scale));
			}
		}
		else
		{
			colors.push_back(Color(PixelDataf[ss]));

		}

	}
	int sampleScale  = pow(2,sampleInter);

	glBindTexture(GL_TEXTURE_2D,texEnvCol);  
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F_ARB,sampleScale,sampleScale,0,GL_RGB,GL_FLOAT,colors[0]); 
	
}



void EvnLight:: inite(GLuint texEvn, GLuint texEnvCol, GLuint texEnvDir, char* filename, int in)
{
	if (inited)
		return;


	interNum = in;
	getSamplePoints(texEnvDir, directions);
	readBmp(texEvn, filename, PixelDataf);
	getSampleColor(texEnvCol, PixelDataf, directions);
	inited = true;
	return;
}



ShadingProgram:: ~ShadingProgram(void)
{

	
	//glDetachShader(programId, vsId);
	//glDetachShader(programId, fsId);
	//glDeleteShader(vsId);
	//glDeleteShader(fsId);
	//glDeleteProgram(programId);
	
	
}

void ShadingProgram:: inite(char* v, char* f)
{

	if (inited)
		return;


	std::ifstream shaderFileV(v);
	std::ifstream shaderFileF(f);

	std::stringstream shaderDataV;
	std::stringstream shaderDataF;

	shaderDataV << shaderFileV.rdbuf();  //Loads the entire string into a string stream.
	shaderDataF << shaderFileF.rdbuf();  //Loads the entire string into a string stream.

	shaderFileV.close();
	shaderFileF.close();

	const std::string &shaderStringV = shaderDataV.str(); //Get the string stream as a std::string.
	const std::string &shaderStringF = shaderDataF.str(); //Get the string stream as a std::string.



	/*
	ifstream ifileV(v);
	ifstream ifileF(f);

	string temp, tempV, tempF;
	string word;
	//read vertex shader
	while(!ifileV.eof())
	{
		getline(ifileV,temp);

		int location = temp.find("//");
		int location2 = temp.find("//");
		if (location <0 )
			tempV.append(temp,0,location2);
		else if (location2 <0 )
			tempV.append(temp,0,location);
		else
		{
			location = location<location2?location:location2;
			tempV.append(temp,0,location);
		}

	}
	//read fregment shader;
	while(!ifileF.eof())
	{
		getline(ifileF,temp);
		int location = temp.find("//");
		int location2 = temp.find("//");
		if (location <0 )
			tempF.append(temp,0,location2);
		else if (location2 <0 )
			tempF.append(temp,0,location);
		else
		{
			location = location<location2?location:location2;
			tempF.append(temp,0,location);
		}
	}

	const GLchar *bufferV = tempV.c_str();
	const GLchar *bufferf = tempF.c_str();

	*/


	GLint vertCompiled, fragCompiled; // status values
	// Create a vertex shader object and a fragment shader object
	vsId = glCreateShader(GL_VERTEX_SHADER);
	fsId = glCreateShader(GL_FRAGMENT_SHADER);


	const char *strShaderVarV = shaderStringV.c_str();
	const char *strShaderVarF = shaderStringF.c_str();
	GLint iShaderLenV = shaderStringV.size();
	GLint iShaderLenF = shaderStringF.size();
	glShaderSource(vsId, 1, (const GLchar**)&strShaderVarV, (GLint*)&iShaderLenV);
	glShaderSource(fsId, 1, (const GLchar**)&strShaderVarF, (GLint*)&iShaderLenF);


	// Load source code strings into shaders;
//	glShaderSource(vsId, 1, &bufferV, NULL);
//	glShaderSource(fsId, 1, &bufferf, NULL);
	// Compile the brick vertex shader and print out
	// the compiler log file.
	glCompileShader(vsId);

	//	printOpenGLError(); // Check for OpenGL errors
	glGetShaderiv(vsId, GL_COMPILE_STATUS, &vertCompiled);
	printShaderInfoLog(vsId);
	// Compile the brick fragment shader and print out
	// the compiler log file.
	glCompileShader(fsId);
	//	printOpenGLError(); // Check for OpenGL errors
	glGetShaderiv(fsId, GL_COMPILE_STATUS, &fragCompiled);
	printShaderInfoLog(fsId);
	if (!vertCompiled || !fragCompiled)
		cout<<"shader compiled error";
	// Create a program object and attach the two compiled shaders;
	programId = glCreateProgram();
	glAttachShader(programId, vsId);
	glAttachShader(programId, fsId);

	GLint linked;
	// Link the program object and print out the info log
	glLinkProgram(programId);

	//	printOpenGLError(); // Check for OpenGL errors
	glGetProgramiv(programId, GL_LINK_STATUS, &linked);
	printProgramInfoLog(programId);
	if (!linked)
		cout<<"shader linked error";

	inited = true;
	return;

}

void ShadingProgram::printShaderInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
		if (infoLog[0] != '\0')
			printf("%s\n",infoLog);
		free(infoLog);
	}
} 

void ShadingProgram::printProgramInfoLog(GLuint obj)
{
	//use to shading on texture
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
		if (infoLog[0] != '\0')
			printf("%s\n",infoLog);
		free(infoLog);
	}
}

GLint RadianceScaling::getUniLoc(GLuint program, const GLchar *name)
{
	GLint loc;
	loc = glGetUniformLocation(program, name);
	if (loc == -1)
		printf("No such uniform named \"%s\"\n", name);
	return loc;
}

void initeTexture(int num, GLuint* textId, GLuint weidth, GLuint Height)
{
	//generate texture 
	//space
	glGenTextures(num,textId); 
	for (int i = 0; i < num; i++)
	{
		glBindTexture(GL_TEXTURE_2D,textId[i]);  
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		//glTexImage2D(GL_TEXTURE_2D,0, GL_FLOAT_RGBA32_NV ,weidth,Height,0,GL_RGBA,GL_FLOAT,NULL); 
		glTexImage2D(GL_TEXTURE_2D,0, GL_RGBA32F_ARB ,weidth,Height,0,GL_RGBA,GL_FLOAT,NULL); 
	}
}

void initeLodTexture(int level, GLuint* textId, GLuint weidth, GLuint Height)
{
	//generate texture 
	//space
	glGenTextures(1,textId); 
	glBindTexture(GL_TEXTURE_2D,textId[0]);  
	glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, level-1);

	for (int i = 0; i < level; i++)
	{
		glTexImage2D(GL_TEXTURE_2D,i, GL_RGBA32F_ARB ,weidth,Height,0,GL_RGBA,GL_FLOAT,NULL);
		weidth *=0.5;
		Height *=0.5;
	}
}

void RadianceScaling::initeFrameBuffer(void)
{
	normalTexture = new GLfloat[viewPort[2]*viewPort[3]*4];

	glEnable(GL_FRAMEBUFFER_EXT); 

	glGenFramebuffersEXT(1,&farameBuffers);  
	glGenRenderbuffersEXT(1, &depthbuffer);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depthbuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, viewPort[2],viewPort[3]);
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depthbuffer);

	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);  

	layerNum = 2;
	texSmoothedScreenNormal = new GLuint[layerNum+1];
	initeTexture(layerNum, texSmoothedScreenNormal, viewPort[2], viewPort[3]);


	int sampleScale  = pow(2,sampleInter);
	//inite textures
	initeTexture(1, &texScreenNormal, viewPort[2], viewPort[3]);
	initeTexture(1, &texCurvAndDepth, viewPort[2], viewPort[3]);
	initeTexture(1, &texCurvDir, viewPort[2], viewPort[3]);
	initeTexture(1, &texColor, viewPort[2], viewPort[3]);
	initeTexture(1, &texScreenCurv, viewPort[2], viewPort[3]);
	initeTexture(1, &texRidgeLines, viewPort[2], viewPort[3]);
	initeTexture(1, &texSilhouette, viewPort[2], viewPort[3]);
	initeTexture(1, &texResult, viewPort[2], viewPort[3]);
	
	

	initeLodTexture(5, &texBackScreenNormal,viewPort[2],viewPort[3]);
	initeTexture(1, &texWordNormal, viewPort[2], viewPort[3]);
	initeTexture(1, &texWordView, viewPort[2], viewPort[3]);
	initeTexture(1, &texWordPos, viewPort[2], viewPort[3]);

	initeTexture(1, &texView,viewPort[2],viewPort[3]);
	initeLodTexture(5, &texBackView,viewPort[2],viewPort[3]);
	initeTexture(1, &texEvn,1,1);
	initeTexture(1, &texTempEvn,viewPort[2],viewPort[3]);
	initeTexture(1, &texMapedEvn,viewPort[2],viewPort[3]);
	initeTexture(1, &texSum,viewPort[2],viewPort[3]);
	initeTexture(1, &texTempSum,viewPort[2],viewPort[3]);

	initeTexture(1, &texResultEvn, viewPort[2], viewPort[3]);

//	initeTexture(1, &texResultEvn, sampleScale, sampleScale);


	initeTexture(1, &texEnvCol,1,1);
	initeTexture(1, &texEnvDir, 1, 1);
	initeTexture(1, &texAeraLightDir, 1, 1);
	
	initeLodTexture(5, &texTSM,viewPort[2],viewPort[3]);
	

	
}

void ComputeBackgroundQuad(float g_bgVertices[][3], GLint* viewPort)
{  
	GLint viewport[4];  
	GLdouble projMatrix[16];  
	GLdouble modelMatrix[16];  


	modelMatrix[0] = modelMatrix[5] = modelMatrix[10] = modelMatrix[15] = 1.0;
	modelMatrix[1] = modelMatrix[2] = modelMatrix[3] = modelMatrix[4] =
		modelMatrix[6] = modelMatrix[7] = modelMatrix[8] = modelMatrix[9] =
		modelMatrix[11] = modelMatrix[12] = modelMatrix[13] = modelMatrix[14] = 0;


	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);  
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);  

	//   
	double objX, objY, objZ;  
	const double winZ = 0.5;  // 0 / 1  
	gluUnProject(0, 0, winZ, modelMatrix, projMatrix, viewport, &objX, &objY, &objZ);  
	g_bgVertices[0][0] = objX;  
	g_bgVertices[0][1] = objY;  
	g_bgVertices[0][2] = objZ;  

	gluUnProject(viewPort[2], 0, winZ, modelMatrix, projMatrix, viewport, &objX, &objY, &objZ);
	g_bgVertices[1][0] = objX;  
	g_bgVertices[1][1] = objY;  
	g_bgVertices[1][2] = objZ;  

	gluUnProject(viewPort[2], viewPort[3], winZ, modelMatrix, projMatrix, viewport, &objX, &objY, &objZ);
	g_bgVertices[2][0] = objX;  
	g_bgVertices[2][1] = objY;  
	g_bgVertices[2][2] = objZ;  

	gluUnProject(0, viewPort[3], winZ, modelMatrix, projMatrix, viewport, &objX, &objY, &objZ);
	g_bgVertices[3][0] = objX;  
	g_bgVertices[3][1] = objY; 
	g_bgVertices[3][2] = objZ;  
}  

void RadianceScaling::inite(TriMesh* m)
{

	GLenum err=glewInit();//初始化glew;
	if(err!=GLEW_OK)
		cout<<"failed to init glew!";

	if (glewIsSupported("GL_VERSION_2_0"))
		printf("Ready for OpenGL 2.0\n");
	else
	{
		printf("OpenGL 2.0 not supported\n");
		exit(1);
	}

	glGetIntegerv(GL_VIEWPORT, viewPort);
	mesh = m;
	initeFrameBuffer();
	evnlight.inite(texEvn, texEnvCol, texEnvDir, "..\\Image\\doge2.bmp", sampleInter);
}

void RadianceScaling::getAdjacentLightDir(vec3 lightDir, int sizeX, int sizeY, vec3 xDirection, float resolution, bool isPointLight, TriMesh::BSphere &global_bsph)
{

	::std::vector<Color> directions;



	vec3 yDirection = lightDir.cross(xDirection);
	normalize(yDirection);
	xDirection = lightDir.cross(yDirection);
	normalize(xDirection);




	int totalX = sizeX * 2 + 1;
	int totalY = sizeY * 2 + 1;


	float startAngleX = -resolution * (float)sizeX;

	for (int x = 0; x < totalX; x++)
	{
		vec3 tempXDir = lightDir * cos(startAngleX) + xDirection * sin(startAngleX);

		normalize(tempXDir);


		float startAngleY = -resolution * (float)sizeY;


		for (int y = 0; y < totalY; y++)
		{
			vec3 tempDirection = tempXDir * cos(startAngleY) + yDirection * sin(startAngleY);
			normalize(tempDirection);
			if (isPointLight)
			{
				tempDirection = global_bsph.center + tempDirection *global_bsph.r;
			}
			directions.push_back(tempDirection);

			startAngleY += resolution;
		}
		startAngleX += resolution;
	}

	glBindTexture(GL_TEXTURE_2D, texAeraLightDir);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, totalX, totalY, 0, GL_RGB, GL_FLOAT, directions[0]);

}


void RadianceScaling::filtEvn(float field_of_view)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	evnFilter.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\evn_filter.fs");
	glUseProgram(evnFilter.getID());

	//get view port size
	GLint V[4]; 
	glGetIntegerv(GL_VIEWPORT, V);  
	int width = V[2], height = V[3];
	float diag = sqrt(float(sqr(width) + sqr(height)));
	//get pixel scale on texture
	float pixScale = field_of_view/diag ;
	float ySpan = -(float) height * pixScale;
	float xSpan = -(float) width * pixScale;
	
	// shading the model using samples on the environment light

		//add shading result of some sample on the final result
		//the shading result is also the input
		//we use two texturesz
		//one is input other is result
		//switch them every turn

		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texMapedEvn,0); 

		GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
		glDrawBuffers(1, buffers);

		glClearColor(0.0, 0.0, 0.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D,texView);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D,texScreenNormal);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D,texEvn);
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_2D,texEnvCol);
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D,texEnvDir);
		glUniform1i(getUniLoc(evnFilter.getID(), "view"), 0);
		glUniform1i(getUniLoc(evnFilter.getID(), "screenNormal"), 1);
		glUniform1i(getUniLoc(evnFilter.getID(), "texEvn"), 2);
		glUniform1i(getUniLoc(evnFilter.getID(), "col"), 3);
		glUniform1i(getUniLoc(evnFilter.getID(), "dir"), 4);

		glUniform1f(getUniLoc(evnFilter.getID(), "increasement"), 1.0/pow(2,sampleInter));
		glUniform1i(getUniLoc(evnFilter.getID(), "bufferSize"), evnlight.getDrectionNum()-2);
		glUniform1f(getUniLoc(evnFilter.getID(), "envsize"), (float)evnlight.getDrectionNum());
		glUniform1f(getUniLoc(evnFilter.getID(), "xSpan"), xSpan);
		glUniform1f(getUniLoc(evnFilter.getID(), "ySpan"), ySpan);



		float g_bgVertices[4][3];
		ComputeBackgroundQuad(g_bgVertices, viewPort);

		glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0);
		glVertex3fv(g_bgVertices[0]);
		glTexCoord2f(1.0, 0.0);
		glVertex3fv(g_bgVertices[1]);
		glTexCoord2f(1.0, 1.0);
		glVertex3fv(g_bgVertices[2]);
		glTexCoord2f(0.0, 1.0);
		glVertex3fv(g_bgVertices[3]);
		glEnd();



}

void RadianceScaling::expectDiffuse(float field_of_view)
{


	float mainDirection[3];
	PCATex(texScreenNormal, field_of_view, mainDirection);


	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	evnExpectation.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\evn_expectDiffuse.fs");
	glUseProgram(evnExpectation.getID());

	GLint V[4]; 
	glGetIntegerv(GL_VIEWPORT, V);  
	int width = V[2], height = V[3];
	float diag = sqrt(float(sqr(width) + sqr(height)));

	float pixScale = field_of_view/diag ;
	float ySpan = -(float) height * pixScale;
	float xSpan = -(float) width * pixScale;



	// shading the model using samples on the environment light
		//add shading result of some sample on the final result
		//the shading result is also the input
		//we use two textures
		//one is input other is result
		//switch them every turn

		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texMapedEvn,0); 

		GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
		glDrawBuffers(1, buffers);

		glClearColor(0.0, 0.0, 0.0, 0.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D,texScreenNormal);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_2D,texEvn);
		glUniform1i(getUniLoc(evnExpectation.getID(), "screenNormal"), 1);
		glUniform1i(getUniLoc(evnExpectation.getID(), "texEvn"), 2);

		glUniform1f(getUniLoc(evnExpectation.getID(), "xSpan"), xSpan);
		glUniform1f(getUniLoc(evnExpectation.getID(), "ySpan"), ySpan);
		glUniform3f(getUniLoc(evnExpectation.getID(), "mainDirection"),mainDirection[0], mainDirection[1], mainDirection[2]);

		float g_bgVertices[4][3];
		ComputeBackgroundQuad(g_bgVertices, viewPort);

		glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0);
		glVertex3fv(g_bgVertices[0]);
		glTexCoord2f(1.0, 0.0);
		glVertex3fv(g_bgVertices[1]);
		glTexCoord2f(1.0, 1.0);
		glVertex3fv(g_bgVertices[2]);
		glTexCoord2f(0.0, 1.0);
		glVertex3fv(g_bgVertices[3]);
		glEnd();

	

}

void RadianceScaling::expectGlossy(float field_of_view)
{

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	evnExpectation.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\evn_expectGlossy.fs");
	glUseProgram(evnExpectation.getID());


	GLint V[4]; 
	glGetIntegerv(GL_VIEWPORT, V);  
	int width = V[2], height = V[3];
	float diag = sqrt(float(sqr(width) + sqr(height)));

	float pixScale = field_of_view/diag ;
	float ySpan = -(float) height * pixScale;
	float xSpan = -(float) width * pixScale;



	// shading the model using samples on the environment light
	//add shading result of some sample on the final result
	//the shading result is also the input
	//we use two textures
	//one is input other is result
	//switch them every turn
	GLuint temp = texTempEvn;
	texTempEvn = texMapedEvn;
	texMapedEvn = temp;

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texMapedEvn,0); 

	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D,texScreenNormal);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D,texEvn);
	glUniform1i(getUniLoc(evnExpectation.getID(), "screenNormal"), 1);
	glUniform1i(getUniLoc(evnExpectation.getID(), "texEvn"), 2);

	glUniform1f(getUniLoc(evnExpectation.getID(), "xSpan"), xSpan);
	glUniform1f(getUniLoc(evnExpectation.getID(), "ySpan"), ySpan);

	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();



}

void RadianceScaling::getNormal(bool needWorldNormal, bool needView, bool needWorldView, int* layerIndex, bool needCurv, bool needDepth, 
	bool needFeature, bool needBack, bool needWorldPos, bool needCurvDir, bool needColor, 
	bool needLightDirection, vec3 lightDirection, TriMesh::BSphere &global_bsph,
	GLint* sn, int zmin, int zmax, void(*draw_tstrips)(const TriMesh *, bool), bool useFaceNormal)
{



	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	int increasement = 0;

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texScreenNormal, 0);
	if (needWorldNormal)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texWordNormal, 0);
	if (needView)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texView, 0);
	if (needWorldView)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texWordView, 0);
	if (needWorldPos)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texWordPos, 0);
	if (layerIndex[0] >= 0)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texSmoothedScreenNormal[0], 0);
	if (layerIndex[1] >= 0)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texSmoothedScreenNormal[1], 0);
	if (needCurv || needDepth || needFeature)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texCurvAndDepth, 0);
	if (needCurvDir)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texCurvDir, 0);
	if (needColor)
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT + increasement++, GL_TEXTURE_2D, texColor, 0);
	
	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT,
		GL_COLOR_ATTACHMENT2_EXT, GL_COLOR_ATTACHMENT3_EXT,
		GL_COLOR_ATTACHMENT4_EXT, GL_COLOR_ATTACHMENT5_EXT,
		GL_COLOR_ATTACHMENT6_EXT, GL_COLOR_ATTACHMENT7_EXT,
		GL_COLOR_ATTACHMENT8_EXT, GL_COLOR_ATTACHMENT9_EXT
	};

	glDrawBuffers(increasement, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	evnNormal.inite("..\\RS_shader_ml\\evn_normal.vs", "..\\RS_shader_ml\\evn_normal.fs");
	glUseProgram(evnNormal.getID());

	glUniform1i(getUniLoc(evnNormal.getID(), "needWorldNormal"), needWorldNormal);
	glUniform1i(getUniLoc(evnNormal.getID(), "needView"), needView);
	glUniform1i(getUniLoc(evnNormal.getID(), "needWorldView"), needWorldView);
	glUniform1i(getUniLoc(evnNormal.getID(), "needWorldPos"), needWorldPos);
	glUniform1i(getUniLoc(evnNormal.getID(), "needLayer0"), (layerIndex[0] >= 0));
	glUniform1i(getUniLoc(evnNormal.getID(), "needLayer1"), (layerIndex[1] >= 0));
	glUniform1i(getUniLoc(evnNormal.getID(), "needCurv"), needCurv);
	glUniform1i(getUniLoc(evnNormal.getID(), "needDepth"), needDepth);
	glUniform1i(getUniLoc(evnNormal.getID(), "needFeature"), needFeature);
	glUniform1i(getUniLoc(evnNormal.getID(), "needCurvDir"), needCurvDir);
	glUniform1i(getUniLoc(evnNormal.getID(), "needColor"), needColor);

	glUniform3f(getUniLoc(evnNormal.getID(), "campos"), mesh->eyePosition[0], mesh->eyePosition[1], mesh->eyePosition[2]);
	glUniform1f(getUniLoc(evnNormal.getID(), "zmin"), zmin);
	glUniform1f(getUniLoc(evnNormal.getID(), "zmax"), zmax);



	for (int lIndex = 0; lIndex<layerNum; lIndex++)
	{
		if (layerIndex[lIndex] >= 0)
		{
			char cIndex[30];
			sprintf_s(cIndex, 30, "smoothedNormal%d", lIndex);
			sn[lIndex] = glGetAttribLocation(evnNormal.getID(), cIndex);
			glEnableVertexAttribArray(sn[lIndex]);
			mesh->smoothNormal(layerIndex[lIndex], mesh->smoothRadius[layerIndex[lIndex]], true, false);
			glVertexAttribPointer(sn[lIndex], 3, GL_FLOAT, GL_TRUE, 0, &mesh->smoothedNormal[layerIndex[lIndex]][0][0]);
		}
	}



	if (needCurv)
	{
		sn[2] = glGetAttribLocation(evnNormal.getID(), "weightedCurv");
		glEnableVertexAttribArray(sn[2]);
		glVertexAttribPointer(sn[2], 1, GL_FLOAT, GL_TRUE, 0, &mesh->weightedCurvature[0]);
	}
	if (needFeature)
	{
		sn[3] = glGetAttribLocation(evnNormal.getID(), "weightedFeature");
		glEnableVertexAttribArray(sn[3]);
		glVertexAttribPointer(sn[3], 1, GL_FLOAT, GL_TRUE, 0, &mesh->halfEdge.reliefFeature[0]);
	}


	if (needCurvDir)
	{
		sn[4] = glGetAttribLocation(evnNormal.getID(), "curDir");
		glEnableVertexAttribArray(sn[4]);
		glVertexAttribPointer(sn[4], 3, GL_FLOAT, GL_TRUE, 0, &mesh->pdir1[0][0]);
	}


	if (needLightDirection)
	{
		glBegin(GL_LINES);
		glVertex3fv(global_bsph.center);
		glVertex3fv(global_bsph.center + lightDirection * global_bsph.r);
		glEnd();
	}

	draw_tstrips(mesh, useFaceNormal);
	glFinish();

	for (int lIndex = 0; lIndex < layerNum; lIndex++)
		glDisableVertexAttribArray(sn[lIndex]);
	if (needCurv)
		glDisableVertexAttribArray(sn[2]);
	if (needFeature)
		glDisableVertexAttribArray(sn[3]);
	if (needCurvDir)
		glDisableVertexAttribArray(sn[4]);

	
}


void RadianceScaling::anlyNormal(vec3 glossyDirection, vec3& resultNormal, vec3& resultView)
{

	//calculate the sum of pixels on a texture
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	analyNormal.inite("..\\RS_shader_ml\\evn_null.vs",
		"..\\RS_shader_ml\\AnalyNormal.fs");
	glUseProgram(analyNormal.getID());





	bool start = true;

	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	int xShift = V[2], yShift = V[3];
	for (int phase = 0; phase <2; phase++)
	{
		//get view port size
		int end = 1;


		while (xShift>end || yShift>end)
		{
			//calculate the sum of all normal.
			//xyz is the sum
			//w is the total number

			GLuint temp = texTempSum;
			texTempSum = texSum;
			texSum = temp;

			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texSum, 0);

			GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
			glDrawBuffers(1, buffers);


			glClearColor(0.0, 0.0, 0.0, 0.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			//first loop

			if (start)
			{
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texWordNormal);
				glUniform1i(getUniLoc(analyNormal.getID(), "normal"), 0);
				glActiveTexture(GL_TEXTURE1);
				glBindTexture(GL_TEXTURE_2D, texWordView);
				glUniform1i(getUniLoc(analyNormal.getID(), "view"), 1);
			}
			else
			{
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texTempSum);
				glUniform1i(getUniLoc(analyNormal.getID(), "tempResult"), 0);
			}

			if (xShift > 1)
				xShift = (xShift + 1)*0.5;

			if (yShift > 1)
				yShift = (yShift + 1)*0.5;

			glUniform1f(getUniLoc(analyNormal.getID(), "xShift"), (float)xShift / (float)V[2]);
			glUniform1f(getUniLoc(analyNormal.getID(), "yShift"), (float)yShift / (float)V[3]);
			glUniform1i(getUniLoc(analyNormal.getID(), "phase"), phase);
			glUniform3f(getUniLoc(analyNormal.getID(), "glossyDirection"), glossyDirection[0], glossyDirection[1], glossyDirection[2]);

			float g_bgVertices[4][3];
			ComputeBackgroundQuad(g_bgVertices, viewPort);

			glBegin(GL_QUADS);
			glTexCoord2f(0.0, 0.0);
			glVertex3fv(g_bgVertices[0]);
			glTexCoord2f(1.0, 0.0);
			glVertex3fv(g_bgVertices[1]);
			glTexCoord2f(1.0, 1.0);
			glVertex3fv(g_bgVertices[2]);
			glTexCoord2f(0.0, 1.0);
			glVertex3fv(g_bgVertices[3]);
			glEnd();

			glFinish();
			
			if (start)
				start = false;

			if (phase == 0)
				break;
		}

		if (phase == 1)
		{
			vec4 averageNormal;
			vec4 averageView;
			glReadPixels(0, 0, 1, 1, GL_RGBA, GL_FLOAT, averageNormal);

			glReadPixels(511, 0, 1, 1, GL_RGB, GL_FLOAT, averageView);
			resultNormal = vec3(averageNormal[0], averageNormal[1], averageNormal[2]);
			resultView = vec3(averageView[0], averageView[1], averageView[2]);
			normalize(resultNormal);
			normalize(resultView);
		}
	}


}


void RadianceScaling::getRidgeInfo(void(*draw_tstrips)(const TriMesh *, const bool), void(*DrawRidgeLines)(void), bool useFaceNormal)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texRidgeLines, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	evnRidgeLines.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\evn_ridge.fs");
	glUseProgram(evnRidgeLines.getID());

	glPolygonOffset(1.f, 1.f);
	glEnable(GL_POLYGON_OFFSET_FILL);
	draw_tstrips(mesh, useFaceNormal);
	glDisable(GL_POLYGON_OFFSET_FILL);

	DrawRidgeLines();

	glFinish();


	glUseProgram(0);

}

void RadianceScaling::calSilhouette(float silThre, float silHard)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texSilhouette, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	evnSilhouette.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\evn_silhouette.fs");
	glUseProgram(evnSilhouette.getID());



	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texCurvAndDepth);
	glUniform1i(getUniLoc(evnSilhouette.getID(), "depth"), 0);
	
	glUniform1f(getUniLoc(evnSilhouette.getID(), "sw"), 1.f / viewPort[2]);
	glUniform1f(getUniLoc(evnSilhouette.getID(), "sh"), 1.f / viewPort[3]);

	glUniform1f(getUniLoc(evnSilhouette.getID(), "threshold"), silThre);
	glUniform1f(getUniLoc(evnSilhouette.getID(), "hardness"), silHard);

	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();


	glFinish();
	glUseProgram(0);

}

void RadianceScaling::addEdgeGlossy(int glossySize, float shiness, float enhanceScale)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texResult, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	evnEdgeGlossy.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\evn_edgeGlossy.fs");
	glUseProgram(evnEdgeGlossy.getID());

	if (enhanceScale < 0)
		enhanceScale = -1.0 / enhanceScale;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texSilhouette);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texRidgeLines);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, texScreenNormal);

	glUniform1i(getUniLoc(evnEdgeGlossy.getID(), "silhouette"), 0);
	glUniform1i(getUniLoc(evnEdgeGlossy.getID(), "ridgeLines"), 1);
	glUniform1i(getUniLoc(evnEdgeGlossy.getID(), "normal"), 2);


	glUniform1f(getUniLoc(evnEdgeGlossy.getID(), "shineness"), shiness);
	glUniform1f(getUniLoc(evnEdgeGlossy.getID(), "gamma"), enhanceScale);
	
	glUniform1i(getUniLoc(evnEdgeGlossy.getID(), "glossySize"), glossySize);
	glUniform1f(getUniLoc(evnEdgeGlossy.getID(), "sw"), 1.f / viewPort[2]);
	glUniform1f(getUniLoc(evnEdgeGlossy.getID(), "sh"), 1.f / viewPort[3]);

	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();


	glFinish();
	glUseProgram(0);

}

void RadianceScaling::renderMaterial(int materialType, bool isPointLight, vec3 lightSource, float roughness, float metalness, float edginess, float backscatter, 
	float fresnel, float alphaX, float alphaY, vec3 IosDir, float transparency, 
	int lightSizeX, int lightSizeY,
	float diffuseValue, float specularValue, vec3 diffColor, vec3 specColor, bool needColor, Color averageColor, bool costumColor)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texResult, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	showMaterial.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\showMaterial.fs");
	glUseProgram(showMaterial.getID());

	if (enhanceScale < 0)
		enhanceScale = -1.0 / enhanceScale;

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texWordNormal);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texWordView);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, texWordPos);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, texAeraLightDir);



	glUniform1i(getUniLoc(showMaterial.getID(), "worldNormal"), 0);
	glUniform1i(getUniLoc(showMaterial.getID(), "worldView"), 1);
	glUniform1i(getUniLoc(showMaterial.getID(), "worldPos"), 2);
	glUniform1i(getUniLoc(showMaterial.getID(), "areaLight"), 3);

	if (needColor)
	{
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_2D, texColor);
		glUniform1i(getUniLoc(showMaterial.getID(), "modelColor"), 4);
		glUniform1i(getUniLoc(showMaterial.getID(), "needColor"), needColor);
	}


	lightSizeX = lightSizeX * 2 + 1;
	lightSizeY = lightSizeY * 2 + 1;
	glUniform1f(getUniLoc(showMaterial.getID(), "xIncrement"), 1.f / (float)lightSizeX);
	glUniform1f(getUniLoc(showMaterial.getID(), "yIncrement"), 1.f / (float)lightSizeY);


	glUniform1i(getUniLoc(showMaterial.getID(), "isPointLight"), isPointLight);
	glUniform1i(getUniLoc(showMaterial.getID(), "materialType"), materialType);

//	glUniform3f(getUniLoc(showMaterial.getID(), "lightSource"), lightSource[0], lightSource[1], lightSource[2]);


	glUniform1f(getUniLoc(showMaterial.getID(), "roughness"), roughness);
	glUniform1f(getUniLoc(showMaterial.getID(), "fresnel"), fresnel);
	
	glUniform1f(getUniLoc(showMaterial.getID(), "diffuseValue"), diffuseValue);
	glUniform1f(getUniLoc(showMaterial.getID(), "specularValue"), specularValue);

	if (costumColor)
	{
		glUniform3f(getUniLoc(showMaterial.getID(), "diffColor"), diffColor[0], diffColor[1], diffColor[2]);
		glUniform3f(getUniLoc(showMaterial.getID(), "specColor"), specColor[0], specColor[1], specColor[2]);
	}
	else
	{
		glUniform3f(getUniLoc(showMaterial.getID(), "diffColor"), averageColor[0], averageColor[1], averageColor[2]);
		glUniform3f(getUniLoc(showMaterial.getID(), "specColor"), averageColor[0], averageColor[1], averageColor[2]);
	}

	glUniform1f(getUniLoc(showMaterial.getID(), "matelness"), metalness);
	glUniform1f(getUniLoc(showMaterial.getID(), "edginess"), edginess);
	glUniform1f(getUniLoc(showMaterial.getID(), "backscatter"), backscatter);
	
	glUniform1f(getUniLoc(showMaterial.getID(), "fTransparency"), transparency);

	glUniform1f(getUniLoc(showMaterial.getID(), "alphaX"), alphaX);
	glUniform1f(getUniLoc(showMaterial.getID(), "alphaY"), alphaY);
	glUniform3f(getUniLoc(showMaterial.getID(), "iosDir"), IosDir[0], IosDir[1], IosDir[2]);
	
	


	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();


	glFinish();
	glUseProgram(0);

}


void RadianceScaling::calScreenCurv(float* weight, float enhanceScale)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texScreenCurv, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);



	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	texSmoothedScreenNormal[2] = texScreenNormal;

	evnScreenCurv.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\evn_screenCurve.fs");
	glUseProgram(evnScreenCurv.getID());




	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texSilhouette);
	glUniform1i(getUniLoc(evnScreenCurv.getID(), "silhouette"), 0);

	
	char cIndex[20];
	int tempLayerNum = 0;
	for (int lIndex = 0; lIndex<3; lIndex++)
	{
		if (weight[lIndex] > 0.f)
		{
			glActiveTexture(GL_TEXTURE2 + tempLayerNum);
			glBindTexture(GL_TEXTURE_2D, texSmoothedScreenNormal[lIndex]);

			sprintf_s(cIndex, 20, "norms[%d]", lIndex);
			glUniform1i(getUniLoc(evnScreenCurv.getID(), cIndex), tempLayerNum+2);
			tempLayerNum++;
		}
		sprintf_s(cIndex, 20, "weights[%d]", lIndex);
		glUniform1f(getUniLoc(evnScreenCurv.getID(), cIndex), weight[lIndex]);
	}
	float scale = enhanceScale > 0.f ? enhanceScale : -1.0 / enhanceScale;

	glUniform1f(getUniLoc(evnScreenCurv.getID(), "sw"), 1.f / viewPort[2]);
	glUniform1f(getUniLoc(evnScreenCurv.getID(), "sh"), 1.f / viewPort[3]);
	glUniform1f(getUniLoc(evnScreenCurv.getID(), "enhanceScale"), scale);
	
	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();


	glFinish();
	glUseProgram(0);

}

void RadianceScaling::EvnRS(float field_of_view, bool useObjectCurv, bool movingObject, float enhanceScale, float diffuseValue, float specularValue, float shiness, vec3 diffColor, vec3 specColor, float highlightTher)
{






	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	evnFilter.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\evn_filter.fs");
	glUseProgram(evnFilter.getID());

	//get view port size
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	int width = V[2], height = V[3];
	float diag = sqrt(float(sqr(width) + sqr(height)));
	//get pixel scale on texture
	float pixScale = field_of_view / diag;
	float ySpan = -(float)height * pixScale;
	float xSpan = -(float)width * pixScale;

	// shading the model using samples on the environment light

	//add shading result of some sample on the final result
	//the shading result is also the input
	//we use two texturesz
	//one is input other is result
	//switch them every turn

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texResultEvn, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	if (movingObject)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texView);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, texScreenNormal);
	}
	else
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texWordView);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, texWordNormal);
	}
	



	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, texEvn);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, texEnvCol);
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D, texEnvDir);

	if (enhanceScale > 0.f)
	{
		glActiveTexture(GL_TEXTURE5);
		if (useObjectCurv)
			glBindTexture(GL_TEXTURE_2D, texCurvAndDepth);
		else
			glBindTexture(GL_TEXTURE_2D, texScreenCurv);
		glUniform1i(getUniLoc(evnFilter.getID(), "curvature"), 5);
	}


	glUniform1i(getUniLoc(evnFilter.getID(), "view"), 0);
	glUniform1i(getUniLoc(evnFilter.getID(), "normal"), 1);
	glUniform1i(getUniLoc(evnFilter.getID(), "environmentMap"), 2);
	glUniform1i(getUniLoc(evnFilter.getID(), "evnSampleCol"), 3);
	glUniform1i(getUniLoc(evnFilter.getID(), "evnSampleDir"), 4);

	glUniform1f(getUniLoc(evnFilter.getID(), "increasement"), 1.0 / pow(2, sampleInter));
	glUniform1i(getUniLoc(evnFilter.getID(), "bufferSize"), evnlight.getDrectionNum() - 2);
	glUniform1f(getUniLoc(evnFilter.getID(), "envsize"), (float)evnlight.getDrectionNum());
	glUniform1f(getUniLoc(evnFilter.getID(), "xSpan"), xSpan);
	glUniform1f(getUniLoc(evnFilter.getID(), "ySpan"), ySpan);
	glUniform1f(getUniLoc(evnFilter.getID(), "enhanceScale"), enhanceScale);


	glUniform1f(getUniLoc(evnFilter.getID(), "diffuseValue"), diffuseValue);
	glUniform1f(getUniLoc(evnFilter.getID(), "specularValue"), specularValue);
	glUniform1f(getUniLoc(evnFilter.getID(), "shiness"), shiness);
	glUniform1f(getUniLoc(evnFilter.getID(), "highlightTher"), highlightTher);
	glUniform3f(getUniLoc(evnFilter.getID(), "diffColor"), diffColor[0], diffColor[1], diffColor[2]);
	glUniform3f(getUniLoc(evnFilter.getID(), "specColor"), specColor[0], specColor[1], specColor[2]);
	



	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();


	glFinish();
	glUseProgram(0);





}

void RadianceScaling::showText(GLuint texId)
{

	//show the shading result in texture
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	showTexture.inite("..\\RS_shader_ml\\evn_null.vs", "..\\RS_shader_ml\\showTex.fs");
	glUseProgram(showTexture.getID());

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texId);

	glUniform1i(getUniLoc(showTexture.getID(), "tex"), 0);

	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);


	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();



	glFinish();

	glUseProgram(0);

}

void RadianceScaling::showResult(bool GLSL_showNormal, bool GLSL_showCurv, bool GLSL_useObjectCurv, bool GLSL_radianceScaling, bool GLSL_enhancePRT, bool GLSL_silhouette, bool showMaterial, bool GLSL_ExaggeratedShading)
{
	
	/*
	showText(texWordNormal);
	return;
	*/

	if (GLSL_radianceScaling)
		showText(texResultEvn);
	else if (GLSL_showNormal)
	{
		showText(texWordNormal);
	}
	else if (GLSL_showCurv)
	{
	if (GLSL_useObjectCurv)
		showText(texCurvAndDepth);
	else
		showText(texScreenCurv);
	}
	else if (GLSL_silhouette)
		showText(texSilhouette);
	else if (GLSL_enhancePRT || showMaterial || GLSL_ExaggeratedShading)
		showText(texResult);
	

}

void RadianceScaling::getFrontNormal(float curv_thre_pos, float curv_thre_neg, bool needDepth)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texScreenNormal,0); 
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT1_EXT,GL_TEXTURE_2D,texView,0); 
	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
	glDrawBuffers(2, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	buffProgramEvn.inite("..\\RS_shader_ml\\evn_buffer.vs", "..\\RS_shader_ml\\evn_buffer.fs");
	glUseProgram(buffProgramEvn.getID());

	glUniform1f(getUniLoc(buffProgramEvn.getID(), "curv_thre_pos"), curv_thre_pos);
	glUniform1f(getUniLoc(buffProgramEvn.getID(), "curv_thre_neg"), curv_thre_neg);
	glUniform1i(getUniLoc(buffProgramEvn.getID(), "needDepth"), needDepth);
}

void RadianceScaling::getBackNormal(float curv_thre_pos, float curv_thre_neg, bool needDepth)
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texBackScreenNormal,0); 
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT1_EXT,GL_TEXTURE_2D,texBackView,0); 
	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
	glDrawBuffers(2, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	buffProgramEvn.inite("..\\RS_shader_ml\\evn_buffer.vs", "..\\RS_shader_ml\\evn_buffer.fs");
	glUseProgram(buffProgramEvn.getID());
	glUniform1f(getUniLoc(buffProgramEvn.getID(), "curv_thre_pos"), curv_thre_pos);
	glUniform1f(getUniLoc(buffProgramEvn.getID(), "curv_thre_neg"), curv_thre_neg);
	glUniform1i(getUniLoc(buffProgramEvn.getID(), "needDepth"), needDepth);
}

void RadianceScaling::transparent(float field_of_view)
{

	float ni = 1.0;
	float nt = 1.15;


	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	refraction.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\refraction.fs");
	glUseProgram(refraction.getID());
	
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texMapedEvn,0); 

	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D,texEvn);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D,texScreenNormal);

	glUniform1i(getUniLoc(refraction.getID(), "texEnv"), 0);
	glUniform1i(getUniLoc(refraction.getID(), "screenNormal"), 1);

	// local2 = { n_i/n_t, (n_i/n_t)^2, <unused>, <unused> }
	//    n_i = index of refraction of external material (air)
	//    n_t = index of refraction of internal material
	glUniform2f(getUniLoc(refraction.getID(), "local2"), ni/nt, ni*ni/(nt*nt));


	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();



}

void RadianceScaling::translucent(float field_of_view, float m_refindex, float m_sigma_s[], float m_sigma_a[])
{


	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	subScattering.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\translucent.fs");
	glUseProgram(subScattering.getID());

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texMapedEvn,0); 

	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D,texEvn);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D,texScreenNormal);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D,texBackScreenNormal);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D,texView);
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D,texBackView);
	glActiveTexture(GL_TEXTURE5);
	glBindTexture(GL_TEXTURE_2D,texTSM);

	glUniform1i(getUniLoc(subScattering.getID(), "texEnv"), 0);
	glUniform1i(getUniLoc(subScattering.getID(), "screenNormal"), 1);
	glUniform1i(getUniLoc(subScattering.getID(), "backScreenNormal"), 2);
	glUniform1i(getUniLoc(subScattering.getID(), "texView"), 3);
	glUniform1i(getUniLoc(subScattering.getID(), "texBackView"), 4);
	glUniform1i(getUniLoc(subScattering.getID(), "texTSM"), 5);

	float Fdr, A, D;
	float sigma_t[3], sigma_tr[3], zr3[3], zv3[3], Rd0[3];

	Fdr = -(1.440f/(m_refindex*m_refindex)) + (0.710f/m_refindex) + 0.668f + (0.0636f*m_refindex);
	A = (1.0f+Fdr)/(1.0f-Fdr);



	for (int d=0; d<3; d++)
	{
		sigma_t[d] = m_sigma_s[d] + m_sigma_a[d];
		D = 1.0f / (3.0f * sigma_t[d]);
		sigma_tr[d] = sqrt(3.0f*m_sigma_a[d]*sigma_t[d]);
		zr3[d] = 1.0/sigma_t[d];
		zv3[d] = zr3[d] + (4.0f*A*D);
		Rd0[d] = m_sigma_s[d] / (sigma_t[d] * 4.0f * 3.141592653589f);
	}
	
	glUniform3fv(getUniLoc(subScattering.getID(), "zr3"), 1, zr3);
	glUniform3fv(getUniLoc(subScattering.getID(), "zv3"), 1, zv3);
	glUniform3fv(getUniLoc(subScattering.getID(), "sigma_t"), 1, sigma_t);
	glUniform3fv(getUniLoc(subScattering.getID(), "sigma_tr"), 1, sigma_tr);
	glUniform3fv(getUniLoc(subScattering.getID(), "Rd0"), 1, Rd0);



	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();



}

void RadianceScaling::createTSM(float field_of_view, float enhanceScale)
{

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	tralTSM.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\TSM.fs");
	glUseProgram(tralTSM.getID());

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texTSM,0); 

	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
	glDrawBuffers(1, buffers);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D,texEvn);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D,texBackScreenNormal);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D,texView);

	glUniform1i(getUniLoc(tralTSM.getID(), "texEnv"), 0);
	glUniform1i(getUniLoc(tralTSM.getID(), "backScreenNormal"), 1);
	glUniform1i(getUniLoc(tralTSM.getID(), "texView"), 2);
	glUniform1f(getUniLoc(tralTSM.getID(), "enhanceScale"), enhanceScale);
	

	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();



}

void RadianceScaling::showResultEnv(GLuint texId)
{

	//show the shading result in texture
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);  

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	showTexResultEnv.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\showResultEnv.fs");
	glUseProgram(showTexResultEnv.getID());


	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D,texEnvCol);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D,texEnvDir);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D,texScreenNormal);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D,texView);
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_2D,texId);

	float a = 1.0;
	glUniform1i(getUniLoc(showTexResultEnv.getID(), "col"), 0);
	glUniform1i(getUniLoc(showTexResultEnv.getID(), "dir"), 1);
	glUniform1i(getUniLoc(showTexResultEnv.getID(), "screenNormal"), 2);
	glUniform1i(getUniLoc(showTexResultEnv.getID(), "view"), 3);
	glUniform1i(getUniLoc(showTexResultEnv.getID(), "scale"), 4);
	glUniform1f(getUniLoc(showTexResultEnv.getID(), "enhance"), a);

	glUniform1f(getUniLoc(showTexResultEnv.getID(), "increasement"), 1.0/pow(2,sampleInter));
	glUniform1i(getUniLoc(showTexResultEnv.getID(), "bufferSize"), evnlight.getDrectionNum()-2);




	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);


	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();




	glUseProgram(0);

}

int Eigs(float a1[][3], double d[], double v[][3], int& nrot)
{
	double a[3][3];
	int n=3;

	double b[3], z[3], g, s, c, h, tau, sss, ddd, t, tresh, sm, theta;
	int ip, iq, i, j;
	for (ip=0; ip<n; ip++)
	{
		for (iq=0; iq<n; iq++)
		{
			v[ip][iq]=0.0;
			a[ip][iq] = a1[ip][iq];
		}
		v[ip][ip]=1.0;
	}
	for (ip=0; ip<n; ip++)
	{
		b[ip]=a[ip][ip];
		d[ip]=b[ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=0; i<50; i++)
	{
		sm=0.0;
		for (ip=0; ip<n-1; ip++)
		{
			for (iq=ip+1; iq<n; iq++)
			{
				sm=sm+fabs(a[ip][iq]);
			}
		}
		if (sm==0.0)
		{
			return 0;
		}
		if (i<3)
		{
			tresh=0.2*sm/double(n*n);
		}
		else
		{
			tresh=0.0;
		}
		for (ip=0; ip<n-1; ip++)
		{
			for (iq=ip+1; iq<n; iq++)
			{
				g=100.0*fabs(a[ip][iq]);
				sss=fabs(d[ip])+g;
				ddd=fabs(d[iq])+g;
				if ((i>3)&&(sss==fabs(d[ip]))&&(ddd==fabs(d[iq])))
				{
					a[ip][iq]=0.0;
				}
				else
				{
					if (fabs(a[ip][iq])>tresh)
					{
						h=d[iq]-d[ip];
						if ((fabs(h)+g)==fabs(h))
						{
							t=a[ip][iq]/h;
						}
						else
						{
							theta=0.5*h/a[ip][iq];
							t=1.0/(fabs(theta)+sqrt(1.0+pow(theta, 2)));
							if (theta<0.0)
								t=-t;
						}
						c=1.0/sqrt(1.0+t*t);
						s=t*c;
						tau=s/(1.0+c);
						h=t*a[ip][iq];
						z[ip]=z[ip]-h;
						z[iq]=z[iq]+h;
						d[ip]=d[ip]-h;
						d[iq]=d[iq]+h;
						a[ip][iq]=0.0;
						for (j=0; j<=ip-1; j++)
						{
							g=a[j][ip];
							h=a[j][iq];
							a[j][ip]=g-s*(h+g*tau);
							a[j][iq]=h+s*(g-h*tau);
						}
						for (j=ip+1; j<=iq-1; j++)
						{
							g=a[ip][j];
							h=a[j][iq];
							a[ip][j]=g-s*(h+g*tau);
							a[j][iq]=h+s*(g-h*tau);
						}
						for (j=iq+1; j<n; j++)
						{
							g=a[ip][j];
							h=a[iq][j];
							a[ip][j]=g-s*(h+g*tau);
							a[iq][j]=h+s*(g-h*tau);
						}
						for (j=0; j<n; j++)
						{
							g=v[j][ip];
							h=v[j][iq];
							v[j][ip]=g-s*(h+g*tau);
							v[j][iq]=h+s*(g-h*tau);
						}
						nrot=nrot+1;
					}
				}
			}
		}
		for (ip=0; ip<n; ip++)
		{
			b[ip]=b[ip]+z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	return 1;
}

void RadianceScaling::PCATex(GLuint texId, float field_of_view, float mainDirection[])
{

	//calculate the sum of pixels on a texture
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	PCATexture.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\PCATex.fs");
	glUseProgram(PCATexture.getID());


	

	vec4 averageNormal;

	for (int phase=0; phase <3; phase++)
	{
		//get view port size
		GLint V[4]; 
		glGetIntegerv(GL_VIEWPORT, V);  
		int xShift = V[2], yShift = V[3];
		int end = 1;
		bool start = true;

		if (phase >= 1)
			end = 1;

		while(xShift>end || yShift>end)
		{
			//calculate the sum of all normal.
			//xyz is the sum
			//w is the total number

			GLuint temp = texTempSum;
			texTempSum = texSum;
			texSum = temp;

			glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texSum,0); 

			GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
			glDrawBuffers(1, buffers);


			glClearColor(0.0, 0.0, 0.0, 0.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			//first loop

			glActiveTexture(GL_TEXTURE0);
			if (start)
				glBindTexture(GL_TEXTURE_2D,texId);
			else
				glBindTexture(GL_TEXTURE_2D,texTempSum);

			if(xShift > 1)
				xShift=(xShift+1)*0.5;

			if(yShift > 1)
				yShift=(yShift+1)*0.5;

			glUniform1i(getUniLoc(PCATexture.getID(), "tempResult"), 0);
			glUniform1f(getUniLoc(PCATexture.getID(), "xShift"), (float)xShift/(float)V[2]);
			glUniform1f(getUniLoc(PCATexture.getID(), "yShift"), (float)yShift/(float)V[3]);
			glUniform1i(getUniLoc(PCATexture.getID(), "phase"), phase);
			glUniform3f(getUniLoc(PCATexture.getID(), "averageNormal"),averageNormal[0], averageNormal[1], averageNormal[2]);

			float g_bgVertices[4][3];
			ComputeBackgroundQuad(g_bgVertices, viewPort);

			glBegin(GL_QUADS);
			glTexCoord2f(0.0, 0.0);
			glVertex3fv(g_bgVertices[0]);
			glTexCoord2f(1.0, 0.0);
			glVertex3fv(g_bgVertices[1]);
			glTexCoord2f(1.0, 1.0);
			glVertex3fv(g_bgVertices[2]);
			glTexCoord2f(0.0, 1.0);
			glVertex3fv(g_bgVertices[3]);
			glEnd();

			glFinish();

			if (start)
				start = false;

			if (phase == 1)
				phase++;
		}

		if (phase == 0)
		{
			glReadPixels(0, 0, 1, 1, GL_RGBA, GL_FLOAT, averageNormal);
			averageNormal/=averageNormal[3];
		}
	}



	float a[3][3];
	double v[3][3];
	double d[3];
	int jt;

	glReadPixels(0, 0, 1, 1, GL_RGB, GL_FLOAT, a[0]);
	glReadPixels(511, 0, 1, 1, GL_RGB, GL_FLOAT, a[1]);
	glReadPixels(0, 511, 1, 1, GL_RGB, GL_FLOAT, a[2]);

	Eigs(a,d,v,jt);

	int largestD = 0;
	if (d[0]>d[1])
		d[0]>d[2] ? largestD=0 : largestD=2;
	else
		d[1]>d[2] ? largestD=1 : largestD=2;


	for (int i=0; i<3; i++)
	{
		if (v[0][largestD] < 0)
			mainDirection[i] = -v[i][largestD];
		else
			mainDirection[i] = v[i][largestD];
	}
}

void RadianceScaling::remapGlossy()
{

//	int texSize1 = pow(2,sampleInter);
//	
//	glBindTexture(GL_TEXTURE_2D,texResultEvn);
//	float (*pixelss)[3] = new float[texSize1*texSize1][3];
//	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pixelss[0]);
//
//	ifstream ridgeLines_if("D:\\d\\EnhEbnShape\\T.txt");
//
//	for (int i=0; i<texSize1*texSize1; i++)
//	{
//		ridgeLines_if>>pixelss[i][0]>>pixelss[i][1]>>pixelss[i][2];
//	}
//
//	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F_ARB,texSize1,texSize1,0,GL_RGB,GL_FLOAT,pixelss[0]); 
//
//
//	ridgeLines_if.close();
//	return;



	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,farameBuffers);  

	evnRemapGlossy.inite("..\\RS_shader_ml\\evn_null.vs", 
		"..\\RS_shader_ml\\evn_expectGlossy.fs");
	glUseProgram(evnRemapGlossy.getID());


	
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT,GL_COLOR_ATTACHMENT0_EXT,GL_TEXTURE_2D,texTempEvn,0); 

	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT};
	glDrawBuffers(1, buffers);

	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	GLint viewport[4];  
	glGetIntegerv(GL_VIEWPORT, viewport); 
	int texSize = pow(2,sampleInter);
	float evnSize  = (float)texSize/(float)viewport[2];


	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D,texScreenNormal);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D,texView);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D,texEnvDir);

	glUniform1i(getUniLoc(evnRemapGlossy.getID(), "screenNormal"), 0);
	glUniform1i(getUniLoc(evnRemapGlossy.getID(), "view"), 1);
	glUniform1i(getUniLoc(evnRemapGlossy.getID(), "texEnvDir"), 2);
	glUniform1f(getUniLoc(evnRemapGlossy.getID(), "sampleScale"), evnSize);
	glUniform1f(getUniLoc(evnRemapGlossy.getID(), "sampleScaleReverse"), 1.0/evnSize);
	glUniform1f(getUniLoc(evnRemapGlossy.getID(), "increasement"), 1.0/(float)viewport[2]);



	float g_bgVertices[4][3];
	ComputeBackgroundQuad(g_bgVertices, viewPort);

	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3fv(g_bgVertices[0]);
	glTexCoord2f(1.0, 0.0);
	glVertex3fv(g_bgVertices[1]);
	glTexCoord2f(1.0, 1.0);
	glVertex3fv(g_bgVertices[2]);
	glTexCoord2f(0.0, 1.0);
	glVertex3fv(g_bgVertices[3]);
	glEnd();



	glBindTexture(GL_TEXTURE_2D,texResultEvn);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, texSize, texSize);

//	glBindTexture(GL_TEXTURE_2D,texResultEvn);
//	float (*pixelss)[3] = new float[texSize*texSize][3];
//	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pixelss[0]);
//
//	ofstream ridgeLines_if("D:\\d\\EnhEbnShape\\T.txt");
//
//	for (int i=0; i<texSize*texSize; i++)
//	{
//		ridgeLines_if<<pixelss[i][0]<<" "<<pixelss[i][1]<<" "<<pixelss[i][2]<<endl;
//	}
//
//
//	ridgeLines_if.close();
		


}

void RadianceScaling::afterWriteMesh(float field_of_view, bool isGlossy, bool isRefract, bool isTranslucent, float m_refindex, float m_sigma_s[], float m_sigma_a[], float enhanceScale11)
{
	if (isRefract)
	{
		transparent(field_of_view);

		showText(texMapedEvn);
		return;
	}
	isGlossy = false;
	if (isGlossy)
	{
		remapGlossy();
		showResultEnv(texResultEvn);
		//		showResultEnv(texEnvCol);
		//showText(texView);
		return;
	}


	if (isTranslucent)
	{
		createTSM(field_of_view, enhanceScale11);

		glBindTexture(GL_TEXTURE_2D, texBackScreenNormal); 
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texBackView); 
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texTSM); 
		glGenerateMipmap(GL_TEXTURE_2D);

		translucent(field_of_view, m_refindex, m_sigma_s, m_sigma_a);
		showText(texMapedEvn);
		return;
	}



	expectDiffuse(field_of_view);
	showText(texMapedEvn);



	//get normal in screen space

	//calculate curvature in screen space
//	comCurvature();
	//shading
//	filtEvn(field_of_view);


	//show result
//	showResultEnv(texEnvCol);
	

}

void RadianceScaling::exaggeratedshading(int layerNum, float layerWeight[], bool isOur, float enhanceScale, float averageIlluminate, 
	void(*draw_tstrips)(const TriMesh *, const bool), vec3 color, int layerIndex[], vec3 lGlobal, bool useFaceNormal)
{
	
	int vertexNum = mesh[0].vertices.size();



	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, farameBuffers);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, texResult, 0);

	GLenum buffers[] = { GL_COLOR_ATTACHMENT0_EXT };
	glDrawBuffers(1, buffers);


	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (exaggeratedShading.unInited())
		exaggeratedShading.inite("..\\RS_shader_ml\\exaggerated.vs", "..\\RS_shader_ml\\exaggerated.fs");
	glUseProgram(exaggeratedShading.getID());

	GLint sn[11];

	float weightSum = 0.0;
	//weight is proportional to smooth radius

	float weight[12];

	if (!isOur)
	{
		layerNum = 11;

		int mode = 1;
		float edgeL = 0.4;
		for (int lIndex = 0; lIndex<layerNum; lIndex++)
		{
			if (mode == 0)
				weight[lIndex] = 1.0 / edgeL;
			if (mode == 1)
				weight[lIndex] = 1.0;
			if (mode == 2)
				weight[lIndex] = edgeL;

			weightSum += weight[lIndex];


			exSmoothedNormals[lIndex].getNormal(edgeL, mesh, false, lIndex);


			edgeL *= sqrt(2.f);
		}

		if (mode == 0)
			weight[layerNum] = 1.0 / edgeL;
		if (mode == 1)
			weight[layerNum] = 1.0;
		if (mode == 2)
			weight[layerNum] = edgeL;
		weightSum += weight[layerNum];



		weightSum = 1.0 / weightSum;
		for (int lIndex = 0; lIndex <= layerNum; lIndex++)
			weight[lIndex] *= weightSum;


	}
	else
	{
		for (int lIndex = 0; lIndex<layerNum; lIndex++)
		{
			exSmoothedNormals[lIndex].getNormal(1.0, mesh, true, layerIndex[lIndex]);
		}
	}

	if (isOur)
		for (int lIndex = 0; lIndex<=layerNum; lIndex++)
		{
			char cIndex[30];
			sprintf_s(cIndex, 20, "weights[%d]", lIndex);
			glUniform1f(getUniLoc(exaggeratedShading.getID(), cIndex), layerWeight[lIndex]);
		}
	else
		for (int lIndex = 0; lIndex <= layerNum; lIndex++)
		{
			char cIndex[30];
			sprintf_s(cIndex, 20, "weights[%d]", lIndex);
			glUniform1f(getUniLoc(exaggeratedShading.getID(), cIndex), weight[lIndex]);
		}

	glFinish();
	cstart = clock();


	for (int lIndex = 0; lIndex<layerNum; lIndex++)
	{
		char cIndex[30];
		sprintf_s(cIndex, 30, "smoothedNormal%d", lIndex);
		sn[lIndex] = glGetAttribLocation(exaggeratedShading.getID(), cIndex);
		glEnableVertexAttribArray(sn[lIndex]);
		glVertexAttribPointer(sn[lIndex], 3, GL_FLOAT, GL_TRUE, 0, &exSmoothedNormals[lIndex].sommthedNormal[0][0]);
	}



	glUniform1f(getUniLoc(exaggeratedShading.getID(), "enhanceScale"), enhanceScale);
	glUniform1f(getUniLoc(exaggeratedShading.getID(), "averageIlluminate"), averageIlluminate);


	glUniform1i(getUniLoc(exaggeratedShading.getID(), "layerNum"), layerNum);



	glUniform3f(getUniLoc(exaggeratedShading.getID(), "color"), color[0], color[1], color[2]);
	glUniform3f(getUniLoc(exaggeratedShading.getID(), "lGlobal"), lGlobal[0], lGlobal[1], lGlobal[2]);

	

	draw_tstrips(mesh, useFaceNormal);

	glFinish();
	

	for (int lIndex = 0; lIndex<layerNum; lIndex++)
		glDisableVertexAttribArray(sn[lIndex]);


}
