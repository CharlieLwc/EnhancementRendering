// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"



#include "../TriMesh2/TriMesh.h"


#ifdef _DEBUG
#pragma comment(lib, "TriMesh2.lib")
#pragma comment(lib, "lbfgs.lib")
#else
#pragma comment(lib, "../Release/TriMesh2.lib")
#pragma comment(lib, "../Release/lbfgs.lib")
#endif




using namespace trimesh;
using namespace std;


#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers


#ifndef EPSILON
#define EPSILON 0.01f
#endif



// TODO: reference additional headers your program requires here
