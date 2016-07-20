// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>


#include "../PrecomputedRadiaceTransfer/SHLighting.h"
#include "../PrecomputedRadiaceTransfer/Transport.h"
#include "../TriMesh2/TriMesh.h"


#ifdef _DEBUG

#pragma comment(lib, "PrecomputedRadiaceTransfer.lib")
#pragma comment(lib, "TriMesh2.lib")
//#pragma comment(lib, "color2gray.lib")
#else
#pragma comment(lib, "../Release/PrecomputedRadiaceTransfer.lib")
#pragma comment(lib, "../Release/TriMesh2.lib")
//#pragma comment(lib, "../Release/color2gray.lib")
#endif



using namespace std;
using namespace trimesh;


#define GLEW_STATIC 1


// TODO: reference additional headers your program requires here
