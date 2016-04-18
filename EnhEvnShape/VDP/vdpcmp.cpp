/**
 * @brief main procedure of the HDR VDP, command line interface
 * 
 * This file is a part of HDR VDP.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2003,2004 Rafal Mantiuk and Grzegorz Krawczyk
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ---------------------------------------------------------------------- 
 *
 * @author Rafal Mantiuk, <mantiuk@mpi-inf.mpg.de>
 *
 * $Id: vdpcmp.cpp,v 1.19 2006/04/07 09:42:32 mantiuk Exp $
 */

#include "stdafx.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "array2d_algorithm.h"

#include "pfs.h"
#include "getopt.h"
#include <math.h>

#include "otf.h"
#include "nonlinearity.h"
#include "fftw_array2d.h"
#include "csf_filter.h"
#include "cortex_transform.h"
#include "fftutils.h"

#include "dump_image_aspect.h"


#define PROG_NAME "vdpcmp"

class QuietException 
{
};

#include <assert.h>

void fixPredictionWithRMS( const pfs::Array2D *inTarget, const pfs::Array2D *inMask, pfs::Array2D *probabilityMap );

void newProcessVDP( const pfs::Array2D *target, const pfs::Array2D *mask,
  pfs::Array2D *probabilityMap,
  OTF *otf,
  PhotoreceptorNonlin *nonlinearity, CSFFilter *csfFilter, DetectionMechanism *dm,
  bool fixRMS );

/*
#define PFSEOL "\x0a"
static void dumpPFS( const char *fileName, const pfs::Array2D *data, const char *channelName )
{
  FILE *fh = fopen( fileName, "wb" );
  assert( fh != NULL );

  int width = data->getCols();
  int height = data->getRows();

  fprintf( fh, "PFS1" PFSEOL "%d %d" PFSEOL "1" PFSEOL "0" PFSEOL
    "%s" PFSEOL "0" PFSEOL "ENDH", width, height, channelName );

  for( int y = 0; y < height; y++ )
    for( int x = 0; x < width; x++ ) {
      fwrite( &((*data)(x,y)), sizeof( float ), 1, fh );
    }
  
  fclose( fh );
}
static void dumpPFS( const char *fileName, FFTWComplexArray *data, const char *channelName )
{
  FILE *fh = fopen( fileName, "wb" );
  assert( fh != NULL );

  int width = data->getFFTCols();
  int height = data->getFFTRows();

  fprintf( fh, "PFS1" PFSEOL "%d %d" PFSEOL "1" PFSEOL "0" PFSEOL
    "%s" PFSEOL "0" PFSEOL "ENDH", width, height, channelName );

  const fftwf_complex *rawData = data->getData();
  
  for( int y = 0; y < height; y++ )
    for( int x = 0; x < width; x++ ) {
      float mod = sqrt( rawData[x+y*width][0]*rawData[x+y*width][0] +
        rawData[x+y*width][1]*rawData[x+y*width][1] );
      fwrite( &mod, sizeof( float ), 1, fh );
    }
  
  fclose( fh );
}
*/

static void errorCheck( bool condition, char *string )
{
  if( !condition ) {
    throw pfs::Exception( string );
  }
}

void printHelp()
{
  fprintf( stderr, PROG_NAME " <mask_file> <target_file> [--ldr] [--mask <val>] [--psycho <val>] [--otf <otfID>] [--peak-contrast <contrast>] [--detection-mechanism <dmID>] [--no-phase-uncertainty] [--dump <pattern>] [--display-x-resolution] [--display-y-resolution] [--display-width] [--display-height] [--min-distance] [--max-distance] [--multiply-lum] [--no-abs-fix] [--verbose] [--help]\n"
    "See man page for more information.\n" );
}

void processVDP(int argc, char **argv)
{
  const char *targetFile = NULL, *maskFile = NULL;
  float maskingSlope = 1;
  float psychometricFunctionSlope = 3.5; // The value from the patent
  float multiplyLum = 1;
  const char *detectionMechID = NULL;
  const char *nonlinearityID = NULL;
  const char *otfID = "DEELEY";
  bool hdr = true;
  bool phaseUncertainty = true;
  float minDistance = 0.5, maxDistance = 0.5;
  float displayWidth = 0.375f, displayHeight = 0.300f;
  int displayXResolution = 1280, displayYResolution = 1024;
  float peakContrast = 0.006f;
  bool absFix = true;

  bool verbose = false;

  static struct option cmdLineOptions[] = {
    { (const wchar_t *)"help", no_argument, NULL, 'h' },
    { (const wchar_t *)"verbose", no_argument, NULL, 'v' },
    { (const wchar_t *)"mask", required_argument, NULL, 'm' },
    { (const wchar_t *)"psycho", required_argument, NULL, 'p' },
    { (const wchar_t *)"detection-mechanism", required_argument, NULL, 't' },
    { (const wchar_t *)"otf", required_argument, NULL, 'o' },    
    { (const wchar_t *)"nonlinearity", required_argument, NULL, 'n' }, // TODO
    { (const wchar_t *)"ldr", no_argument, NULL, 'l' }, 
    { (const wchar_t *)"no-phase-uncertainty", no_argument, NULL, 'u' }, 
    { (const wchar_t *)"no-pu", no_argument, NULL, 'u' }, 
    { (const wchar_t *)"dump", required_argument, NULL, 'd' },
    { (const wchar_t *)"display-x-resolution", required_argument, NULL, 'x' },
    { (const wchar_t *)"display-y-resolution", required_argument, NULL, 'y' },
    { (const wchar_t *)"display-width", required_argument, NULL, 'w' },
    { (const wchar_t *)"display-height", required_argument, NULL, 'e' },
    { (const wchar_t *)"min-distance", required_argument, NULL, 'i' },
    { (const wchar_t *)"max-distance", required_argument, NULL, 'a' },
    { (const wchar_t *)"peak-contrast", required_argument, NULL, 'c' },
    { (const wchar_t *)"multiply-lum", required_argument, NULL, '*' },
    { (const wchar_t *)"no-abs-fix", no_argument, NULL, 'R' },
    { NULL, 0, NULL, 0 }
  };

  initializeDumpImageAspect();
  
  int optionIndex = 0;
  while( 1 ) {
    int c = getopt_long (argc, (wchar_t *const *)argv, (wchar_t *)"m:p:", cmdLineOptions, &optionIndex);
    if( c == -1 ) break;
    switch( c ) {
    case 'h':
      printHelp();
      throw QuietException();
    case 'v':
      verbose = true;
      break;
    case 'm':
      maskingSlope = (float)strtod( (const char *)optarg, NULL );
      break;
    case 'p':
		psychometricFunctionSlope = (float)strtod((const char *)optarg, NULL);
      break;
    case 'd':
		dumpImageAspect->addIncludePattern(( char *)optarg);
      break;
    case 't':
		detectionMechID = (const char *)optarg;
      break;
    case 'o':
		otfID = (const char *)optarg;
      break;
    case 'u':
      phaseUncertainty = false;
      break;
    case 'l':
      hdr = false;
      break;
    case 'x':
		displayXResolution = strtol((const char *)optarg, NULL, 10);
      break;
    case 'y':
		displayYResolution = strtol((const char *)optarg, NULL, 10);
      break;
    case 'w':
		displayWidth = (float)strtod((const char *)optarg, NULL);
      break;
    case 'e':
		displayHeight = (float)strtod((const char *)optarg, NULL);
      break;
    case 'i':
		minDistance = (float)strtod((const char *)optarg, NULL);
      break;      
    case 'a':
		maxDistance = (float)strtod((const char *)optarg, NULL);
      break;      
    case 'c':
		peakContrast = (float)strtod((const char *)optarg, NULL);
      break;
    case 'R':
      absFix = false;
      break;
    case '*':
		multiplyLum = (float)strtod((const char *)optarg, NULL);
      break;
    case '?':
      throw QuietException();
    case ':':
      throw QuietException();
    }
  } 

  
//    int oi = optind;
//    while( oi < argc )
//      cerr << "option: " << argv[oi++] << "\n";
  
  // Exactly two non-option arguments: target, and mask
  errorCheck( optind == argc - 2, "Both target and mask files must be specified" );

  targetFile = argv[optind++];
  maskFile = argv[optind++];

  if( verbose ) {
    std::cerr << "Psychometric function slope: " << psychometricFunctionSlope << "\n";
    std::cerr << "Threshold elevation slope: " << maskingSlope << "\n";
  }
  
  FILE *targetFh;
  FILE *maskFh;

  bool targetIsStdin = !strcmp( "-", targetFile );
  bool maskIsStdin = !strcmp( "-", maskFile );

  errorCheck( !(targetIsStdin && maskIsStdin),
    "Only one of the target and mask pair of images can be loaded from standard input" );  
  


  if (targetIsStdin)
	  targetFh = stdin;
  else
	  fopen_s(&targetFh, targetFile, "rb");





  errorCheck( targetFh != NULL, "Can not open target image file" );
  

  if (maskIsStdin)
	  maskFh = stdin;
  else
	  fopen_s(&maskFh, maskFile, "rb");






  errorCheck( maskFh != NULL, "Can not open mask image file" );

  pfs::DOMIO pfsio;
    
  // Load target image (with artifacts)
  pfs::Frame *targetFrame = pfsio.readFrame( targetFh );
  pfs::Array2D *targetY = targetFrame->getChannel( "Y" );
  errorCheck( targetY != NULL, "Luminance channel missing in the target image" );
  clampToPositive( targetY, targetY ); // Remove <=0 and nan(s)
    
  // Load mask image (original)
  pfs::Frame *maskFrame = pfsio.readFrame( maskFh );
  pfs::Array2D *maskY = maskFrame->getChannel( "Y" );
  errorCheck( maskY != NULL, "Luminance channel missing in the mask image" );
  clampToPositive( maskY, maskY ); // Remove <=0 and nan(s)
  
  // Check if sizes match
  errorCheck( maskY->getRows() == targetY->getRows() &&
    maskY->getCols() == targetY->getCols(), "Can compare only images of the same size" );

  // Multiply luminance if necessary
  if( multiplyLum != 1 ) {
    const int pixels = maskY->getRows()*maskY->getCols();
    for( int i = 0; i < pixels; i++ ) {
      (*maskY)(i) *= multiplyLum;
      (*targetY)(i) *= multiplyLum;      
    }
  }
  
  // VDP processing blocks
  OTF *otf = NULL;
  PhotoreceptorNonlin *nonlinearity = NULL;
  CSFFilter *csfFilter = NULL;
  DetectionMechanism *dm = NULL;
  
  {                           // Configure processing blocks of the VDP
    float contrastNormalizationFactor = 0; // Used if nonlinearity == photoreceptor
    ViewingConditions viewCond( displayXResolution, displayYResolution,
      displayWidth, displayHeight, minDistance, maxDistance );

    if( verbose )
      fprintf( stderr, "Cycles per degree: %g\n", viewCond.getPixelsPerDegree( (viewCond.minDistance + viewCond.maxDistance)/2 ) );
    
    
    if( hdr == false ) {        // Use original Daly's VDP
      nonlinearity = new VDPNonlinearity( );
      //Normalization factor used to compute contrast between mask and
      //target Sugested values are display device mean and logAverage(
      //maskY ), but those do not give stable results. This
      //implementation follows the patent and set normalization factor
      //to 1 (no normalization)
      contrastNormalizationFactor = 1; // Should be 
      fprintf( stderr, "conrast norm: %g\n", contrastNormalizationFactor );
      csfFilter = new VDPCSF( maskY->getCols(), maskY->getRows(), viewCond, 30 );
      otfID="NONE";
    }    

	if (!_stricmp(otfID, "SPENCER")) {
      otf = new SpencerOTF( viewCond, maskY->getCols(), maskY->getRows(),
        logAverage( maskY ) );
	}
	else if (!_stricmp(otfID, "NORMANN")) {
      otf = new NormannBaxterOTF( viewCond, maskY->getCols(), maskY->getRows() );
	}
	else if (!_stricmp(otfID, "WESTHEIMER")) {
      otf = new WestheimerOTF( viewCond, maskY->getCols(), maskY->getRows() );
	}
	else if (!_stricmp(otfID, "MARIMONT")) {
      otf = new MarimontOTF( viewCond, maskY->getCols(), maskY->getRows(), logAverage( maskY ) );
	}
	else if (!_stricmp(otfID, "DEELEY")) {
      otf = new DeeleyOTF( viewCond, maskY->getCols(), maskY->getRows(), logAverage( maskY ) );
	}
	else if (!_stricmp(otfID, "NONE")) {
      otf = NULL;
    } else {
      throw pfs::Exception( "Unknown optical transfer function (OTF). Possible values: MARIMONT, DEELEY, SPENCER, NORMANN, WESTHEIMER, NONE (default)" );
    }    
    
    if( nonlinearity == NULL )
      nonlinearity = new JNDScaledNonlinearity( peakContrast );
    if( csfFilter == NULL ) {
      csfFilter = new MultiAdaptationCSF( maskY->getCols(), maskY->getRows(), viewCond );
    }
    if( detectionMechID == NULL || !_stricmp( detectionMechID, "VDP" ) ||
      !_stricmp( detectionMechID, "WATSON" ) ) { // default dm
      dm = new WatsonDM( psychometricFunctionSlope, maskingSlope, contrastNormalizationFactor != 0, contrastNormalizationFactor, phaseUncertainty );
    } else if( !_stricmp( detectionMechID, "SIMPLE" ) ) {
      dm = new SingleChannelDM( psychometricFunctionSlope );
    } else {
      throw pfs::Exception( "Unknown detection mechanism. Possible values: VDP, WATSON, SIMPLE" );
    }
  }
  
  // Prepare output for saving
  const int cols = targetY->getCols(), rows = targetY->getRows();
  pfs::Channel *resultMap = targetFrame->createChannel( "VDP" );
  
  newProcessVDP( targetY, maskY, resultMap, otf, nonlinearity, csfFilter, dm, absFix );

  // Write output
  pfsio.writeFrame( targetFrame, stdout );
  pfsio.freeFrame( targetFrame );
  pfsio.freeFrame( maskFrame );
  
  delete nonlinearity;
  delete csfFilter;
  delete dm;
  
  if( targetFh != stdin ) fclose( targetFh );
  if( maskFh != stdin ) fclose( maskFh );  
  
  
}

int main( int argc, char* argv[] )
{
  try {
    processVDP( argc, argv );
  }
  catch( pfs::Exception ex ) {
    fprintf( stderr, PROG_NAME " error: %s\n", ex.getMessage() );
    return EXIT_FAILURE;
  }        
  catch( QuietException  ex ) {
    return EXIT_FAILURE;
  }        
  return EXIT_SUCCESS;
}

void newProcessVDP( const pfs::Array2D *inTarget, const pfs::Array2D *inMask,
  pfs::Array2D *probabilityMap,
  OTF *otf,
  PhotoreceptorNonlin *nonlinearity, CSFFilter *csfFilter, DetectionMechanism *dm,
  bool absFix )
{
  const int rows = inTarget->getRows(), cols = inTarget->getCols();

  BidomainArray2D target( cols, rows ), mask( cols, rows );

  // Copy images to structures that can handle fft
  pfs::copyArray( inTarget, target.setSpatial() );
  pfs::copyArray( inMask, mask.setSpatial() );
  
  if( otf != NULL ) {
    std::cerr << "OTF: " << otf->getName() << "...\n";
    otf->process( &target, &target );
    otf->process( &mask, &mask );
  }
  
  dumpImageAspect->dump( "target_otf.pfs", &target, "Y" );
  dumpImageAspect->dump( "mask_otf.pfs", &mask, "Y" );

  BidomainArray2D adaptationMap( cols, rows );
  adaptationMap = mask;
  dumpImageAspect->dump( "adaptation_map.pfs", &adaptationMap, "Y" );
  
  std::cerr << "Nonlinearity: " << nonlinearity->getName() << "...\n";
  
  nonlinearity->process( &target, &target );
  nonlinearity->process( &mask, &mask );

  dumpImageAspect->dump( "target_photoreceptor_nonlin.pfs", &target, "Y" );
  dumpImageAspect->dump( "mask_photoreceptor_nonlin.pfs", &mask, "Y" );

//   FFTWComplexArray freqTarget( cols, rows );
//   FFTWComplexArray freqMask( cols, rows );

  std::cerr << "CSF filtering: " << csfFilter->getName() << "...\n";
  
  csfFilter->process( &target, &target, &adaptationMap );
  csfFilter->process( &mask, &mask, &adaptationMap );

  dumpImageAspect->dump( "target_csf.pfs", &target, "Y" );
  dumpImageAspect->dump( "mask_csf.pfs", &mask, "Y" );
  
  dm->process( &target, &mask, probabilityMap );

  // Compare target and mask and zero predictions if pixels are
  // almost identical
  if( absFix ) {
    std::cerr << "Applying ABS fix\n";
    fixPredictionWithRMS( inTarget, inMask, probabilityMap );
  }
  
}

/**
 * Compare target and mask and zero predictions if pixels are almost
 * identical
 */
void fixPredictionWithRMS( const pfs::Array2D *inTarget, const pfs::Array2D *inMask, pfs::Array2D *probabilityMap ) 
{
  const int pixels = inTarget->getCols()*inTarget->getRows();

  for( int i = 0; i < pixels; i++ ) {
    float d = fabsf( (*inTarget)(i) - (*inMask)(i) );
    if( d < 0.0001 ) (*probabilityMap)(i) = 0;
  }
  
}
