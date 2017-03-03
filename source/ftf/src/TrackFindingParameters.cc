#include "TrackFindingParameters.h"
#include "TrackUtil.h"
#include <algorithm>
#include <iostream>
#include <cstring>
using namespace ftf;
using std::cout;
using std::endl;

void TrackFindingParameters::read ( char* inputFile ) 
{

	FILE* dataFile = fopen( inputFile, "r");
	if (dataFile == NULL) {
		printf ( "TrackFindingParameters::write: Error opening input file %s \n", inputFile ) ;
		return ;
	}

	char* name = new char[100] ;
	while ( 1 ) {
		if ( fscanf ( dataFile, "%s", name ) == EOF ) break ;
		if ( !strncmp(name,"infoLevel"       ,8) ) { 
			fscanf ( dataFile, "%d", &infoLevel ) ;
			continue ;
		}
		if ( !strncmp(name,"segmentRowSearch",8) ) {
			fscanf ( dataFile, "%d", &segmentRowSearchRange ) ;
			continue ;
		}
		if ( !strncmp(name,"trackRowSearch",  8) ) {
			fscanf ( dataFile, "%d", &trackRowSearchRange ) ;
			continue ;
		}
		if ( !strncmp(name,"getErrors     ",  8) ) {
			fscanf ( dataFile, "%d", &getErrors      ) ;
			continue ;
		}
		if ( !strncmp(name,"fillTracks    ",  8) ) {
			fscanf ( dataFile, "%d", &fillTracks     ) ;
			continue ;
		}
		if ( !strncmp(name,"ghostFlag     ",  8) ) {
			fscanf ( dataFile, "%d", &ghostFlag      ) ;
			continue ;
		}
		if ( !strncmp(name,"goBackwards   ",  8) ) {
			fscanf ( dataFile, "%d", &goBackwards    ) ;
			continue ;
		}
		if ( !strncmp(name,"mergePrimaries",  8) ) {
			fscanf ( dataFile, "%d", &mergePrimaries ) ;
			continue ;
		}
		if ( !strncmp(name,"minHitsPerTrack", 8) ) {
			fscanf ( dataFile, "%d", &minHitsPerTrack ) ;
			continue ;
		}
		if ( !strncmp(name,"modRow",          6) ) {
			fscanf ( dataFile, "%d", &modRow          ) ;
			continue ;
		}
		if ( !strncmp(name,"nHitsForSegment", 8) ) {
			fscanf ( dataFile, "%d", &nHitsForSegment ) ;
			continue ;
		}
		if ( !strncmp(name,"minHitsForFit",   8) ) {
			fscanf ( dataFile, "%d", &minHitsForFit   ) ;
			continue ;
		}
		if ( !strncmp(name,"nEtaTrack",   8) ) {
			fscanf ( dataFile, "%i", &nEtaTrack     ) ;
			continue;
		}
		if ( !strncmp(name,"nEta",   4) ) {
			fscanf ( dataFile, "%d", &nEta          ) ;
			continue ;
		}
		if ( !strncmp(name,"nPhiTrack",   8) ) {
			fscanf ( dataFile, "%d", &nPhiTrack     ) ;
			continue ;
		}
		if ( !strncmp(name,"nPhi",   4) ) {
			fscanf ( dataFile, "%d", &nPhi          ) ;
			continue ;
		}
		if ( !strncmp(name,"detaMerge    ",   8) ) {
			fscanf ( dataFile, "%le", &detaMerge     ) ;
			continue ;
		}
		if ( !strncmp(name,"deta         ",   4) ) {
			fscanf ( dataFile, "%le", &deta          ) ;
			continue ;
		}  
		if ( !strncmp(name,"dphiMerge    ",   8) ) {
			fscanf ( dataFile, "%le", &dphiMerge     ) ;
			continue ;
		}
		if ( !strncmp(name,"dphi         ",   4) ) {
			fscanf ( dataFile, "%le", &dphi          ) ;
			continue ;
		}  
		if ( !strncmp(name,"etaMinTrack  ",   8) ) {
			fscanf ( dataFile, "%le", &etaMinTrack   ) ;
			continue ;
		}
		if ( !strncmp(name,"etaMin",   6) ) {
			fscanf ( dataFile, "%le", &etaMin        ) ;
			continue ;
		}  
		if ( !strncmp(name,"etaMaxTrack  ",   8) ) {
			fscanf ( dataFile, "%le", &etaMaxTrack   ) ;
			continue ;
		}
		if ( !strncmp(name,"etaMax       ",   6) ) {
			fscanf ( dataFile, "%le", &etaMax        ) ;
			continue ;
		}  
		if ( !strncmp(name,"phiMinTrack  ",   8) ) {
			fscanf ( dataFile, "%le", &phiMinTrack   ) ;
			continue ;
		}
		if ( !strncmp(name,"phiMin       ",   6) ) {
			fscanf ( dataFile, "%le", &phiMin        ) ;
			continue ;
		}  
		if ( !strncmp(name,"phiMaxTrack  ",   8) ) {
			fscanf ( dataFile, "%le", &phiMaxTrack   ) ;
			continue ;
		}
		if ( !strncmp(name,"phiMax       ",   6) ) {
			fscanf ( dataFile, "%le", &phiMax        ) ;
			continue ;
		}  
		if ( !strncmp(name,"phiShift     ",   8) ) {
			fscanf ( dataFile, "%le", &phiShift      ) ;
			continue ;
		}  
		if ( !strncmp(name,"distanceMerge",   8) ) {
			fscanf ( dataFile, "%le", &distanceMerge ) ;
			continue ;
		}  
		if ( !strncmp(name,"nPrimaryPasses",  8) ) {
			fscanf ( dataFile, "%d", &nPrimaryPasses ) ;
			continue ;
		}  
		if ( !strncmp(name,"nSecondary",  8) ) {
			fscanf ( dataFile, "%d", &nSecondaryPasses ) ;
			continue ;
		}  
		if ( !strncmp(name,"vertexConstrainedFit",  8) ) {
			fscanf ( dataFile, "%d", &vertexConstrainedFit ) ;
			continue ;
		}  
		if ( !strncmp(name,"parameterLocation",  8) ) {
			fscanf ( dataFile, "%d", &parameterLocation ) ;
			continue ;
		}  
		if ( !strncmp(name,"rowInnerMost",  8) ) {
			fscanf ( dataFile, "%d", &rowInnerMost ) ;
			continue ;
		}  
		if ( !strncmp(name,"rowOuterMost",  8) ) {
			fscanf ( dataFile, "%d", &rowOuterMost ) ;
			continue ;
		}  
		if ( !strncmp(name,"rowStart",         8) ) {
			fscanf ( dataFile, "%d", &rowStart ) ;
			continue ;
		}  
		if ( !strncmp(name,"rowEnd",           6) ) {
			fscanf ( dataFile, "%d", &rowEnd ) ;
			continue ;
		}  
		if ( !strncmp(name,"szFitFlag",        8) ) {
			fscanf ( dataFile, "%d", &szFitFlag ) ;
			continue ;
		}  
		if ( !strncmp(name,"bField",       6) ) {
			fscanf ( dataFile, "%le", &bField ) ;
			continue ;
		}  
		if ( !strncmp(name,"maxChi2Primary",       12) ) {
			fscanf ( dataFile, "%le", &maxChi2Primary ) ;
			continue ;
		}  
		if ( !strncmp(name,"hitChi2Cut",           8) ) {
			fscanf ( dataFile, "%le", &hitChi2Cut ) ;
			continue ;
		}  
		if ( !strncmp(name,"goodHitChi2",           8) ) {
			fscanf ( dataFile, "%le", &goodHitChi2 ) ;
			continue ;
		}  
		if ( !strncmp(name,"trackChi2Cut",           8) ) {
			fscanf ( dataFile, "%le", &trackChi2Cut ) ;
			continue ;
		} 
		if ( !strncmp(name,"goodDistance",           8) ) {
			fscanf ( dataFile, "%le", &goodDistance ) ;
			continue ;
		} 
		if ( !strncmp(name,"ptMinHelixFit",           8) ) {
			fscanf ( dataFile, "%le", &ptMinHelixFit ) ;
			continue ;
		} 
		if ( !strncmp(name,"maxDistanceSegment",   15) ) {
			fscanf ( dataFile, "%le", &maxDistanceSegment ) ;
			continue ;
		} 
		if ( !strncmp(name,"xyErrorScale",   10) ) {
			fscanf ( dataFile, "%le", &xyErrorScale ) ;
			continue ;
		} 
		if ( !strncmp(name,"szErrorScale",   10) ) {
			fscanf ( dataFile, "%le", &szErrorScale ) ;
			continue ;
		} 
		if ( !strncmp(name,"xVertex",   7) ) {
			fscanf ( dataFile, "%le", &xVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"yVertex",   7) ) {
			fscanf ( dataFile, "%le", &yVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"zVertex",   7) ) {
			fscanf ( dataFile, "%le", &zVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"dxVertex",   8) ) {
			fscanf ( dataFile, "%le", &dxVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"dyVertex",   8) ) {
			fscanf ( dataFile, "%le", &dyVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"xyWeightVertex",   8) ) {
			fscanf ( dataFile, "%le", &xyWeightVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"phiVertex",   8) ) {
			fscanf ( dataFile, "%le", &phiVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"rVertex",   7) ) {
			fscanf ( dataFile, "%le", &rVertex ) ;
			continue ;
		} 
		if ( !strncmp(name,"maxTime",   7) ) {
			fscanf ( dataFile, "%le", &maxTime ) ;
			continue ;
		} 
		printf ( "TrackFindingParameters::read: parameter %s not found \n", name ) ;
		float variable ; 
		fscanf ( dataFile, "%e", &variable ) ;

	}

	fclose ( dataFile ) ;

}

void TrackFindingParameters::write ( char* outputFile ) 
{
	FILE* dataFile = fopen( outputFile, "w");
	if (dataFile == NULL) {
		printf ( "TrackFindingParameters::write: Error opening output file %s \n ", outputFile ) ;
		return ;
	}
	write ( dataFile ) ;
	fclose ( dataFile ) ;
}

void TrackFindingParameters::write ( FILE* dataFile ) 
{

	fprintf ( dataFile, "infoLevel            %10d \n", infoLevel ) ;
	fprintf ( dataFile, "segmentRowSearch     %10d  \n", segmentRowSearchRange ) ;
	fprintf ( dataFile, "trackRowSearch       %10d  \n", trackRowSearchRange ) ;
	fprintf ( dataFile, "getErrors            %10d  \n", getErrors ) ;
	fprintf ( dataFile, "fillTracks           %10d  \n", fillTracks ) ;
	fprintf ( dataFile, "ghostFlag            %10d  \n", ghostFlag ) ;
	fprintf ( dataFile, "goBackwards          %10d  \n", goBackwards ) ;
	fprintf ( dataFile, "mergePrimaries       %10d  \n", mergePrimaries ) ;
	fprintf ( dataFile, "minHitsPerTrack      %10d  \n", minHitsPerTrack ) ;
	fprintf ( dataFile, "modRow               %10d  \n", modRow ) ;
	fprintf ( dataFile, "nHitsForSegment      %10d  \n", nHitsForSegment ) ;
	fprintf ( dataFile, "minHitsForFit        %10d  \n", minHitsForFit   ) ;
	fprintf ( dataFile, "nEta                 %10d  \n", nEta            ) ;
	fprintf ( dataFile, "nPhi                 %10d  \n", nPhi            ) ;
	fprintf ( dataFile, "deta                 %10.2e\n", deta         ) ;
	fprintf ( dataFile, "dphi                 %10.2e\n", dphi         ) ;
	fprintf ( dataFile, "etaMin               %10.2e\n", etaMin        ) ;
	fprintf ( dataFile, "etaMax               %10.2e\n", etaMax        ) ;
	fprintf ( dataFile, "phiMin               %10.2e\n", phiMin        ) ;
	fprintf ( dataFile, "phiMax               %10.2e\n", phiMax        ) ;
	fprintf ( dataFile, "phiShift             %10.2e\n", phiShift   ) ;
	fprintf ( dataFile, "detaMerge            %10.2e\n", detaMerge    ) ;
	fprintf ( dataFile, "distanceMerge        %10.2e\n", distanceMerge ) ;
	fprintf ( dataFile, "nEtaTrack            %10d  \n", nEtaTrack       ) ;
	fprintf ( dataFile, "nPhiTrack            %10d  \n", nPhiTrack       ) ;
	fprintf ( dataFile, "etaMinTrack          %10.2e\n", etaMinTrack   ) ;
	fprintf ( dataFile, "etaMaxTrack          %10.2e\n", etaMaxTrack   ) ;
	fprintf ( dataFile, "phiMinTrack          %10.2e\n", phiMinTrack   ) ;
	fprintf ( dataFile, "phiMaxTrack          %10.2e\n", phiMaxTrack   ) ;
	fprintf ( dataFile, "nPrimaryPasses       %10d  \n", nPrimaryPasses  ) ;
	fprintf ( dataFile, "nSecondaryPasses     %10d  \n", nSecondaryPasses  ) ;
	fprintf ( dataFile, "vertexConstrainedFit %10d  \n", vertexConstrainedFit ) ;
	fprintf ( dataFile, "parameterLocation    %10d  \n", parameterLocation ) ;
	fprintf ( dataFile, "rowInnerMost         %10d  \n", rowInnerMost      ) ;
	fprintf ( dataFile, "rowOuterMost         %10d  \n", rowOuterMost      ) ;
	fprintf ( dataFile, "rowStart             %10d  \n", rowStart          ) ;
	fprintf ( dataFile, "rowEnd               %10d  \n", rowEnd            ) ;
	fprintf ( dataFile, "szFitFlag            %10d  \n", szFitFlag         ) ;
	fprintf ( dataFile, "maxChi2Primary       %10.2e\n", maxChi2Primary    ) ;
	fprintf ( dataFile, "bField               %10.2e\n", bField            ) ;
	fprintf ( dataFile, "hitChi2Cut           %10.2e\n", hitChi2Cut ) ;
	fprintf ( dataFile, "goodHitChi2          %10.2e\n", goodHitChi2 ) ;
	fprintf ( dataFile, "trackChi2Cut         %10.2e\n", trackChi2Cut ) ;
	fprintf ( dataFile, "goodDistance         %10.2e\n", goodDistance ) ;
	fprintf ( dataFile, "ptMinHelixFit        %10.2e\n", ptMinHelixFit ) ;
	fprintf ( dataFile, "maxDistanceSegment   %10.2e\n", maxDistanceSegment ) ;
	fprintf ( dataFile, "xyErrorScale         %10.2e\n", xyErrorScale       ) ;
	fprintf ( dataFile, "szErrorScale         %10.2e\n", szErrorScale       ) ;
	fprintf ( dataFile, "xVertex              %10.2e\n", xVertex            ) ;
	fprintf ( dataFile, "yVertex              %10.2e\n", yVertex            ) ;
	fprintf ( dataFile, "dxVertex             %10.2e\n", dxVertex           ) ;
	fprintf ( dataFile, "dyVertex             %10.2e\n", dyVertex           ) ;
	fprintf ( dataFile, "zVertex              %10.2e\n", zVertex            ) ;
	fprintf ( dataFile, "xyWeightVertex       %10.2e\n", xyWeightVertex     ) ;
	fprintf ( dataFile, "phiVertex            %10.2e\n", phiVertex          ) ;
	fprintf ( dataFile, "rVertex              %10.2e\n", rVertex            ) ;
	fprintf ( dataFile, "maxTime              %10.2e\n", maxTime            ) ;
}

void TrackFindingParameters::setDefaults (void)
{
	/*  Define cuts - this should be obsolete */

	modRow          = 1    ;
	infoLevel       = 0 ;
	hitChi2Cut      = 500.F  ;
	goodHitChi2     = 100.F ;
	trackChi2Cut    = 250.F ;
	maxChi2Primary  = 0. ;
	segmentRowSearchRange = 1 ;
	trackRowSearchRange   = 3 ;
	dEdx              = 0     ;
	dEdxNTruncate     = 20    ;
	dphi              = 0.10F * modRow ;
	deta              = 0.10F * modRow ;
	dphiMerge         = 0.02F  ;
	detaMerge         = 0.02F  ;
	distanceMerge     = 2. ;
	etaMin            = -2.5F  ;
	etaMinTrack       = -2.2F  ;
	etaMax            =  2.5F  ;
	etaMaxTrack       =  2.2F  ;
	eventReset        =  1     ;
	getErrors         =  0     ;
	fillTracks        =  1     ;
	ghostFlag         =  0     ;
	goBackwards       =  0     ;
	goodDistance      =  1.F * modRow ;
	init              =  0 ;
	mergePrimaries    =  1    ;
	parameterLocation =  1    ;
	phiMin            =  -0.000001/toDeg  ;
	phiMinTrack       =  -0.000001/toDeg  ;
	phiMax            =  360.2/toDeg  ;
	phiMaxTrack       =  360.2/toDeg  ;
	maxDistanceSegment = 100.F * modRow ;
	minHitsPerTrack   = 5      ;
	nHitsForSegment   = 2      ;
	nEta              = 60     ;
	nEtaTrack         = 60     ;
	nPhi              = 20     ;
	nPhiTrack         = 60     ;
	nPrimaryPasses    = 1      ;
	nSecondaryPasses  = 0      ;
	vertexConstrainedFit = 0 ;
	rowInnerMost      = 0      ;
	rowOuterMost      = 500     ;
	rowStart          = 500     ;
	rowEnd            =  0     ;
	segmentMaxAngle   = 10.F/toDeg ;
	szFitFlag         = 1      ;
	xyErrorScale      = 1.0F   ;
	szErrorScale      = 1.0F   ;
	bField            = 0.5F   ;
	phiShift          = 0.0    ;

	ptMinHelixFit     = 0.F  ;
	rVertex           = 0.F    ;
	xVertex           = 0.F    ;
	yVertex           = 0.F    ;
	zVertex           = 0.F    ;
	dxVertex          = 0.005F ;
	dyVertex          = 0.005F ;
	phiVertex         = 0.F    ;
	maxTime           = 1.e18F ; // by default tracker can run as long as the age of the Universe

	return  ;
}

void TrackFindingParameters::print()
{
	cout << "Track Finding Parameters: " << endl;

	cout << "infoLevel        " << infoLevel<< endl;
	cout << "segmentRowSearch  " << segmentRowSearchRange<< endl;
	cout << "trackRowSearch    " << trackRowSearchRange << endl;
	cout << "getErrors         " << getErrors << endl;
	cout << "fillTracks        " << fillTracks << endl;
	cout << "ghostFlag         " << ghostFlag << endl;
	cout << "goBackwards       " << goBackwards << endl;
	cout << "mergePrimaries    " << mergePrimaries << endl;
	cout << "minHitsPerTrack   " << minHitsPerTrack << endl;
	cout << "modRow            " << modRow << endl;
	cout << "nHitsForSegment   " << nHitsForSegment << endl;
    cout << "minHitsForFit     " << minHitsForFit   << endl;
	cout << "nEta              " << nEta            << endl;
	cout << "nPhi              " << nPhi            << endl;
	cout << "deta              " <<  deta         << endl;
	cout << "dphi              " <<  dphi         << endl;
	cout << "etaMin            " <<  etaMin        << endl;
	cout << "etaMax            " <<  etaMax        << endl;
	cout << "phiMin            " <<  phiMin        << endl;
	cout << "phiMax            " <<  phiMax        << endl;
	cout << "phiShift          " <<  phiShift   << endl;
	cout << "detaMerge         " <<  detaMerge    << endl;
	cout << "distanceMerge     " <<  distanceMerge << endl;
	cout << "nEtaTrack         " << nEtaTrack       << endl;
	cout << "nPhiTrack         " << nPhiTrack       << endl;
	cout << "etaMinTrack       " <<  etaMinTrack   << endl;
	cout << "etaMaxTrack       " <<  etaMaxTrack   << endl;
	cout << "phiMinTrack       " <<  phiMinTrack   << endl;
	cout << "phiMaxTrack       " <<  phiMaxTrack   << endl;
	cout << "nPrimaryPasses    " << nPrimaryPasses  << endl;
	cout << "nSecondaryPasses  " << nSecondaryPasses  << endl;
	cout << "vertexConstrainedFit "<< vertexConstrainedFit << endl;
	cout << "parameterLocation " << parameterLocation << endl;
	cout << "rowInnerMost      " << rowInnerMost      << endl;
	cout << "rowOuterMost      " << rowOuterMost      << endl;
	cout << "rowStart          " << rowStart          << endl;
	cout << "rowEnd            " << rowEnd            << endl;
	cout << "szFitFlag         " << szFitFlag         << endl;
	cout << "maxChi2Primary    " <<  maxChi2Primary    << endl;
	cout << "bField            " <<  bField            << endl;
	cout << "hitChi2Cut        " <<  hitChi2Cut << endl;
	cout << "goodHitChi2       " <<  goodHitChi2 << endl;
	cout << "trackChi2Cut      " <<  trackChi2Cut << endl;
	cout << "goodDistance      " <<  goodDistance << endl;
	cout << "ptMinHelixFit     " <<  ptMinHelixFit << endl;
	cout << "maxDistanceSegment" <<  maxDistanceSegment << endl;
	cout << "xyErrorScale      " <<  xyErrorScale       << endl;
	cout << "szErrorScale      " <<  szErrorScale       << endl;
	cout << "xVertex           " <<  xVertex            << endl;
	cout << "yVertex           " <<  yVertex            << endl;
	cout << "dxVertex          " <<  dxVertex           << endl;
	cout << "dyVertex          " <<  dyVertex           << endl;
	cout << "zVertex           " <<  zVertex            << endl;
	cout << "xyWeightVertex    " <<  xyWeightVertex     << endl;
	cout << "phiVertex         " <<  phiVertex          << endl;
	cout << "rVertex           " <<  rVertex            << endl;
	cout << "maxTime           " <<  maxTime            << endl;

}
