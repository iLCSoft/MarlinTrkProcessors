#ifndef TRACKFINDINGPARAMETERS_H
#define TRACKFINDINGPARAMETERS_H
//
// Based on the FTF code written by Pablo Yepes. 
//
// P. Yepes, “A fast track pattern recognition,” 
// Nuclear Instruments & Methods in Physics Research A 
// 380(1996) pp. 582-585.
//

#include <stdio.h>

namespace ftf
{
	class TrackFindingParameters 
	{          
	public:
		TrackFindingParameters( ) { setDefaults() ; } ;
		void       setDefaults( ) ;
		void       read( char* inputFile ) ;
		void       write( char* outputFile ) ;
		void       write( FILE* outFile ) ;
		void       print();

		int        infoLevel;       // Level of information printed about progress
		int        segmentRowSearchRange;       // Row search range for segments 
		int        trackRowSearchRange;         // Row search range for tracks 
		int        dEdx  ;          // dEdx switch
		int        dEdxNTruncate ;  // # points to truncate in dEdx
		int        eventReset   ;   // Flag to reset event in fft 
		int        getErrors    ;   // Flag to switch error calculation
		int        fillTracks   ;   // Flag to switch FtfTrack class filling
		int        ghostFlag    ;   // =1 when there are ghost hits
		int        goBackwards  ;   // Flag to go backwards at the end of track reco
		int        init;            // Control initialization 
		int        mergePrimaries ; // Switch to control primary merging 
		int        minHitsPerTrack; // Minimum # hits per track 
		int        modRow;          // Modulo pad row number to use 
		int        nHitsForSegment; // # hits in initial segments 
		int        minHitsForFit;
		int        nEta;            // # volumes in eta 
		int        nEtaTrack;       // # Track areas in eta 
		int        nPhi;            // # volumes in nphi 
		int        nPhiTrack;       // # Track areas in nphi 
		int        nPrimaryPasses;  // # iterations looking for primaries
		int        nSecondaryPasses;// # iterations looking for secondaries
		int        vertexConstrainedFit; // 
		int        parameterLocation; // 0=inner most point, 1=closest approach
		double     maxChi2Primary ; // maximum chi2 to be considered primary 
		int        rowInnerMost;    // Row where end track search 
		int        rowOuterMost;    // Outer most row to consider tin tracking
		int        rowStart;        // Row where start track search
		int        rowEnd  ;        // Row where end   track search
		int        szFitFlag;       // Switch for sz fit 
		double     bField      ;    // Magnetic field  
		double     hitChi2Cut;      // Maximum hit chi2 
		double     goodHitChi2;     // Chi2 to stop looking for next hit 
		double     trackChi2Cut;    // Maximum track chi2 
		double     deta;            // Eta search range 
		double     dphi;            // Phi search range 
		double     detaMerge ;      // Eta difference for track merge 
		double     dphiMerge ;      // Phi difference for track merge
		double     distanceMerge ;  // Maximum distance for reference point to merge secondaries
		double     etaMin;          // Min eta to consider 
		double     etaMinTrack ;    // Track min eta to consider 
		double     etaMax;          // Max eta to consider 
		double     etaMaxTrack ;    // Track max eta to consider 
		double     goodDistance ;   // In segment building
		// distance consider good enough 
		double     phiMin;          // Min phi to consider 
		double     phiMinTrack ;    // Track min phi to consider 
		double     phiMax;          // Max phi to consider 
		double     phiMaxTrack ;    // Track max phi to consider 
		double     phiShift      ;  // Shift in phi when calculating phi
		double     ptMinHelixFit ;  // Minimum pt to apply helix fit
		double     maxDistanceSegment; // Maximum distance for segments 
		double     segmentMaxAngle; // Maximum angle between to consecutive track pieces 
		// when forming segments. A piece is the connection 
		// two hits
		double     szErrorScale;    // sz error scale 
		double     xyErrorScale;    // xy error scale 
		double     xVertex      ;   // x position primary vertex 
		double     yVertex      ;   // y position primary vertex 
		double     dxVertex     ;
		double     dyVertex     ;
		double     zVertex      ;
		double     xyWeightVertex;  // Weight vertex in x-y
		double     phiVertex      ;
		double     rVertex        ;
		double     maxTime        ; // maxTime tracker can run
		int        phiClosed ;
		int        primaries  ;
		int        nRowsPlusOne, nPhiPlusOne   ; // Number volumes + 1
		int        nEtaPlusOne, nPhiEtaPlusOne ; // Number volumes + 1 
		int        nPhiTrackPlusOne, nEtaTrackPlusOne ;                
		double     phiSlice, etaSlice ;
		double     phiSliceTrack, etaSliceTrack ;

    //#define TRDEBUG 1
#ifdef TRDEBUG
		int       trackDebug ;
		int       hitDebug ;
		int       debugLevel ;
#endif
	} ;

} // end namespace ftf
#endif
