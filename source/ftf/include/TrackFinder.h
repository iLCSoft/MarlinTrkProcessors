#ifndef TRACKFINDER_H
#define TRACKFINDER_H
//
// Based on the FTF code written by Pablo Yepes. 
//
// P. Yepes, “A fast track pattern recognition,” 
// Nuclear Instruments & Methods in Physics Research A 
// 380(1996) pp. 582-585.
//

#include <string.h>

#include "TrackUtil.h"
#include "TrackFindingParameters.h"
#include "Hit.h"
#include "Track.h"

namespace ftf
{
	class TrackFinder {

	public:
		TrackFinder( ) ;
		~TrackFinder( ) ;

		void    dEdx                    ( );
		int     getTracks               ( );
		void    mergePrimaryTracks      ( );
		double  process                 ( );
		int     reset                   ( );
		int     setConformalCoordinates ( );
		int     setPointers             ( );
		double  CpuTime                 ( );

		int                    nHits;  
		int                    nHitsOutOfRange;
		int                    maxHits;  
		Hit*                   hit;  
		int                    nTracks; 
		Track*                 track;  
		TrackFindingParameters para;
		int                    maxTracks;
		Container*             volumeC;
		Container*             rowC;
		TrackContainer*        trackC;
		double                 initialCpuTime;
		double                 cpuTime;

	private: 
		Track*     currentTrack;

	} ;
} // end namespace ftf
#endif
