#include "TrackFinder.h"
using namespace ftf;
using std::max;

//*********************************************************************
//      Initializes the package
//*********************************************************************
TrackFinder::TrackFinder ( ) 
{
	hit        = 0 ;  
	track      = 0 ;
	trackC     = 0 ;
	volumeC    = 0 ;
	rowC       = 0 ;
	nHitsOutOfRange = 0 ;
}
//*********************************************************************
//      Initializes the package
//*********************************************************************
TrackFinder::~TrackFinder ( ) 
{
	if ( volumeC != 0 ) delete[] volumeC  ;
	if ( rowC  != 0 ) delete[] rowC ;
	if ( trackC != 0 ) delete[] trackC ;
}
//*********************************************************************
//      Steers the tracking 
//*********************************************************************
double TrackFinder::process (  ) 
{  
	//-----------------------------------------------------------------
	//        Make sure there is something to work with
	//------------------------------------------------------------------ 
	if ( nHits <= 0 ) {
		if ( para.infoLevel > 2 )         
			fprintf ( stderr, "fft: Hit structure is empty \n " ) ;
		return 1 ;
	}
	//
	initialCpuTime  = CpuTime ( );
	//
	//        General initialization 
	//
	if ( para.init == 0 ) {
		if ( reset ( ) ) return 1 ;
	}
	//
	//      Event reset and set pointers
	//
	if ( para.eventReset  && setPointers ( ) ) return 1 ;   
	//
	//      Build primary tracks now
	//
	short i ;
	para.primaries = 1 ;
	for ( i = 0 ; i < para.nPrimaryPasses ; i++ )
		if ( getTracks() != 0 ) break ;
	//
	//      Look for secondaries    
	//
	para.primaries = 0 ;
	for ( i = 0 ; i < para.nSecondaryPasses ; i++ )
		if ( getTracks() != 0 ) break ;

	//   if ( para.dEdx ) dEdx ( ) ;

	cpuTime  = CpuTime  ( ) - initialCpuTime  ;
#ifdef DEBUG
	if ( para.infoLevel > 0 )
		fprintf ( stderr, "TrackFinder::process: cpu %7.3f \n", cpuTime ) ;
#endif
	return cpuTime ;
} 
//********************************************************************
//     Calculates deposited Energy
//********************************************************************
void TrackFinder::dEdx ( ) 
{
	for ( int i = 0 ; i<nTracks ; i++ ){
		track[i].dEdx( ) ;
	}
}

//********************************************************************** 
//	Reconstruct primary tracks 
//**********************************************************************
int TrackFinder::getTracks ( ) 
{
	//
	//     Set conformal coordinates if we are working with primaries
	//
	int nHitsSegment   = (short)para.nHitsForSegment;  
	if ( para.primaries != 0 ) {
		setConformalCoordinates ( ) ;
		para.minHitsForFit = 1 ;
		para.nHitsForSegment = max(2,nHitsSegment);
	}
	else {
		para.minHitsForFit = 2 ;
		para.nHitsForSegment = max(3,nHitsSegment);
	}
	//
	//               Loop over rows   
	//
	for ( int ir = para.nRowsPlusOne - 1 ; ir>=para.minHitsPerTrack ; ir--) {
		//
		//           Loop over hits in this particular row
		//
		if ( rowC[ir].first != 0 &&  ((rowC[ir].first)->row) < para.rowEnd ) break ;
		//    if ( ((rowC[ir].first)->row) < para.rowEnd ) break ;
		for ( Hit* firstHit = rowC[ir].first ;
			firstHit != 0 ;
			firstHit =(firstHit->nextRowHit) ) {
				//
				//     Check hit was not used before
				//
				if ( firstHit->track != 0  ) continue ;
				//
				//     One more track 
				//
				nTracks++ ;
				//
				//
				if ( nTracks > maxTracks ){
					fprintf(stderr,"\n TrackFinder::getTracks: Max nr tracks reached !") ;
					nTracks = maxTracks  ;
					return 1 ;
				}
				//
				//     Initialize variables before going into track hit loop
				//
				Track* thisTrack    = &track[nTracks-1];
				thisTrack->para     = &para ;
				thisTrack->id       = nTracks ;
				thisTrack->firstHit = thisTrack->lastHit = firstHit ;
				thisTrack->innerMostRow = thisTrack->outerMostRow = firstHit->row ;
				thisTrack->xRefHit  = firstHit->x ;
				thisTrack->yRefHit  = firstHit->y ;
				thisTrack->xLastHit = firstHit->x ;
				thisTrack->yLastHit = firstHit->y ;
#ifdef TRDEBUG
				thisTrack->debugNew ( ) ;
#endif
				//
				//              Set fit parameters to zero
				//
				thisTrack->reset ( ) ;
				//
				//      Go into hit looking loop
				//
				if ( thisTrack->buildTrack ( firstHit, volumeC ) != 0 ) {
					//
					//    Merge Tracks if requested
					//
					if ( para.primaries &&
						para.mergePrimaries == 1 &&
						para.fillTracks &&
						thisTrack->mergePrimary( trackC )  ) {
							nTracks-- ;
							thisTrack->deleteCandidate ( ) ;
					}
				}
				else{
					//
					//      If track was not built delete candidate
					//
					thisTrack->deleteCandidate ( ) ;
					nTracks-- ;
				}
				//    
				//       End loop over hits inside row               
				//
		}
		//       End loop over rows                           
		//
		//    Check time
		//
		if ( CpuTime() - initialCpuTime > para.maxTime ) {
			fprintf ( stderr, "TrackFinder::getTracks: tracker time out after %f\n", 
				CpuTime() - initialCpuTime ) ;
			break ;
		}
	}
	//
	para.nHitsForSegment = nHitsSegment ;  
	//
	return 0 ;
}
//********************************************************************
//
void TrackFinder::mergePrimaryTracks ( ) 
{
	//
	//   Reset area keeping track pointers
	// 
	memset ( trackC, 0, para.nPhiTrackPlusOne*para.nEtaTrackPlusOne*sizeof(Container) ) ;  
	//
	//    Loop over tracks
	//

	for ( int i = 0 ; i < nTracks ; i++ ) {
		currentTrack = &(track[i]);
		if ( currentTrack->flag < 0 ) continue ;
		//
		//  reset link to following track
		//
		currentTrack->nxatrk = 0 ;
		//
		//    Try to merge this track 
		//    if track is not merged is added
		//    to the track volume (area)
		//
		if ( currentTrack->mergePrimary ( trackC ) ) {
			currentTrack->flag = -1 ;
		}
	}
}
//********************************************************************
//      Resets program
//*********************************************************************
int TrackFinder::reset (void)
{
	double phiDiff ;
	//
	//   Initialization flag in principle assume failure
	//
	para.init = 0 ;
	//----------------------------------------------------------------------------
	//     Allocate volumes 
	//---------------------------------------------------------------------------*/
	para.nRowsPlusOne = ( para.rowOuterMost - para.rowInnerMost ) / para.modRow + 2 ;
	if ( para.nRowsPlusOne < 1 ) {
		fprintf ( stderr, " =====> Error <===== \n" ) ;
		fprintf ( stderr, " Rows: Outer Most Inner Most %d % d \n", 
			para.rowOuterMost,  para.rowInnerMost ) ;
		return 1 ;
	}
	para.nPhiPlusOne    = para.nPhi + 1 ;
	para.nEtaPlusOne    = para.nEta + 1 ;
	para.nPhiEtaPlusOne = para.nPhiPlusOne * para.nEtaPlusOne ;
	if ( para.mergePrimaries ) {
		para.nPhiTrackPlusOne = para.nPhiTrack + 1 ;
		para.nEtaTrackPlusOne = para.nEtaTrack + 1 ;
	}
	//
	//-->    Allocate volume memory
	//
	if (volumeC != NULL) delete[] volumeC; 
	int nVolumes = para.nRowsPlusOne*para.nPhiPlusOne *
		para.nEtaPlusOne ;
	volumeC = new Container[nVolumes];
	if(volumeC == NULL) {
		fprintf ( stderr, "Problem with memory allocation... exiting\n" ) ;
		return 1 ;
	}
	// 
	//      Allocate row memory
	//
	if ( rowC != NULL ) delete[] rowC ;
	rowC = new Container[para.nRowsPlusOne];
	if ( rowC == NULL) {
		fprintf ( stderr, "Problem with memory allocation... exiting\n" ) ;
		exit(0);
	}
	//
	//       Allocate track area memory
	//
	if ( para.mergePrimaries ) {
		if (trackC    != NULL) delete []trackC     ;
		int nTrackVolumes = para.nPhiTrackPlusOne*
			para.nEtaTrackPlusOne ;
		trackC    = new TrackContainer[nTrackVolumes];
		if(trackC == NULL) {
			fprintf ( stderr, "Problem with memory allocation... exiting\n" ) ;
			return 1 ;
		}
		else{
			//
			//   Check there is some memory allocated
			//
			if ( trackC == 0 ){
				fprintf ( stderr, "TrackFinder::reset: Merging option not available \n " ) ; 
				printf ( " when option was not used the first time         \n " ) ; 
				return 1 ;
			}
		}
	}
	/*--------------------------------------------------------------------------
	Look whether the phi range is closed (< 5 degrees )
	-------------------------------------------------------------------------- */
	phiDiff = para.phiMax - para.phiMin ;
	if ( phiDiff > 2. * pi + 0.1 ) 
	{
		fprintf ( stderr, "TrackFinder::reset: Wrong phi range %f, %f ", 
			para.phiMin*toDeg, para.phiMax*toDeg ) ;
		return 1 ;
	}
	if ( fabs(phiDiff-twoPi ) < pi / 36. )
	{
		para.phiClosed = 1 ;
	}
	else
	{
		para.phiClosed = 0 ;
	}
	/*--------------------------------------------------------------------------
	Calculate volume dimensions
	-------------------------------------------------------------------------- */
	para.phiSlice   = (para.phiMax - para.phiMin)/para.nPhi ;
	para.etaSlice   = (para.etaMax - para.etaMin)/para.nEta ;
	/*--------------------------------------------------------------------------
	If needed calculate track area dimensions
	-------------------------------------------------------------------------- */
	para.phiSliceTrack   = (para.phiMaxTrack - para.phiMinTrack)/para.nPhiTrack ;
	para.etaSliceTrack   = (para.etaMaxTrack - para.etaMinTrack)/para.nEtaTrack ;
	//
	//    Set vertex parameters
	//
	if ( para.xVertex != 0 || para.yVertex != 0 ) 
	{ 
		para.rVertex   = (double)sqrt (para.xVertex*para.xVertex +
			para.yVertex*para.yVertex) ;
		para.phiVertex = (double)atan2(para.yVertex,para.xVertex);
	}
	else 
	{
		para.rVertex   = 0.F ;
		para.phiVertex = 0.F ;
	}

	if ( para.dxVertex != 0 || para.dyVertex != 0 )
	{
		para.xyWeightVertex = 1.F / ((double)sqrt(para.dxVertex*para.dxVertex + para.dyVertex*para.dyVertex) ) ;
	}
	else 
	{
		para.xyWeightVertex = 1.0F ;
	}
	//
	//   Set # hits & tracks to zero
	//
	// nHits   = 0 ;
	// nTracks = 0 ;
	//
	//    Set initialization flag to true
	//
	para.init = 1 ;
	return 0 ;
}

//*********************************************************************
//	Set hit pointers
//*********************************************************************
int TrackFinder::setConformalCoordinates ( )
{
	/*-------------------------------------------------------------------------
	Loop over hits 
	-------------------------------------------------------------------------*/
	Hit* thisHit ;
	double x, y, r2, invR2 ;
	for ( int ihit = 0 ; ihit<nHits ; ihit++ )
	{
		/*-------------------------------------------------------------------------
		Transform coordinates
		-------------------------------------------------------------------------*/
		thisHit = &(hit[ihit]) ;

		x            = thisHit->x - para.xVertex ;
		y            = thisHit->y - para.yVertex ;
		r2           = x * x + y * y ;
		invR2        = 1.F / r2 ;

		thisHit->xp    =     x * invR2 ;
		thisHit->yp    =   - y * invR2 ;
		thisHit->wxy   =   r2 * r2 /  ( square(para.xyErrorScale) * ( square(thisHit->dx) + square(thisHit->dy) ) ) ;
	} 

	return 0 ;
} 
//********************************************************************
//	Set hit pointers
//********************************************************************
int TrackFinder::setPointers ( )
{
	int ihit, localRow ;
	register int volumeIndex;
	double r, r2, phi, eta ;
	Hit* thisHit ;
	//
	nHitsOutOfRange = 0 ;
	//
	//   Set volumes to zero
	//
	memset ( rowC,   0, para.nRowsPlusOne*sizeof(Container) ) ;
	int n = para.nRowsPlusOne*para.nEtaPlusOne*para.nPhiPlusOne ;
	memset ( volumeC, 0, n*sizeof(Container) ) ;
	if ( para.mergePrimaries )
	{ 
		memset ( trackC, 0, para.nPhiTrackPlusOne*para.nEtaTrackPlusOne*sizeof(Container) ) ;  
	}
	/*-------------------------------------------------------------------------
	Loop over hits 
	-------------------------------------------------------------------------*/
	for ( ihit = 0 ; ihit<nHits ; ihit++ )
	{
		/*-------------------------------------------------------------------------
		Check whether row is to be considered
		-------------------------------------------------------------------------*/
		thisHit = &(hit[ihit]) ;
		localRow = thisHit->row - para.rowInnerMost ;
		if ( fmod((double)localRow,(double)para.modRow) != 0 ) continue ;

		if ( thisHit->row < para.rowInnerMost ) continue ;
		if ( thisHit->row > para.rowOuterMost ) continue ;
		/*
		*->    Get "renormalized" index
		*/
		localRow = localRow / para.modRow + 1 ;

		/*-------------------------------------------------------------------------
		Transform coordinates
		-------------------------------------------------------------------------*/
		r2            = thisHit->x * thisHit->x + thisHit->y * thisHit->y ;
		r             = sqrt ( r2 ) ;
		phi           = atan2(thisHit->y,thisHit->x) + para.phiShift ;
		if ( phi < 0 ) phi = phi + twoPi ;
		eta           = seta(r,thisHit->z) ;

		if ( para.szFitFlag ) {
			thisHit->s  = 0.F ;
			thisHit->wz = (1./ square ( para.szErrorScale * thisHit->dz ));
		}

		thisHit->r   = r   ;
		thisHit->phi = phi ;
		thisHit->eta = eta ;
		/*-------------------------------------------------------------------------
		Set pointers
		-------------------------------------------------------------------------*/
		thisHit->nextVolumeHit  = 
			thisHit->nextRowHit     = 0 ;
		/*-------------------------------------------------------------------------
		Get phi index for hit
		-------------------------------------------------------------------------*/

		thisHit->phiIndex = (int)( (thisHit->phi-para.phiMin)/para.phiSlice + 1.);
		if ( thisHit->phiIndex < 1 || thisHit->phiIndex > para.nPhi ) {
			if ( para.infoLevel > 2 ) {
				printf ( " === > Hit %d has Phi = %f \n", (int)thisHit->id,    
					thisHit->phi*toDeg ) ;
				printf ( " Phi index %d out of range \n", thisHit->phiIndex ) ;
			}
			nHitsOutOfRange++ ;
			continue ;
		} 

		/*-------------------------------------------------------------------------
		Get eta index for hit
		-------------------------------------------------------------------------*/

		thisHit->etaIndex = (int)((thisHit->eta - para.etaMin)/para.etaSlice + 1.);
		if ( thisHit->etaIndex < 1 || thisHit->etaIndex > para.nEta ) {
			if ( para.infoLevel > 2 ) {
				printf ( " \n === > Hit %d has Eta = %f  z %f ", (int)thisHit->id, 
					thisHit->eta, thisHit->z ) ;
				printf ( " \n Eta min/max %f %f ", para.etaMin, para.etaMax ) ;
				printf ( " \n Eta slice   %f    ", para.etaSlice ) ;
				fprintf ( stderr, " \n Eta index %d out of range ", thisHit->etaIndex ) ;	  
			}
			nHitsOutOfRange++ ;
			continue ;
		}
		//
		//    Reset track assigment
		//
		thisHit->nextTrackHit  = 0 ;
		thisHit->track         = 0 ;
		/* ------------------------------------------------------------------------- 
		Increase nr. of hits in small volume  WARNING! C-arrays go from 0
		Set Id of first hit in this vol. and link to next hit in previous
		hit in the same volume
		-------------------------------------------------------------------------*/
		volumeIndex = localRow  * para.nPhiEtaPlusOne + 
			thisHit->phiIndex * para.nEtaPlusOne + thisHit->etaIndex ;


		if (volumeC[volumeIndex].first == 0 ) 
			volumeC[volumeIndex].first = thisHit ;
		else
			((volumeC[volumeIndex].last))->nextVolumeHit = thisHit ;
		volumeC[volumeIndex].last = thisHit ;

		/*-------------------------------------------------------------------------
		Set row pointers
		-------------------------------------------------------------------------*/

		if ( rowC[localRow].first == NULL ){
			rowC [localRow].first = thisHit ;
		}
		else
			((rowC[localRow].last))->nextRowHit = thisHit ;
		rowC[localRow].last = thisHit ;
	}
	return 0 ;
} 
//***********************************************************************
//     For timing
//***********************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double TrackFinder::CpuTime( void )
{
	return (double)(clock()) / CLOCKS_PER_SEC;
}

