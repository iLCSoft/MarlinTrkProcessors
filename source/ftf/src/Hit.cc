#include "Hit.h"
using namespace ftf;
void Hit::print (  ) { print (11) ; }  
void Hit::print ( int point_level ) 
{
	//--
	//--     print hit info
	//  

	if ( point_level > 9 ) 
		printf ( "hit Row    x      y     z\n" ) ;

	if ( fmod((double)point_level,10) > 0 )
		printf ( "%3d %2d %6.2f %6.2f %6.2f \n", 
		id, row, x, y, z ) ;
}


void Hit::printLinks (  ) { print (11) ; } 
void Hit::printLinks ( int point_level ) 
{
	//--
	//--     print hit info
	//  

	if ( point_level > 9 ) 
		printf ( "hit ir iphi ieta   phi   eta      x      y     z\n" ) ;

	if ( fmod((double)point_level,10) > 0 )
		printf ( "%3d %3d %3d  %3d  %6.2f %5.2f  %6.2f %6.2f %6.2f \n", 
		(int)id, (int)row, (int)phiIndex, (int)etaIndex, 
		phi*toDeg, eta, x, y, z ) ;
	int vhit ;
	if ( nextVolumeHit != 0 ) vhit = (nextVolumeHit)->id ;
	else vhit = -1 ;
	int rhit ;
	if ( nextRowHit != 0 ) rhit = (nextRowHit)->id ;
	else rhit = -1 ;
	int thit ;
	if ( nextTrackHit != 0 ) thit = (nextTrackHit)->id ;
	else thit = -1 ;

	if ( fmod((double)point_level,10) > 1 ) 
		printf ( "pointers:vol,row,tr (%4d,%4d,%4d)\n ",
		vhit, rhit, thit) ; 
	//int tid ;
	//if ( track != 0 ) tid = track->id ;
	//else tid = -1 ;
	//if ( fmod((double)point_level,10) > 2 )
	//   printf ( "\n Tracks  :reco            (%4d) ", tid ) ;

	/*   if ( fmod((double)point_level,10)  printf ( "\n  " ) ; */
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    This function assigns this hit to a given track
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hit::setTrack ( Track* this_track ) 
{
	track = this_track ;
}
