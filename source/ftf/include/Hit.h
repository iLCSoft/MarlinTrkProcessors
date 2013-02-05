#ifndef HIT_H
#define HIT_H
//
// Based on the FTF code written by Pablo Yepes. 
//
// P. Yepes, “A fast track pattern recognition,” 
// Nuclear Instruments & Methods in Physics Research A 
// 380(1996) pp. 582-585.
//

#include "TrackUtil.h"

namespace ftf
{
	class Track;

	class Hit
	{         
	public:
		// attributes
		int          id ;
		short        row;            // Row #     
		Track*       track;          // track to which the pnt was assigned
		Hit*         nextTrackHit;   // Next track hit             

		double       xyChi2 ;        // Chi2 in xy                 
		double       szChi2 ;        // Chi2 in sz                 
		double       x    ; 
		double       y    ;
		double       z    ;
		double       dx   ;          // error on the x coordinate  
		double       dy   ;          // error on the y coordinate  
		double       dz   ;          // error on the z coordinate  
		double       q    ;          // total charge assigned to this point 
		double       wxy  ;          // x-y weight x-y
		double       wz   ;          // z weight on z
		double       s    ;          // Track trajectory
		short        phiIndex ;        // Phi index    
		short        etaIndex ;        // Eta index      
		//
		Hit*         nextVolumeHit ;  // Next volume hit            
		Hit*         nextRowHit    ;  // Next row hit               
		double       r    ;            // radius                     
		double       phi  ;            // azimuthal angle            
		double       dphi ;            // Error in phi               
		double       eta  ;            // hit pseudorapidity         
		double       xp   ;            // x conformal coordinate 
		double       yp   ;            // y conformal coordinate 

		// methods
		void         print ( ) ;
		void         print ( int level ) ;
		void         printLinks ( ) ;
		void         printLinks ( int level ) ;
		void         setTrack ( Track* this_track ) ;
	} ;
} // end namespace ftf;
#endif
