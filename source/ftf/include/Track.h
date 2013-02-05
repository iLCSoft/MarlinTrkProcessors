#ifndef  TRACK_H
#define  TRACK_H
//
// Based on the FTF code written by Pablo Yepes. 
//
// P. Yepes, “A fast track pattern recognition,” 
// Nuclear Instruments & Methods in Physics Research A 
// 380(1996) pp. 582-585.
//

#include "TrackUtil.h"
#include "Hit.h"
#include "ThreeDPoint.h"
#include "TrackFindingParameters.h"

int const USE_SEGMENT= 1 ;
int const USE_FOLLOW = 2 ;
int const GO_DOWN    =-1 ;
int const GO_UP      = 1 ;

namespace ftf
{
	class Track 
	{ 

	public:

		//attributes

		Hit* firstHit;    // First hit belonging to track
		Hit* lastHit ;    // Last  hit belonging to track
		Hit* currentHit ;
		Hit* getCurrentHit ( ) { return currentHit ; } ;

		double    bField ; 

		int       id     ;  // primary key 
		short     flag   ;  // Primaries flag=1, Secondaries flag=0      
		int       innerMostRow ;
		int       outerMostRow ;
		short     nHits  ;  // Number of points assigned to that track
		short     nDedx  ;  // Number of points used for dEdx
		short     q  ;      // charge 
		double    chi2[2];  // chi squared of the momentum fit 
		double    dedx;     // dE/dx information 
		double    pt  ;     // pt (transverse momentum) at (r,phi,z) 
		double    phi0;     // azimuthal angle of point where parameters are given 
		double    psi ;     // azimuthal angle of the momentum at (r,.. 
		double    r0  ;     // r (in cyl. coord.) for point where parameters given 
		double    tanl;     // tg of the dip angle at (r,phi,z) 
		double    z0  ;     // z coordinate of point where parameters are given 
		double    length ;
		double    dpt ;
		double    dpsi;
		double    dz0 ;
		double    eta ;
		double    dtanl ;

		TrackFindingParameters*  para  ;    // Parameters pointer

		double    lastXyAngle ;    // Angle in the xy plane of line connecting to last hits        

		double    xRefHit ;
		double    yRefHit ;
		double    xLastHit ;
		double    yLastHit ;

		double    s11Xy  ;       // Fit Parameters
		double    s12Xy  ;
		double    s22Xy  ;
		double    g1Xy   ;
		double    g2Xy   ;       
		double    s11Sz  ;
		double    s12Sz  ;
		double    s22Sz  ;
		double    g1Sz   ;
		double    g2Sz   ; 

		double    ddXy, a1Xy, a2Xy ;    /*fit par in xy */
		double    ddSz, a1Sz, a2Sz ;    /*fit par in sz */

		Track*    nxatrk  ;

		//methods

		Track ( ) ;
		virtual ~Track(){};

		int         fitHelix        (  ) ;
		int         fitCircle       (  ) ;
		int         fitLine         (  ) ;
		TrackFindingParameters* getPara()  { return para ; } ;
		int         getErrorsCircleFit ( double a, double b, double r ) ;


		double      arcLength       ( double x1, double y1, double x2, double y2 ) ;
		ThreeDPoint closestApproach ( double xBeam, double yBeam ) ;
		ThreeDPoint extraRadius     ( double r ) ;
		int         extraRCyl       ( double& r, double& phi, double& z, double& rc, double& xc, double& yc ) ;
		int         intersectorZLine( double a, double b, ThreeDPoint& cross ) ;
		ThreeDPoint getClosest      ( double xBeam, double yBeam, double& rc, double& xc, double& yc ) ;
		int         getClosest      ( double xBeam, double yBeam, double rc, double xc, double yc, double& xClosest, double& yClosest ) ;

		void     updateToRadius         ( double r ) ;
		void     updateToClosestApproach ( double xBeam, double yBeam ) ;
		int      phiRotate              ( double deltaPhi ) ;

		inline virtual   void startLoop ( ){ currentHit = firstHit ; } ;
		inline virtual   int  done      ( ) { return currentHit != 0 ; } ;
		void       Print                ( int level ) ;   

		void      add                   ( Hit* thisHit, int way ) ;
		void      add                   ( Track* thisTrack ) ;
		int       buildTrack            ( Hit* firstHit, Container* volume ) ;
		void      dEdx                  ( ) ;
		void      deleteCandidate       ( ) ;
		void      fill                  ( ) ;
		void      fillPrimary           ( double& xc, double& yc, double& rc,double xPar, double yPar ) ;
		void      fillSecondary         ( double& xc, double& yc, double xPar, double yPar ) ;
		int       follow                ( Container* volume, int way, int rowToStop ) ;
		int       followHitSelection    ( Hit* baseHit, Hit* candidateHit ) ;
		Track*    getNextTrack          ( )      { return nxatrk ; } ; 
		int       mergePrimary          ( TrackContainer* trackArea ) ;
		void      reset                 ( ) ;
		Hit*      seekNextHit           ( Container* volume, Hit* baseHit, int nradiusSteps, int whichFunction ) ;
		int       segment               ( Container* volume, int way ) ;
		int       segmentHitSelection   ( Hit* baseHit, Hit* candidateHit ) ;

    //#define TRDEBUG 1
#ifdef TRDEBUG
		void debugAsk                 ( ) ;
		void debugDeleteCandidate     ( ) ;
		void debugFill                ( ) ;
		void debugFollowCandidate     ( Hit* candidate_hit ) ;
		void debugFollowSuccess       ( double dxy, double dsz, double lchi2_xy,
			double lchi2_sz, double chi2_min,
			Hit* candidate_hit ) ;
		void debugInVolume            ( Hit* base_hit, Hit* current_hit ) ;
		void debugNew                 ( ) ;
#endif


		//private:
		inline virtual   void nextHit (){ currentHit = (currentHit)->nextTrackHit ; } ;

	} ;
} // end namepsace ftf
#endif
