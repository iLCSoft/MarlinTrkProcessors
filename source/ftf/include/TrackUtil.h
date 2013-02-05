#ifndef TRACKUTIL_H
#define TRACKUTIL_H
//
// Based on the FTF code written by Pablo Yepes. 
//
// P. Yepes, “A fast track pattern recognition,” 
// Nuclear Instruments & Methods in Physics Research A 
// 380(1996) pp. 582-585.
//

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

namespace ftf
{
	class Hit;
	class Track;
}
//
//    Some constants 
//
const double pi      = acos(-1.);
const double twoPi   = 2. * pi ; 
const double piHalf  = 0.5 * pi ;
const double toDeg   = 360./twoPi ;
const double bFactor = 0.00299792458 ;
//
//-->   Functions
//
//#ifndef min
//#define min(a,b)    ( ( (a) < (b) ) ? (a) : (b) )
//#endif
//#ifndef max
//#define max(a,b)    ( ( (a) > (b) ) ? (a) : (b) )
//#endif
#define seta(r,z)   (3.0F * (z) / (fabs(z)+2.0F*(r)))
//#define reta(eta,r) ((2.F*(r)*eta / ( 3 - fabs(eta)) )) 
//#define sgn(a)      ( ( (a) > 0   ) ? (1) :(-1) )
#define square(a)   ( (a) * (a) )


extern double fmod(double,double);
extern double sqrt(double);
extern double fabs(double);
extern double atan2(double,double);

namespace ftf
{
	class Container
	{
	public:
		Hit* first; 
		Hit* last;  
	};

	class TrackContainer
	{
	public:
		Track* first; 
		Track* last;  
	} ;

} // end namespace ftf
#endif
