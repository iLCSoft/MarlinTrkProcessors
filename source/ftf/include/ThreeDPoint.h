#ifndef THREEDPOINT_H
#define THREEDPOINT_H
//
// Based on the FTF code written by Pablo Yepes. 
//
// P. Yepes, “A fast track pattern recognition,” 
// Nuclear Instruments & Methods in Physics Research A 
// 380(1996) pp. 582-585.
//

namespace ftf
{
	class ThreeDPoint 
	{
	public:
		ThreeDPoint (  ) : x(0), y(0), z(0){ } ;
		ThreeDPoint ( double _x, double _y, double _z ) : x(_x), y(_y), z(_z){ } ; 
		void set ( double _x, double _y, double _z ){ x = _x ; y = _y ; z = _z ; } ; 
		double x ;
		double y ;
		double z ;
	};
} // end namespace ftf
#endif
