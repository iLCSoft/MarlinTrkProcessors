#include "Track.h"
#include <cstring>
using namespace ftf;
using std::min;
using std::max;

void ftfMatrixDiagonal ( double *h, double &h11, double &h22, double &h33 ) ;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//      Track Initialization
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Track::Track ( ){
	firstHit = 0 ;
	lastHit  = 0 ;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//       Fit a circle
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::fitHelix (  ) 
{
	if ( fitCircle ( ) != 0 ){
		printf ( " Problem in Fit_Circle " ) ;
		return 1 ;
	}
	//
	//     Fit line in s-z plane now
	//
	if ( fitLine ( ) != 0) {
		printf ( " Problem fitting a line " ) ;
		return 1 ;
	}
	return 0 ;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//       End of Fit Helix
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    
//  Fits circle parameters using algorithm
//  described by ChErnov and Oskov in Computer Physics
//  Communications.
// 
//  Written in FORTRAN by Jawluen Tang, Physics department , UT-Austin 
//  Moved to C by Pablo Yepes
//---------------------------------------------------------------
int Track::fitCircle (  ) 
{
	double wsum  = 0.0 ;
	double xav   = 0.0 ;
	double yav   = 0.0 ;
	//
	//     Loop over hits calculating average
	//
	for ( startLoop() ; done() ; nextHit() ) { 

		Hit* cHit = currentHit ;
		cHit->wxy = 1.0F/ (double)(cHit->dx*cHit->dx + cHit->dy*cHit->dy) ;
		wsum      += cHit->wxy ;
		xav       += cHit->wxy * cHit->x ;
		yav       += cHit->wxy * cHit->y ;
	}
  
	if ( getPara()->vertexConstrainedFit != 0 ) {
		wsum += getPara()->xyWeightVertex ;
		xav  += getPara()->xVertex ;
		yav  += getPara()->yVertex ;
	}
  
	xav = xav / wsum ;
	yav = yav / wsum ;

	//
	//  CALCULATE <X**2>, <XY>, AND <Y**2> WITH <X> = 0, & <Y> = 0
	//
	double xxav  = 0.0 ;
	double xyav  = 0.0 ; 
	double yyav  = 0.0 ;
	double xi, yi ;

	for ( startLoop() ; done() ; nextHit() ) { 
		Hit* cHit = currentHit ;
		xi        = cHit->x - xav ;
		yi        = cHit->y - yav ;
		xxav     += xi * xi * cHit->wxy ;
		xyav     += xi * yi * cHit->wxy ;
		yyav     += yi * yi * cHit->wxy ;
	}
  
	if ( getPara()->vertexConstrainedFit != 0 ) {
		xi        = getPara()->xVertex - xav ;
		yi        = getPara()->yVertex - yav ;
		xxav     += xi * xi * getPara()->xyWeightVertex ;
		xyav     += xi * yi * getPara()->xyWeightVertex ;
		yyav     += yi * yi * getPara()->xyWeightVertex ; 
	}
  
	xxav = xxav / wsum ;
	xyav = xyav / wsum ;
	yyav = yyav / wsum ;
	//
	//-->  ROTATE COORDINATES SO THAT <XY> = 0
	//
	//-->  SIGN(C**2 - S**2) = SIGN(XXAV - YYAV) >
	//-->  &                                     > ==> NEW : (XXAV-YYAV) > 0
	//-->  SIGN(S) = SIGN(XYAV)                  >

	double a = fabs( xxav - yyav ) ;
	double b = 4.0 * xyav * xyav ;

	double asqpb  = a * a + b  ;
	double rasqpb = sqrt ( asqpb) ;

	double splus  = 1.0 + a / rasqpb ;
	double sminus = b / (asqpb * splus) ;

	splus  = sqrt (0.5 * splus ) ;
	sminus = sqrt (0.5 * sminus) ;
	//
	//->  FIRST REQUIRE : SIGN(C**2 - S**2) = SIGN(XXAV - YYAV)
	//
	double sinrot, cosrot ;
	if ( xxav <= yyav ) {
		cosrot = sminus ;
		sinrot = splus  ;
	}
	else {
		cosrot = splus ;
		sinrot = sminus ;
	}
	//
	//->  REQUIRE : SIGN(S) = SIGN(XYAV) * SIGN(C) (ASSUMING SIGN(C) > 0)
	//
	if ( xyav < 0.0 ) sinrot = - sinrot ;
	//
	//-->  WE NOW HAVE THE SMALLEST ANGLE THAT GUARANTEES <X**2> > <Y**2>
	//-->  TO GET THE SIGN OF THE CHARGE RIGHT, THE NEW X-AXIS MUST POINT
	//-->  OUTWARD FROM THE ORGIN.  WE ARE FREE TO CHANGE SIGNS OF BOTH
	//-->  COSROT AND SINROT SIMULTANEOUSLY TO ACCOMPLISH THIS.
	//
	//-->  CHOOSE SIGN OF C WISELY TO BE ABLE TO GET THE SIGN OF THE CHARGE
	//
	if ( cosrot*xav+sinrot*yav < 0.0 ) {
		cosrot = -cosrot ;
		sinrot = -sinrot ;
	}
	//
	//->  NOW GET <R**2> AND RSCALE= SQRT(<R**2>)
	//
	double rrav   = xxav + yyav ;
	double rscale = sqrt(rrav) ;

	xxav   = 0.0 ;
	yyav   = 0.0 ;
	xyav   = 0.0 ;
	double xrrav	 = 0.0 ;
	double yrrav	 = 0.0 ;
	double rrrrav  = 0.0 ;

	double xixi, yiyi, riri, wiriri, xold, yold ;
	for ( startLoop() ; done() ; nextHit() ) { 
		Hit* cHit = currentHit ;
		xold = cHit->x - xav ;
		yold = cHit->y - yav ;
		//
		//-->  ROTATE SO THAT <XY> = 0 & DIVIDE BY RSCALE SO THAT <R**2> = 1
		//
		xi = (  cosrot * xold + sinrot * yold ) / rscale ;
		yi = ( -sinrot * xold + cosrot * yold ) / rscale ;

		xixi   = xi * xi ;
		yiyi   = yi * yi ;
		riri   = xixi + yiyi ;
		wiriri = cHit->wxy * riri ;

		xyav   += cHit->wxy * xi * yi ;
		xxav   += cHit->wxy * xixi ;
		yyav   += cHit->wxy * yiyi ;

		xrrav  += wiriri * xi ;
		yrrav  += wiriri * yi ;
		rrrrav += wiriri * riri ;
	}
	//
	//   Include vertex if required
	//
	if ( getPara()->vertexConstrainedFit != 0 ) {
		xold = getPara()->xVertex - xav ;
		yold = getPara()->yVertex - yav ;
		//
		//-->  ROTATE SO THAT <XY> = 0 & DIVIDE BY RSCALE SO THAT <R**2> = 1
		//
		xi = (  cosrot * xold + sinrot * yold ) / rscale ;
		yi = ( -sinrot * xold + cosrot * yold ) / rscale ;

		xixi   = xi * xi ;
		yiyi   = yi * yi ;
		riri   = xixi + yiyi ;
		wiriri = getPara()->xyWeightVertex * riri ;

		xyav   += getPara()->xyWeightVertex * xi * yi ;
		xxav   += getPara()->xyWeightVertex * xixi ;
		yyav   += getPara()->xyWeightVertex * yiyi ;

		xrrav  += wiriri * xi ;
		yrrav  += wiriri * yi ;
		rrrrav += wiriri * riri ;
	}
	//
	//    
	//
	//-->  DIVIDE BY WSUM TO MAKE AVERAGES
	//
	xxav    = xxav   / wsum ;
	yyav    = yyav   / wsum ;
	xrrav   = xrrav  / wsum ;
	yrrav   = yrrav  / wsum ;
	rrrrav  = rrrrav / wsum ;
	xyav    = xyav   / wsum ;

	int const ntry = 5 ;
	//
	//-->  USE THESE TO GET THE COEFFICIENTS OF THE 4-TH ORDER POLYNIMIAL
	//-->  DON'T PANIC - THE THIRD ORDER TERM IS ZERO !
	//
	double xrrxrr = xrrav * xrrav ;
	double yrryrr = yrrav * yrrav ;
	double rrrrm1 = rrrrav - 1.0  ;
	double xxyy   = xxav  * yyav  ;        

	double c0  =          rrrrm1*xxyy - xrrxrr*yyav - yrryrr*xxav ;
	double c1  =        - rrrrm1      + xrrxrr      + yrryrr   - 4.0*xxyy ;        
	double c2  =   4.0  + rrrrm1                               - 4.0*xxyy ;           
	double c4  = - 4.0  ;                
	//
	//-->  COEFFICIENTS OF THE DERIVATIVE - USED IN NEWTON-RAPHSON ITERATIONS
	//
	double c2d =   2.0 * c2 ;
	double c4d =   4.0 * c4 ;
	//
	//-->  0'TH VALUE OF LAMDA - LINEAR INTERPOLATION BETWEEN P(0) & P(YYAV)
	//
	//   LAMDA = YYAV * C0 / (C0 + YRRSQ * (XXAV-YYAV))
	double lamda  = 0.0 ;
	double dlamda = 0.0 ;
	//
	double chiscl = wsum * rscale * rscale ;
	double dlamax = 0.001 / chiscl ;   

	double p, pd ;
	for ( int itry = 1 ; itry <= ntry ; itry++ ) {
		p      = c0 + lamda * (c1 + lamda * (c2 + lamda * lamda * c4 )) ;
		pd     = (c1 + lamda * (c2d + lamda * lamda * c4d)) ;
		dlamda = -p / pd ;
		lamda  = lamda + dlamda ;
		if (fabs(dlamda)<   dlamax) break ;
	}

	chi2[0]  = (double)(chiscl * lamda) ;
	// double dchisq = chiscl * dlamda ;	     
	//
	//-->  NOW CALCULATE THE MATRIX ELEMENTS FOR ALPHA, BETA & KAPPA
	//
	double h11   = xxav  -     lamda ;
	double h14   = xrrav ;
	double h22   = yyav  -     lamda ; 
	double h24   = yrrav ;
	double h34   = 1.0   + 2.0*lamda ;
	if ( h11 == 0.0 || h22 == 0.0 ){
		printf ( " Problems fitting a circle " ) ;
		return 1 ;
	}
	double rootsq = (h14*h14)/(h11*h11) + 4.0*h34 ;

	double ratio, kappa, beta ;
	if ( fabs(h22) > fabs(h24) ) {
		ratio  = h24 / h22 ;
		rootsq = ratio * ratio + rootsq ;
		kappa = 1.0 / sqrt(rootsq) ;
		beta  = - ratio * kappa ;
	}
	else {
		ratio  = h22 / h24 ;
		rootsq = 1.0 + ratio * ratio * rootsq ;
		beta  = 1.0 / sqrt(rootsq) ;
		if ( h24 > 0 ) beta = - beta ;
		kappa = -ratio * beta ;
	}            
	double alpha = - (h14/h11) * kappa ;
	//
	//-->  transform these into the lab coordinate system
	//-->  first get kappa and back to real dimensions
	//
	double kappa1 = kappa / rscale ;
	double dbro   = 0.5   / kappa1 ;
	//
	//-->  next rotate alpha and beta and scale
	//
	double alphar = (cosrot * alpha - sinrot * beta)* dbro ;
	double betar  = (sinrot * alpha + cosrot * beta)* dbro ;
	//
	//-->  then translate by (xav,yav)
	//
	double acent  = (double)(xav - alphar) ;
	double bcent  = (double)(yav - betar ) ;
	double radius = (double)dbro ;
	//
	//   Get charge
	//
	q = ( ( yrrav < 0 ) ? 1 : -1 ) ;
	//
	//    Get other track parameters
	//
	double x0, y0 ;
	if ( getPara()->vertexConstrainedFit ) {
		flag = 1 ; // primary track flag
		x0   = getPara()->xVertex ;
		y0   = getPara()->yVertex ;
		phi0 = getPara()->phiVertex ;
		r0   = getPara()->rVertex ;
	} 
	else {
		Hit* lHit = lastHit ;
		flag =  0 ; // primary track flag
		x0   =  lHit->x  ;
		y0   =  lHit->y  ;
		phi0 =  atan2(lHit->y,lHit->x);
		if ( phi0 < 0 ) phi0 += twoPi ;
		r0   =  sqrt ( lHit->x * lHit->x + lHit->y * lHit->y )  ;
	}
	//
	psi  = (double)atan2(bcent-y0,acent-x0) ;
	psi  = psi + q * 0.5F * pi ;
	if ( psi < 0 ) psi = psi + twoPi ;
	pt   = (double)(2.9979e-3 * getPara()->bField * radius ) ;
	//
	//    Get errors from fast fit
	//
	if ( getPara()->getErrors ) getErrorsCircleFit ( acent, bcent, radius ) ;
	//
	return 0 ;
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//    End Fit Circle
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Fit Line in s-z plane
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::fitLine ( )
{
	//
	//     initialization 
	//
	double sum = 0.F ;
	double ss  = 0.F ;
	double sz  = 0.F ;
	double sss = 0.F ;
	double ssz = 0.F ;
	//
	//     find sum , sums ,sumz, sumss 
	// 
	double dx, dy ;
	double radius = (double)(pt / ( 2.9979e-3 * getPara()->bField ) ) ;
	if ( getPara()->vertexConstrainedFit != 0 ) {
		dx   = (firstHit)->x - getPara()->xVertex ;
		dy   = (firstHit)->y - getPara()->yVertex ;
	}
	else {
		dx   = (firstHit)->x - (lastHit)->x ;
		dy   = (firstHit)->y - (lastHit)->y ;
	}
	double localPsi = 0.5F * sqrt ( dx*dx + dy*dy ) / radius ;
	double total_s ;
	if ( fabs(localPsi) < 1. ) {
		total_s = 2.0F * radius * asin ( localPsi ) ;
	} 
	else { 
		total_s = 2.0F * radius * pi ;
	} 

	//
	Hit* previousHit = 0  ;

	for ( startLoop() ; done() ; nextHit() ) {
		Hit* cHit = currentHit ;
		if ( currentHit != firstHit ) {
			dx   = cHit->x - previousHit->x ;
			dy   = cHit->y - previousHit->y ;
			dpsi = 0.5F * (double)sqrt ( dx*dx + dy*dy ) / radius ;
			if ( dpsi > 1.) {
				fprintf (stderr, "Hit::fitLine(): dpsi=%f\n", dpsi);
				dpsi = 1.;
			}
			cHit->s = previousHit->s - 2.0F * radius * (double)asin ( dpsi ) ;
		}
		else
			cHit->s = total_s ;

		sum += cHit->wz ;
		ss  += cHit->wz * cHit->s ;
		sz  += cHit->wz * cHit->z ;
		sss += cHit->wz * cHit->s * cHit->s ;
		ssz += cHit->wz * cHit->s * cHit->z ;
		previousHit = cHit ;
	}

	double det = sum * sss - ss * ss;
	if ( fabs(det) < 1e-20){ 
		chi2[1] = 99999.F ;
		return 0 ;
	}
	//
	//     compute the best fitted parameters A,B
	//
	tanl = (double)((sum * ssz - ss * sz ) / det );
	z0   = (double)((sz * sss - ssz * ss ) / det );
	//
	//     calculate chi-square 
	//
	chi2[1] = 0.F ;
	double r1 ;
	for ( startLoop() ; done() ; nextHit() ) {
		Hit* cHit = currentHit ;
		r1   = cHit->z - tanl * cHit->s - z0 ;
		chi2[1] += (double) ( (double)cHit->wz * (r1 * r1) );
	}
	//
	//     calculate estimated variance
	//      varsq=chi/(double(n)-2.) 
	//     calculate covariance matrix 
	//      siga=sqrt(varsq*sxx/det) 
	//      sigb=sqrt(varsq*sum/det) 
	//
	dtanl = (double) ( sum / det );
	dz0   = (double) ( sss / det );

	return 0 ;
} 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    End Fit Line
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// CIRCOV - a covariance matrix calculation program for circle fitting 
// DESCRIPTION: 
// Compute the covariance matrix of an effective circle fitting algorithm 
// The circle equation is (X(I)-A)**2 + (Y(I)-B)**2 = R**2. 
// The functional minimum is W(I)*[(X(I)-A)**2+(Y(I)-B)**2-R*R]**2/(R*R) 
// For details about the algorithm, see 
// N.I. CHERNOV, G.A. OSOSKOV, COMPUT. PHYS. COMMUN. 33(1984) 329-333 
// INPUT ARGUMENTS: */
//      A              - Best fitted circle center in X axis, REAL 
//      B              - Best fitted circle center in Y axis, REAL 
//      R              - Best fitted radius                   REAL 
//
// From a routine written in Fortran by  AUTHOR:
//  Jawluen Tang, Physics department , UT-Austin 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::getErrorsCircleFit ( double a, double b, double r ) {

	double h[9] = { 0. };
	double dx, dy ;
	double h11, h22, h33 ;
	static int j ;
	static double ratio, c1, s1;
	static double hyp;


	for (j = 0; j < 9; j++ ) {
		h[j] = 0.;
	}
	//
	//    If circle fit was not used the
	//    errors in the real space need to
	//    be calculated
	//
	if ( pt < getPara()->ptMinHelixFit ) {
		for ( startLoop() ; done() ; nextHit() ) { 

			Hit* cHit = currentHit ;
			cHit->wxy = 1.0F/ (double)(cHit->dx*cHit->dx + cHit->dy*cHit->dy) ;
		}
	}
	//
	//    Loop over points in fit
	//
	for ( startLoop() ; done() ; nextHit() ) {
		Hit* cHit = currentHit ;
		dx = cHit->x - a;
		dy = cHit->y - b;
		hyp = (double)sqrt( dx * dx + dy * dy );
		s1 = dx / hyp;
		c1 = dy / hyp;
		ratio = r / hyp;
		h[0] += cHit->wxy * (ratio * (s1 * s1 - 1) + 1);
		h[1] += cHit->wxy * ratio * s1 * c1;
		h[2] += cHit->wxy * s1;
		h[4] += cHit->wxy * (ratio * (c1 * c1 - 1) + 1);
		h[5] += cHit->wxy * c1;
		h[8] += cHit->wxy ;
	}
	h[3]  = h[1];
	h[6]  = h[2];
	h[7]  = h[5];

	ftfMatrixDiagonal  ( h, h11, h22, h33 ) ;
	//
	//   Calculate pt error now
	//
	dpt          = (double)(2.9979e-3 * getPara()->bField * h33 );
	//
	//     Get error in psi now
	//
	if ( getPara()->vertexConstrainedFit != 0 ) {
		dx = a ;
		dy = b ;
	}
	else {
		dx = (lastHit)->x + a - (firstHit)->x ;
		dy = (lastHit)->y + b + (firstHit)->y ;
	}
	double w   = dy / dx ;
	dpsi  = (double)(( 1. / ( 1. + w*w ) ) * ( h22 / dx - dy * h11 / ( dx*dx ) )) ;

	return 0 ;
}

//*************************************************************************
//   Prints one track information
//*************************************************************************
void Track::Print ( int level )
{
	double pmom, pz;   
	/*
	----->   Print info
	*/
	if ( level > 9 ) {
		pz   = pt * tanl ;
		pmom = (double)sqrt ( pz * pz + pt * pt  ) ;
		printf ( " \n =======> Track      %d  <======== ", id ) ;
		printf ( " \n p,  pt, q         %7.2f  %7.2f  %2d ", pmom, pt, q ) ;
	}   
	if ( level > 19 ) {
		printf ( " \n r0,   z0,  nHits  %7.2f  %7.2f %d    ", r0, z0, nHits ) ;
		printf ( " \n phi0, psi, tanl   %7.2f  %7.2f %7.2f ", phi0, psi, tanl ) ;
	}  
	else printf ( "\n " ) ; 

	if ( level > 29 ) {
		printf ( " \n chi2 (s,z)        %6.2e  %6.2e ", chi2[0],
			chi2[1] ) ;
	}
	else printf ( "\n " ) ;


	if ( fmod((double)level,10.) > 0 ) {
		printf ( " \n *** Clusters in this track *** \n" ) ;
		(firstHit)->print ( 10 ) ;
		for ( startLoop() ; done() ; nextHit()  ) { 
			(currentHit)->print ( 1 ) ;
		}
	}
	printf ( "\n " ) ; 
}
/*:>--------------------------------------------------------------------
**: METHOD:   Finds point of closest approach
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:    xBeam, yBeam: beam position
**:
**: RETURNS:    
**:             tHit            - Point closest approach to center
*:>------------------------------------------------------------------*/
ThreeDPoint Track::closestApproach ( double xBeam, double yBeam ) {
	double rc, xc, yc ;
	return getClosest ( xBeam, yBeam, rc, xc, yc ) ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Finds point of closest approach
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:    xBeam, yBeam: beam position
**:         OUT:
**:                  rc, xc, yc  track circle radius and center
**:
**: RETURNS:    
**:             tHit            - Point closest approach to center
*:>------------------------------------------------------------------*/
ThreeDPoint Track::getClosest ( double xBeam, double yBeam,  double &rc, double &xc, double &yc ) 
{
	double xp, yp, zp ;
	xp = yp = 0. ;
	//--------------------------------------------------------
	//     Get track parameters
	//--------------------------------------------------------
	double tPhi0 = psi + double(q) * 0.5 * pi / fabs((double)q) ;

	double x0   = r0 * cos(phi0) ;
	double y0   = r0 * sin(phi0) ;
	rc   = pt / ( bFactor * getPara()->bField )  ;
	xc   = x0 - rc * cos(tPhi0) ;
	yc   = y0 - rc * sin(tPhi0) ;

	getClosest ( xBeam, yBeam, rc, xc, yc, xp, yp ) ;

	//-----------------------------------------------------------------
	//     Get the z coordinate
	//-----------------------------------------------------------------
	double angle  = atan2 ( (yp-yc), (xp-xc) ) ;
	if ( angle < 0. ) angle = angle + 2.0 * pi ;

	double dangle = angle  - tPhi0  ;
	dangle = fmod ( dangle, 2.0 * pi ) ;
	if ( fabs(dangle) < 1.e-10 ) dangle = 0 ; // Problems with -0.000 values
	if ( (q * dangle) < 0 )
		dangle = dangle + q * 2. * pi ;

	double stot   = fabs(dangle) * rc ;
	zp   = z0 - stot * tanl ;

	ThreeDPoint tHit(xp,yp,zp) ;
	return tHit ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Finds point of closest approach
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:    xBeam, yBeam: beam position
**:                 rc, xc, yc  : track circle radius and center
**:         OUT:
**:                 double xClosest, yClosest
**:
**: RETURNS: 0 if Ok
*:>------------------------------------------------------------------*/
int Track::getClosest ( double xBeam, double yBeam, double  rc, double xc, double yc, double& xClosest, double& yClosest ) 
{
	//----------------------------------------------------------
	//     Shift center to respect beam axis
	//----------------------------------------------------------
	double dx = xc - xBeam ;
	double dy = yc - yBeam ;
	//----------------------------------------------------------
	//     Solve the equations 
	//----------------------------------------------------------
	double fact = rc / sqrt ( dx * dx + dy * dy ) ;
	double f1   = 1. + fact ;
	double f2   = 1. - fact ;

	double dx1 = dx * f1 ;
	double dy1 = dy * f1 ;
	double d1 = sqrt ( dx1 * dx1 + dy1 * dy1 ) ;

	double dx2 = dx * f2 ;
	double dy2 = dy * f2 ;
	double d2 = sqrt ( dx2 * dx2 + dy2 * dy2 ) ;
	//---------------------------------------------------------------
	//     Choose the closest
	//---------------------------------------------------------------
	if ( d1 < d2 ) {
		xClosest = dx1 + xBeam ;
		yClosest = dy1 + yBeam ;
	}
	else {
		xClosest = dx2 + xBeam ;
		yClosest = dy2 + yBeam ;
	}
	return 0 ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Extrapolates track to cylinder with radius r
**:
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             track           - Global track pointer
**:             r               - Cylinder radius     
**:         OUT:
**:             x,y,z           - Extrapolated point
**:             xc,yc,rr        - Center and radius track circle in x-y plane
**:
**: RETURNS:    0=ok, <>0 error
**:>------------------------------------------------------------------*/
ThreeDPoint Track::extraRadius ( double r ) { 
	double phi ;
	//
	// Default values 
	//
	double x, y, z, rc, xc, yc ;
	x = y = z = 0.F ;
	//
	//    If error return with error 
	//
	ThreeDPoint tHit(0,0,0) ;
	if ( extraRCyl ( r, phi, z, rc, xc, yc ) ) return tHit ;
	//
	//   Otherwise get point in cartesian coordinates
	//
	x = r * cos(phi) ;
	y = r * sin(phi) ; 
	tHit.x = x ;
	tHit.y = y ;
	tHit.z = z ;

	return tHit ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Extrapolates track to cylinder with radius r
**:
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             r               - Cylinder radius
**:         OUT:
**:             phi,z           - Extrapolated point
**:             xc,yc,rc        - Center and radius track circle in x-y plane
**:
**: RETURNS:    0=ok, <>0 error
**:>------------------------------------------------------------------*/
int Track::extraRCyl ( double &r,  double &phi, double &z,
					  double &rc, double &xc,  double &yc )  
{

	double td  ;
	double fac1,sfac, fac2,deltat ;
	//--------------------------------------------------------
	//     Get track parameters
	//--------------------------------------------------------
	double tPhi0 = psi + double(q) * 0.5 * pi / fabs((double)q) ;
	double x0    = r0 * cos(phi0) ;
	double y0    = r0 * sin(phi0) ;
	rc    = fabs(pt) / ( bFactor * getPara()->bField )  ;
	xc    = x0 - rc * cos(tPhi0) ;
	yc    = y0 - rc * sin(tPhi0) ;
	//  
	//     Check helix and cylinder intersect
	// 
	fac1 = xc*xc + yc*yc ;
	sfac = sqrt( fac1 ) ;
	//  
	//  If they don't intersect return
	//  Trick to solve equation of intersection of two circles
	//  rotate coordinates to have both circles with centers on x axis
	//  pretty simple system of equations, then rotate back
	//  
	if ( fabs(sfac-rc) > r || fabs(sfac+rc) < r ) {
		//    printf ( "particle does not intersect \n" ) ;
		return  1 ;
	}
	//  
	//     Find intersection
	//  
	fac2   = ( r*r + fac1 - rc*rc) / (2.00 * r * sfac ) ;
	phi    = atan2(yc,xc) + q*acos(fac2) ;
	td     = atan2(r*sin(phi) - yc,r*cos(phi) - xc) ;
	/*
	double xx = x0 + rc * cos(phi);
	double yy = x0 + rc * cos(phi);
	double ttphi = atan2((yy-yc),(xx-xc));
	double ppsi = ttphi - double(q) * 0.5 * pi / fabs((double)q) ;
	*/
	//    Intersection in z

	if ( td < 0 ) td = td + 2. * pi ;
	deltat = fmod((-q*td + q*tPhi0),2*pi) ;

	// if ( deltat < 0.      ) deltat += 2. * pi ;
	// if ( deltat > 2.*pi ) deltat -= 2. * pi ;
	z = z0 + rc * tanl * deltat ;
	// 
	//    That's it
	// 
	return 0 ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Calculates intersection of track with plane define by line
**:           y = a x + b and the z 
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             a, b            - Line parameters
**:         OUT:
**:             crossPoint      - intersection point
**:
**: RETURNS:    0=ok, <>0 track does not cross the plane
**:>------------------------------------------------------------------*/
int Track::intersectorZLine ( double a, double b, ThreeDPoint& cross ) {
	//
	//   Calculate circle center and radius
	//
	double x0    = r0 * cos(phi0) ;
	double y0    = r0 * sin(phi0) ;
	double trackPhi0 = psi + q * 0.5 * pi / fabs((double)q) ;
	double rc   = pt  / ( bFactor * bField )  ;
	double xc   = x0 - rc * cos(trackPhi0) ;
	double yc   = y0 - rc * sin(trackPhi0) ;

	double ycPrime = yc - b ;
	double aa = ( 1. + a * a ) ;
	double bb = -2. * ( xc + a * ycPrime ) ;
	double cc = ( xc * xc + ycPrime * ycPrime - rc * rc ) ;

	double racine = bb * bb - 4. * aa * cc ;
	if ( racine < 0 ) return 1 ;
	double rootRacine = sqrt(racine) ;

	double oneOverA = 1./aa;
	//
	//   First solution
	//
	double x1 = 0.5 * oneOverA * ( -1. * bb + rootRacine ) ; 
	double y1 = a * x1 + b ;
	double r1 = sqrt(x1*x1+y1*y1);
	//
	//   Second solution
	//
	double x2 = 0.5 * oneOverA * ( -1. * bb - rootRacine ) ; 
	double y2 = a * x2 + b ;
	double r2 = sqrt(x2*x2+y2*y2);
	//
	//    Choose close to (0,0) 
	//
	double xHit ;
	double yHit ;
	if ( r1 < r2 ) {
		xHit = x1 ;
		yHit = y1 ;
	}
	else {
		xHit = x2 ;
		yHit = y2 ;
	}
	//-------------------------------------------------------------------
	//     Get the z coordinate
	//-------------------------------------------------------------------
	double angle  = atan2 ( (yHit-yc), (xHit-xc) ) ;
	if ( angle < 0. ) angle = angle + 2.0 * pi ;
	//   printf ( " angle %f trackPhi0 %f \n ", angle, trackPhi0 ) ;
	double dangle = angle  - trackPhi0  ;
	dangle = fmod ( dangle, 2.0 * pi ) ;
	if ( (q * dangle) > 0 ) dangle = dangle - q * 2. * pi  ;

	double stot   = fabs(dangle) * rc ;
	double zHit   = z0 + stot * tanl ;
	//
	cross.set ( xHit, yHit, zHit ) ;
	//
	return 0 ;
}
//
/*:>--------------------------------------------------------------------
**: METHOD:   Calculates trajectory length between two points on a track
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             track           - Track pointer
**:             x1, y1          - Point 1
**:             x2, y2          - Point 2
**:         OUT:
**:
**: RETURNS:    0=ok, <>0 error
**:>------------------------------------------------------------------*/
double Track::arcLength ( double x1, double y1, double x2, double y2 )
{
	double x0, y0, xc, yc, rc ;
	double angle_1, angle_2, d_angle, sleng_xy, sleng ;
	/*----------------------------------------------------------
	Get track parameters
	----------------------------------------------------------*/

	x0   = r0 * cos(phi0) ;
	y0   = r0 * sin(phi0) ;
	rc   = pt / ( bFactor * getPara()->bField )  ;
	double tPhi0 = psi + double(q) * 0.5 * pi / fabs((double)q) ;
	xc   = x0 - rc * cos(tPhi0) ;
	yc   = y0 - rc * sin(tPhi0) ;
	/*
	Get angle difference 
	*/
	angle_1  = atan2 ( (y1-yc), (x1-xc) ) ;
	if ( angle_1 < 0 ) angle_1 = angle_1 + 2. * pi ;
	angle_2  = atan2 ( (y2-yc), (x2-xc) ) ;
	if ( angle_2 < 0 ) angle_2 = angle_2 + 2. * pi ;
	d_angle  = double(q) * ( angle_1 - angle_2 ) ;
	d_angle  = fmod ( d_angle, 2. * pi ) ;
	if ( d_angle < 0 ) d_angle = d_angle + 2. * pi ;
	/*----------------------------------------------------------
	Get total angle and total trajectory
	----------------------------------------------------------*/
	sleng_xy = fabs ( rc ) * d_angle ;
	sleng    = sleng_xy * sqrt ( 1.0 + tanl * tanl )  ;
	return sleng ;

}
/*:>--------------------------------------------------------------------
**: METHOD:   Phi rotates the track
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             deltaPhi        - Angle to rotate in phi
**:
**: RETURNS:    0=ok, <>0 error
**:>------------------------------------------------------------------*/
int Track::phiRotate ( double deltaPhi ) {
	phi0 += deltaPhi ;
	if ( phi0 > 2. * pi ) phi0 -= 2. * pi ;
	if ( phi0 <         0 ) phi0 += 2. * pi ;
	psi  += deltaPhi ;
	if ( psi > 2. * pi ) psi -= 2. * pi ;
	if ( psi <         0 ) psi += 2. * pi ;

	return 0 ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Updates track parameters to point of intersection with 
**:           cylinder of radius r
**:
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             radius         - Cylinder radius to extrapolate track
**:         OUT:
**:
**:>------------------------------------------------------------------*/
void Track::updateToRadius  ( double radius ) {

	double phiExtra, zExtra, rCircle, xCircleCenter, yCircleCenter ;

	int ok = extraRCyl ( radius,  phiExtra, zExtra, rCircle, xCircleCenter, yCircleCenter ) ;
	if ( ok ) {
		//    printf ( "Track::updateToRadius: track %d does not intersect radius %f\n", 
		//              id, radius ) ;
		return ;
	}

	double xExtra = radius * cos(phiExtra) ;
	double yExtra = radius * sin(phiExtra) ;

	double tPhi = atan2(yExtra-yCircleCenter,xExtra-xCircleCenter);

	// if ( tPhi < 0 ) tPhi += 2. * pi ;

	double tPsi = tPhi - double(q) * 0.5 * pi / fabs((double)q) ;
	if ( tPsi > 2. * pi ) tPsi -= 2. * pi ;
	if ( tPsi < 0.        ) tPsi += 2. * pi ;
	//
	//    Update track parameters
	//
	r0   = radius ;
	phi0 = phiExtra ;
	z0   = zExtra ;
	psi  = tPsi ;
}
/*:>--------------------------------------------------------------------
**: METHOD:   Updates track parameters to point of closest approach
**:
**:
**: AUTHOR:     ppy - P.P. Yepes,  yepes@rice.edu
**: ARGUMENTS:
**:          IN:
**:             xBeam          - x Beam axis
**:             yBeam          - y Beam axis
**:
**:>------------------------------------------------------------------*/
void Track::updateToClosestApproach ( double xBeam, double yBeam ) {
	double rc, xc, yc ;
	ThreeDPoint closest = getClosest ( xBeam, yBeam, rc, xc, yc ) ;
	//
	double tPhi = atan2(closest.y-yc,closest.x-xc);

	// if ( tPhi < 0 ) tPhi += 2. * pi ;

	double tPsi = tPhi - double(q) * 0.5 * pi / fabs((double)q) ;
	if ( tPsi > 2. * pi ) tPsi -= 2. * pi ;
	if ( tPsi < 0.        ) tPsi += 2. * pi ;
	//
	//   Update track parameters
	//
	r0   = sqrt(closest.x*closest.x+closest.y*closest.y) ;
	phi0 = atan2(closest.y,closest.x) ;
	if ( phi0 < 0 ) phi0 += 2. * pi ;
	z0   = closest.z ;
	psi  = tPsi ;
}


// part 2
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Add hits to track
// Arguments:
//        thisHit:  hit pointer
//        way     :  >0 add at beginning, <0 at end
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Track::add ( Hit* thisHit, int way )
{
	//
	//      Increment # hits in this track
	//
	nHits++ ; 
	//
	//         Update pointers
	//
	if ( way < 0 || nHits == 1 ) {
		if ( nHits > 1 ) (lastHit)->nextTrackHit = thisHit ;
		lastHit = thisHit ;
		innerMostRow = (lastHit)->row ;
		xLastHit = (lastHit)->x ;
		yLastHit = (lastHit)->y ;
	}
	else {
		(thisHit)->nextTrackHit = firstHit ; 
		firstHit = thisHit ;
		outerMostRow = (firstHit)->row ;
	}
	//
	//        Declare hit as used and fill chi2
	//
	thisHit->setTrack ( this ) ;
	//
	//    Check whether a fit update is needed
	//
	if ( nHits < getPara()->minHitsForFit ) return ;
	//
	//    Include hit in xy fit parameter calculation
	//
	if ( nHits > getPara()->minHitsForFit  )
	{
	s11Xy = s11Xy + thisHit->wxy ;
	s12Xy = s12Xy + thisHit->wxy * thisHit->xp ;
	s22Xy = s22Xy + thisHit->wxy * square(thisHit->xp) ;
	g1Xy  = g1Xy  + thisHit->wxy * thisHit->yp ;
	g2Xy  = g2Xy  + thisHit->wxy * thisHit->xp * thisHit->yp ;



		ddXy  = s11Xy * s22Xy - square ( s12Xy ) ;
		if ( ddXy != 0 ) {
			a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;
			a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;
		}
		else {
			if ( getPara()->infoLevel > 0 ) {
				fprintf ( stderr, "Track:add: ddXy = 0 \n" ) ;
			}
		}
	}
	//
	//     Now in the sz plane
	//
	if ( getPara()->szFitFlag ) 
	{
		if ( nHits > getPara()->minHitsForFit ) 
		{
		s11Sz = s11Sz + thisHit->wz ;
		s12Sz = s12Sz + thisHit->wz * thisHit->s ;
		s22Sz = s22Sz + thisHit->wz * thisHit->s * thisHit->s ;
		g1Sz  = g1Sz  + thisHit->wz * thisHit->z ;
		g2Sz  = g2Sz  + thisHit->wz * thisHit->s * thisHit->z ;



			ddSz  = s11Sz * s22Sz -  s12Sz * s12Sz ;
			if ( ddSz != 0 ) {
				a1Sz  = ( g1Sz * s22Sz - g2Sz * s12Sz ) / ddSz ;
				a2Sz  = ( g2Sz * s11Sz - g1Sz * s12Sz ) / ddSz ;
			}
			else
			{
				if ( getPara()->infoLevel > 0 ) {
					fprintf ( stderr, "Track:add: ddSz = 0 \n" ) ;
				}
			}
		}
	}
}
//****************************************************************************
//   Fill track information tables
//****************************************************************************
void Track::add ( Track* piece ) 
{
	//
	//   Get circle parameters
	//
	s11Xy += piece->s11Xy  ;
	s12Xy += piece->s12Xy  ;
	s22Xy += piece->s22Xy  ;
	g1Xy  += piece->g1Xy   ;
	g2Xy  += piece->g2Xy   ;

	ddXy  =   s11Xy * s22Xy - square ( s12Xy ) ;
	a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;
	a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;
	//
	//     Now in the sz plane
	//
	if ( getPara()->szFitFlag != 0 ) {
		double det1 = s11Sz * s22Sz - s12Sz * s12Sz ;
		dtanl = (double) ( s11Sz / det1 );
		dz0   = (double) ( s22Sz / det1 );

		double det2 = piece->s11Sz * piece->s22Sz - piece->s12Sz * piece->s12Sz ;
		piece->dtanl = (double) ( piece->s11Sz / det2 );
		piece->dz0   = (double) ( piece->s22Sz / det2 );

		double weight1 = 1./(dtanl*dtanl);
		double weight2 = 1./(piece->dtanl*piece->dtanl);
		double weight  = (weight1+weight2);
		tanl = ( weight1 * tanl + weight2 * piece->tanl ) / weight ; 

		weight1 = 1./(dz0*dz0);
		weight2 = 1./(piece->dz0*piece->dz0);
		weight  = (weight1+weight2);
		z0   = ( weight1 * z0 + weight2 * piece->z0 ) / weight ; 
	}

	//
	//  Add space points to first track
	//
	int counter ;
	if ( piece->outerMostRow < outerMostRow ){
		if ( lastHit != NULL ) {
			counter = 0 ;
			for ( currentHit   = piece->firstHit ; 
				currentHit != 0 && counter < piece->nHits ;
				currentHit  = (currentHit)->nextTrackHit  ) {
					(currentHit)->track = this   ;
					counter++ ;
			}
			(lastHit)->nextTrackHit = piece->firstHit ;
			lastHit         = piece->lastHit ;
		}
		piece->firstHit = 0 ;
		innerMostRow = piece->innerMostRow ;
		xLastHit     = piece->xLastHit ;
		yLastHit     = piece->yLastHit ;
	}
	else {
		if ( piece->lastHit != NULL ) {
			counter = 0 ;
			for ( currentHit   = piece->firstHit ; 
				currentHit != 0 && counter < piece->nHits ;
				currentHit  = (currentHit)->nextTrackHit  ) {
					(currentHit)->track = this;
					counter++;
			}
			(piece->lastHit)->nextTrackHit = firstHit ;
			firstHit               = piece->firstHit ;
		}
		outerMostRow = piece->outerMostRow ;
		piece->firstHit = 0 ;
	}
	//
	//
	nHits  += piece->nHits ;
	chi2[0] += piece->chi2[0] ;
	chi2[1] += piece->chi2[1] ;
	//
	//   Update track parameters
	//
	//
	getPara()->szFitFlag = 0 ;
	if ( getPara()->fillTracks != 0) fill ( ) ;
	getPara()->szFitFlag = 1 ;
	//
	//
	//   Declare track 2 not to be used
	//
	piece->flag    = -1 ;
}

//****************************************************************************
//   Control how the track gets built
//****************************************************************************
int Track::buildTrack ( Hit* frstHit, Container* volume ) {
	//
	//   Add first hit to track
	//
	add ( frstHit, GO_DOWN ) ;
	//
	//    Try to build a segment first
	//
	if ( !segment ( volume, GO_DOWN ) ) return 0 ;
	//
	//    If segment build go for a real track with a fit
	//
	int rowToStop = getPara()->rowInnerMost ;
	if ( !follow ( volume, GO_DOWN, rowToStop ) ) return 0 ;
	//
	//    Now to extent track the other direction if requested
	//
	if ( getPara()->goBackwards != 0 ) follow ( volume, GO_UP, getPara()->rowOuterMost ) ;
	//
	//  Fill tracks
	//
	if ( getPara()->fillTracks != 0 ) fill ( ) ;
#ifdef TRDEBUG
	debugFill ( ) ;
#endif

	return 1 ;
}
//***************************************************************************
//   Calculates dEdx
//***************************************************************************
void Track::dEdx (  ){
	int i, j ;
	Hit* nextHit ;
	int nTruncate = max(1,
		getPara()->dEdxNTruncate*nHits/100) ;
	nTruncate = min(nHits/2,nTruncate) ;
	//
	//   Define array to keep largest de's
	//
	double *de = new double[nTruncate] ;
	//
	//    Reset
	//
	dedx = 0.F ;
	memset ( de, 0, nTruncate*sizeof(double) ) ;
	//
	//
	//
	for  ( nextHit = firstHit ; 
		nextHit != 0 ;
		nextHit = nextHit->nextTrackHit) { 

			dedx += nextHit->q ;

			if ( nextHit->q < de[0] ) continue ;

			for ( i = nTruncate-1 ; i>=0 ; i-- ){
				if ( nextHit->q > de[i] ){
					for ( j=0 ; j<i ; j++ ) de[j] = de[j+1] ;
					de[i] = nextHit->q ;
					break ;
				}
			}
	}
	//
	//    Subtract largest de
	//
	for ( i=0 ; i<nTruncate ; i++ ) dedx -= de[i] ;
	dedx = dedx / length ;
	/*   End track in required volume condition */

}
//*********************************************************************** 
//   Delete track candidate 
//***********************************************************************
void Track::deleteCandidate(void)
{
	Hit* curentHit = firstHit ;
	Hit* nextHit ;
#ifdef TRDEBUG
	debugDeleteCandidate ( ) ;
#endif
	while ( curentHit != 0 )
	{
		nextHit            = curentHit->nextTrackHit;
		curentHit->nextTrackHit     =  0 ;
		curentHit->xyChi2   =
			curentHit->szChi2   =  
			curentHit->s        =  0.F ;

		curentHit->setTrack ( 0 ) ;
		curentHit = nextHit;
	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Fills track variables with or without fit
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Track::fill (  ) {
	//
	//   Get circle parameters
	//
	double xc, yc ;
	double rc   = sqrt ( a2Xy * a2Xy + 1 ) / ( 2 * fabs(a1Xy) ) ;
	pt          = (double)(bFactor * getPara()->bField * rc );
	double xParameters = 0.; 
	double yParameters = 0.;
	//
	if ( pt > getPara()->ptMinHelixFit ) {
		double combinedChi2 = 0.5*(chi2[0]+chi2[1])/nHits ;
		if ( getPara()->primaries && combinedChi2 < getPara()->maxChi2Primary ) 
			getPara()->vertexConstrainedFit = 1 ;
		else 
			getPara()->vertexConstrainedFit = 0 ;

		fitHelix ( ) ;

		if ( getPara()->vertexConstrainedFit && getPara()->parameterLocation ){ 
			updateToRadius ( sqrt(xLastHit*xLastHit+yLastHit*yLastHit) ) ;
		}
		else if ( !getPara()->vertexConstrainedFit && !getPara()->parameterLocation ) {
			updateToClosestApproach ( getPara()->xVertex, getPara()->yVertex ) ;
		}
	}
	else{
		if ( getPara()->primaries ){ 
			fillPrimary ( xc, yc, rc, getPara()->xVertex, getPara()->yVertex ) ;
			if ( getPara()->parameterLocation ) {// give track parameters at inner most point
				updateToRadius ( sqrt(xLastHit*xLastHit+yLastHit*yLastHit) ) ;
			}
		}
		else { // Secondaries now
			xc = - a2Xy / ( 2. * a1Xy ) + xRefHit ;
			yc = - 1.   / ( 2. * a1Xy ) + yRefHit ;
			if ( getPara()->parameterLocation ) { // give track parameters at inner most point
				xParameters = xLastHit ; 
				yParameters = yLastHit ; 
			}
			else { // give parameters at point of closest approach
				getClosest ( getPara()->xVertex, getPara()->yVertex,
					rc, xc, yc, xParameters, yParameters ) ;
			}
			fillSecondary ( xc, yc, xParameters, yParameters ) ;
		}
		//
		//    Get Errors
		//
		if ( getPara()->getErrors ) {
			getErrorsCircleFit (  (double)xc, (double)yc, (double)rc ) ;
			double det = s11Sz * s22Sz - s12Sz * s12Sz ;
			dtanl = (double) ( s11Sz / det );
			dz0   = (double) ( s22Sz / det );
		}
	}
}
//****************************************************************************     
//     Fill track information variables
//****************************************************************************
void Track::fillPrimary (  double &xc, double &yc, double &rc,
						 double xPar, double yPar ) {
							 //
							 //   Get circle parameters
							 //
							 xc = getPara()->xVertex - a2Xy / ( 2. * a1Xy ) ;
							 yc = getPara()->yVertex - 1.   / ( 2. * a1Xy ) ;
							 //
							 //   Get track parameters
							 //
							 double angle_vertex  = atan2 ( yPar-yc, xPar-xc ) ;
							 if ( angle_vertex < 0. ) angle_vertex = angle_vertex + twoPi ;

							 double dx_last    = xLastHit - xc ;
							 double dy_last    = yLastHit - yc ;
							 double angle_last = atan2 ( dy_last, dx_last ) ;
							 if ( angle_last < 0. ) angle_last = angle_last + twoPi ;
							 //
							 //       Get the rotation
							 //
							 double d_angle = angle_last - angle_vertex ;
							 // double d_angle = angle_vertex - angle_last ;

							 // if ( d_angle >  pi ) d_angle -= twoPi  ;
							 if ( d_angle < -pi ) d_angle += twoPi  ;

							 q = ( ( d_angle < 0 ) ? 1 : -1 ) ;
							 r0   = sqrt(xPar*xPar+yPar*yPar) ;
							 phi0 = atan2(yPar,xPar) ;
							 if ( phi0 < 0 ) phi0 += 2. * pi ;
							 psi  = (double)(angle_vertex - q * 0.5F * pi) ;
							 if ( psi < 0     )  psi = (double)(psi + twoPi );
							 if ( psi > twoPi )  psi = (double)(psi - twoPi );
							 //
							 //      Get z parameters if needed       
							 //
							 if ( getPara()->szFitFlag == 1 ){
								 tanl = -(double)a2Sz ;
								 z0   =  (double)(a1Sz + a2Sz * ( length - rc * d_angle * q ) );
							 }
							 else if ( getPara()->szFitFlag == 2 ) {
								 tanl = (firstHit)->z /
									 (double)sqrt ( (firstHit)->x*(firstHit)->x + 
									 (firstHit)->y*(firstHit)->y ) ;
								 z0      = 0.F ;
							 }
							 //
							 //    Store some more track info
							 //
							 eta     = seta(1.,tanl )   ;
							 //
							 //   Set primary track
							 //
							 flag = 1 ;

}
//****************************************************************************
//   
//   Fill track information tables
//
//****************************************************************************
void Track::fillSecondary ( double &xc, double &yc,
						   double xPar, double yPar )
{
	/*--------------------------------------------------------------------------
	Get angles for initial and final points
	------------------------------------------------------------------------------*/
	double dx1    = (firstHit)->x - xc ;
	double dy1    = (firstHit)->y - yc ;
	double angle1 = atan2 ( dy1, dx1 ) ;
	if ( angle1 < 0. ) angle1 = angle1 + twoPi ;

	double dx2    = xLastHit - xc ;
	double dy2    = yLastHit - yc ;
	double angle2 = atan2 ( dy2, dx2 ) ;
	if ( angle2 < 0. ) angle2 = angle2 + twoPi ;
	/*--------------------------------------------------------------------------
	Get the rotation
	------------------------------------------------------------------------------*/
	double dangle = angle2 - angle1 ;
	//  if ( dangle >  pi ) dangle =   dangle - twoPi  ;
	if ( dangle < -pi ) dangle =   dangle + twoPi  ;

	q    = ( ( dangle > 0 ) ? 1 : -1 ) ;
	r0   = (lastHit)->r   ;
	phi0 = (lastHit)->phi ;
	psi  = (double)(angle2 - q * piHalf );
	if ( psi < 0     ) psi = (double)(psi + twoPi );
	//
	//      Get z parameters if needed       
	//
	if ( getPara()->szFitFlag ){
		tanl = -(double)a2Sz ;
		z0   =  (double)(a1Sz + a2Sz * length  );
	}
	else{
		tanl = (firstHit)->z /
			(double)sqrt ( (firstHit)->x*(firstHit)->x + 
			(firstHit)->y*(firstHit)->y ) ;
		z0      = 0.F ;
	}
	//
	//-->    Store some more track info
	//
	eta     = seta(1., tanl )   ;
	//
	//    Set primary track flag
	//
	flag = 0 ;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//        Adds hits to a track chosing the closest to fit
// Arguments:
//              volume:	      volume pointer
//              way   :       which way to procede in r (negative or positive)
//              row_to_stop:  row index where to stop
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::follow ( Container* volume, int way, int ir_stop ) {

	Hit* nextHit ;

	if ( way < 0 )
		nextHit = lastHit ;
	else
		nextHit = firstHit ; 
#ifdef TRDEBUG
	if ( getPara()->trackDebug && getPara()->debugLevel >= 2 )
		printf ( "Track::follow: ===> Going into Track extension <===\n" );
#endif
	//
	//     Define variables to keep total chi2
	//
	double xyChi2 = chi2[0] ;
	double szChi2 = chi2[1] ;

	//
	//    Loop as long as a hit is found and the segment
	//    is shorter than n_hit_segm
	//
	while ( way * nextHit->row < way * ir_stop ) {
		//
		//      Select next hit
		//
		chi2[0] = getPara()->hitChi2Cut ;

		nextHit = seekNextHit ( volume, nextHit, way*getPara()->trackRowSearchRange, USE_FOLLOW ) ;

#ifdef TRDEBUG
		if ( getPara()->trackDebug && getPara()->debugLevel >= 1 ){
			if ( nextHit != 0 ){
				printf ( "Track::follow: Search succesful, hit selected %d\n", 
					nextHit->id );
				//		    nextHit->Show ( getPara()->color_track ) ;
			}
			else{
				printf ( "Track::follow: Search unsuccesful\n" );
				if ( chi2[0]+chi2[1] > getPara()->hitChi2Cut )
				{
					printf ( " hit chi2 %f larger than cut %f ", chi2[0]+chi2[1], getPara()->hitChi2Cut ) ;
					printf ( " xyChi2 %f : szChi2 %f ", chi2[0], chi2[1]);
				}
			}
//			debugAsk () ;
		}
#endif
		//
		//    Stop if nothing found
		//
		if ( nextHit == 0 ) break ;
		//
		//   Keep total chi2
		//
		double lxyChi2 = chi2[0]-chi2[1] ;
		xyChi2 += lxyChi2 ;
		nextHit->xyChi2 = lxyChi2 ;
		//
		//   if sz fit update track length
		//
		if ( getPara()->szFitFlag  ) {
			length = nextHit->s ;
			szChi2 += chi2[1]  ;
			nextHit->szChi2 = chi2[1] ;
		}
		//
		//     Add hit to track
		//
		add ( nextHit, way ) ;

	} // End while
	//
	//    Check # hits
	//
	if ( nHits < getPara()->minHitsPerTrack ) return 0 ; 
	//
	//   Store track chi2
	//
	chi2[0] = xyChi2 ;
	chi2[1] = szChi2 ;
	//
	//        Check total chi2
	//
	double normalized_chi2 = (chi2[0]+chi2[1])/nHits ;
	if ( normalized_chi2 > getPara()->trackChi2Cut ) return 0 ;
	//
	return 1 ;
}
/*******************************************************************************
Reconstructs tracks
*********************************************************************************/
int Track::followHitSelection ( Hit* baseHit, Hit* candidateHit ){
	//
	double lszChi2 = 0 ;
	double lchi2 ;
	double slocal, deta, dphi ;
	double dx, dy, dxy, dsz, temp ;
	//
	//           Check delta eta 
	//
	//   if ( baseHit->dz < 1000. && candidateHit->dz < 1000 ){
	deta = fabs((baseHit->eta)-(candidateHit->eta)) ;
	if ( deta > getPara()->deta ) return 0 ; 
	//   }
	//   else deta = 0.F ;
	//
	//           Check delta phi
	//
	dphi = fabs((baseHit->phi)-(candidateHit->phi)) ;
	if ( dphi > getPara()->dphi && dphi < twoPi-getPara()->dphi ) return 0 ;
	//
	//      If looking for secondaries calculate conformal coordinates
	//
	if ( getPara()->primaries == 0 ){
		double xx = candidateHit->x - xRefHit ;
		double yy = candidateHit->y - yRefHit ;
		double rr = xx * xx + yy * yy ;
		candidateHit->xp =   xx / rr ;
		candidateHit->yp = - yy / rr ;

		candidateHit->wxy  = rr * rr /
			( square(getPara()->xyErrorScale)  *
			( square(candidateHit->dx) + square(candidateHit->dy) ) ) ;
	}
	//
	//      Calculate distance in x and y
	//
	temp = (a2Xy * candidateHit->xp - candidateHit->yp + a1Xy) ;
	dxy  = temp * temp / ( a2Xy * a2Xy + 1.F ) ;
	//
	//    Calculate chi2
	//
	lchi2    = (dxy * candidateHit->wxy) ;

	if ( lchi2 > chi2[0] ) return 0 ;
	//
	//      Now in the sz plane
	//
	if ( getPara()->szFitFlag != 0 ){
		//
		//        Get "s" and calculate distance hit-line
		//
		dx     = baseHit->x - candidateHit->x ;
		dy     = baseHit->y - candidateHit->y ;
		slocal = length + sqrt ( dx * dx + dy * dy ) ;

		temp = (a2Sz * slocal - candidateHit->z + a1Sz) ;
		dsz  = temp * temp / ( a2Sz * a2Sz + 1 ) ;
		//
		//              Calculate chi2
		//
		lszChi2 = dsz * candidateHit->wz ;
		lchi2 += lszChi2 ;
	} 
	else {
		lszChi2 = 0.F ;
	}
	//
	//         Check whether the chi2 square is better than previous one
	//
	if ( lchi2 < chi2[0] ) {
		chi2[0]       = (double)lchi2    ;
		chi2[1]       = (double)lszChi2 ;

		if ( getPara()->szFitFlag != 0 ) candidateHit->s = (double)slocal ;
		//
		//       if a good chi2 is found let's stop here
		//
		if ( lchi2 < getPara()->goodHitChi2 ) return 2 ;

		return 1 ;
	}
	//
	//     Return the selected hit
	//
	return 0 ;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Merges tracks
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::mergePrimary ( TrackContainer* trackArea ){
	short  track_merged ;
	register int  areaIndex ;
	int    i_phi, i_eta ;
	Track *i_track = 0 ;
	int    ip, ie ;
	double  delta_psi ;
	//
	//   Check Track is primary
	//
	if ( flag != 1 ) return 0 ;
	//-
	//   Get track area       
	//
	i_phi = (int)(( psi - getPara()->phiMinTrack ) / getPara()->phiSliceTrack + 1 );
	if ( i_phi < 0 ) {
		printf ( " Track phi index too low  %d \n", i_phi ) ;
		i_phi = 1 ;
	}
	if ( i_phi >= getPara()->nPhiTrackPlusOne ) {
		printf ( " Track phi index too high %d \n", i_phi ) ;
		i_phi = getPara()->nPhiTrack ;
	}
	//
	//     Now eta
	//
	i_eta = (int)(( eta - getPara()->etaMinTrack ) / getPara()->etaSliceTrack + 1 );
	if ( i_eta <= 0 ) {
		printf ( " Track eta index too low  %d \n", i_eta ) ;
		i_eta = 1 ;
	}
	if ( i_eta >= getPara()->nEtaTrackPlusOne ) {
		printf ( " Track eta index too high %d \n", i_eta ) ;
		i_eta = getPara()->nEtaTrack ;
	}
	//
	//     Loop around selected area
	//
	track_merged = 0 ;
	for ( ip = max(i_phi-1,1) ; ip < min(i_phi+2,getPara()->nPhiTrackPlusOne) ; ip++ ) {
		for ( ie = max(i_eta-1,1) ; ie < min(i_eta+2,getPara()->nEtaTrackPlusOne) ; ie++ ) {
			areaIndex = ip * getPara()->nEtaTrackPlusOne + ie ;
			//
			//    Loop over tracks
			//
			for ( i_track = trackArea[areaIndex].first ; 
				i_track != 0 ;
				i_track = i_track->getNextTrack()  ) {
					//
					//    Reject track if it is not good
					//
					if ( i_track->flag < 0 ) continue ; 
					//
					// Compare both tracks
					//
					//   No overlapping tracks
					short delta1 = i_track->outerMostRow - outerMostRow ;
					short delta2 = i_track->innerMostRow - innerMostRow ;
					if ( delta1 * delta2 <= 0 ) continue ;
					//
					//    Tracks close enough
					//
					if ( fabs(eta-i_track->eta) > getPara()->detaMerge ) continue ;
					delta_psi = (double)fabs(psi - i_track->psi) ;
					if ( delta_psi > getPara()->dphiMerge && delta_psi < twoPi - getPara()->dphiMerge ) continue ;

					i_track->add ( this ) ;
#ifdef TRDEBUG
					if ( getPara()->debugLevel > 1 )
						printf ( " \n Track %d merge into %d ", this->id, i_track->id ) ;
#endif
					track_merged = 1 ;
					break ;
			}
		}
	}
	//
	//->  If track not matched add it
	//
	if ( track_merged == 0 ) {
		areaIndex = i_phi * getPara()->nEtaTrackPlusOne + i_eta ;
		if ( trackArea[areaIndex].first == 0 )
			trackArea[areaIndex].first = 
			trackArea[areaIndex].last = this  ;
		else {
			(trackArea[areaIndex].last)->nxatrk = this ; 
			trackArea[areaIndex].last = this ;
		}
	}
	return track_merged ;
}
/************************************************************************* 
Recontruct primary tracks 
*************************************************************************/
void Track::reset (void)
{
	/*----------------------------------------------------------------------
	Set fit parameters to zero
	----------------------------------------------------------------------*/

	flag     = getPara()->primaries ;
	nHits    = 0 ;
	s11Xy   = 
		s12Xy   = 
		s22Xy   = 
		g1Xy    = 
		g2Xy    = 
		chi2[0]  = 0.F ;
	nxatrk   = 0 ;
	if ( getPara()->szFitFlag ) 
	{
		s11Sz =
			s12Sz =
			s22Sz =
			g1Sz  =
			g2Sz  =
			chi2[1]  = 
			length         = 0.F ;
	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     Function to look for next hit
// Input:	volume:         Volume pointer
//          baseHit:       Last point in track
//          n_r_steps:      How many rows search and which way (up or down)
//		    which_function: Function to be used to decide whether the hit is good
// Returns:	Selected hit
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Hit* Track::seekNextHit ( Container  *volume, 
						 Hit* baseHit,
						 int     n_r_steps,
						 int which_function ) {
#define N_LOOP 9 
							 int loop_eta[N_LOOP] = { 0, 0, 0,-1,-1,-1, 1, 1, 1 } ;
							 int loop_phi[N_LOOP] = { 0,-1, 1, 0,-1, 1, 0,-1, 1 };


							 int ir, irp, ipp, itp, k;
							 register int areaIndex ; 
							 int result ;

							 //-------------------------------------------------------------------------------
							 //     Calculate limits on the volume loop
							 //-----------------------------------------------------------------------------*/
							 int initialRow, way ;
							 if ( n_r_steps < 0 ) {
								 initialRow = max(1, (baseHit->row - getPara()->rowInnerMost)/getPara()->modRow);
								 n_r_steps  = min(initialRow,-n_r_steps ) ;
								 way        = -1 ;
							 }
							 else {
								 initialRow = max(1, (baseHit->row - getPara()->rowInnerMost + 2)/getPara()->modRow);
								 n_r_steps  = min((getPara()->rowOuterMost-initialRow+1),n_r_steps) ;
								 way = 1 ;
							 }


							 Hit* selected_hit  = 0 ;
							 //
							 //      Loop over modules
							 //
							 for ( ir = 0 ; ir < n_r_steps ; ir++ ){
								 irp = initialRow + way * ir ;
								 for ( k=0; k< N_LOOP; k++){ 
									 ipp = baseHit->phiIndex + loop_phi[k];
									 //
									 //--   Gymnastics if phi is closed
									 //
									 if ( ipp < 1 ) {
										 if ( getPara()->phiClosed )
											 ipp = getPara()->nPhi + ipp ;
										 else
											 continue ;
									 }
									 else if ( ipp > getPara()->nPhi ) {
										 if ( getPara()->phiClosed )
											 ipp = ipp - getPara()->nPhi ;
										 else
											 continue ;
									 }
									 //
									 //     Now get eta index
									 //
									 itp = baseHit->etaIndex + loop_eta[k];
									 if ( itp <     1      ) continue ;
									 if ( itp > getPara()->nEta ) continue ;
									 //
#ifdef TRDEBUG
									 if ( getPara()->trackDebug && getPara()->debugLevel >= 4 )
									 {
										 printf ( "Track::seekNextHit: search in row %d \n",irp ) ;
										 printf ( "base iphi %d, loop ipp %d, base ieta %d, loop itp %d \n",baseHit->phiIndex, ipp, baseHit->etaIndex, itp); 
									 }
#endif
									 //
									 //       Now loop over hits in each volume 
									 //
									 areaIndex = irp   * getPara()->nPhiEtaPlusOne + ipp * getPara()->nEtaPlusOne + itp ;
									 for ( Hit* candidateHit = volume[areaIndex].first ; 
										 candidateHit != 0 ;
										 candidateHit = candidateHit->nextVolumeHit ){
#ifdef TRDEBUG
											 debugInVolume ( baseHit, candidateHit ) ;
#endif
											 //----------------------------------------------------------------------------
											 //         Check whether the hit was used before
											 //--------------------------------------------------------------------------*/
											 if ( candidateHit->track != 0 ) continue ;
											 //--------------------------------------------------------------------------
											 //         If first points, just choose the closest hit
											 //-------------------------------------------------------------------------- */
											 if ( which_function == USE_SEGMENT ) 
												 result = segmentHitSelection ( baseHit, candidateHit ) ;
											 else 
												 result = followHitSelection  ( baseHit, candidateHit ) ;
											 //
											 //     Check result
											 //
											 if ( result > 0 ) {
												 selected_hit = candidateHit ;
												 if ( result ==2  ) goto found ; 
											 }
											 //
											 //       End hit loop  
											 //
									 }
									 //
									 //     End row loop      
									 //
								 }
								 //
								 //   End volume loop inside cone      
								 //
							 }
found: ;

							 return selected_hit ;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   Forms segments
//   Arguments:
//             volume     :    volume pointer
//             way        :    whether to go to negative or positive ir
//             row_to_stop:    row index where to stop
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::segment( Container* volume, int way ){
	//
	//   Define some variables
	//
	double dx, dy, rr ;
	Hit* nextHit ;
	//
	//   Check which way to go
	//
	if ( way < 0 ) 
		nextHit = lastHit ;
	else
		nextHit = firstHit ;
#ifdef TRDEBUG
	if ( getPara()->trackDebug && getPara()->debugLevel >= 4 )
		printf ( "Track:segment: **** Trying to form segment ****\n" );
#endif
	//
	//    Loop as long a a hit is found and the segment
	//    is shorter than n_hit_segm
	//
	while ( nextHit != 0 && nHits < getPara()->nHitsForSegment ) {
		chi2[0] = getPara()->maxDistanceSegment ; ;
		nextHit = seekNextHit ( volume, nextHit, way*getPara()->segmentRowSearchRange, 
			USE_SEGMENT ) ;
#ifdef TRDEBUG
		if ( getPara()->trackDebug && getPara()->debugLevel > 0 ) {
			if ( nextHit != 0 ) {
				printf ( "Track::segment: Search succesful, hit %d selected\n",
					nextHit->id );
				//       nextHit->Show ( getPara()->color_track ) ;
			}
			else
				printf ( "Track::segment: Search unsuccesful\n" );
//			debugAsk () ;
		}
#endif
		//
		//     If sz fit update s
		//
		if ( nextHit != 0 ){
			//
			//   Calculate track length if sz plane considered
			//
			if ( getPara()->szFitFlag  ){
				dx = (nextHit)->x - (lastHit)->x ;
				dy = (nextHit)->y - (lastHit)->y ;
				length    += (double)sqrt ( dx * dx + dy * dy ) ;
				nextHit->s      = length ;
			}
			//
			//   Calculate conformal coordinates
			//
			if ( getPara()->primaries == 0 ){
				rr = square ( xRefHit - nextHit->x ) +
					square ( yRefHit - nextHit->y ) ;


				nextHit->xp    =   ( nextHit->x - xRefHit ) / rr ;
				nextHit->yp    = - ( nextHit->y - yRefHit ) / rr ;
				nextHit->wxy   = rr * rr / ( square(getPara()->xyErrorScale)  *
					square(nextHit->dx) + square(nextHit->dy) ) ;
			}
			//
			//     Add hit to track
			//
			add ( nextHit, way ) ;
		}
	} // End while ( lastHit ...
	//
	//    If number of hits is as expected return 1 
	//
	if ( nHits == getPara()->nHitsForSegment )
		return 1 ;
	else
		return 0 ;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     Routine to look for segments.
//	 Arguments:
//	 baseHit:       Hit from which track is being extrapolated
//   candidateHit:  Hit being examined as a candidate to which extend track
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int Track::segmentHitSelection ( Hit* baseHit, Hit* candidateHit ){

	double dx, dy, dr, d3, dangle ;
	double dphi, deta ;
	double   angle ;
	//
	//   select hit with the
	//   the smallest value of d3 (defined below)
	//
	dphi  = (double)fabs((baseHit->phi) - (candidateHit->phi)) ; 
	if ( dphi > pi ) dphi = (double)fabs( twoPi - dphi ) ;
	if ( dphi > getPara()->dphi && dphi < twoPi -getPara()->dphi ) return 0 ;
	//
	//    Make sure we want to look at the difference in eta
	//
	if ( baseHit->dz < 1000. && candidateHit->dz < 1000. ){
		deta  = (double)fabs((baseHit->eta) - (candidateHit->eta)) ; 
		if ( deta > getPara()->deta ) return 0 ;
	}
	else deta = 0.F ;

	dr    = (double)fabs((double)(baseHit->row - candidateHit->row));
	d3    = (double)(toDeg * dr * ( dphi  + deta ) ) ;
	//
	//     If initial segment is longer than 2 store angle info in 
	//     a1Xy and a1_sz
	//
	if ( getPara()->nHitsForSegment > 2 && nHits-1 < getPara()->nHitsForSegment ) {
		dx = candidateHit->x - baseHit->x ;
		dy = candidateHit->y - baseHit->y ;
		angle = (double)atan2 ( dy, dx ) ;
		if ( angle < 0  ) angle = angle + twoPi ;
		lastXyAngle = angle ;
	}
#ifdef TRDEBUG
	if ( getPara()->trackDebug && getPara()->debugLevel >= 3 ) {
		printf ( "Track::segmentHitSelection:\n");
		printf ( "dr,dphi,deta,distance, Min distance  %7.2f %7.2f %7.2f %7.2f %7.2f\n",
			dr,dphi,deta,d3,chi2[0] ) ;
		if ( d3 < chi2[0] )
			printf ( "Best point, keep it !!!\n" );  
		else{
			printf ( "Worse than previous, reject !!\n" );
			//       candidateHit->Show ( getPara()->color_transparent );
		}
//		debugAsk() ;
	}
#endif
	if ( d3 < chi2[0] ) {
		//
		//   For second hit onwards check the difference in angle 
		//   between the last two track segments
		//
		if ( nHits > 1 ) {
			dx     = candidateHit->x - baseHit->x ;
			dy     = candidateHit->y - baseHit->y ;
			angle  = (double)atan2 ( dy, dx ) ;
			if ( angle < 0  ) angle = angle + twoPi ;
			dangle = (double)fabs ( lastXyAngle - angle );
			lastXyAngle = angle ;
			if ( dangle > getPara()->segmentMaxAngle ) return 0 ;
		}
		//
		//    Check whether this is the "closest" hit
		//
		chi2[0]          = d3 ;
		if ( d3 < getPara()->goodDistance ) return 2 ;
		return 1 ;
	}
	//
	//    If hit does not fulfill criterai return 0
	//
	return 0 ;
}
#ifdef TRDEBUG
//*****************************************************************************  
//    Ask for a character to keep going
//******************************************************************************/
void Track::debugAsk (void) 
{
	char cc;

	printf ( "stop(s), continue (any other key)\n" );
	cc = getchar();
	if ( cc == 's' ) getPara()->trackDebug = 0 ;
	if ( cc == '1' ) getPara()->debugLevel = 1 ;
	if ( cc == '2' ) getPara()->debugLevel = 2 ;
	if ( cc == '3' ) getPara()->debugLevel = 3 ;
	if ( cc == '4' ) getPara()->debugLevel = 4 ;
	if ( cc == '5' ) getPara()->debugLevel = 5 ;
	if ( cc == '6' ) getPara()->debugLevel = 6 ;
	if ( cc == '7' ) getPara()->debugLevel = 7 ;

	printf ( "Track::debugAsk: Debug Level %d\n", getPara()->debugLevel ) ;

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//    Debug Delete Candidate
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Track::debugDeleteCandidate(void)
{
	if ( getPara()->trackDebug == 0 || getPara()->debugLevel < 1 ) return ;

	//  for ( startLoop() ; done() ; nextHit() ) {
	//  currentHit->Show ( getPara()->color_back ) ;
	//  }
	printf ( "Track::debugDeleteCandidate: Track %d has %d hits <==\n"
		,id, nHits );
	printf ( "Track::debugDeleteCandidate: Minimum is %d, delete it \n",
		getPara()->minHitsPerTrack  );
	//  print ( 31 ) ;
//	debugAsk () ;
}
/*****************************************************************************  
Fill track information tables
******************************************************************************/
void Track::debugFill (  )
{
	if ( getPara()->trackDebug && getPara()->debugLevel >= 1 ) {
		printf ( "\n ===> Track %d added  <=== ",id+1 );
		Print ( 31 ) ;
	}
}
/*****************************************************************************
Reconstructs tracks
******************************************************************************/
void Track::debugFollowCandidate ( Hit* candidateHit )
{
	if ( !getPara()->trackDebug || getPara()->debugLevel >= 4 ) return ;
	//
	//    Show the whole track and fit
	//
	//for ( startLoop() ; done() ; nextHit() ) {
	//  currentHit->Show ( getPara()->color_track ) ;
	//}
	//  candidateHit->Show ( getPara()->color_candidate );
	//
	//        Print relevant information
	//
	printf ( "Track::debugFollowCandidate ===> Extension in Follow <===\n" ) ;
	// print ( 31 ) ;

	printf ( "Track::debugFollowCandidate: Try hit %d\n", candidateHit->id ) ;
	candidateHit->print ( 11 ) ;
	//
	//     If the hit already used say it and forget about it
	//
	if ( candidateHit->track != 0 )
	{
		printf ( "Track::debugFollowCandidate: hit %d used in track %d\n",
			candidateHit->id, id );
//		debugAsk () ;
		//    candidateHit->Show ( getPara()->color_candidate ) ;
		candidateHit->print ( 3 ) ;
	}
}
/*******************************************************************************
Reconstructs tracks
*********************************************************************************/
void Track::debugFollowSuccess ( double dxy,      double dsz, double lxyChi2, 
								double  lszChi2, double chi2_min,
								Hit* candidateHit ) {
									//
									//     Check whether track needs to be debugged
									//
									if ( !getPara()->trackDebug     ) return ;
									if (  getPara()->debugLevel < 2 ) return ;
									//
									//      Show first level of info
									//
									double lchi2 = lxyChi2 + lszChi2 ;

									printf ( " \n ------------------------------------- " ) ;
									if ( lchi2 < chi2_min ){
										printf ( "Track::debugFollowSuccess: %f Best Chi2, keep point !!!\n", 
											lchi2 );
										if ( lchi2 < getPara()->goodHitChi2 ){
											printf ( "This Chi2 is better than the good cut %f\n",
												lchi2, getPara()->goodHitChi2 );
											printf ( "Stop search !!! " );
										}
									}
									else{
										printf ( "Track::debugFollowSuccess: Hit %d worse than previous, forget it !! ",
											candidateHit->id );
										//    candidateHit->Show ( getPara()->color_track ) ;
									}


									printf ( " \n ------------------------------------- " ) ;
									//
									//   Show second level of info
									//
									if ( getPara()->debugLevel > 2 ) {
										printf ( "Track::debugFollowSuccess:\n");
										printf ( "dis_xy dis_sz   %7.2e %7.2e\n ", dxy, dsz );
										printf ( "Error xy   sz   %7.2e %7.2e\n ",  
											candidateHit->wxy, candidateHit->wz );
										printf ( "xy:a1,a2;sz:a1,a2  %7.2f %7.2f %7.2f %7.2f\n",
											a1Xy, a2Xy, a1Sz, a2Sz );
										printf ( "ch2:xy sz tot min  %7.2f %7.2f %7.2f %7.2f\n", 
											lxyChi2,lszChi2, lchi2, chi2_min );
									}
//									debugAsk() ;
									// candidateHit->Show ( getPara()->color_transparent ) ;
}
/*********************************************************************************
Routine to look for segments.
Segments are track starting chains
*******************************************************************************/
void Track::debugInVolume ( Hit* baseHit, Hit* candidateHit )
{

	if ( getPara()->trackDebug && getPara()->debugLevel >= 2 ) {
		/*----------------------------------------------------------------------------
		Show the whole segment
		----------------------------------------------------------------------------*/
		//    for ( startLoop() ; done() ; nextHit() ) {
		//       currentHit->Show ( getPara()->color_track ) ;
		//    }

		//    candidateHit->Show ( getPara()->color_candidate ) ;
		/*----------------------------------------------------------------------------
		Print relevant information
		----------------------------------------------------------------------------*/
		if ( nHits > getPara()->nHitsForSegment+1 ) Print ( 31 ) ;

		printf ( "Track:debugInVolume: Try hit %d\n", candidateHit->id ) ; 
		candidateHit->print ( 11 ) ;
		/*----------------------------------------------------------------------------
		If the hit already used say it and forget about it
		----------------------------------------------------------------------------*/
		if ( candidateHit->track != 0 ) {
			printf ( "Track:debugInVolume: hit %d used in track %d\n", 
				candidateHit->id, id+1 );
			//       candidateHit->Show ( 0 );
		}
		else {
			double dphi  = (double)fabs(baseHit->phi - candidateHit->phi) ;
			double deta ;
			if ( baseHit->dz < 1000 && candidateHit->dz < 1000 )
				deta  = (double)fabs(baseHit->eta - candidateHit->eta) ;
			else
				deta  = 0.F ;

			if ( dphi > getPara()->dphi )
				printf ( "Track:debugInVolume: Hits too far apart in phi: %f \n", 
				dphi ) ;
			if ( deta > getPara()->deta )
				printf ( "Track:debugInVolume: Hits too far apart in eta: %f \n", 
				deta ) ;
		}
//		debugAsk () ;
	}
}
/*****************************************************************************
Fill track information tables
******************************************************************************/
void Track::debugNew (  )
{

	if ( firstHit->id == getPara()->hitDebug ) getPara()->trackDebug = 1 ;
	if ( getPara()->trackDebug && getPara()->debugLevel >= 1 )
	{
		printf ( "================================================ \n" );
		printf ( "Track::debugNew:Starting track %d from point %d\n", 
			id, firstHit->id );
		printf ( "================================================ \n" );

		//   firstHit->Show ( getPara()->color_track ) ;
//		debugAsk () ;
	}
}
#endif
