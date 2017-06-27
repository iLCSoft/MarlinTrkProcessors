
#include "TPCModularEndplate.h"

#include <DD4hep/DD4hepUnits.h>

#include <math.h>


TPCModularEndplate::TPCModularEndplate( const dd4hep::rec::FixedPadSizeTPCData* tpc ) : _tpc( tpc) { }


void TPCModularEndplate::addModuleRing( unsigned nModules, double phi0 ){

  if( isInitialized ){
    throw std::runtime_error("TPCModularEndplate: addModuleRing() called after initialize() ") ; 
  }
  
  _moduleRings.push_back( { nModules, phi0, 0 , 0 , 0 } )  ;
}


void TPCModularEndplate::initialize(){

  unsigned nRing = _moduleRings.size() ;

  double deltaR = ( _tpc->rMaxReadout  - _tpc->rMinReadout ) /  nRing ;

  double rMin = _tpc->rMinReadout ;

  double rMax = rMin + deltaR ;

  for( auto modRing : _moduleRings){

    modRing.rMin = rMin ;
    modRing.rMax = rMax ;

    modRing.deltaPhi =  2. * M_PI / modRing.nModules ;  
    
    rMin += deltaR ;
    rMax += deltaR ;
  }

  isInitialized = true ; 
}



double TPCModularEndplate::computeDistanceRPhi(const dd4hep::rec::Vector3D& hit){

  if( !isInitialized ){
    throw std::runtime_error("TPCModularEndplate: computeDistanceRPhi() called before initialize() ") ; 
  }
  
  unsigned indexR = _moduleRings.size() * ( hit.r()  - _tpc->rMinReadout )  / ( _tpc->rMaxReadout  - _tpc->rMinReadout )  ;

  auto modRing = _moduleRings.at( indexR )  ;

  double phi = hit.phi() - modRing.phi0   ;

  double deltaPhi  = std::fabs( std::fmod(  ( phi - modRing.phi0 ) , modRing.deltaPhi ) );

  if( deltaPhi > modRing.deltaPhi/2. ) {

    deltaPhi = modRing.deltaPhi - deltaPhi ;
  }  

  return hit.r() * deltaPhi ;

}
  
  
