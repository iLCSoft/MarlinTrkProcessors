/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SimpleDiscDigiProcessor.h"
#include <iostream>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/ILDConf.h>



#include <marlin/Global.h>
#include "marlin/ProcessorEventSeeder.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <cmath>

#include <gsl/gsl_randist.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;





SimpleDiscDigiProcessor aSimpleDiscDigiProcessor ;


SimpleDiscDigiProcessor::SimpleDiscDigiProcessor() : Processor("SimpleDiscDigiProcessor") {
  
  // processor description
  _description = "SimpleDiscDigiProcessor creates FTD TrackerHits from SimTrackerHits" ;
  
  registerProcessorParameter( "Sub_Detector_ID" , 
                             "ID of Sub-Detector using UTIL/ILDConf.h from lcio." ,
                             _sub_det_id ,
                             int(ILDDetID::FTD));
  
  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "CollectionName" , 
                          "Name of the SimTrackerHit collection"  ,
                          _inColName ,
                          std::string("FTDCollection") ) ;
  
  registerProcessorParameter( "PointResolution" ,
                             "Point Resolution"  ,
                             _pointReso ,
                             (float)0.010) ;
  
  registerOutputCollection( LCIO::TRACKERHIT,
                           "OutputCollectionName" , 
                           "Name of the TrackerHit output collection"  ,
                           _outColName ,
                           std::string("FTDTrackerHits") ) ;
  
  
  registerProcessorParameter( "keepHitsFromDeltas" ,
                             "Whether to put deltas (secondary particles) in the collection"  ,
                             _keepHitsFromDeltas ,
                             false) ;                            
  
  registerProcessorParameter( "PetalsPerDisk" ,
                             "Number of petals per disk"  ,
                             _petalsPerDisk ,
                             1 ) ;    
  
  registerProcessorParameter( "SensorsPerPetal" ,
                             "Number of Sensors Per Petal"  ,
                             _sensorsPerPetal ,
                             1 ) ;    
  
}


void SimpleDiscDigiProcessor::init() { 
  
  // usually a good idea to
  printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;
  
  //intialise random number generator 
  r = gsl_rng_alloc(gsl_rng_ranlxs2);
  Global::EVENTSEEDER->registerProcessor(this);
  
  const gear::GearParameters& pFTD = Global::GEAR->getGearParameters("FTD");
  
  _FTDZCoordinate = pFTD.getDoubleVals( "FTDZCoordinate" ) ;
  _diskInnerRadius = pFTD.getDoubleVals( "FTDInnerRadius" ) ;
  _diskOuterRadius = pFTD.getDoubleVals( "FTDOuterRadius" ) ;
  
  
}

void SimpleDiscDigiProcessor::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
} 

void SimpleDiscDigiProcessor::processEvent( LCEvent * evt ) { 
  
  gsl_rng_set( r, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  
  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _inColName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG) << "Collection " << _inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }
  
  
  if( STHcol != 0 ){    
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( ILDCellID0::encoder_string , trkhitVec ) ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName 
    << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    for(int i=0; i< nSimHits; i++){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      const int celId = SimTHit->getCellID0() ;
      streamlog_out( DEBUG2 ) << "Get Layer Number using celId - 1 : celId : " << celId << std::endl ;
      int layerNumber = abs(celId) - 1 ;
      int side = abs(celId)/celId;
      
      const double *pos ;
      pos =  SimTHit->getPosition() ;  
      gear::Vector3D hitvec(pos[0],pos[1],pos[2]);
      
      if ( ( _keepHitsFromDeltas == true ) || ( hasCorrectZPos (pos[2]) == true ) ){
        
        streamlog_out(DEBUG) << "Hit = "<< i << " has celId " << celId << " layer number = " << layerNumber  << endl;
        
        streamlog_out(DEBUG) <<"Position of hit before smearing = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " r = " << hitvec.rho() << endl;
        
        double xSmear = gsl_ran_gaussian(r,_pointReso);
        double ySmear = gsl_ran_gaussian(r,_pointReso);
        
        double smearedPos[3] ;
        smearedPos[0] = pos[0] + xSmear;
        smearedPos[1] = pos[1] + ySmear;
        // No semaring of Z coordinate
        smearedPos[2] = pos[2] ;
        
        streamlog_out(DEBUG) <<"Position of hit after smearing = "<<smearedPos[0]<<" "<<smearedPos[1]<<" "<<smearedPos[2] << std::endl ;
        
        
        //store hit variables
        TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;        
        
        trkHit->setType( 200+abs(celId));  // needed for FullLDCTracking et al.
        
        
        int petalNumber = getPetalNumber( layerNumber , smearedPos[0] , smearedPos[1] );
        int sensorNumber = getSensorNumber( layerNumber , smearedPos[0] , smearedPos[1] );
        
        
        cellid_encoder[ ILDCellID0::subdet ] = _sub_det_id ;
        cellid_encoder[ ILDCellID0::side   ] = side ;
        cellid_encoder[ ILDCellID0::layer  ] = layerNumber ;
        cellid_encoder[ ILDCellID0::module ] = petalNumber ;
        cellid_encoder[ ILDCellID0::sensor ] = sensorNumber ;
        
        cellid_encoder.setValue( lcio::long64(cellid_encoder.lowWord() ) << 32 );
        
        cellid_encoder[ ILDCellID0::subdet ] = _sub_det_id ;
        cellid_encoder[ ILDCellID0::side   ] = side ;
        cellid_encoder[ ILDCellID0::layer  ] = layerNumber ;
        cellid_encoder[ ILDCellID0::module ] = 0 ;
        cellid_encoder[ ILDCellID0::sensor ] = 0 ;
        
        cellid_encoder.setCellID( trkHit ) ;
        
        streamlog_out(DEBUG2) <<"side = "<< side << std::endl ;
        streamlog_out(DEBUG2) <<"layerNumber = "<< layerNumber << std::endl ;
        streamlog_out(DEBUG2) <<"moduleNumber = "<< petalNumber << std::endl ;
        streamlog_out(DEBUG2) <<"sensorNumber = "<< sensorNumber << std::endl ;
        
        
        trkHit->setPosition(  smearedPos  ) ;
        
        float u_direction[2] ; // x
        u_direction[0] = 0.0 ; 
        u_direction[1] = M_PI/2.0 ;
        
        float v_direction[2] ; // y
        v_direction[0] = M_PI/2.0 ;
        v_direction[1] = M_PI/2.0 ;
        
        trkHit->setU( u_direction ) ;
        trkHit->setV( v_direction ) ;
        
        trkHit->setdU( _pointReso ) ;
        trkHit->setdV( _pointReso ) ;
        
        trkHit->setEDep( SimTHit->getEDep() ) ;
        
        MCParticle *mcp ;
        mcp = SimTHit->getMCParticle() ;
        if( mcp != 0 )  {
          trkHit->rawHits().push_back( SimTHit ) ;
        }
        else{
          streamlog_out( DEBUG0 ) << " ignore simhit pointer as MCParticle pointer is NULL ! " << std::endl ;
        }
        
        trkhitVec->addElement( trkHit ) ; 
        
        streamlog_out(DEBUG) << "-------------------------------------------------------" << std::endl;
        
      }
      else{
        
        streamlog_out(DEBUG) << "Hit "<< i << " is NOT KEPT! The z value is does not exactly correspond to a disk"  << endl;
        
        streamlog_out(DEBUG) <<"Position of hit = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " r = " << hitvec.rho() << std::endl << endl;
        
        streamlog_out(DEBUG) << "-------------------------------------------------------" << std::endl;
        
      }
      
    }
    
    evt->addCollection( trkhitVec ,  _outColName ) ;
  }
  
  _nEvt ++ ;
}



void SimpleDiscDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SimpleDiscDigiProcessor::end(){ 
  
  gsl_rng_free(r);  
  //   std::cout << "SimpleDiscDigiProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
  
}


bool SimpleDiscDigiProcessor::hasCorrectZPos ( double z ){
  
  double zPos = fabs ( z );
  
  for (unsigned int i=0; i < _FTDZCoordinate.size(); i++){
    
    if ( fabs ( _FTDZCoordinate[i] - zPos ) < 0.0001 ) return true;
    
  }
  
  return false;
  
}


int SimpleDiscDigiProcessor::getPetalNumber ( int layer , double x , double y ){
  
  
  
  
  //find out the petal (= modul)
  
  double phi = atan2 ( y , x ); // the phi angle
  if  ( phi < 0. ) phi += 2*M_PI; // to the range of 0-2Pi
  
  double phiRel = phi / 2. /M_PI; // phi in a range from 0 to 1
  
  int petal = int ( phiRel * _petalsPerDisk ); //the number of the corresponding petal
  if ( petal == _petalsPerDisk ) petal--; //just in case, the hit really is at the exact boarder of a petal
  
  
  return petal;
  
}


int SimpleDiscDigiProcessor::getSensorNumber ( int layer , double x , double y ){
  
  
  //find out the sensor
  
  double r = sqrt ( x*x + y*y ); //radius in xy plane
  
  // The relative radial position: from 0 to 1. 
  // 0.0 means at the inner radius of the petal.
  // 1.0 means at the outer radius of the petal.
  double posRel = (r - _diskInnerRadius[layer]) / ( _diskOuterRadius[layer] - _diskInnerRadius[layer] ); 
  
  
  int sensor = int ( posRel * _sensorsPerPetal ); //the number of the sensor
  if ( sensor == _sensorsPerPetal ) sensor--; //just in case, the hit really is at the exact boarder
  if ( sensor < 0) sensor = 0;  //hits might get smeared off disk and therefore get sensor = -1. this is prevented here   
  
  return sensor;   
  
  
}




