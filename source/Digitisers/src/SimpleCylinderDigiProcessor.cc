/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SimpleCylinderDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>

#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitZCylinderImpl.h>

#include <EVENT/MCParticle.h>

#include <UTIL/CellIDEncoder.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "marlin/ProcessorEventSeeder.h"

#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

SimpleCylinderDigiProcessor aSimpleCylinderDigiProcessor ;


SimpleCylinderDigiProcessor::SimpleCylinderDigiProcessor() : Processor("SimpleCylinderDigiProcessor") {
  
  // modify processor description
  _description = "SimpleCylinderDigiProcessor creates TrackerHits from SimTrackerHits, smearing in R-Phi and Z according to the input parameters. The geometry should be supplied via parameters in a Gear File" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  registerProcessorParameter( "PointResolutionRPhi" ,
                             "R-Phi Resolution"  ,
                             _pointResoRPhi ,
                             float(0.0070)) ;
  
  registerProcessorParameter( "PointResolutionZ" , 
                             "Z Resolution" ,
                             _pointResoZ ,
                             float(0.0070));
  
  registerProcessorParameter( "HitsEncodedWithCellID" , 
                             "Mokka has encoded the hits in the cellID0 according to UTIL/ILDConf.h" ,
                             _hits_encoded_with_cellID ,
                             bool(true));
  
  registerProcessorParameter( "Sub_Detector_ID" , 
                             "ID of Sub-Detector from UTIL/ILDConf.h from lcio. Either for VXD, SIT or SET. Only used if Ladder_Number_encoded_in_cellID == false" ,
                             _sub_det_id ,
                             int(lcio::ILDDetID::SIT));
 
  registerProcessorParameter( "GearParametersName" , 
                             "Name of the Gear Parameters to be used" ,
                             _gearParametersName ,
                             std::string("SIT_Simple"));
  
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "SimTrackHitCollectionName" , 
                          "Name of the Input SimTrackerHit collection"  ,
                          _inColName ,
                          std::string("SITCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHITZCYLINDER,
                           "TrackerHitCollectionName" , 
                           "Name of the TrackerHit output collection"  ,
                           _outColName ,
                           std::string("SITTrackerHits") ) ;
  
  registerOutputCollection(LCIO::LCRELATION,
                           "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection",
                           _outRelColName,
                           std::string("SITTrackerHitRelations"));
  
  
  
}


void SimpleCylinderDigiProcessor::init() { 
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  Global::EVENTSEEDER->registerProcessor(this);
  
  // Here we should read in the gear parameters so we can check that the hits are on and within the surface
  
}

void SimpleCylinderDigiProcessor::processRunHeader( LCRunHeader* ) {
  ++_nRun ;
} 

void SimpleCylinderDigiProcessor::processEvent( LCEvent * evt ) { 
  
  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG4 ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  
  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _inColName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG4) << "Collection " << _inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }
  
  if( STHcol != 0 ){    
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITZCYLINDER )  ;
    
    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);

    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;

    
    CellIDEncoder<TrackerHitZCylinderImpl> cellid_encoder( lcio::LCTrackerCellID::encoding_string() , trkhitVec ) ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
      
    //get geometry info
    
    int det_id = 0 ;
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
    
    if ( nSimHits>0 ) {
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( 0 ) ) ;
      if (_hits_encoded_with_cellID) {
        encoder.setValue(SimTHit->getCellID0()) ;
        det_id  = encoder[lcio::LCTrackerCellID::subdet()] ;
      }
      else {
        det_id = _sub_det_id;
      }
    }
    
    //smearing
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName 
    << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    
    for(int i=0; i< nSimHits; ++i){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      const int celId = SimTHit->getCellID0() ;
      
      const double *pos ;
      pos =  SimTHit->getPosition() ;  
      
      double smearedPos[3];
      
    
      int layerNumber = 0 ;
      
      encoder.setValue(celId) ;
      
      if(_hits_encoded_with_cellID) {
        streamlog_out( DEBUG3 ) << "Get Layer Number using Standard ILD Encoding from ILDConf.h : celId = " << celId << std::endl ;
      
        encoder.setValue(celId) ;
        layerNumber  = encoder[lcio::LCTrackerCellID::layer()] ;
        streamlog_out( DEBUG3 ) << "layerNumber = " <<  layerNumber << std::endl ;
        
      }
      else{
        streamlog_out( DEBUG3 ) << "Get Layer Number using celId - 1 : celId : " << celId << std::endl ;
        layerNumber = celId  - 1 ;
      }
      

      float edep = SimTHit->getEDep() ;
             
      streamlog_out(DEBUG3) << "Hit = "<< i << " has celId " << celId << " layer number = " << layerNumber << endl;
      
      CLHEP::Hep3Vector point(pos[0], pos[1], pos[2]);
      
      streamlog_out(DEBUG3) <<"Position of hit before smearing = "<< point.x() <<" "<< point.y() <<" "<< point.z()<< " r = " << point.perp() << endl;
      
      
      streamlog_out(DEBUG3) << ":" 
      << "  layer: " << layerNumber 
      << "  phi: " <<  point.phi()
      << "  R-Phi: " <<  point.phi()* point.perp()
      << std::endl ;
            
      double rPhiSmear(0) ;
      double zSmear(0) ;
            
      rPhiSmear  = gsl_ran_gaussian(_rng, _pointResoRPhi);
      zSmear  = gsl_ran_gaussian(_rng, _pointResoZ);
      
      point.setPhi( point.phi() + rPhiSmear / point.perp() );
      point.setZ( point.z() + zSmear );
      
      smearedPos[0] = point.x();
      smearedPos[1] = point.y();
      smearedPos[2] = point.z();      
      
      //store hit variables
      TrackerHitZCylinderImpl* trkHit = new TrackerHitZCylinderImpl ;
                  
      streamlog_out(DEBUG3) <<"Position of hit after smearing = "<<smearedPos[0]<<" "<<smearedPos[1]<<" "<<smearedPos[2]
      << " :" 
      << "  R-Phi: " <<  point.phi()*point.perp()
      << "  z: " <<  point.z()
      << std::endl ;
      
      
      cellid_encoder[ lcio::LCTrackerCellID::subdet() ] = det_id ;
      cellid_encoder[ lcio::LCTrackerCellID::side()   ] = 0 ;
      cellid_encoder[ lcio::LCTrackerCellID::layer()  ] = layerNumber ;
      cellid_encoder[ lcio::LCTrackerCellID::module() ] = 0 ;
      cellid_encoder[ lcio::LCTrackerCellID::sensor() ] = 0 ;
      
      cellid_encoder.setCellID( trkHit ) ;
      
      trkHit->setPosition( smearedPos ) ;
            
      trkHit->setCenter(0.0, 0.0);
      
      trkHit->setdRPhi( _pointResoRPhi ) ;
      trkHit->setdZ( _pointResoZ ) ;
      
      trkHit->setEDep( edep ) ;
      
       
      LCRelationImpl* rel = new LCRelationImpl;

      rel->setFrom (trkHit);
      rel->setTo (SimTHit);
      rel->setWeight( 1.0 );
      relCol->addElement(rel);

      trkhitVec->addElement( trkHit ) ; 
      
      streamlog_out(DEBUG3) << "-------------------------------------------------------" << std::endl;
      
    }      
    
    
    evt->addCollection( trkhitVec , _outColName ) ;
    evt->addCollection( relCol , _outRelColName ) ;
    
  }
  _nEvt ++ ;
}



void SimpleCylinderDigiProcessor::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SimpleCylinderDigiProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}


double SimpleCylinderDigiProcessor::correctPhiRange( double Phi ) const {
  
  while( (Phi < -1.0*M_PI) || (Phi > 1.0*M_PI) )
    {
    if( Phi > 1.0*M_PI )
      {
      Phi -= 2.0 * M_PI;
      }
    else
      {
      Phi += 2.0 * M_PI;
      }
    }
  
  return Phi ;
  
} // function correctPhiRange


