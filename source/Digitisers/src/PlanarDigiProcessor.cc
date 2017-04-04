/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "PlanarDigiProcessor.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/CellIDEncoder.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>


//FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
#include "MarlinTrk/Factory.h"

#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "marlin/ProcessorEventSeeder.h"


#include "UTIL/BitSet32.h"

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

PlanarDigiProcessor aPlanarDigiProcessor ;


PlanarDigiProcessor::PlanarDigiProcessor() : Processor("PlanarDigiProcessor") {
  
  // modify processor description
  _description = "PlanarDigiProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters." ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  FloatVec resUEx ;
  resUEx.push_back( 0.0040 ) ;
  
  registerProcessorParameter( "ResolutionU" ,
                              "resolution in direction of u - either one per layer or one for all layers "  ,
                              _resU ,
                              resUEx) ;
  
  FloatVec resVEx ;
  resVEx.push_back( 0.0040 ) ;

  registerProcessorParameter( "ResolutionV" , 
                              "resolution in direction of v - either one per layer or one for all layers " ,
                             _resV ,
                              resVEx );
  
  registerProcessorParameter( "IsStrip",
                              "whether hits are 1D strip hits",
                              _isStrip,
                              bool(false) );

  
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "SimTrackHitCollectionName" , 
                          "Name of the Input SimTrackerHit collection"  ,
                          _inColName ,
                          std::string("VXDCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHITPLANE,
                           "TrackerHitCollectionName" , 
                           "Name of the TrackerHit output collection"  ,
                           _outColName ,
                           std::string("VXDTrackerHits") ) ;
  
  registerOutputCollection(LCIO::LCRELATION,
                           "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection",
                           _outRelColName,
                           std::string("VTXTrackerHitRelations"));
  
 
  
}


void PlanarDigiProcessor::init() { 
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;

  
  if( _resU.size() !=  _resV.size() ) {
    
    std::stringstream ss ;
    ss << name() << "::init() - Inconsistent number of resolutions given for U and V coordinate: " 
       << "ResolutionU  :" <<   _resU.size() << " != ResolutionV : " <<  _resV.size() ;

    throw EVENT::Exception( ss.str() ) ;
  }

  
  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  Global::EVENTSEEDER->registerProcessor(this);
  
  //FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
  MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  
  if( trksystem == 0 ) {
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  trksystem->init() ;  
  
  //FIXME:SJA gear surface store has now been filled so we can dispose of the MarlinTrkSystem
  //delete trksystem;

  
}

void PlanarDigiProcessor::processRunHeader( LCRunHeader* run) { 
  ++_nRun ;
} 

void PlanarDigiProcessor::processEvent( LCEvent * evt ) { 
  
  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG4 ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  
  LCCollection* STHcol = 0 ;
  int nSimHits = 0 ;

  try{
    STHcol = evt->getCollection( _inColName ) ;
    nSimHits = STHcol->getNumberOfElements()  ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG4) << "Collection " << _inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }
  
  if( nSimHits != 0){    
    
    unsigned nCreatedHits=0;
    unsigned nDismissedHits=0;
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
    
    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);

    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;

    
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::LCTrackerCellID::encoding_string() , trkhitVec ) ;
    
    
    int det_id = 0 ;
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
    
    if ( nSimHits>0 ) {
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( 0 ) ) ;
      encoder.setValue(SimTHit->getCellID0()) ;
      det_id  = encoder[lcio::LCTrackerCellID::subdet()] ;
      
    }
    
    
    if     ( det_id == lcio::ILDDetID::VXD ){}
    else if( det_id == lcio::ILDDetID::SIT ){}
    else if( det_id == lcio::ILDDetID::SET ){}
    else if( det_id == lcio::ILDDetID::FTD ){}
    else{
      std::stringstream errorMsg;
      errorMsg << "PlanarDigiProcessor::processEvent: unsupported detector ID = " << det_id << " in collection " << _inColName  << ": file " << __FILE__ << " line " << __LINE__ ;
      throw Exception( errorMsg.str() );  
    }
    
    //smearing
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName 
    << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    
    for(int i=0; i< nSimHits; ++i){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      
      const int celId = SimTHit->getCellID0() ;
      
     
      encoder.setValue(celId) ;
      int side   = encoder[lcio::LCTrackerCellID::side()];
      int layer  = encoder[lcio::LCTrackerCellID::layer()];
      int module = encoder[lcio::LCTrackerCellID::module()];
      int sensor = encoder[lcio::LCTrackerCellID::sensor()];
      
      streamlog_out( DEBUG3 ) << "Hit = "<< i << " has celId " << encoder.valueString() << std::endl;

      //      //************************************************************
      //      // Quick check if the MCParticle is in the list of MCParticles
      //      //************************************************************
      //
      //      if (SimTHit->getMCParticle() == 0) {
      //        streamlog_out(ERROR) << " SimHit " << SimTHit << " Created by zero MCParticle which is not in the list of MCParticles: "
      //        << std::endl;
      //        continue;
      //      }
      //
      //      //************************************************************
      //      // Quick check if the MCParticle is of zero charge
      //      //************************************************************
      //
      //      if( abs( SimTHit->getMCParticle()->getCharge()) < 0.01  ){
      //        streamlog_out(ERROR) << " SimHit Created by zero charge particle: "
      //        << " Charge =  " << SimTHit->getMCParticle()->getCharge()
      //        << " EDep =  " << SimTHit->getEDep()
      //        << " x =  " << SimTHit->getPosition()[0]
      //        << " PDG =  " << SimTHit->getMCParticle()->getPDG()
      //        << " ID =  " << SimTHit->getMCParticle()->id()
      //        << std::endl;
      //        continue;
      //      }

      
      float resU = ( _resU.size() > 1 ?   _resU.at(  layer )     : _resU.at(0)   )  ;
      float resV = ( _resV.size() > 1 ?   _resV.at(  layer )     : _resV.at(0)   )  ; 

      streamlog_out( DEBUG3 ) << " --- will smear hit with resU = " << resU << " and resV = " << resV << std::endl ; 

      const double *pos ;
      pos =  SimTHit->getPosition() ;  
      
      double smearedPos[3];
      
      //      GearSurfaces::MeasurementSurface* ms = GearSurfaces::MeasurementSurfaceStore::Instance().GetMeasurementSurface( SimTHit->getCellID0() );
      
      gear::MeasurementSurface const* ms = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( encoder.lowWord() );;
      CLHEP::Hep3Vector globalPoint(pos[0],pos[1],pos[2]);
      CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
      CLHEP::Hep3Vector localPointSmeared = localPoint;
      
      
      streamlog_out(DEBUG3) <<"Position of hit before smearing global: ( "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " ) "
                            << "local: ( " << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )" << endl;
      
      // A small check, if the hit is in the boundaries:
      if ( !ms->isLocalInBoundary( localPoint ) ){
        
        streamlog_out( ERROR )<< "Hit local: ( " << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )"
                              << " is not within boundaries. Hit is skipped.\n";
        
        nDismissedHits++;
        continue;
        
      }
      
      unsigned  tries = 0;              
      
      bool accept_hit = false;
      
      // Now try to smear the hit and make sure it is still on the surface
      while( tries < 100 ) {
        
        if(tries > 0) streamlog_out(DEBUG0) << "retry smearing for side" << side << " layer"<< layer<< " module" << module << " sensor" << sensor << " : retries " << tries << std::endl;
        
        localPointSmeared.setX( localPoint.x() + gsl_ran_gaussian(_rng, resU) );
        localPointSmeared.setY( localPoint.y() + gsl_ran_gaussian(_rng, resV) );
          
        //check if hit is in boundaries
        if ( ms->isLocalInBoundary( localPointSmeared ) ){
          
          accept_hit = true;
          break;
          
        }
        
        tries++;
        
      }
      
      if( accept_hit == false ){
        
        streamlog_out(DEBUG4) << "hit could not be smeared within ladder after 100 tries: hit dropped"  << std::endl;
        continue; 
        
      }  
      
      // for 1D strip measurements: set v to 0! Only the measurement in u counts!
      if( _isStrip ) localPointSmeared[1] = 0. ;
      
      // convert back to global position for TrackerHitPlaneImpl
      CLHEP::Hep3Vector globalPointSmeared = ms->getCoordinateSystem()->getGlobalPoint(localPointSmeared);
      
      streamlog_out(DEBUG3) <<"Position of hit after smearing global: ( "  
                            << globalPointSmeared.x() <<" "<< globalPointSmeared.y() <<" "<< globalPointSmeared.z() << " ) "
                            << "local: ( "
                            << localPointSmeared.x() << " " << localPointSmeared.y() << " " << localPointSmeared.z() << " )" << endl;
      
      smearedPos[0] = globalPointSmeared.x();
      smearedPos[1] = globalPointSmeared.y();
      smearedPos[2] = globalPointSmeared.z();
      
      
      //make the TrackerHitPlaneImpl
      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;
      
      trkHit->setCellID0( encoder.lowWord() );
      
      trkHit->setPosition( smearedPos ) ;
      
      gear::CartesianCoordinateSystem* cartesian = dynamic_cast< gear::CartesianCoordinateSystem* >( ms->getCoordinateSystem() ); 
      CLHEP::Hep3Vector uVec = cartesian->getLocalXAxis();
      CLHEP::Hep3Vector vVec = cartesian->getLocalYAxis();
      
      float u_direction[2] ;
      u_direction[0] = uVec.theta();
      u_direction[1] = uVec.phi();
      
      float v_direction[2] ;
      v_direction[0] = vVec.theta();
      v_direction[1] = vVec.phi();
      
      streamlog_out(DEBUG3) 
      << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1] 
      << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
      << std::endl ;
      
      trkHit->setU( u_direction ) ;
      trkHit->setV( v_direction ) ;
      
      trkHit->setdU( resU ) ;
      
      if( _isStrip ) trkHit->setdV( 0 ); // no error in v direction for strip hits as there is no meesurement information in v direction
      else trkHit->setdV( resV ) ;
      
      if( _isStrip ){
        trkHit->setType( UTIL::set_bit( trkHit->getType() , UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) ) ;
      }
      
      trkHit->setEDep( SimTHit->getEDep() );
      
      trkhitVec->addElement( trkHit ) ; 
      
      
      // make the relation
      LCRelationImpl* rel = new LCRelationImpl;
      
      float weight = 1.0;
      
      streamlog_out(DEBUG3) <<" Set relation between "
      << " sim hit " << SimTHit 
      << " to tracker hit " << trkHit
      << " with a weight of " << weight 
      << std::endl;
      
      rel->setFrom (trkHit);
      rel->setTo (SimTHit);
      rel->setWeight( weight );
      relCol->addElement(rel);
      
      nCreatedHits++;
      streamlog_out(DEBUG3) << "-------------------------------------------------------" << std::endl;
      
    }      
    
    
    evt->addCollection( trkhitVec , _outColName ) ;
    evt->addCollection( relCol , _outRelColName ) ;
    
    streamlog_out(DEBUG4) << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits got dismissed for being out of boundary\n";
    
  }
  _nEvt ++ ;
}



void PlanarDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void PlanarDigiProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}


