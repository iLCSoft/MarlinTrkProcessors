/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SimplePlanarDigiProcessor.h"
#include <iostream>
#include <climits>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>

//FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
#include "MarlinTrk/Factory.h"

#include "gear/gearsurf/MeasurementSurfaceStore.h"

#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "marlin/ProcessorEventSeeder.h"

#include "UTIL/ILDConf.h"

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

SimplePlanarDigiProcessor aSimplePlanarDigiProcessor ;


SimplePlanarDigiProcessor::SimplePlanarDigiProcessor() : Processor("SimplePlanarDigiProcessor") {
  
  // modify processor description
  _description = "SimplePlanarDigiProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. The plannar geometry should be either VXD, SIT or SET described using ZPlannarLayout" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  registerProcessorParameter( "PointResolutionRPhi" ,
                             "R-Phi Resolution"  ,
                             _pointResU ,
                             float(0.0040)) ;
  
  registerProcessorParameter( "PointResolutionZ" , 
                             "Z Resolution" ,
                             _pointResV ,
                             float(0.0040));
  
  registerProcessorParameter( "Ladder_Number_encoded_in_cellID" , 
                             "Mokka has encoded the ladder number in the cellID" ,
                             _ladder_Number_encoded_in_cellID ,
                             bool(true));
  
  registerProcessorParameter( "Sub_Detector_ID" , 
                             "ID of Sub-Detector from UTIL/ILDConf.h from lcio. Either for VXD, SIT or SET. Only used if Ladder_Number_encoded_in_cellID == false" ,
                             _sub_det_id ,
                             int(lcio::ILDDetID::VXD));
  
  
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
                           std::string("VTXTrackerHits") ) ;
  
  registerOutputCollection(LCIO::LCRELATION,
                           "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection",
                           _outRelColName,
                           std::string("VTXTrackerHitRelations"));
  
  
  // setup the list of supported detectors
  
  
}


void SimplePlanarDigiProcessor::init() { 
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
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
  delete trksystem;
  
}

void SimplePlanarDigiProcessor::processRunHeader( LCRunHeader* run) { 
  ++_nRun ;
} 

void SimplePlanarDigiProcessor::processEvent( LCEvent * evt ) { 
  
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
    
    unsigned nCreatedHits=0;
    unsigned nDismissedHits=0;
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
    
    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);

    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;

    
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::ILDCellID0::encoder_string , trkhitVec ) ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
      
    //get geometry info
        
    const gear::ZPlanarParameters* gearDet = NULL ;
    
    int det_id = 0 ;
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
    
    if ( nSimHits>0 ) {
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( 0 ) ) ;
      if (_ladder_Number_encoded_in_cellID) {
        encoder.setValue(SimTHit->getCellID0()) ;
        det_id  = encoder[lcio::ILDCellID0::subdet] ;
      }
      else {
        det_id = _sub_det_id;
      }
    }
    
    
    //***********************************************************
    // Check which 
    //***********************************************************   
    
    if( det_id == lcio::ILDDetID::VXD ) {      
      gearDet = &(Global::GEAR->getVXDParameters()) ;
    }
    else if(det_id == lcio::ILDDetID::SIT) {
      gearDet = &(Global::GEAR->getSITParameters()) ;
    }
    else if(det_id == lcio::ILDDetID::SET) {
      gearDet = &(Global::GEAR->getSETParameters()) ;
    }
    else{
      std::stringstream errorMsg;
      errorMsg << "SimplePlanarDigiProcessor::processEvent: unsupported detector ID = " << det_id << ": file " << __FILE__ << " line " << __LINE__ ;
      throw Exception( errorMsg.str() );  
    }
    
    //smearing
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName 
    << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    
    const gear::ZPlanarLayerLayout& layerLayout = gearDet->getZPlanarLayerLayout() ;
    
    
    for(int i=0; i< nSimHits; ++i){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      //**************************************************************
      // the first issue is to deal with the cellID 
      // whether it is set using the encoding described in ILDConf.h 
      // or alla LoI i.e. ( layer = celId  - 1 )
      //**************************************************************
      
      const int celId = SimTHit->getCellID0() ;
  
      int layerNumber = 0 ;
      int ladderNumber = 0 ;
      int side = 0 ;
      int sensor = 0 ;

      encoder.setValue(celId) ;
      
      if(_ladder_Number_encoded_in_cellID) {

        streamlog_out( DEBUG3 ) << "Get Layer Number using Standard ILD Encoding from ILDConf.h : celId = " << celId << std::endl ;
      
        encoder.setValue(celId) ;
        layerNumber  = encoder[lcio::ILDCellID0::layer] ;
        ladderNumber = encoder[lcio::ILDCellID0::module] ;
        side = encoder[lcio::ILDCellID0::side] ;
        sensor = encoder[lcio::ILDCellID0::sensor] ;
        
      }
      else{
        streamlog_out( DEBUG3 ) << "Get Layer Number using celId - 1 : celId : " << celId << std::endl ;
        layerNumber = celId  - 1 ;
        ladderNumber = this->getLadderNumber(gearDet, SimTHit, ladderNumber); 
        side = 0 ; // cannot be discerned unless set alla ILDConf.h
        sensor = 0 ; // cannot be discerned unless set alla ILDConf.h
      }
      
      cellid_encoder[ lcio::ILDCellID0::subdet ] = det_id ;
      cellid_encoder[ lcio::ILDCellID0::side   ] = side ;
      cellid_encoder[ lcio::ILDCellID0::layer  ] = layerNumber ;
      cellid_encoder[ lcio::ILDCellID0::module ] = ladderNumber ;
      cellid_encoder[ lcio::ILDCellID0::sensor ] = sensor ;

      streamlog_out( DEBUG3 ) << "det_id = " <<  det_id << std::endl ;
//      streamlog_out( DEBUG3 ) << "layerNumber = " <<  layerNumber << std::endl ;
//      streamlog_out( DEBUG3 ) << "ladderNumber = " << ladderNumber << std::endl ;
//      streamlog_out( DEBUG3 ) << "side = " << side << std::endl ;
//      streamlog_out( DEBUG3 ) << "sensor = " << sensor << std::endl ;

      
      if (layerNumber > layerLayout.getNLayers()) {
        streamlog_out( ERROR ) << "Layer Number " << layerNumber << " greater than that in Gear File: exit(1) called from " << __FILE__ << " line " << __LINE__ << std::endl ;
        exit(1);
      }

      
      if (ladderNumber > layerLayout.getNLadders(layerNumber)) {
        streamlog_out( ERROR ) << "Ladder Number " << ladderNumber << " greater than that in Gear File for layer "<< layerNumber << " : exit(1) called from " << __FILE__ << " line " << __LINE__ << std::endl ;
        exit(1);
      }

      
      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************
      
      gear::MeasurementSurface const* ms = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( cellid_encoder.lowWord() );
      
      
      //**********************************************
      // Now get the geometric information from gear 
      //**********************************************

      double deltaPhi = ( 2 * M_PI ) / layerLayout.getNLadders(layerNumber) ;
      double sensitive_length  = layerLayout.getSensitiveLength(layerNumber) * 2.0 ; // note: gear for historical reasons uses the halflength 
      double sensitive_width  = layerLayout.getSensitiveWidth(layerNumber);
      double sensitive_offset = layerLayout.getSensitiveOffset(layerNumber);
      double ladder_r         = layerLayout.getSensitiveDistance(layerNumber);
      
      double ladderPhi = correctPhiRange( layerLayout.getPhi0( layerNumber ) + ladderNumber*deltaPhi ) ;
      double ladder_incline = correctPhiRange( (M_PI/2.0 ) + ladderPhi );

      
      
      //************************************************************
      // Now get the energy and position in global and local coord.
      //************************************************************

      float edep = SimTHit->getEDep() ;
      const double *pos ;
      pos =  SimTHit->getPosition() ;        
      double smearedPos[3];
    
      CLHEP::Hep3Vector globalPoint(pos[0],pos[1],pos[2]);
      CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
      CLHEP::Hep3Vector localPointSmeared = localPoint;

      double PhiInLocal = globalPoint.phi() - ladderPhi;
      
      double u = localPoint[0];
      double v = localPoint[1];
      
      
      //************************************************************
      // Debug Print Out
      //************************************************************

      
      streamlog_out(DEBUG3) << "Hit = "<< i << " has celId " << celId << " layer number = " << layerNumber << " ladderNumber = " << ladderNumber << endl;
      
      streamlog_out(DEBUG3) <<"Position of hit before smearing = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " r = " << globalPoint.rho() << "  phi: " <<  globalPoint.phi()<< endl;
      
      
      streamlog_out(DEBUG2) << ":" 
      << "  layer: " << layerNumber 
      << "  ladderIndex: " << ladderNumber 
      << "  half ladder width " << sensitive_width * 0.5 
      << "  u: " <<  u
      << "  half ladder length " << sensitive_length * 0.5 
      << "  v: " <<  v
      << "  layer sensitive_offset " << sensitive_offset
      << "  layer phi0 " << layerLayout.getPhi0( layerNumber )
      << "  phi: " <<  globalPoint.phi()
      << "  PhiInLocal: " << PhiInLocal
      << "  ladderPhi: " << ladderPhi
      << "  ladder_incline: " << ladder_incline
      << "  nladders: " << layerLayout.getNLadders(layerNumber) 
      << "  ladder r: " << ladder_r
      << std::endl ;


      //************************************************************
      // Check if Hit is inside senstive 
      //************************************************************
      
      if ( ms->isLocalInBoundary(localPoint) == false ) {
        streamlog_out(DEBUG4) << "hit not in sensitive boundary hit dropped" << std::endl;        
        ++nDismissedHits;
        continue; 
      }
      
      
      //**************************************************************************
      // Try to smear the hit but ensure the hit is inside the sensitive region
      //**************************************************************************

      int  tries = 0;              
      bool accept_hit = false;
      
      double uSmear(0) ;
      double vSmear(0) ;
      
      while( tries < 100 ) {
        
        if(tries > 0) streamlog_out(DEBUG0) << "retry smearing for " << layerNumber << " " << ladderNumber << " : retries " << tries << std::endl;
        
        uSmear  = gsl_ran_gaussian(_rng, _pointResU);
        vSmear  = gsl_ran_gaussian(_rng, _pointResV);
        
        CLHEP::Hep3Vector temp_localPoint( u+uSmear, v+vSmear, localPoint[2] );
        
        if (ms->isLocalInBoundary(temp_localPoint)) {
          accept_hit =true;
          localPoint = temp_localPoint;
        
          break;
          
        }
        
        ++tries;
        
      }   
      
      if( accept_hit == false ) {
        streamlog_out(DEBUG4) << "hit could not be smeared within ladder after 100 tries: hit dropped"  << std::endl;
        ++nDismissedHits;
        continue; 
      } 

      //**************************************************************************
      // Convert back to global pos for TrackerHitPlaneImpl
      //**************************************************************************
      
      // note for real strip hits this needs to be done properly. The local v coordinate should be set to the centre of the strip
      globalPoint = ms->getCoordinateSystem()->getGlobalPoint(localPoint);
      

      for (int ipos=0; ipos<3; ++ipos) {
        smearedPos[ipos] = globalPoint[ipos];
      }
      
      streamlog_out(DEBUG3) <<"Position of hit after smearing = "<<smearedPos[0]<<" "<<smearedPos[1]<<" "<<smearedPos[2]
      << " :" 
      << "  u: " <<  localPoint[0]
      << "  v: " <<  localPoint[1]
      << "  w: " <<  localPoint[2]
      << std::endl ;

      
      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************
      

      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;
                  
      cellid_encoder.setCellID( trkHit ) ;
      
      trkHit->setPosition( smearedPos ) ;

      trkHit->setEDep( edep ) ;

      //**************************************************************************
      // Store planar variables to TrackerHitPlaneImpl
      //**************************************************************************
      
      gear::CartesianCoordinateSystem* cartesian = dynamic_cast< gear::CartesianCoordinateSystem* >( ms->getCoordinateSystem() ); 

      if ( ! cartesian ) {
        streamlog_out( ERROR ) << "dynamic_cast of  ICoordinateSystem to CartesianCoordinateSystem failed: exit(1) called from " << __FILE__ << " line " << __LINE__ << std::endl ;
        exit(1);
      }

      
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
      
      trkHit->setdU( _pointResU ) ;
      trkHit->setdV( _pointResV ) ;
      
      
      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************    
         
      LCRelationImpl* rel = new LCRelationImpl;

      rel->setFrom (trkHit);
      rel->setTo (SimTHit);
      rel->setWeight( 1.0 );
      relCol->addElement(rel);

      
      //**************************************************************************
      // Add hit to collection
      //**************************************************************************    
      
      trkhitVec->addElement( trkHit ) ; 
      
      ++nCreatedHits;
      
      streamlog_out(DEBUG3) << "-------------------------------------------------------" << std::endl;
      
    }      
    
    
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
    
    evt->addCollection( trkhitVec , _outColName ) ;
    evt->addCollection( relCol , _outRelColName ) ;
    
    streamlog_out(DEBUG4) << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed as not on sensitive element\n";
    
  }
  _nEvt ++ ;
}



void SimplePlanarDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SimplePlanarDigiProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}


double SimplePlanarDigiProcessor::correctPhiRange( double Phi ) const {
  
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

int SimplePlanarDigiProcessor::getLadderNumber(const gear::ZPlanarParameters* det, SimTrackerHit* sh, int layerNumber){
  
  int ladderNumber = INT_MAX;
  
  const gear::ZPlanarLayerLayout& layerLayout = det->getZPlanarLayerLayout() ;

  gear::Vector3D hitvec(sh->getPosition()[0],sh->getPosition()[1],sh->getPosition()[2]);
  
  //phi between each ladder
  double deltaPhi = ( 2 * M_PI ) / layerLayout.getNLadders(layerNumber) ;

  //  double sensitive_length  = layerLayout.getSensitiveLength(layerNumber) * 2.0 ; // note: gear for historical reasons uses the halflength 
//  double sensitive_width  = layerLayout.getSensitiveWidth(layerNumber);
//  double sensitive_offset = layerLayout.getSensitiveOffset(layerNumber);
//  double ladder_r         = layerLayout.getSensitiveDistance(layerNumber);
  
  if( ! _ladder_Number_encoded_in_cellID) {
    
    for (int ic=0; ic < layerLayout.getNLadders(layerNumber); ++ic) {
      
      double ladderPhi = correctPhiRange( layerLayout.getPhi0( layerNumber ) + ic*deltaPhi ) ;
      
      double PhiInLocal = hitvec.phi() - ladderPhi;
      double RXY = hitvec.rho();
      
      streamlog_out(DEBUG) << "ladderPhi = " << ladderPhi << " PhiInLocal = "<< PhiInLocal << " RXY = " << RXY << " (RXY*cos(PhiInLocal) = " << RXY*cos(PhiInLocal) << " layerLayout.getSensitiveDistance(layerNumber) = " << layerLayout.getSensitiveDistance(layerNumber) << endl;
      
      
      // check if point is in range of ladder
      if (RXY*cos(PhiInLocal) - layerLayout.getSensitiveDistance(layerNumber) > -layerLayout.getSensitiveThickness(layerNumber) && 
          RXY*cos(PhiInLocal) - layerLayout.getSensitiveDistance(layerNumber) <  layerLayout.getSensitiveThickness(layerNumber) ) {
        ladderNumber = ic;
        break;
      }
    }
  }


  return ladderNumber ;
  
}
