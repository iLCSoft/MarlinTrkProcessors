/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SimplePlanarDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
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
                              _pointResoRPhi ,
                              float(0.0040)) ;
	
  registerProcessorParameter( "PointResolutionZ" , 
                              "Z Resolution" ,
                              _pointResoZ ,
                              float(0.0040));

  registerProcessorParameter( "Ladder_Number_encoded_in_cellID" , 
                              "Mokka has encoded the ladder number in the cellID" ,
                              _ladder_Number_encoded_in_cellID ,
                              bool(false));

  registerProcessorParameter( "Sub_Detector_ID" , 
                              "ID of Sub-Detector using UTIL/ILDConf.h from lcio. Either VXD, SIT or SET" ,
                              _sub_det_id ,
                              int(ILDDetID::VXD));


  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "SimTrackHitCollectionName" , 
                           "Name of the Input SimTrackerHit collection"  ,
                           _inColName ,
                           std::string("VXDCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHIT,
                            "TrackerHitCollectionName" , 
                            "Name of the TrackerHit output collection"  ,
                            _outColName ,
                            std::string("VTXTrackerHits") ) ;
  
 
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


}

void SimplePlanarDigiProcessor::processRunHeader( LCRunHeader* run) { 
  ++_nRun ;
} 

void SimplePlanarDigiProcessor::processEvent( LCEvent * evt ) { 

  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
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
    
    
    //get geometry info

    const gear::ZPlanarParameters* gearDet = NULL ;


    int det_id_for_type = 0 ;
    if( _sub_det_id == ILDDetID::VXD ) {      
      gearDet = &(Global::GEAR->getVXDParameters()) ;
      det_id_for_type = 100 ;
    }
    else if(_sub_det_id == ILDDetID::SIT) {
      gearDet = &(Global::GEAR->getSITParameters()) ;
      det_id_for_type = 400 ;
    }
    else if(_sub_det_id == ILDDetID::SET) {
      gearDet = &(Global::GEAR->getSETParameters()) ;
    }
    else{
      streamlog_out( DEBUG ) << "unknown detector ID" << std::endl;
    }

    //smearing
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName 
                           << " with " <<  nSimHits  << " hits ... " << std::endl ;


    const gear::ZPlanarLayerLayout& layerLayout = gearDet->getZPlanarLayerLayout() ;

    
    for(int i=0; i< nSimHits; ++i){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      const int celId = SimTHit->getCellID0() ;
      
      const double *pos ;
      pos =  SimTHit->getPosition() ;  
      
      double smearedPos[3];
      
      gear::Vector3D hitvec(pos[0],pos[1],pos[2]);
      gear::Vector3D smearedhitvec(pos[0],pos[1],pos[2]);


      int layerNumber = 0 ;
      int ladderNumber = 0 ;

      if(_ladder_Number_encoded_in_cellID) {
        streamlog_out( DEBUG2 ) << "Get Layer Number using Standard ILD Encoding from ILDConf.h : celId = " << celId << std::endl ;
        UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
        encoder.setValue(celId) ;
        layerNumber  = encoder[ILDCellID0::layer] ;
        ladderNumber = encoder[ILDCellID0::module] ;
        streamlog_out( DEBUG2 ) << "layerNumber = " <<  layerNumber << std::endl ;
        streamlog_out( DEBUG2 ) << "ladderNumber = " << ladderNumber << std::endl ;

      }
      else{
        streamlog_out( DEBUG2 ) << "Get Layer Number using celId - 1 : celId : " << celId << std::endl ;
        layerNumber = celId  - 1 ;
      }

      float edep = SimTHit->getEDep() ;
      
      //phi between each ladder
      double deltaPhi = ( 2 * M_PI ) / layerLayout.getNLadders(layerNumber) ;
      //      double sensitive_length  = layerLayout.getSensitiveLength(layerNumber);
      double sensitive_width  = layerLayout.getSensitiveWidth(layerNumber);
      double sensitive_offset = layerLayout.getSensitiveOffset(layerNumber);
      double ladder_r         = layerLayout.getSensitiveDistance(layerNumber);

      if( ! _ladder_Number_encoded_in_cellID) {
        
        for (int ic=0; ic < layerLayout.getNLadders(layerNumber); ++ic) {
          
          double ladderPhi = correctPhiRange( layerLayout.getPhi0( layerNumber ) + ic*deltaPhi ) ;

          double PhiInLocal = hitvec.phi() - ladderPhi;
          double RXY = hitvec.rho();
          
          //          streamlog_out(DEBUG) << "ladderPhi = " << ladderPhi << " PhiInLocal = "<< PhiInLocal << " RXY = " << RXY << " (RXY*cos(PhiInLocal) = " << RXY*cos(PhiInLocal) << " layerLayout.getSensitiveDistance(layerNumber) = " << layerLayout.getSensitiveDistance(layerNumber) << endl;
      

          // check if point is in range of ladder
          if (RXY*cos(PhiInLocal) - layerLayout.getSensitiveDistance(layerNumber) > -layerLayout.getSensitiveThickness(layerNumber) && 
              RXY*cos(PhiInLocal) - layerLayout.getSensitiveDistance(layerNumber) <  layerLayout.getSensitiveThickness(layerNumber) )
            {
              ladderNumber = ic;
              break;
            }
        }
        
      }

      double ladderPhi = correctPhiRange( layerLayout.getPhi0( layerNumber ) + ladderNumber*deltaPhi ) ;
      double ladder_incline = correctPhiRange( (M_PI/2.0 ) + ladderPhi );

      double PhiInLocal = hitvec.phi() - ladderPhi;

      double u = (hitvec.rho() * sin(PhiInLocal) - sensitive_offset );
      //double u = (hitvec.rho() * sin(PhiInLocal) );

      streamlog_out(DEBUG) << "Hit = "<< i << " has celId " << celId << " layer number = " << layerNumber << " ladderNumber = " << ladderNumber << endl;
      
      streamlog_out(DEBUG) <<"Position of hit before smearing = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " r = " << hitvec.rho() << endl;


      streamlog_out(DEBUG) << ":" 
                           << "  layer: " << layerNumber 
                           << "  ladderIndex: " << ladderNumber 
                           << "  half ladder width " << sensitive_width * 0.5 
                           << "  u: " <<  u
                           << "  layer sensitive_offset " << sensitive_offset
                           << "  layer phi0 " << layerLayout.getPhi0( layerNumber )
                           << "  phi: " <<  hitvec.phi()
                           << "  PhiInLocal: " << PhiInLocal
                           << "  ladderPhi: " << ladderPhi
                           << "  ladder_incline: " << ladder_incline
                           << "  nladders: " << layerLayout.getNLadders(layerNumber) 
                           << "  ladder r: " << ladder_r
                           << std::endl ;

      if( u > sensitive_width * 0.5 || u < -sensitive_width * 0.5)
        {
          streamlog_out(DEBUG) << "hit not in sensitive: u: " << u << " half ladder width = " << sensitive_width * 0.5 << std::endl;
        }

      int  tries = 0;              
      // try to smear the hit within the ladder
      bool accept_hit = false;

      double rPhiSmear(0) ;
      double zSmear(0) ;

      while( tries < 100 )
        {
          
          if(tries > 0) streamlog_out(DEBUG) << "retry smearing for " << layerNumber << " " << ladderNumber << " : retries " << tries << std::endl;
          
          rPhiSmear  = gsl_ran_gaussian(_rng, _pointResoRPhi);
          
          if( (u+rPhiSmear) < sensitive_width * 0.5 && (u+rPhiSmear) > -sensitive_width * 0.5)
            //          if( true )
            {
              accept_hit =true;
              zSmear  = gsl_ran_gaussian(_rng, _pointResoZ);


              //find smearing for x and y, so that hit is smeared along ladder plane
              smearedPos[0] = hitvec.x() + rPhiSmear * cos(ladder_incline);
              smearedPos[1] = hitvec.y() + rPhiSmear * sin(ladder_incline); 
              smearedPos[2] = hitvec.z() + zSmear;
              
              break;
              
            }
          ++tries;
        }
 
      if( accept_hit == false ) 
        {
          streamlog_out(DEBUG) << "hit could not be smeared within ladder after 100 tries: hit dropped"  << std::endl;
          continue; 
        } // 


      //store hit variables
      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;

      trkHit->setType(det_id_for_type+layerNumber );         // needed for FullLDCTracking et al.

      
      streamlog_out(DEBUG) <<"Position of hit after smearing = "<<smearedPos[0]<<" "<<smearedPos[1]<<" "<<smearedPos[2]
                           << " :" 
                           << "  u: " <<  u+rPhiSmear
                           << "  v: " <<  hitvec.z()+zSmear
                           << std::endl ;


      cellid_encoder[ ILDCellID0::subdet ] = _sub_det_id ;
      cellid_encoder[ ILDCellID0::side   ] = 0 ;
      cellid_encoder[ ILDCellID0::layer  ] = layerNumber ;
      cellid_encoder[ ILDCellID0::module ] = ladderNumber ;
      cellid_encoder[ ILDCellID0::sensor ] = 0 ;

      cellid_encoder.setCellID( trkHit ) ;

      trkHit->setPosition( smearedPos ) ;

      float u_direction[2] ;
      u_direction[0] = ladder_incline ;
      u_direction[1] = M_PI/2.0 ;

      float v_direction[2] ;
      v_direction[0] = 0.0 ;
      v_direction[1] = 0.0 ;

      trkHit->setU( u_direction ) ;
      trkHit->setV( v_direction ) ;
      
      trkHit->setdU( _pointResoRPhi ) ;
      trkHit->setdV( _pointResoZ ) ;

      trkHit->setEDep( edep ) ;
      
//      float covMat[TRKHITNCOVMATRIX]={0.,0.,_pointResoRPhi*_pointResoRPhi,0.,0.,_pointResoZ*_pointResoZ};
//      trkHit->setCovMatrix(covMat);      
        
      // 	  push back the SimTHit for this TrackerHit
      // fg: only if we have a sim hit with proper link to MC truth
    
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
    
    
    evt->addCollection( trkhitVec , _outColName ) ;
    
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


