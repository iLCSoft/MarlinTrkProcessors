/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SimplePlanarDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/CellIDEncoder.h>

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
  _description = "SimplePlanarDigiProcessor should create VTX and SIT TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "PointResolutionRPhi_Inner" ,
                              "R-Phi Resolution in the Inner Layers"  ,
                              _pointResoRPhi_Inner ,
                              float(0.0040)) ;
	
  registerProcessorParameter( "PointResolutionZ_Inner" , 
                              "Z Resolution in the Inner Layers" ,
                              _pointResoZ_Inner ,
                              float(0.0040));

  registerProcessorParameter( "PointResolutionRPhi_Outer" ,
                              "R-Phi Resolution in the Outer Layers"  ,
                              _pointResoRPhi_Outer ,
                              float(0.010)) ;
	
  registerProcessorParameter( "PointResolutionZ_Outer" , 
                              "Z Resolution in the Outer Layers" ,
                              _pointResoZ_Outer ,
                              float(0.010));


  registerProcessorParameter( "Last_Inner_Layer" , 
                              "Layer Number counting from 0 of the last Inner Layer" ,
                              _last_inner_layer ,
                              int(6));

  registerProcessorParameter( "Ladder_Number_encoded_in_cellID" , 
                              "Mokka has encoded the ladder number in the cellID" ,
                              _ladder_Number_encoded_in_cellID ,
                              bool(false));



  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VXDCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHIT,
                            "VTXHitCollection" , 
                            "Name of the vxd TrackerHit output collection"  ,
                            _outColNameVTX ,
                            std::string("VTXTrackerHits") ) ;
  
  
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
  _nRun++ ;
} 

void SimplePlanarDigiProcessor::processEvent( LCEvent * evt ) { 

  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;

  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _colNameVTX ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG) << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }
  
  if( STHcol != 0 ){    
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;


    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( ILDCellID0::encoder_string , trkhitVec ) ;

    
    int nSimHits = STHcol->getNumberOfElements()  ;
    

    //VXD smearing
    streamlog_out( DEBUG ) << " processing collection " << _colNameVTX 
                           << " with " <<  nSimHits  << " hits ... " << std::endl ;

    _pointResoRPhi = _pointResoRPhi_Inner;
    _pointResoZ = _pointResoZ_Inner;               
    

    //get VXD geometry info
    const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
    const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout() ; 
    
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
        layerNumber  = ( celId / 10000 ) - 1 ;
      }
      else{
        layerNumber = celId  - 1 ;
      }

      float edep = SimTHit->getEDep() ;
      
      MCParticle *mcp ;
      mcp = SimTHit->getMCParticle() ;
      
      //phi between each ladder
      double deltaPhi = ( 2 * M_PI ) / layerVXD.getNLadders(layerNumber) ;
      //      double sensitive_length  = layerVXD.getSensitiveLength(layerNumber);
      double sensitive_width  = layerVXD.getSensitiveWidth(layerNumber);
      double sensitive_offset = layerVXD.getSensitiveOffset(layerNumber);
      double ladder_r         = layerVXD.getSensitiveDistance(layerNumber);

      if(_ladder_Number_encoded_in_cellID) {
        ladderNumber = ( celId % ( 10000 * (layerNumber + 1) ) ) -1 ;
      }
      else{
        
        for (int ic=0; ic < layerVXD.getNLadders(layerNumber); ++ic) {
          
          double ladderPhi = correctPhiRange( layerVXD.getPhi0( layerNumber ) + ic*deltaPhi ) ;

          double PhiInLocal = hitvec.phi() - ladderPhi;
          double RXY = hitvec.rho();
          
          //          streamlog_out(DEBUG) << "ladderPhi = " << ladderPhi << " PhiInLocal = "<< PhiInLocal << " RXY = " << RXY << " (RXY*cos(PhiInLocal) = " << RXY*cos(PhiInLocal) << " layerVXD.getSensitiveDistance(layerNumber) = " << layerVXD.getSensitiveDistance(layerNumber) << endl;
      

          // check if point is in range of ladder
          if (RXY*cos(PhiInLocal) - layerVXD.getSensitiveDistance(layerNumber) > -layerVXD.getSensitiveThickness(layerNumber) && 
              RXY*cos(PhiInLocal) - layerVXD.getSensitiveDistance(layerNumber) <  layerVXD.getSensitiveThickness(layerNumber) )
            {
              ladderNumber = ic;
              break;
            }
        }
        
      }

      double ladderPhi = correctPhiRange( layerVXD.getPhi0( layerNumber ) + ladderNumber*deltaPhi ) ;
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
                           << "  layer phi0 " << layerVXD.getPhi0( layerNumber )
                           << "  phi: " <<  hitvec.phi()
                           << "  PhiInLocal: " << PhiInLocal
                           << "  ladderPhi: " << ladderPhi
                           << "  ladder_incline: " << ladder_incline
                           << "  nladders: " << layerVXD.getNLadders(layerNumber) 
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
          
          if( layerNumber <= _last_inner_layer ){
            _pointResoRPhi = _pointResoRPhi_Inner;
            _pointResoZ    = _pointResoZ_Inner;
          }
          else{
            _pointResoRPhi = _pointResoRPhi_Outer;
            _pointResoZ    = _pointResoZ_Outer;
          }

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
            
      trkHit->setType(100+layerNumber );         // needed for FullLDCTracking et al.

      
      streamlog_out(DEBUG) <<"Position of hit after smearing = "<<smearedPos[0]<<" "<<smearedPos[1]<<" "<<smearedPos[2]
                           << " :" 
                           << "  u: " <<  u+rPhiSmear
                           << "  v: " <<  hitvec.z()+zSmear
                           << std::endl ;

      // SJA:FIXME: here you can use the value 2 but not 3 which is odd as the width of the field is 1, only 0 and 1 should be allowed?
      int side = 1 ;

      if( smearedPos[2] < 0.0 ) side = 1 ;

      cellid_encoder[ ILDCellID0::subdet ] = ILDDetID::VXD ;
      cellid_encoder[ ILDCellID0::layer  ] = layerNumber ;
      cellid_encoder[ ILDCellID0::module ] = ladderNumber ;

      //SJA:FIXME: for now don't use side
      //  (*_cellid_encoder)[ ILDCellID0::Fields::side   ] = side ;
      cellid_encoder[ ILDCellID0::side   ] = 0 ;

      cellid_encoder.setCellID( trkHit ) ;

      trkHit->setPosition( smearedPos ) ;

      float u_direction[2] ;
      u_direction[0] = ladder_incline ;
      u_direction[1] = 0.0 ;

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
      if( mcp != 0 )  {
        trkHit->rawHits().push_back( SimTHit ) ;
      }
      else{
        streamlog_out( DEBUG0 ) << " ignore simhit pointer as MCParticle pointer is NULL ! " << std::endl ;
      }
      

      trkhitVec->addElement( trkHit ) ; 

      streamlog_out(DEBUG) << "-------------------------------------------------------" << std::endl;

    }      
    
    
    evt->addCollection( trkhitVec , _outColNameVTX ) ;
    
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


