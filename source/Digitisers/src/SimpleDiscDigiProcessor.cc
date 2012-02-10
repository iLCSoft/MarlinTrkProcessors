/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SimpleDiscDigiProcessor.h"
#include <iostream>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/ILDConf.h>



#include <marlin/Global.h>
#include "marlin/ProcessorEventSeeder.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"
#include <cmath>

#include <gsl/gsl_randist.h>

//FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
#include "MarlinTrk/Factory.h"

#include "gear/gearsurf/MeasurementSurfaceStore.h"

#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"


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
                             int(lcio::ILDDetID::FTD));
                           
  
  registerProcessorParameter( "SimHits_encoded_with_cellID" , 
                             "Mokka has encoded the hits with the cellID as specified in UTIL/ILDConf.h from lcio" ,
                             _SimHits_encoded_with_cellID ,
                             bool(false));

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
  
  registerOutputCollection( LCIO::TRACKERHITPLANE,
                           "OutputCollectionName" , 
                           "Name of the TrackerHit output collection"  ,
                           _outColName ,
                           std::string("FTDTrackerHits") ) ;
  
 
  registerOutputCollection(LCIO::LCRELATION,
                           "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection",
                           _outRelColName,
                           std::string("FTDTrackerHitRelations"));
  
  
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
  
  _use_FTDLayerLayout_from_GEAR = false ;
  
  try {
    
    const gear::FTDParameters& ftdParams = Global::GEAR->getFTDParameters() ;
    const gear::FTDLayerLayout& ftdlayers = ftdParams.getFTDLayerLayout() ;
    streamlog_out( MESSAGE ) << "  SimpleDiscDigiProcessor - Use FTDLayerLayout with " << ftdlayers.getNLayers() << " layers" << std::endl ;
    _use_FTDLayerLayout_from_GEAR = true ;
    
  } catch (gear::UnknownParameterException& e) {
    
    _use_FTDLayerLayout_from_GEAR = false ;
  }
  
  if ( _use_FTDLayerLayout_from_GEAR == false ) {

    streamlog_out( MESSAGE ) << "  SimpleDiscDigiProcessor - Use Loi style FTDParameters" << std::endl ;
    
    const gear::GearParameters& pFTD = Global::GEAR->getGearParameters("FTD");
    
    _FTDZCoordinate = pFTD.getDoubleVals( "FTDZCoordinate" ) ;
    _diskInnerRadius = pFTD.getDoubleVals( "FTDInnerRadius" ) ;
    _diskOuterRadius = pFTD.getDoubleVals( "FTDOuterRadius" ) ;
    
  }

  //FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
  MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  
  if( trksystem == 0 ) {
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  trksystem->init() ;  
  
  //FIXME:SJA gear surface store has now been filled so we can dispose of the MarlinTrkSystem
  delete trksystem;

  
}

void SimpleDiscDigiProcessor::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void SimpleDiscDigiProcessor::process_hits_loi( LCEvent * evt, LCCollection* STHcol ) { 

  if( STHcol != 0 ){    
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;

    LCCollectionVec * relCol = new LCCollectionVec(LCIO::LCRELATION);
    
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;
    
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::ILDCellID0::encoder_string , trkhitVec ) ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << " process hits alla loi " << _inColName 
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
      
      if ( ( _keepHitsFromDeltas == true ) || ( hasCorrectZPos (SimTHit) == true ) ){
        
        streamlog_out(DEBUG) << "Hit = "<< i << " has celId " << celId << " layer number = " << layerNumber  << endl;
        
        streamlog_out(DEBUG) <<"Position of hit before smearing = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " r = " << hitvec.rho() << endl;
        
        double xSmear = gsl_ran_gaussian(r,_pointReso);
        double ySmear = gsl_ran_gaussian(r,_pointReso);
        
        double smearedPos[3] ;
        smearedPos[0] = pos[0] + xSmear;
        smearedPos[1] = pos[1] + ySmear;

        // No semaring of Z coordinate
        //smearedPos[2] = pos[2] ;

        smearedPos[2] = side*_FTDZCoordinate[layerNumber] ;
        
        streamlog_out(DEBUG) <<"Position of hit after smearing with z set to disk z-coord. = "<<smearedPos[0]<<" "<<smearedPos[1]<<" "<<smearedPos[2] << std::endl ;
        
        //store hit variables
        TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;        
        
        //        trkHit->setType( 200+layerNumber );  // needed for FullLDCTracking et al.
        
        int petalNumber = getPetalNumber( layerNumber , smearedPos[0] , smearedPos[1] );
        int sensorNumber = getSensorNumber( layerNumber , smearedPos[0] , smearedPos[1] );
        
        
        cellid_encoder[ lcio::ILDCellID0::subdet ] = _sub_det_id ;
        cellid_encoder[ lcio::ILDCellID0::side   ] = side ;
        cellid_encoder[ lcio::ILDCellID0::layer  ] = layerNumber ;
        cellid_encoder[ lcio::ILDCellID0::module ] = petalNumber ;
        cellid_encoder[ lcio::ILDCellID0::sensor ] = sensorNumber ;
        
        cellid_encoder.setValue( lcio::long64(cellid_encoder.lowWord() ) << 32 );
        
        cellid_encoder[ lcio::ILDCellID0::subdet ] = _sub_det_id ;
        cellid_encoder[ lcio::ILDCellID0::side   ] = side ;
        cellid_encoder[ lcio::ILDCellID0::layer  ] = layerNumber ;
        cellid_encoder[ lcio::ILDCellID0::module ] = 0 ;
        cellid_encoder[ lcio::ILDCellID0::sensor ] = 0 ;
        
        cellid_encoder.setCellID( trkHit ) ;
        
        streamlog_out(DEBUG2) <<"side = "<< side << std::endl ;
        streamlog_out(DEBUG2) <<"layerNumber = "<< layerNumber << std::endl ;
        streamlog_out(DEBUG2) <<"moduleNumber = "<< petalNumber << std::endl ;
        streamlog_out(DEBUG2) <<"sensorNumber = "<< sensorNumber << std::endl ;
        
        trkHit->setPosition(  smearedPos  ) ;
        
        float u_direction[2] ; // x
        u_direction[0] = M_PI/2.0 ;
        u_direction[1] = 0.0 ; 
        
        float v_direction[2] ; // y
        v_direction[0] = M_PI/2.0 ;
        v_direction[1] = M_PI/2.0 ;
        
        trkHit->setU( u_direction ) ;
        trkHit->setV( v_direction ) ;
        
        trkHit->setdU( _pointReso ) ;
        trkHit->setdV( _pointReso ) ;
        
        trkHit->setEDep( SimTHit->getEDep() ) ;
        trkHit->setTime( SimTHit->getTime() ) ;
        
        LCRelationImpl* rel = new LCRelationImpl;
        
        rel->setFrom (trkHit);
        rel->setTo (SimTHit);
        rel->setWeight( 1.0 );
        relCol->addElement(rel);
        
//        MCParticle *mcp ;
//        mcp = SimTHit->getMCParticle() ;
//        if( mcp != 0 )  {
//          trkHit->rawHits().push_back( SimTHit ) ;
//        }
//        else{
//          streamlog_out( DEBUG0 ) << " ignore simhit pointer as MCParticle pointer is NULL ! " << std::endl ;
//        }
        
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
    evt->addCollection( relCol , _outRelColName ) ;
    
  }
  
}

void SimpleDiscDigiProcessor::process_hits_new( LCEvent * evt, LCCollection* STHcol ) { 

  if( STHcol != 0 ){    
    
    unsigned nCreatedHits=0;
    unsigned nDismissedHits=0;
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
    LCCollectionVec * relCol = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;
    
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::ILDCellID0::encoder_string , trkhitVec ) ;
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << " processing hit using new FTDLayerLayout " << _inColName 
    << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    for(int i=0; i< nSimHits; i++){
      
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      
      //*********************************************
      // the first issue is to deal with the cellID 
      //*********************************************
      
      const int celId = SimTHit->getCellID0() ;
      
      int layerNumber(0);

      UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;       
      
      encoder.setValue(celId) ;  
      
      cellid_encoder.setValue(encoder.getValue());
      
      layerNumber = encoder[lcio::ILDCellID0::layer];
      
      streamlog_out( DEBUG2 ) << "CelId : " << celId <<
      " subdet = " << encoder[lcio::ILDCellID0::subdet] <<
      " side = " << encoder[lcio::ILDCellID0::side] <<
      " layer = " << encoder[lcio::ILDCellID0::layer] <<
      " module = " << encoder[lcio::ILDCellID0::module] <<
      " sensor = " << encoder[lcio::ILDCellID0::sensor] <<
      std::endl ;
        
      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************
      
      gear::MeasurementSurface const* ms = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( cellid_encoder.lowWord() );

      
      //************************************************************
      // Now get the energy and position in global and local coord.
      //************************************************************
      
      const double *pos ;
      pos =  SimTHit->getPosition() ;
      CLHEP::Hep3Vector globalPoint(pos[0],pos[1],pos[2]);
      CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
      
      double u = localPoint[0];
      double v = localPoint[1];
      
      streamlog_out(DEBUG3) << ":" 
      << "  layer: " << layerNumber 
      << "  u: " <<  u
      << "  v: " <<  v
      << "  phi: " <<  globalPoint.phi()
      << std::endl ;

      
      if ( ( _keepHitsFromDeltas == true ) || ( hasCorrectZPos (SimTHit) == true ) ){
        
        streamlog_out(DEBUG) << "Hit = "<< i << " has celId " << celId << " layer number = " << layerNumber  << endl;
        
        streamlog_out(DEBUG) <<"Position of hit before smearing = "<<globalPoint[0]<<" "<<globalPoint[1]<<" "<<pos[2]<< " r = " << globalPoint.rho() << "  phi: " <<  globalPoint.phi()<< endl;

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
          
          if(tries > 0) streamlog_out(DEBUG0) << "retry smearing for " << encoder[lcio::ILDCellID0::layer] << " " << encoder[lcio::ILDCellID0::module] << " : retries " << tries << std::endl;
          
          uSmear  = gsl_ran_gaussian(r, _pointReso);
          vSmear  = gsl_ran_gaussian(r, _pointReso);
          
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
        
        double smearedPos[3] ;
        
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
              
        cellid_encoder.setValue(encoder.getValue());
        
        cellid_encoder.setCellID( trkHit ) ;
        
        trkHit->setPosition(  smearedPos  ) ;
        
        trkHit->setEDep( SimTHit->getEDep() ) ;
        
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
        
        trkHit->setdU( _pointReso ) ;
        trkHit->setdV( _pointReso ) ;
        
        
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
        
        streamlog_out(DEBUG) << "-------------------------------------------------------" << std::endl;
        
      }
      else{
        
        ++nDismissedHits;
        streamlog_out(DEBUG) << "Hit "<< i << " is NOT KEPT! The z value is does not exactly correspond to a disk"  << endl;
        
        streamlog_out(DEBUG) <<"Position of hit = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<< " r = " << globalPoint.rho() << std::endl << endl;
        
        streamlog_out(DEBUG) << "-------------------------------------------------------" << std::endl;
        
      }
      
    }
    
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
    
    evt->addCollection( trkhitVec ,  _outColName ) ;
    evt->addCollection( relCol , _outRelColName ) ;
    
  
    streamlog_out(DEBUG4) << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed as not on sensitive element\n";

    
  }

  
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
  
  
  if (_use_FTDLayerLayout_from_GEAR) {
    this->process_hits_new(evt,STHcol);
  }
  else{
    this->process_hits_loi(evt, STHcol);
  }
  
  _nEvt ++ ;
}



void SimpleDiscDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SimpleDiscDigiProcessor::end(){ 
  
  gsl_rng_free(r);  
  //   std::cout << "SimpleDiscDigiProcessor::end()  " << name() 
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;
  
}

bool SimpleDiscDigiProcessor::hasCorrectZPos ( SimTrackerHit* hit ){
  
  double zPos = fabs ( hit->getPosition()[2] );
  
  if(_use_FTDLayerLayout_from_GEAR){
    
    streamlog_out(DEBUG) << "SimpleDiscDigiProcessor::hasCorrectZPos using FTDLayerLayout_from_GEAR" << std::endl;
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.setValue(hit->getCellID0()) ;  
   
    if( encoder[lcio::ILDCellID0::sensor] < 1 || encoder[lcio::ILDCellID0::sensor] > 4 ) {
      
      streamlog_out(ERROR) << "FTD sensor value is not in the range 1-4 : value = " << encoder[lcio::ILDCellID0::sensor] << std::endl;
      return false;
      
    }
    
    if( encoder[lcio::ILDCellID0::sensor] == 3 || encoder[lcio::ILDCellID0::sensor] == 4 ) {
      
      streamlog_out(DEBUG4) << "FTD sensor value is not in the range 1-2 : Hit Skipped: value = " << encoder[lcio::ILDCellID0::sensor] << " Currently only the sensitive faces facing the IP are taken" << std::endl;
      return false;
      
    }

    
    int layer = encoder[lcio::ILDCellID0::layer];
    int petal = encoder[lcio::ILDCellID0::module];
    int sensor = encoder[lcio::ILDCellID0::sensor];
    
    double zSensitive = Global::GEAR->getFTDParameters().getFTDLayerLayout().getSensitiveZposition(layer, petal, sensor);
    
    streamlog_out(DEBUG) << " zSensitive = " << zSensitive;
    streamlog_out(DEBUG) << " zPos = " << zPos << std::endl;
    
    if ( fabs ( zSensitive - zPos ) < 0.0001 ) return true;
    
  }
  else{

    streamlog_out(DEBUG) << "SimpleDiscDigiProcessor::hasCorrectZPos using FTD ZCoordinate only " << std::endl;
    
    for (unsigned int i=0; i < _FTDZCoordinate.size(); i++){
      
      streamlog_out(DEBUG) << " zSensitive = " << _FTDZCoordinate[i];
      streamlog_out(DEBUG) << " zPos = " << zPos << " dz = " << fabs ( _FTDZCoordinate[i] - zPos ) << std::endl;
      
      //      if ( (fabs ( _FTDZCoordinate[i] - zPos ) < 0.0001) || (fabs ( _FTDZCoordinate[i] - 0.275 - zPos ) < 0.0001) ) return true;
      //if ((fabs ( _FTDZCoordinate[i] - 0.275 - zPos ) < 0.0001) ) return true;
      if ((fabs ( _FTDZCoordinate[i] - zPos ) < 0.01) ) return true;
      
    }
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




