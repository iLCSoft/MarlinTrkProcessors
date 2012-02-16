#include "TruthTracker.h"
#include <iostream>

#include <vector>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <IMPL/TrackImpl.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include <gear/BField.h>

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

using namespace lcio ;
using namespace marlin ;

using namespace MarlinTrk ;



TruthTracker aTruthTracker ;


TruthTracker::TruthTracker() : Processor("TruthTracker") {
  
  // modify processor description
  _description = "Creates Track Collection from MC Truth. Can handle composite spacepoints as long as they consist of two TrackerHits" ;
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  // register steering parameters: name, description, class-variable, default value
  
  
  std::vector< std::string > trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back( "FTDTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "SITTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "TPCTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "VXDTrackerHitRelations" );
//   trackerHitsRelInputColNamesDefault.push_back( "SETTrackerHitRelations" );
//   trackerHitsRelInputColNamesDefault.push_back( "ETDTrackerHitRelations" );
  
  registerProcessorParameter( "TrackerHitsRelInputCollections",
                              "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits. Have to be in same order as TrackerHitsInputCollections!!!",
                              _colNamesTrackerHitRelations,
                              trackerHitsRelInputColNamesDefault );
  
  
  std::vector< std::string > trackerHitsInputColNamesDefault;
  trackerHitsInputColNamesDefault.push_back( "FTDTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "SITTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "TPCTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "VXDTrackerHits" );

  
  registerProcessorParameter( "TrackerHitsInputCollections",
                              "Name of the tracker hit input collections",
                              _colNamesTrackerHits,
                              trackerHitsInputColNamesDefault );

  
  registerOutputCollection( LCIO::TRACK,
                           "OutputTrackCollectionName" , 
                           "Name of the output track collection"  ,
                           _output_track_col_name ,
                           std::string("TruthTracks") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                           "OutputTrackRelCollection" , 
                           "Name of the MCParticle-Track Relations collection for output tracks"  ,
                           _output_track_rel_name ,
                           std::string("TruthTracksMCP") ) ;
  
//   registerProcessorParameter( "nEventPrintout",
//                              "Print out progress every N events ",
//                              _nEventPrintout,
//                              int(1000));
  
  registerProcessorParameter( "MCpThreshold",
                             "Transverse Momentum Threshold MC particles which will produce tracks GeV",
                             _MCpThreshold,
                             float(0.1));
  
  registerProcessorParameter( "FitTracksWithMarlinTrk",
                             "Fit the Tracks with MarlinTrk, otherwise take track parameters from MCParticle",
                             _FitTracksWithMarlinTrk,
                             bool(true));
  
  registerProcessorParameter("MultipleScatteringOn",
                             "Use MultipleScattering in Fit",
                             _MSOn,
                             bool(true));
  
  registerProcessorParameter("EnergyLossOn",
                             "Use Energy Loss in Fit",
                             _ElossOn,
                             bool(true));
  
  registerProcessorParameter("SmoothOn",
                             "Smooth All Mesurement Sites in Fit",
                             _SmoothOn,
                             bool(false));
  
  
  
  _n_run = 0 ;
  _n_evt = 0 ;
  

  

  
}


void TruthTracker::init() { 
  
  streamlog_out(DEBUG) << "   init called  " 
  << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  
  if( _trksystem == 0 ) {
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
  
}

void TruthTracker::processRunHeader( LCRunHeader* run) { 
  
  ++_n_run ;
} 

void TruthTracker::processEvent( LCEvent * evt ) { 
  
  

  streamlog_out(DEBUG3) << "   processing event: " << _n_evt << std::endl ;
  
  _nMCP = 0 ;
  
  _colTrackerHits.clear();
  _navTrackerHitRel.clear();
  _nCreatedTracks = 0;
  
  /**********************************************************************************************/
  /*                Prepare the collections                                                     */
  /**********************************************************************************************/
  
  // get the input collections and fill the vectors
  this->SetupInputCollections(evt) ;
  
  // establish the track collection that will be created 
  _trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
  
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  _trackVec->setFlag( trkFlag.getFlag()  ) ;
  
  // establish the track relations collection that will be created 
  _trackRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;
  
  // create the encoder to decode cellID0
  UTIL::BitField64 cellID_encoder( ILDCellID0::encoder_string ) ;
  
  
  
  /**********************************************************************************************/
  /*                Store the TrackerHits and SimTrackerHits as pairs                           */
  /**********************************************************************************************/
  
  std::vector< std::pair<SimTrackerHit*, TrackerHit* > > simHitTrkHit;
  
  
  for( unsigned iCol=0; iCol<_colTrackerHits.size(); iCol++){
    
    LCCollection* trackerHitCol = _colTrackerHits[iCol];
    int nHits = trackerHitCol->getNumberOfElements();
    
    LCRelationNavigator* nav = _navTrackerHitRel[iCol];
    
    for( int j=0; j<nHits; j++ ){
      
      TrackerHit * trkhit = dynamic_cast<TrackerHit*>( trackerHitCol->getElementAt( j ));      
      
     
      const LCObjectVec& to = nav->getRelatedToObjects( trkhit );
      
      if( trkhit->getType() == UTIL::ILDTrkHitType::COMPOSITE_SPACEPOINT ){ //it is a composite spacepoint
        
        if( to.size() == 2 ){
          
          SimTrackerHit* simhitA = dynamic_cast<SimTrackerHit*>(to.at(0));
          SimTrackerHit* simhitB = dynamic_cast<SimTrackerHit*>(to.at(1));
          
          // Check if the simHits are from the same particel in order to avoid ghost hits
          if( simhitA->getMCParticle() == simhitB->getMCParticle() ) simHitTrkHit.push_back(std::make_pair(simhitA, trkhit));
          else streamlog_out( DEBUG0 ) << "spacepoint discarded, because simHits are not equal " << simhitA->getMCParticle() << " != " 
            << simhitB->getMCParticle() << "\n";
          
        }        
        else{ streamlog_out( DEBUG0 ) << "spacepoint discarded, because it is related to " << to.size() << "SimTrackerHits. It should be 2!\n"; } 
        
      }
      else{  // no composite spacepoint
        
        if( to.size() == 1){ // only take trackerHits, that have only one related SimHit
          
          SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(to.at(0));
          simHitTrkHit.push_back(std::make_pair(simhit, trkhit));
          
        }
        else{ streamlog_out( DEBUG0 ) << "TrackerHit discarded, because it is related to " << to.size() << "SimTrackerHits. It should be 1!\n"; }
      
      }
      
      
    }
    
  }
  streamlog_out( DEBUG4 ) << "Number of Tracker hits = " << simHitTrkHit.size() << std::endl;     
  
  // now order the hits by MCParticle and then in time  
  std::sort( simHitTrkHit.begin(),simHitTrkHit.end(), SimTrackerHitSortPredicate() );
  

  
 
  // some information output
  for(unsigned int i=0; i< simHitTrkHit.size() ; ++i ) {
    
    SimTrackerHit * simHit = simHitTrkHit[i].first;
    TrackerHit* trackerHit = simHitTrkHit[i].second;
    streamlog_out( DEBUG2 ) << "Tracker hit: [" << i << "] = " << trackerHit << "  mcp = " << simHit->getMCParticle() << " time = " << simHit->getTime() << " cellid (trkHit) = " << trackerHit->getCellID0() 
    << " (de" << getDetectorID( trackerHit ) 
    << ",si" << getSideID( trackerHit ) 
    << ",la" <<  getLayerID( trackerHit ) 
    << ",mo"<< getModuleID( trackerHit ) 
    << ",se"<< getSensorID( trackerHit ) 
    << ")" << std::endl;  
    
    
  }

  /**********************************************************************************************/
  /*                Take hits with same MCP and create tracks from them                         */
  /**********************************************************************************************/
  
  streamlog_out( DEBUG4 ) << "Add Tracker hits to Tracks" << std::endl;
  
  std::vector<TrackerHit*> hit_list;
    
  if( simHitTrkHit.size() > 0) {
    
    MCParticle* mcplast = NULL;
      
      
    for(unsigned int i=0; i< simHitTrkHit.size() ; ++i) {
      
      
      SimTrackerHit* simhit = simHitTrkHit[i].first;
      
      MCParticle* mcp = simhit->getMCParticle();
      double const* p    = mcp->getMomentum() ;
      float  const pt2   = p[0]*p[0] + p[1]*p[1] ;
//         float  const pmag2 =  p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ; 
          
      if ( i == 0 ) {
        mcplast = mcp ;
      }
    
      if( mcp != mcplast ) { 
      // new MCParticle
        
        if ( hit_list.size() >= 3) {
          // create track from vector of hits                           
          streamlog_out( DEBUG2 ) << "Create New Track for MCParticle " << mcplast << std::endl;
          this->createTrack(mcplast, cellID_encoder, hit_list );
          
        }
        
        hit_list.clear();     // clear the list for the new mcp
        
      }
      
      if( pt2  > (_MCpThreshold*_MCpThreshold) ) { // if momentum is greater than cut add hit to list
        
        streamlog_out( DEBUG2 ) << "Add hit from det " <<  simhit->getCellID0()  << " to track from MCParticle " << mcp << " : current number of hits = " << hit_list.size() << std::endl;
        hit_list.push_back( simHitTrkHit[i].second ) ;
      }
      
      // set last mcparticle 
      mcplast = mcp;
    
    } // end of loop over hits
    
    // check if there is still a track to be created 
    if( hit_list.size() >= 3 ) { 
      // then create a new track
      streamlog_out( DEBUG3 ) << "Create New Track for Last MCParticle " << mcplast << std::endl;
      this->createTrack(mcplast, cellID_encoder, hit_list );
      
      
      hit_list.clear();  
      
    }
    
  }    
  
  evt->addCollection( _trackVec , _output_track_col_name) ;
  evt->addCollection( _trackRelVec , _output_track_rel_name) ;
  
  
  streamlog_out( DEBUG4 ) << "Created " << _nCreatedTracks << " truth tracks\n";
  
  ++_n_evt ;
  
}



void TruthTracker::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TruthTracker::end() { 
  
  streamlog_out(DEBUG4) << "TruthTracker::end()  " << name() 
  << " processed " << _n_evt << " events in " << _n_run << " runs "
  << std::endl ;

  delete _encoder ;
  
}





void TruthTracker::SetupInputCollections( LCEvent * evt ) {
  
  

  // Check if there are as many tracker hit input collections as relation collections
  if(  _colNamesTrackerHits.size() !=  _colNamesTrackerHitRelations.size() ){
    
    streamlog_out( ERROR ) << "There must be as many input collections of tracker Hits as of relations. At the moment, there are "
    << _colNamesTrackerHits.size() << " tracker hit collections and " << _colNamesTrackerHitRelations.size() << " relation collections passed as steering paremeters!\n";
    
    exit(1);
  }
  
  for( unsigned i=0; i< _colNamesTrackerHits.size(); i++ ){
    
    
    // the tracker hits
    LCCollection* colTrkHits = GetCollection( evt, _colNamesTrackerHits[i] );
    if( colTrkHits == NULL ) continue;
    
    // the relations of them
    LCRelationNavigator* nav = GetRelations( evt, _colNamesTrackerHitRelations[i] );
    if( nav == NULL ) continue;
    
    
    _colTrackerHits.push_back( colTrkHits );
    _navTrackerHitRel.push_back( nav );
    
    
  }
 
  
}


void TruthTracker::createTrack( MCParticle* mcp, UTIL::BitField64& cellID_encoder, std::vector<TrackerHit*>& hit_list ) {
  
  
  TrackImpl* Track = new TrackImpl ; 
  
  streamlog_out( DEBUG3 ) << "Create track with " << hit_list.size() << " hits" << std::endl;  
  
  double d0;
  double phi0;
  double omega;
  double z0;
  double tanL;
  
  double chi2 = 0 ;
  int ndf = 0 ;
  
  float ref[3];

  TrackerHitVec added_hits;
  
  if(_FitTracksWithMarlinTrk) {
    
    MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
    
    
    
    // sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to in. 
    //sort(_hit_list.begin(), _hit_list.end(), TruthTracker::compare_time() );
    sort( hit_list.begin(), hit_list.end(), TruthTracker::compare_r() );

    unsigned counter=0; // degrees of freedom for the fit (this is actually not the number of degrees of freedom. This number minus 5 is it. But how to name this?)
    
    for(unsigned int j=0; j< hit_list.size(); ++j) {
      
      TrackerHit* trkHit = hit_list[j];
      
      // Check for spacepoints
      if( trkHit->getType() == UTIL::ILDTrkHitType::COMPOSITE_SPACEPOINT ){
        
        //Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();
        
        bool isSuccessful = false; 
        
        for( unsigned k=0; k< rawObjects.size(); k++ ){
          
          TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
          if( marlin_trk->addHit( rawHit ) == IMarlinTrack::success ){
            
            isSuccessful = true; //if at least one hit from the spacepoint gets added
            counter++; // one strip hit gives 1 degree of freedom
            
          }
        }
        
        if ( isSuccessful ) added_hits.push_back( trkHit );
        
      }
      else{ // normal non composite tracker hit
        
        if( marlin_trk->addHit( trkHit ) == IMarlinTrack::success ){
          
          added_hits.push_back( trkHit );
          counter += 2; // a spacepoint gives 2 degrees of freedom
        }
      }
    }  
    
    // only fit tracks, that have enough degrees of freedom!
    if ( counter < 8) {
      streamlog_out( DEBUG3 ) << "Less than 8 degrees of freedom: only " << counter << "\n";  
      delete Track;
      delete marlin_trk;
      return;
    }
    
    
    marlin_trk->initialise( IMarlinTrack::backward ) ;
    int fit_status = marlin_trk->fit() ; 
        
    if( fit_status == 0 ) { 
      
      const gear::Vector3D point(0.,0.,0.); // nominal IP
      
      TrackStateImpl* trkState = new TrackStateImpl() ;
      int return_code = marlin_trk->propagate(point, *trkState, chi2, ndf ) ;
      if ( return_code == 0 ) {
        
        Track->addTrackState(trkState);
        Track->setChi2(chi2) ;
        Track->setNdf(ndf) ;
        
      }
    }
    
    delete marlin_trk;
    
  }
  else {
    
    float bZ = float(marlin::Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z()) ;
    // use mcp pos and mom to set the track parameters
    HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), bZ );
    d0    = hel.getD0();
    phi0  = hel.getPhi0();
    omega = hel.getOmega();
    z0    = hel.getZ0();
    tanL  = hel.getTanLambda();
    
    ref[0] = hel.getRefPointX() ;
    ref[1] = hel.getRefPointY() ;
    ref[2] = hel.getRefPointZ() ;
    
    Track->setD0(d0);
    Track->setPhi(phi0);
    Track->setOmega(omega);
    Track->setZ0(z0);
    Track->setTanLambda(tanL);
    
    Track->setReferencePoint(ref) ;
    
    
    added_hits = hit_list;
    
  }
  
  
  
  std::map<int, int> hitNumbers; 
  
  hitNumbers[lcio::ILDDetID::VXD] = 0;
  hitNumbers[lcio::ILDDetID::SIT] = 0;
  hitNumbers[lcio::ILDDetID::FTD] = 0;
  hitNumbers[lcio::ILDDetID::TPC] = 0;
  hitNumbers[lcio::ILDDetID::SET] = 0;
  hitNumbers[lcio::ILDDetID::ETD] = 0;
  
  for(unsigned int j=0; j<added_hits.size(); ++j) {
    
    Track->addHit(added_hits.at(j)) ;
    
    cellID_encoder.setValue(added_hits.at(j)->getCellID0()) ;
    int detID = cellID_encoder[ILDCellID0::subdet];
    ++hitNumbers[detID];
    //    streamlog_out( DEBUG1 ) << "Hit from Detector " << detID << std::endl;     
  }
  
  //SJA:FIXME no distiction made for hits in fit or not
  Track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ] = hitNumbers[lcio::ILDDetID::VXD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ] = hitNumbers[lcio::ILDDetID::FTD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 2 ] = hitNumbers[lcio::ILDDetID::SIT];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ] = hitNumbers[lcio::ILDDetID::TPC];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 2 ] = hitNumbers[lcio::ILDDetID::SET];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 2 ] = hitNumbers[lcio::ILDDetID::ETD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 1 ] = hitNumbers[lcio::ILDDetID::VXD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 1 ] = hitNumbers[lcio::ILDDetID::FTD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 1 ] = hitNumbers[lcio::ILDDetID::SIT];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ] = hitNumbers[lcio::ILDDetID::TPC];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 1 ] = hitNumbers[lcio::ILDDetID::SET];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 1 ] = hitNumbers[lcio::ILDDetID::ETD];

  _trackVec->addElement(Track);

  LCRelationImpl* rel = new LCRelationImpl;
  rel->setFrom (Track);
  rel->setTo (mcp);
  rel->setWeight(1.0);
  _trackRelVec->addElement(rel);

  _nCreatedTracks++;
  
}

LCCollection* TruthTracker::GetCollection(  LCEvent * evt, std::string colName ){
  
  LCCollection* col = NULL;
  
  int nElements = 0;
  
  try {
    col = evt->getCollection( colName.c_str() ) ;
    nElements = col->getNumberOfElements()  ;
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() << " collection found, number of elements = " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e) {
    streamlog_out( ERROR ) << " --> " << colName.c_str() <<  " collection absent" << std::endl;     
  }
  
  return col; 
  
}

LCRelationNavigator* TruthTracker::GetRelations(LCEvent * evt , std::string RelName ) {
  
  LCRelationNavigator* nav = NULL ;
  LCCollection* col = NULL;
  
  try{
    
    col = evt->getCollection( RelName.c_str() );
    nav = new LCRelationNavigator( col );
    streamlog_out( DEBUG4 ) << " --> " << RelName << " track relation collection found, number of elements = " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( ERROR ) << " --> " << RelName.c_str() << " track relation collection absent" << std::endl;     
  }
  
  return nav;
  
}
