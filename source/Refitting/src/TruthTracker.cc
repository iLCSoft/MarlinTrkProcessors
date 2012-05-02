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


#include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif


#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/LCIOTrackPropagators.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

using namespace lcio ;
using namespace marlin ;

using namespace MarlinTrk ;



TruthTracker aTruthTracker ;

TruthTracker::TruthTracker() : Processor("TruthTracker") {
    
  // modify processor description
  _description = "Creates Track Collection from MC Truth. Can handle composite spacepoints as long as they consist of two TrackerHits" ;
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  // register steering parameters: name, description, class-variable, default value
  
  
  StringVec trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back( "VXDTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "SITTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "FTDPixelTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "FTDSpacePointRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "TPCTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "SETTrackerHitRelations" );
  
  
  registerInputCollections("LCRelation",
                           "TrackerHitsRelInputCollections",
                           "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits. Have to be in same order as TrackerHitsInputCollections!!!",
                           _colNamesTrackerHitRelations,
                           trackerHitsRelInputColNamesDefault );
  
  
  StringVec trackerHitsInputColNamesDefault;
  
  trackerHitsInputColNamesDefault.push_back( "VXDTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "SITTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "FTDPixelTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "FTDSpacePointRelations" );
  trackerHitsInputColNamesDefault.push_back( "TPCTrackerHits" );
  trackerHitsInputColNamesDefault.push_back( "SETTrackerHits" );
  
  registerInputCollections("TrackerHit",
                           "TrackerHitsInputCollections", 
                           "Name of the tracker hit input collections",
                           _colNamesTrackerHits,
                           trackerHitsInputColNamesDefault);
  
  
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
  
  registerProcessorParameter( "UseMCParticleParametersFotInitOfFit",
                             "When fitting take track parameters from MCParticle for the initialisation of the Track Fit",
                             _useMCParticleParametersFotInitOfFit,
                             bool(false));
  
  
  registerProcessorParameter( "CreatePrefitUsingMarlinTrk",
                             "When fitting take track parameters a full pre-fit for the initialisation of the Track Fit",
                             _create_prefit_using_MarlinTrk,
                             bool(false));
  
  
  registerProcessorParameter( "InitialTrackErrorD0",
                             "Value used for the initial d0 variance of the trackfit",
                             _initialTrackError_d0,
                             float(1.e6));
  
  registerProcessorParameter( "InitialTrackErrorPhi0",
                             "Value used for the initial phi0 variance of the trackfit",
                             _initialTrackError_phi0,
                             float(1.e2));
  
  registerProcessorParameter( "InitialTrackErrorOmega",
                             "Value used for the initial omega variance of the trackfit",
                             _initialTrackError_omega,
                             float(1.e-4));
  
  registerProcessorParameter( "InitialTrackErrorZ0",
                             "Value used for the initial z0 variance of the trackfit",
                             _initialTrackError_z0,
                             float(1.e6));
  
  registerProcessorParameter( "InitialTrackErrorTanL",
                             "Value used for the initial tanL variance of the trackfit",
                             _initialTrackError_tanL,
                             float(1.e2));
  
#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  registerOptionalParameter("RunMarlinTrkDiagnostics", "Run MarlinTrk Diagnostics. MarlinTrk must be compiled with MARLINTRK_DIAGNOSTICS_ON defined", _runMarlinTrkDiagnostics, bool(true));
  
  registerOptionalParameter("DiagnosticsName", "Name of the root file and root tree if running Diagnostics", _MarlinTrkDiagnosticsName, std::string("TruthTrackerDiagnostics"));    
  
#endif
  
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
  
  
#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  void * dcv = _trksystem->getDiagnositicsPointer();
  DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
  dc->init(_MarlinTrkDiagnosticsName,_MarlinTrkDiagnosticsName, _runMarlinTrkDiagnostics);
  
#endif
  
  _Bz = Global::GEAR->getBField().at( gear::Vector3D(0., 0., 0.) ).z();    //The B field in z direction
  
}

void TruthTracker::processRunHeader( LCRunHeader* run) { 
  
  ++_n_run ;
} 

void TruthTracker::processEvent( LCEvent * evt ) { 
  
  
  
  streamlog_out(DEBUG3) << "   processing event: " << _n_evt << std::endl ;
  
  _current_evt_number = evt->getEventNumber();
  
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
      
      if ( trackerHitCol->getElementAt( j ) == 0 ) {
        streamlog_out( DEBUG0 ) << "Pointer to TrackerHit" << j << " is NULL " << std::endl; 
      }
      
      TrackerHit * trkhit = dynamic_cast<TrackerHit*>( trackerHitCol->getElementAt( j ));      
      
      if ( trkhit == 0 ) {
        
        std::stringstream errorMsg;                
        errorMsg << "dynamic_cast to TrackerHit for hit " << j << " failed. Pointer = " <<  trackerHitCol->getElementAt( j ) << std::endl; 
        
        throw lcio::Exception(errorMsg.str());
        
      }
      
      const LCObjectVec& to = nav->getRelatedToObjects( trkhit );
      
      if( BitSet32( trkhit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
        
        if( to.size() == 2 ){
          
          SimTrackerHit* simhitA = dynamic_cast<SimTrackerHit*>(to.at(0));
          SimTrackerHit* simhitB = dynamic_cast<SimTrackerHit*>(to.at(1));
          
          // Check if the simHits are from the same particel in order to avoid ghost hits
          if( simhitA->getMCParticle() == simhitB->getMCParticle() ) simHitTrkHit.push_back(std::make_pair(simhitA, trkhit));
          else streamlog_out( DEBUG0 ) << "spacepoint discarded, because simHits are not equal " << simhitA->getMCParticle() << " != " 
            << simhitB->getMCParticle() << "\n";
          
          
#ifdef MARLINTRK_DIAGNOSTICS_ON
          
          // set the pointer to the simhit via lcio extention MCTruth4HitExt
          
          //Split it up and add both hits to the MarlinTrk
          const LCObjectVec rawObjects = trkhit->getRawHits();
          
          for( unsigned k=0; k< rawObjects.size(); k++ ){
            
            TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
            if( rawHit ){
              
              if( rawHit->getCellID0() == simhitA->getCellID0() ) {
                streamlog_out( DEBUG4 ) << "link simhit = " << simhitA << " Cell ID = " << simhitA->getCellID0() << " with trkhit = " << rawHit << " Cell ID = " <<  rawHit->getCellID0() << std::endl;     
                rawHit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;    
                rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhitA;                 
              }
              if( rawHit->getCellID0() == simhitB->getCellID0() ) {
                streamlog_out( DEBUG4 ) << "link simhit = " << simhitB << " Cell ID = " << simhitB->getCellID0() << " with trkhit = " << rawHit << " Cell ID = " <<  rawHit->getCellID0() << std::endl;     
                rawHit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;    
                rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhitB;                 
              }
              
            } 
          }    
#endif  
          
          
          
        }        
        else{ streamlog_out( DEBUG0 ) << "spacepoint discarded, because it is related to " << to.size() << "SimTrackerHits. It should be 2!\n"; } 
        
      }
      else{  // no composite spacepoint
        
        if( to.size() == 1){ // only take trackerHits, that have only one related SimHit
          
          SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(to.at(0));
          simHitTrkHit.push_back(std::make_pair(simhit, trkhit));
          
#ifdef MARLINTRK_DIAGNOSTICS_ON
          trkhit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;    
          trkhit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhit;  
#endif       
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
  /*                Take hits with the same MCP and create tracks from them                         */
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
  
  delete _trksystem ;
  
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
  
  ///////////////////////////////////////////////////////
  // check inputs 
  ///////////////////////////////////////////////////////
  
  if( mcp == 0 ){
    throw EVENT::Exception( std::string("TruthTracker::createTrack: MCParticle == NULL ")  ) ;
  }
  
  if ( hit_list.empty() ) return;
    
  HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz );
  streamlog_out( DEBUG3 ) << "\n MCParticle paramters: " 
  << " d0 " <<  hel.getD0()
  << " phi0 " << hel.getPhi0()
  << " omega "<< hel.getOmega()
  << " z0 "<< hel.getZ0()
  << " tanl "<< hel.getTanLambda()
  << " total number of hits = " << hit_list.size() 
  << "\n" << std::endl;

  
  //////////////////////////////////////////////////////////////////////////////////////
  // Only fit the hits before looping over SJA:FIXME: FTD may need special treatment
  //////////////////////////////////////////////////////////////////////////////////////
  
  std::vector<TrackerHit*> hit_list_inner_r;
  
  hit_list_inner_r.reserve(300);
  
  float delta_phi = 0.0;
  
  HelixTrack hel_copy(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz );
  
  for( unsigned ihit = 0; ihit < hit_list.size(); ++ihit){
    
    float x = hit_list[ihit]->getPosition()[0];
    float y = hit_list[ihit]->getPosition()[1];
    float z = hit_list[ihit]->getPosition()[2];
    
    delta_phi += fabsf(hel_copy.moveRefPoint(x, y, z));      
    
    if ( delta_phi < M_PI ) {
      hit_list_inner_r.push_back(hit_list[ihit]);
    }
    
  }
  
  if (hit_list_inner_r.size() < 3) {
    streamlog_out( MESSAGE ) << " Reject Track as the number of hits before the turn " << hit_list_inner_r.size() << " is less than 3 hits. Total delta phi before last hit = " << delta_phi << " . Total number of hits = " << hit_list.size()  << std::endl;
    return;
  }

  
  TrackImpl* Track = new TrackImpl ; 
  
  
  if( _FitTracksWithMarlinTrk ) {
    
    // setup initial dummy covariance matrix
    EVENT::FloatVec covMatrix;
    covMatrix.resize(15);
    
    for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
    }
    
    covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
    covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
    covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
    covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
    covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
    
    
      
    streamlog_out( DEBUG3 ) << "Create track with " << hit_list_inner_r.size() << " hits" << std::endl;
    
    // First sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to in. 
    sort( hit_list_inner_r.begin(), hit_list_inner_r.end(), TruthTracker::compare_r() );
    
    TrackStateImpl* prefit_trackState = 0;
    
    bool fit_backwards = IMarlinTrack::backward;
    
    MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
    
    try {
      
      int error = 0;
      
      if( _useMCParticleParametersFotInitOfFit ){
        
        hel.moveRefPoint(hit_list_inner_r.front()->getPosition()[0], hit_list_inner_r.front()->getPosition()[1], hit_list_inner_r.front()->getPosition()[2]);
        
        const float referencePoint[3] = { hel.getRefPointX() , hel.getRefPointY() , hel.getRefPointZ() };
        
        prefit_trackState = new TrackStateImpl( lcio::TrackState::AtIP, 
                                               hel.getD0(), 
                                               hel.getPhi0(), 
                                               hel.getOmega(), 
                                               hel.getZ0(), 
                                               hel.getTanLambda(), 
                                               covMatrix, 
                                               referencePoint) ;
        
        
        
        error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, hit_list_inner_r, Track, fit_backwards, prefit_trackState, _Bz);
                
      } else {
        
        error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, hit_list_inner_r, Track, fit_backwards, covMatrix, _Bz);
                
      }    
        
      if( error != IMarlinTrack::success || Track->getNdf() < 0 ) {       
        streamlog_out(MESSAGE) << "TruthTracker::createTrack: EVENT: << " << _current_evt_number << " >> Track fit returns error code " << error << " NDF = " << Track->getNdf() <<  ". Number of hits = "<< hit_list_inner_r.size() << std::endl;       
        return ;
      }
      
#ifdef MARLINTRK_DIAGNOSTICS_ON
      if (error != 0) {        
        void * dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();
      }        
#endif

      
      
    } catch (...) {
      
      streamlog_out(ERROR) << "TruthTracker::createTrack: EVENT: << " << _current_evt_number << " >> exception caught and rethown. MCParticle = " << mcp << std::endl;       
      
//      delete Track;
//      delete marlinTrk;
      
      throw ;
      
    }

      
    std::vector<std::pair<EVENT::TrackerHit* , double> > hits_in_fit ;  
    std::vector<std::pair<EVENT::TrackerHit* , double> > outliers ;
    std::vector<TrackerHit*> all_hits;    
    all_hits.reserve(300);

    marlinTrk->getHitsInFit(hits_in_fit);
            
    for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }
        
    MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);
    
    marlinTrk->getOutliers(outliers);
    
    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }
        
    MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);
    
    delete marlinTrk;
  
    
  } else {
    
    // use mcp pos and mom to set the track parameters
    HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz );
    
    hel.moveRefPoint(0.0, 0.0, 0.0);
    
    TrackStateImpl* ts = new TrackStateImpl();
    
    ts->setD0(hel.getD0());
    ts->setPhi(hel.getPhi0());
    ts->setOmega(hel.getOmega());
    ts->setZ0(hel.getZ0());
    ts->setTanLambda(hel.getTanLambda());
    
    float ref[3];
    
    ref[0] = hel.getRefPointX() ;
    ref[1] = hel.getRefPointY() ;
    ref[2] = hel.getRefPointZ() ;
    
    ts->setReferencePoint(ref);    
    
    ts->setLocation(lcio::TrackState::AtIP);
    
    Track->addTrackState(ts);
    
    MarlinTrk::addHitNumbersToTrack(Track, hit_list, true,  cellID_encoder);
    MarlinTrk::addHitNumbersToTrack(Track, hit_list, false, cellID_encoder);
    
    for( unsigned i=0; i<hit_list.size(); i++ ){
      
      Track->addHit( hit_list[i] );
      
    }
    
    
  }
  
  streamlog_out( DEBUG3 ) << "Add Track " << Track << " to collection related to mcp -> " << mcp << std::endl;  
  
  _trackVec->addElement(Track);
  
  LCRelationImpl* rel = new LCRelationImpl;
  rel->setFrom (Track);
  rel->setTo (mcp);
  rel->setWeight(1.0);
  _trackRelVec->addElement(rel);
  
  _nCreatedTracks++;
  
  
}


void TruthTracker::createTrack_old( MCParticle* mcp, UTIL::BitField64& cellID_encoder, std::vector<TrackerHit*>& hit_list ) {
  
  
  TrackImpl* Track = new TrackImpl ; 
  
  streamlog_out( DEBUG3 ) << "Create track with " << hit_list.size() << " hits" << std::endl;  
  
  // check that there are actually hits to fit
  if ( hit_list.empty() ) return;
  
  HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz );
  streamlog_out( DEBUG3 ) << "\n MCParticle paramters: " 
  << " d0 " <<  hel.getD0()
  << " phi0 " << hel.getPhi0()
  << " omega "<< hel.getOmega()
  << " z0 "<< hel.getZ0()
  << " tanl "<< hel.getTanLambda()
  << "\n" << std::endl;
  
  double d0;
  double phi0;
  double omega;
  double z0;
  double tanL;
  
  double chi2 = 0 ;
  int ndf = 0 ;
  
  float ref[3];
  
  TrackerHitVec added_hits;
  
  if( _FitTracksWithMarlinTrk ) {
    
    TrackStateImpl* prefit_trackState = 0;
    
    
    // setup initial dummy covariance matrix
    EVENT::FloatVec covMatrix;
    covMatrix.resize(15);
    
    for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
    }
    
    covMatrix[0]  = ( _initialTrackError_d0    ); //sigma_d0^2
    covMatrix[2]  = ( _initialTrackError_phi0  ); //sigma_phi0^2
    covMatrix[5]  = ( _initialTrackError_omega ); //sigma_omega^2
    covMatrix[9]  = ( _initialTrackError_z0    ); //sigma_z0^2
    covMatrix[14] = ( _initialTrackError_tanL  ); //sigma_tanl^2
    
    
    // First sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to in. 
    sort( hit_list.begin(), hit_list.end(), TruthTracker::compare_r() );
    
    
    // loop over all the hits and create a list consisting only 2D hits 
    
    TrackerHitVec twoD_hits;
    
    for (unsigned ihit=0; ihit < hit_list.size(); ++ihit) {
      
      // check if this a space point or 2D hit 
      if(BitSet32( hit_list[ihit]->getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] == false ){
        // then add to the list 
        twoD_hits.push_back(hit_list[ihit]);
        
      }
    }
    
    
    // now either take the MCParticle parameters directly, use a full prefit or create a helix from 3 2-D hits
    
    if( _useMCParticleParametersFotInitOfFit ){
      
      hel.moveRefPoint(hit_list.front()->getPosition()[0], hit_list.front()->getPosition()[1], hit_list.front()->getPosition()[2]);
      
      const float referencePoint[3] = { hel.getRefPointX() , hel.getRefPointY() , hel.getRefPointZ() };
      
      prefit_trackState = new TrackStateImpl( lcio::TrackState::AtIP, 
                                             hel.getD0(), 
                                             hel.getPhi0(), 
                                             hel.getOmega(), 
                                             hel.getZ0(), 
                                             hel.getTanLambda(), 
                                             covMatrix, 
                                             referencePoint) ;
      
    }         
    else { // use 2-D hits to start fit either simple helix construction or full prefit
      
      // check that there are enough 2-D hits to create a helix 
      if (twoD_hits.size() < 3) { // no chance to initialise print warning and return
        streamlog_out(WARNING) << "TruthTracker::createTrack Cannot create helix from less than 3 2-D hits" << std::endl;
        delete Track;
        return;
      }
      
      // make a helix from 3 hits to get a trackstate
      // SJA:FIXME: this may not be the optimal 3 hits to take in certain cases where the 3 hits are not well spread over the track length 
      const double* x1 = twoD_hits[0]->getPosition();
      const double* x2 = twoD_hits[ twoD_hits.size()/2 ]->getPosition();
      const double* x3 = twoD_hits.back()->getPosition();
      
      HelixTrack helixTrack( x1, x2, x3, _Bz, IMarlinTrack::backward );
      
      helixTrack.moveRefPoint(hit_list.back()->getPosition()[0], hit_list.back()->getPosition()[1], hit_list.back()->getPosition()[2]);
      
      const float referencePoint[3] = { helixTrack.getRefPointX() , helixTrack.getRefPointY() , helixTrack.getRefPointZ() };
      
      prefit_trackState = new TrackStateImpl( lcio::TrackState::AtIP, 
                                             helixTrack.getD0(), 
                                             helixTrack.getPhi0(), 
                                             helixTrack.getOmega(), 
                                             helixTrack.getZ0(), 
                                             helixTrack.getTanLambda(), 
                                             covMatrix, 
                                             referencePoint) ;
      
      
      // if performing full pre fit use prefit_trackState and improve with full fit
      
      if ( _create_prefit_using_MarlinTrk){
        
        
        
        MarlinTrk::IMarlinTrack* prefit_trk = _trksystem->createTrack();
        
#ifdef MARLINTRK_DIAGNOSTICS_ON
        
        void * dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();
        
#endif
        
        
        EVENT::TrackerHitVec::iterator it = hit_list.begin();
        
        streamlog_out(DEBUG2) << "Start Pre Fit: AddHits: number of hits to fit " << hit_list.size() << std::endl;
        
        int pre_fit_ndof_added = 0;
        TrackerHitVec pre_fit_added_hits;
        
        for( it = hit_list.begin() ; it != hit_list.end() ; ++it ) {
          
          TrackerHit* trkHit = *it;
          bool isSuccessful = false; 
          
          if( BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
            
            //Split it up and add both hits to the MarlinTrk
            const LCObjectVec rawObjects = trkHit->getRawHits();                    
            
            for( unsigned k=0; k< rawObjects.size(); k++ ){
              
              TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
              
              if( prefit_trk->addHit( rawHit ) == IMarlinTrack::success ){
                
                isSuccessful = true; //if at least one hit from the spacepoint gets added
                ++pre_fit_ndof_added;
              }
            }
          }
          else { // normal non composite hit
            
            if (prefit_trk->addHit( trkHit ) == IMarlinTrack::success ) {
              isSuccessful = true;
              pre_fit_ndof_added += 2;
            }
          }
          
          if (isSuccessful) {
            pre_fit_added_hits.push_back(trkHit);
          }
          
          
          else{
            streamlog_out(DEBUG4) << "Hit " << it - hit_list.begin() << " Dropped " << std::endl;          
          }
          
        }
        
        if( pre_fit_ndof_added < 8 ) {
          streamlog_out(DEBUG3) << "TruthTracker::Prefit failed: Cannot fit less with less than 8 degrees of freedom. Number of hits =  " << pre_fit_added_hits.size() << " ndof = " << pre_fit_ndof_added << std::endl;
          delete prefit_trk ;
          delete prefit_trackState;
          return;
        }
        
        prefit_trk->initialise( *prefit_trackState, _Bz, IMarlinTrack::backward ) ;
        
        delete prefit_trackState;
        
        int fit_status = prefit_trk->fit() ; 
        
        if( fit_status != 0 ){ 
          streamlog_out(DEBUG3) << "TruthTracker::Prefit failed: fit_status = " << fit_status << std::endl; 
          delete prefit_trk ;
          return;
        } else {
          
          prefit_trackState = new TrackStateImpl() ;
          
          const gear::Vector3D point(hit_list.back()->getPosition()[0],hit_list.back()->getPosition()[1],hit_list.back()->getPosition()[2]); // position of outermost hit
          
          int return_code = prefit_trk->propagate(point, *prefit_trackState, chi2, ndf ) ;
          
          // make sure that the track state can be propagated to the IP and that the NDF is not less than 0
          if ( return_code != 0 || ndf < 0 ) { 
            streamlog_out(WARNING) << "TruthTracker::createTrack Prefit: Track droped return_code for prefit propagation = " << return_code << " NDF = " << ndf << std::endl;
            delete prefit_trackState;
            delete prefit_trk;
            return;
          }
          
          prefit_trackState->setLocation(  lcio::TrackState::AtLastHit ) ;
          
          // blow up covariance matrix
          
          prefit_trackState->setCovMatrix(covMatrix);
          
        }
        
        delete prefit_trk;
        
      }
      
    }
    
    MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
    
    
    // used to keep track of whether hits are accepted by the fitter before fitting
    bool isSuccessful = false; 
    
    // needed to keep track of which hit should used for the first addandfit call
    int index_of_next_hit_to_fit = hit_list.size() - 1; 
    
    
    // to initialise we need to add at least one hit. 
    // given that we need to also work with COMPOSITE_SPACEPOINT hits we will simply added both strip hits in this case.
    // fit will then be called on the hits added and then addandFit will be used for the remaining hits to provide more feedback and control
    
    
    // start by trying to add the last hit and break at the first accepted hit
    for(int j=index_of_next_hit_to_fit; j >= 0; --j) {
      
      TrackerHit* trkHit = hit_list[j];
      
      --index_of_next_hit_to_fit; // first decrement the index of the next hit to fit
      
      // Check for spacepoints
      if( BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
        
        //Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();
        
        for( unsigned k=0; k< rawObjects.size(); k++ ){
          
          TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
          if( marlin_trk->addHit( rawHit ) == IMarlinTrack::success ){
            
            isSuccessful = true; //if at least one hit from the spacepoint gets added
            
          }          
        }
      } 
      else { // normal non-composite tracker hit
        
        isSuccessful = marlin_trk->addHit( trkHit ) == IMarlinTrack::success;
        
      }
      
      if ( isSuccessful ) { 
        added_hits.push_back( trkHit ); 
        // at this point one hit added successfully so break the for loop
        break; 
      }        
      
    } // end of reverse loop over hits 
    
    // We now have either one 2-D hit or two 1-D hits added. 
    // So proceed with initialisation.  
    
    
    streamlog_out( DEBUG3 ) << "\n Helix for prefit: " 
    <<  " d0 =  " << prefit_trackState->getD0()
    <<  " phi0 =  " << prefit_trackState->getPhi() 
    <<  " omega =  " << prefit_trackState->getOmega() 
    <<  " z0 =  " << prefit_trackState->getZ0() 
    <<  " tanl =  " << prefit_trackState->getTanLambda() 
    <<  " ref =  " <<  prefit_trackState->getReferencePoint()[0] << " " << prefit_trackState->getReferencePoint()[1] << " " << prefit_trackState->getReferencePoint()[2]
    << "\n" << std::endl;
    
    
    // set the initial track state for the track    
    marlin_trk->initialise( *prefit_trackState, _Bz, IMarlinTrack::backward ) ;
    
    delete prefit_trackState;
    
    // filter the first 1 or 2 hits 
    int fit_status = marlin_trk->fit() ; 
    
    streamlog_out(DEBUG4) << "fit_status = " << fit_status << std::endl ;
    
    // check that first hit is accepted by the fit, if this fails we bail here as there is little chance of recovering. 
    
    if ( fit_status != IMarlinTrack::success ) { // no chance to initialise print warning and return
      streamlog_out(WARNING) << "TruthTracker::createTrack Initial Hit not accepted by the fit, track droped." << std::endl;
      delete Track;
      delete marlin_trk;
      return;
    }
    
    
    // now used addAndFit to add the remaining hits 
    for(int j=index_of_next_hit_to_fit; j >= 0; --j) {
      
      TrackerHit* trkHit = hit_list[j];
      
      isSuccessful = false;
      
      // Check for spacepoints
      if( BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
        
        //Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();        
        
        for( int k=rawObjects.size()-1 ; k >= 0; --k ){
          
          TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
          
          double chi_incr=DBL_MAX;
          
          if( marlin_trk->addAndFit(rawHit,chi_incr) == IMarlinTrack::success ){
            
            isSuccessful = true; //if at least one hit from the spacepoint gets added
            
          }
        }      
      }
      else { // normal non composite tracker hit
        
        double chi_incr=DBL_MAX;
        
        isSuccessful =  marlin_trk->addAndFit(trkHit,chi_incr) == IMarlinTrack::success ;                    
        
      }
      if ( isSuccessful ) added_hits.push_back( trkHit );
    }  
    
    
    const gear::Vector3D point(0.,0.,0.); // nominal IP
    
    TrackStateImpl* trkStateIP = new TrackStateImpl() ;
    
    int return_code = marlin_trk->propagate(point, *trkStateIP, chi2, ndf ) ;
    
    // make sure that the track state can be propagated to the IP and that the NDF is not less than 0
    if ( return_code != 0 || ndf < 0 ) { 
      streamlog_out(WARNING) << "TruthTracker::createTrack: Track droped return_code for propagation = " << return_code << " NDF = " << ndf << std::endl;
      delete trkStateIP;
      delete marlin_trk;
      return;
    }
    
    trkStateIP->setLocation(  lcio::TrackState::AtIP ) ;
    Track->trackStates().push_back(trkStateIP);
    Track->setChi2(chi2);
    Track->setNdf(ndf);
    
    std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
    
    marlin_trk->getHitsInFit(hits_in_fit);
    
    
    EVENT::TrackerHit* first_trkhit = hits_in_fit.front().first;
    EVENT::TrackerHit* last_trkhit = hits_in_fit.back().first;
    
    // now loop over the hits and get the track State there 
    
    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      
      EVENT::TrackerHit* trkhit = hits_in_fit[ihit].first;
      
      if ( trkhit == 0 ) {
        throw EVENT::Exception( std::string("TruthTracker::createTrack: TrackerHit pointer == NULL ")  ) ;
      }
      
      TrackStateImpl* ts = new TrackStateImpl;
      return_code = marlin_trk->getTrackState(trkhit, *ts, chi2, ndf ) ;
      
      if(return_code !=MarlinTrk::IMarlinTrack::success){
        streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> TruthTracker::createTrack:  could not get TrackState at Hit index " << ihit << std::endl ;
      }
      
      // SJA:FIXME: if the hit is 1D then the pivot will be set to 0,0,0 : what should we do? Probably propagate to the Space Point if there is one ....      
      // set the location for the first hit 
      if ( trkhit == first_trkhit ) {                
        ts->setLocation(  lcio::TrackState::AtFirstHit ) ;
      } 
      // set the location for the last hit 
      else if( trkhit == last_trkhit ) {
        ts->setLocation(  lcio::TrackState::AtLastHit ) ;      
      } 
      // set the location for the ihit'th  hit with an offset of lcio::TrackState::AtFirstHit + 1000
      else {
        ts->setLocation( lcio::TrackState::AtFirstHit + 1000 + ihit  ) ;        
      }
      
      Track->trackStates().push_back(ts);
      
    }
    
    
    // now get the state at the calo face
    
    TrackStateImpl* trkStateCalo = new TrackStateImpl;
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.reset() ;  // reset to 0
    
    encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::ECAL ;
    encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::barrel;
    encoder[lcio::ILDCellID0::layer]  = 0 ;
    
    int detElementID = 0;
    return_code = marlin_trk->propagateToLayer(encoder.lowWord(), last_trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    
    if (return_code == MarlinTrk::IMarlinTrack::no_intersection ) { // try forward or backward
      if (trkStateIP->getTanLambda()>0) {
        encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::fwd;
      }
      else{
        encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::bwd;
      }
      return_code = marlin_trk->propagateToLayer(encoder.lowWord(), last_trkhit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
    }
    
    if (return_code !=MarlinTrk::IMarlinTrack::success ) {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at Calo Face: return_code = " << return_code << std::endl ;
    }
    
    
    trkStateCalo->setLocation(  lcio::TrackState::AtCalorimeter ) ;      
    Track->trackStates().push_back(trkStateCalo);
    
    delete marlin_trk;    
    
    
  }
  else {
    
    
    // use mcp pos and mom to set the track parameters
    HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz );
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
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() <<  " collection absent" << std::endl;     
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







