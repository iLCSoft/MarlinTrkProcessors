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
  _description = "Creates Track Collection from MC Truth" ;
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  // register steering parameters: name, description, class-variable, default value
  
  registerInputCollection( LCIO::TRACKERHIT,
                          "VTXCollectionName" , 
                          "Name of the collection of VTX tracker hits"  ,
                          _colNameTrkHitsVTX ,
                          std::string("VTXTrackerHits") ) ;
  
  registerInputCollection( LCIO::TRACKERHIT,
                          "FTDCollectionName" , 
                          "Name of the collection of FTD tracker hits"  ,
                          _colNameTrkHitsFTD ,
                          std::string("FTDTrackerHits") ) ;
  
  registerInputCollection( LCIO::TRACKERHIT,
                          "SITCollectionName" , 
                          "Name of the collection of SIT tracker hits"  ,
                          _colNameTrkHitsSIT ,
                          std::string("SITTrackerHits") ) ;
  
  registerInputCollection( LCIO::TRACKERHIT,
                          "TPCCollectionName" , 
                          "Name of the collection of TPC tracker hits"  ,
                          _colNameTrkHitsTPC ,
                          std::string("TPCTrackerHits") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "VTXSimTrkhitCollectionName" ,
                          "Name of the VTX sim trk hit collection" ,
                          _colNameSimTrkHitsVTX ,
                          std::string("VXDCollection") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "SITSimTrkhitCollectionName" ,
                          "Name of the SIT sim trk hit collection" ,
                          _colNameSimTrkHitsSIT ,
                          std::string("SITCollection") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "FTDSimTrkhitCollectionName" ,
                          "Name of the FTD sim trk hit collection" ,
                          _colNameSimTrkHitsFTD ,
                          std::string("FTDCollection") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "TPCSimTrkhitCollectionName" , 
                          "Name of the TPC sim trk hit collection"  ,
                          _colNameSimTrkHitsTPC,
                          std::string("TPCCollection") ) ;
  
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
  
  registerProcessorParameter( "nEventPrintout",
                             "Print out progress every N events ",
                             _nEventPrintout,
                             int(1000));
  
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
  
  
  registerInputCollection( LCIO::LCRELATION,"VXDTrackerHitRelInputCollection" , 
                          "Name of the rel collection for VXD TrackerHit collection"  ,
                          _vxdTrackerHitRelInputColName ,
                          std::string("VTXTrackerHitRelations") ) ;
  
  registerInputCollection( LCIO::LCRELATION,"FTDTrackerHitRelInputCollection" , 
                          "Name of the rel collection for FTD TrackerHit collection"  ,
                          _ftdTrackerHitRelInputColName ,
                          std::string("FTDTrackerHitRelations") ) ;
  
  registerInputCollection( LCIO::LCRELATION,"SITTrackerHitRelInputCollection" , 
                          "Name of the rel collection for SIT TrackerHit collection"  ,
                          _sitTrackerHitRelInputColName ,
                          std::string("SITTrackerHitRelations") ) ;
  
  registerInputCollection( LCIO::LCRELATION,"TPCTrackerHitRelInputCollection" , 
                          "Name of the rel collection for TPC TrackerHit collection"  ,
                          _tpcTrackerHitRelInputColName ,
                          std::string("TPCTrackerHitRelations") ) ;
  
  registerInputCollection( LCIO::LCRELATION,"SETTrackerHitRelInputCollection" , 
                          "Name of the rel collection for SET TrackerHit collection"  ,
                          _setTrackerHitRelInputColName ,
                          std::string("SETTrackerHitRelations") ) ;
  
  registerInputCollection( LCIO::LCRELATION,"ETDTrackerHitRelInputCollection" , 
                          "Name of the rel collection for ETD TrackerHit collection"  ,
                          _etdTrackerHitRelInputColName ,
                          std::string("ETDTrackerHitRelations") ) ;

  
  
  _n_run = 0 ;
  _n_evt = 0 ;
  
  _colMCP = NULL ;
  _nMCP = 0 ;
  
  _VTXSimHits = NULL ;
  _SITSimHits = NULL ;
  _FTDSimHits = NULL ;
  _TPCSimHits = NULL ;
  
  _nVTXSimHits = 0 ;
  _nSITSimHits = 0 ;
  _nFTDSimHits = 0 ;
  _nTPCSimHits = 0 ;
  
  _VTXTrkHits = NULL ;
  _SITTrkHits = NULL ;
  _FTDTrkHits = NULL ;
  _TPCTrkHits = NULL ;
  
  _nVTXTrkHits = 0 ;
  _nSITTrkHits = 0 ;
  _nFTDTrkHits = 0 ;
  _nTPCTrkHits = 0 ;
  
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
  
  
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  
  //  if( _n_evt % _nEventPrintout == 0 ) streamlog_out(DEBUG) << "   processing event: " << _n_evt 
  //              << std::endl ;
  
  streamlog_out(DEBUG) << "   processing event: " << _n_evt << std::endl ;
  
  _colMCP = NULL ;
  _nMCP = 0 ;
  
  _VTXSimHits = NULL ;
  _SITSimHits = NULL ;
  _FTDSimHits = NULL ;
  _TPCSimHits = NULL ;
  
  _nVTXSimHits = 0 ;
  _nSITSimHits = 0 ;
  _nFTDSimHits = 0 ;
  _nTPCSimHits = 0 ;
  
  _VTXTrkHits = NULL ;
  _SITTrkHits = NULL ;
  _FTDTrkHits = NULL ;
  _TPCTrkHits = NULL ;
  
  _nVTXTrkHits = 0 ;
  _nSITTrkHits = 0 ;
  _nFTDTrkHits = 0 ;
  _nTPCTrkHits = 0 ;
  
  // get the input collections and set the memeber variable pointers to them
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
  
  std::vector< std::pair<SimTrackerHit*, TrackerHit* > > simTrackerHitAndTrackerHitList;
  
  for(int i=0; i< _nVTXTrkHits ; ++i) {
    TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_VTXTrkHits->getElementAt( i ));      
    if(  this->getSimHits(trkhit)->size() == 1 )  {
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(this->getSimHits(trkhit)->at(0));
      simTrackerHitAndTrackerHitList.push_back(std::make_pair(simhit, trkhit));
    }
  }
  
  for(int i=0; i< _nSITTrkHits ; ++i) {
    TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_SITTrkHits->getElementAt( i ));      
    if(  this->getSimHits(trkhit)->size() == 1 )  {
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(this->getSimHits(trkhit)->at(0));
      simTrackerHitAndTrackerHitList.push_back(std::make_pair(simhit, trkhit));
    }
  }
  
  for(int i=0; i< _nFTDTrkHits ; ++i) {
    TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_FTDTrkHits->getElementAt( i ));      
    if(  this->getSimHits(trkhit)->size() == 1 )  {
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(this->getSimHits(trkhit)->at(0));
      simTrackerHitAndTrackerHitList.push_back(std::make_pair(simhit, trkhit));
    }
  }
  
  for(int i=0; i< _nTPCTrkHits ; ++i) {
    TrackerHit * trkhit = dynamic_cast<TrackerHit*>(_TPCTrkHits->getElementAt( i ));      
    if(  this->getSimHits(trkhit)->size() == 1 )  {
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(this->getSimHits(trkhit)->at(0));
      simTrackerHitAndTrackerHitList.push_back(std::make_pair(simhit, trkhit));
    }
  }
  
  
  
  // now order the hits by MCParticle and then in time  
//  std::sort(trkhits.begin(), trkhits.end(), TruthTracker::SimTrackerHitSortPredicate() );
  
  std::sort( simTrackerHitAndTrackerHitList.begin(),simTrackerHitAndTrackerHitList.end(), SimTrackerHitSortPredicate() );
  

  // loop over all tracker hits and check which MCPartilce created them 
  // check each TrackerHit to see that it has only one SimTrackerHit
  std::vector<TrackerHit*> trkhits; 

  std::vector< std::pair<SimTrackerHit*, TrackerHit* > >::iterator it;
  for ( it = simTrackerHitAndTrackerHitList.begin(); it != simTrackerHitAndTrackerHitList.end(); ++it) {
    trkhits.push_back((*it).second);
  }

  streamlog_out( DEBUG4 ) << "Number of Tracker hits = " << trkhits.size() << std::endl;     
  
  for(unsigned int i=0; i< trkhits.size() ; ++i) {
    SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*>(this->getSimHits(trkhits[i])->at(0)) ;
    streamlog_out( DEBUG1 ) << "Tracker hit: [" << i << "] = " << trkhits[i] << "  mcp = " << simHit->getMCParticle() << " time = " << simHit->getTime() << " cellid = " << simHit->getCellID0() << std::endl;     
  }

  streamlog_out( DEBUG4 ) << "Add Tracker hits to Tracks" << std::endl;
  
  _hit_list.clear();
  
  if( trkhits.size() > 0) {
    MCParticle* mcplast = NULL;
    
    //    for(unsigned int i=0; i< trkhits.size() ; ++i) {
      for(unsigned int i=0; i< trkhits.size() ; ++i) {

        SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(this->getSimHits(trkhits[i])->at(0));
        
        MCParticle* mcp = simhit->getMCParticle();
      double const* p    = mcp->getMomentum() ;
      float  const pt2   = p[0]*p[0] + p[1]*p[1] ;
      //float  const pmag2 =  p[0]*p[0] + p[1]*p[1] + p[2]*p[2] ; 
            
      if ( i == 0 ) {
        mcplast = mcp ;
      }
      
      if( mcp != mcplast ) { 
        // new MCParticle

        if (_hit_list.size() >= 3) {
          // create track from vector of hits                           
          streamlog_out( DEBUG2 ) << "Create New Track for MCParticle " << mcplast << std::endl;
          this->createTrack(mcplast, cellID_encoder);
          _hit_list.clear();      
        }
        else { 
          // clear the hit list to start a new list of hits from the same mcparticle
          _hit_list.clear();   
        }
      }
      
      if( pt2  > (_MCpThreshold*_MCpThreshold) ) {
        // if momentum is greater than cut add hit to list
        streamlog_out( DEBUG3 ) << "Add hit from det " <<  simhit->getCellID0()  << " to track from MCParticle " << mcp << " : current number of hits = " << _hit_list.size() << std::endl;
        _hit_list.push_back(trkhits.at(i)) ;
      }
      
      // set last mcparticle 
      mcplast = mcp;
      
    } // end of loop over hits
    
    // check if there is still a track to be created 
    if( _hit_list.size() >= 3 ) { 
      // then create a new track
      streamlog_out( DEBUG4 ) << "Create New Track for Last MCParticle " << mcplast << std::endl;
      this->createTrack(mcplast, cellID_encoder);

      
      _hit_list.clear();      
    }
    
    
  }    
  
  evt->addCollection( _trackVec , _output_track_col_name) ;
  evt->addCollection( _trackRelVec , _output_track_rel_name) ;
  
  
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

LCCollection* TruthTracker::GetCollection(  LCEvent * evt, std::string colName ){
  
  LCCollection* col = NULL;
  
  int nElements = 0;
  
  try {
    col = evt->getCollection( colName.c_str() ) ;
    nElements = col->getNumberOfElements()  ;
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() << " collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e) {
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;     
  }
  
  return col; 
  
}

//LCRelationNavigator* TruthTracker::GetRelations( EVENT::LCEvent * evt , std::string RelName ) {
//  
//  LCRelationNavigator* nav = NULL ;
//  
//  try {
//    nav = new LCRelationNavigator(evt->getCollection( RelName.c_str() ));
//    streamlog_out( DEBUG2 ) << "TruthTracker --> " << RelName << " track relation collection in event = " << nav << std::endl;
//  }
//  catch(DataNotAvailableException &e) {
//    streamlog_out( DEBUG2 ) << "TruthTracker --> " << RelName.c_str() << " track relation collection absent in event" << std::endl;     
//  }
//  
//  return nav;
//  
//}


void TruthTracker::SetupInputCollections( LCEvent * evt ) {
  
  _colMCP = GetCollection( evt, _colNameMC) ;
  streamlog_out( DEBUG2 ) << "number of MCParticles: " << _nMCP << std::endl ;
  
  _VTXSimHits = GetCollection( evt, _colNameSimTrkHitsVTX ) ;
  _SITSimHits = GetCollection( evt, _colNameSimTrkHitsSIT ) ;
  _FTDSimHits = GetCollection( evt, _colNameSimTrkHitsFTD ) ;
  _TPCSimHits = GetCollection( evt, _colNameSimTrkHitsTPC ) ;
  
  if(_VTXSimHits) _nVTXSimHits = _VTXSimHits->getNumberOfElements() ; 
  if(_SITSimHits) _nSITSimHits = _SITSimHits->getNumberOfElements() ;
  if(_FTDSimHits) _nFTDSimHits = _FTDSimHits->getNumberOfElements() ;
  if(_TPCSimHits) _nTPCSimHits = _TPCSimHits->getNumberOfElements() ;
  
  streamlog_out( DEBUG2 ) << "number of VTXSimHits: " << _nVTXSimHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of SITSimHits: " << _nSITSimHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of FTDSimHits: " << _nFTDSimHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of TPCSimHits: " << _nTPCSimHits << std::endl ;
  
  _VTXTrkHits = GetCollection( evt, _colNameTrkHitsVTX ) ;
  _SITTrkHits = GetCollection( evt, _colNameTrkHitsSIT ) ;
  _FTDTrkHits = GetCollection( evt, _colNameTrkHitsFTD ) ;
  _TPCTrkHits = GetCollection( evt, _colNameTrkHitsTPC ) ;
  
  if(_VTXTrkHits) _nVTXTrkHits = _VTXTrkHits->getNumberOfElements() ; 
  if(_SITTrkHits) _nSITTrkHits = _SITTrkHits->getNumberOfElements() ;
  if(_FTDTrkHits) _nFTDTrkHits = _FTDTrkHits->getNumberOfElements() ;
  if(_TPCTrkHits) _nTPCTrkHits = _TPCTrkHits->getNumberOfElements() ;
  
  streamlog_out( DEBUG2 ) << "number of VTXTrkHits: " << _nVTXTrkHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of SITTrkHits: " << _nSITTrkHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of FTDTrkHits: " << _nFTDTrkHits << std::endl ;
  streamlog_out( DEBUG2 ) << "number of TPCTrkHits: " << _nTPCTrkHits << std::endl ;
  
  _navVXDTrackerHitRel = GetRelations( evt, _vxdTrackerHitRelInputColName );
  _navSITTrackerHitRel = GetRelations( evt, _sitTrackerHitRelInputColName );
  _navFTDTrackerHitRel = GetRelations( evt, _ftdTrackerHitRelInputColName );
  _navTPCTrackerHitRel = GetRelations( evt, _tpcTrackerHitRelInputColName );
  _navSETTrackerHitRel = GetRelations( evt, _setTrackerHitRelInputColName );
  _navETDTrackerHitRel = GetRelations( evt, _etdTrackerHitRelInputColName );
  
  _hit_rels_map.clear();
  
  if (_navVXDTrackerHitRel) _hit_rels_map[ILDDetID::VXD]=_navVXDTrackerHitRel;
  if (_navSITTrackerHitRel) _hit_rels_map[ILDDetID::SIT]=_navSITTrackerHitRel;
  if (_navFTDTrackerHitRel) _hit_rels_map[ILDDetID::FTD]=_navFTDTrackerHitRel;
  if (_navTPCTrackerHitRel) _hit_rels_map[ILDDetID::TPC]=_navTPCTrackerHitRel;
  if (_navSETTrackerHitRel) _hit_rels_map[ILDDetID::SET]=_navSETTrackerHitRel;
  if (_navETDTrackerHitRel) _hit_rels_map[ILDDetID::ETD]=_navETDTrackerHitRel;
  
}


void TruthTracker::createTrack( MCParticle* mcp, UTIL::BitField64& cellID_encoder ) {
  
  TrackImpl* Track = new TrackImpl ; 
  
  streamlog_out( DEBUG3 ) << "Create track with " << _hit_list.size() << " hits" << std::endl;  
  
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
    sort(_hit_list.begin(), _hit_list.end(), TruthTracker::compare_r() );

    for(unsigned int j=0; j<_hit_list.size(); ++j) {
      if( IMarlinTrack::success == marlin_trk->addHit( _hit_list[j] ) ) added_hits.push_back( _hit_list[j] );
    }  

    if (added_hits.size() < 4) {
      streamlog_out( DEBUG3 ) << "Less than 4 hits: only " << _hit_list.size() << " hits skip" << std::endl;  
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
    
  }
  
  
  
  std::map<int, int> hitNumbers; 
  
  hitNumbers[lcio::ILDDetID::VXD] = 0;
  hitNumbers[lcio::ILDDetID::SIT] = 0;
  hitNumbers[lcio::ILDDetID::FTD] = 0;
  hitNumbers[lcio::ILDDetID::TPC] = 0;
  
  
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
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 2 ] = 0;
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 2 ] = 0;
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 1 ] = hitNumbers[lcio::ILDDetID::VXD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 1 ] = hitNumbers[lcio::ILDDetID::FTD];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 1 ] = hitNumbers[lcio::ILDDetID::SIT];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ] = hitNumbers[lcio::ILDDetID::TPC];
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 1 ] = 0;
  Track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 1 ] = 0;

  _trackVec->addElement(Track);

  LCRelationImpl* rel = new LCRelationImpl;
  rel->setFrom (Track);
  rel->setTo (mcp);
  rel->setWeight(1.0);
  _trackRelVec->addElement(rel);

  
}

const LCObjectVec* TruthTracker::getSimHits( TrackerHit* trkhit, const FloatVec* weights ){
  
  std::map<int, LCRelationNavigator*>::iterator it;
  
  it = _hit_rels_map.find( getDetectorID(trkhit));
  
  if (it != _hit_rels_map.end()) {
    if(weights) weights = &(it->second->getRelatedToWeights(trkhit));
    return &(it->second->getRelatedToObjects(trkhit));
  }
  else{
    streamlog_out(ERROR) << "Relations not present for Hit Collection from Detector ID " << this->getDetectorID(trkhit) << "  : exit(1) "<< std::endl;
    exit(1);
  }
  
}

LCRelationNavigator* TruthTracker::GetRelations(LCEvent * evt , std::string RelName ) {
  
  LCRelationNavigator* nav = NULL ;
  
  try{
    nav = new LCRelationNavigator(evt->getCollection( RelName.c_str() ));
    streamlog_out( DEBUG2 ) << "TruthTracker --> " << RelName << " track relation collection in event = " << nav << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG2 ) << "TruthTracker --> " << RelName.c_str() << " track relation collection absent in event" << std::endl;     
  }
  
  return nav;
  
}
