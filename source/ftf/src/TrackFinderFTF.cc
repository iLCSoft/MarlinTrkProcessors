#include "TrackFinderFTF.h"
#include <iostream>
#include <algorithm>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <EVENT/Track.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"

#include <cmath>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include <gear/BField.h>
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>

#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif


#include "TrackFinder.h"

using namespace lcio ;
//using namespace marlin ;

using namespace MarlinTrk ;

using namespace ftf;

TrackFinderFTF aTrackFinderFTF ;


TrackFinderFTF::TrackFinderFTF() : Processor("TrackFinderFTF") {
  
  // modify processor description
  _description = "TrackFinderFTF uses FTF to find tracks" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  
  // Input Collections
  // ^^^^^^^^^^^^^^^^^
  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "InputVXDTrackerHitCollectionName",
                          "VTX Hit Collection Name",
                          _input_vxd_hits_col_name,
                          std::string("VTXTrackerHits"));
  
  
  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "InputFTDPXLTrackerHitCollectionName",
                          "FTD Pixel Hit Collection Name",
                          _input_ftd_pixel_hits_col_name,
                          std::string("FTDPixelTrackerHits"));  
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "InputFTDSSPTrackerHitCollectionName",
                          "FTD FTDSpacePoint Collection Name",
                          _input_ftd_spacepoint_hits_col_name,
                          std::string("FTDSpacePoints"));  
  
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "InputSITTrackerHitCollectionName",
                          "SIT Hit Collection Name",
                          _input_sit_hits_col_name,
                          std::string("SITTrackerHits"));  

  
//  registerInputCollection( LCIO::TRACKERHITPLANE,
//                          "InputVXDTrackerHitCollectionName" , 
//                          "Name of the input tracker hit collection"  ,
//                          _input_vxd_hits_col_name ,
//                          std::string("VTXTrackerHits") ) ;
//  
//  registerInputCollection( LCIO::TRACKERHIT,
//                          "InputTrackerHitCollectionName" , 
//                          "Name of the input tracker hit collection"  ,
//                          _input_sit_hits_col_name ,
//                          std::string("SITTrackerHits") ) ;
//  
//  registerInputCollection( LCIO::TRACKERHITPLANE,
//                          "InputTrackerHitCollectionName" , 
//                          "Name of the input tracker hit collection"  ,
//                          _input_ftd_hits_col_name ,
//                          std::string("FTDTrackerHits") ) ;
//
//  registerInputCollection( LCIO::TRACKERHITPLANE,
//                          "InputTrackerHitCollectionName" , 
//                          "Name of the input tracker hit collection"  ,
//                          _input_ftd_hits_col_name ,
//                          std::string("FTDTrackerHits") ) ;
//
//  
//  registerInputCollection( LCIO::TRACKERHIT,
//                          "InputTPCTrackerHitCollectionName" , 
//                          "Name of the input tracker hit collection"  ,
//                          _input_tpc_hits_col_name ,
//                          std::string("TPCTrackerHits") ) ;
//  
//  registerInputCollection( LCIO::TRACKERHITPLANE,
//                          "InputTrackerHitCollectionName" , 
//                          "Name of the input tracker hit collection"  ,
//                          _input_set_hits_col_name ,
//                          std::string("SETTrackerHits") ) ;
  
  registerOutputCollection( LCIO::TRACK,
                           "OutputTrackCollectionName" , 
                           "Name of the output track collection"  ,
                           _output_track_col_name ,
                           std::string("FTFTracks") ) ;
  
    
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
  
  registerProcessorParameter( "MaxChi2PerHit",
                             "Maximum Chi-squared value allowed when assigning a hit to a track",
                             _maxChi2PerHit,
                             double(1.e2));
  
  
#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  registerOptionalParameter("RunMarlinTrkDiagnostics", "Run MarlinTrk Diagnostics. MarlinTrk must be compiled with MARLINTRK_DIAGNOSTICS_ON defined", _runMarlinTrkDiagnostics, bool(false));
  
  registerOptionalParameter("DiagnosticsName", "Name of the root file and root tree if running Diagnostics", _MarlinTrkDiagnosticsName, std::string("TruthTrackerDiagnostics"));    
  
#endif

  
}


void TrackFinderFTF::init() { 
  
  streamlog_out(DEBUG) << "   init called  " 
  << std::endl ;
   
  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

  // usually a good idea to
  printParameters() ;

  _Bz = marlin::Global::GEAR->getBField().at( gear::Vector3D(0., 0., 0.) ).z();    //The 
  
  try {
    
    const gear::ZPlanarParameters& pVXDDetMain = marlin::Global::GEAR->getVXDParameters();
    const gear::ZPlanarLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();
    _nlayers_vxd = pVXDLayerLayout.getNLayers();; 
    
  }
  catch (gear::UnknownParameterException& e) {
    streamlog_out( MESSAGE ) << "  TrackFinderFTF - VXD missing in gear file exit(1) " << std::endl ;
    exit(1);
  }
  
  try {
    
    const gear::ZPlanarParameters& pSITDetMain = marlin::Global::GEAR->getSITParameters();
    const gear::ZPlanarLayerLayout& pSITLayerLayout = pSITDetMain.getZPlanarLayerLayout();
    _nlayers_sit = pSITLayerLayout.getNLayers();; 
    
  }
  catch (gear::UnknownParameterException& e) {
    streamlog_out( MESSAGE ) << "  TrackFinderFTF - SIT missing in gear file exit(1) " << std::endl ;
    exit(1);
  }
  
//  try {
//    
//    const gear::TPCParameters& gearTPC = marlin::Global::GEAR->getTPCParameters() ;
//    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;
//    _nlayers_tpc = padLayout.getNRows();
//    
//  } catch (gear::UnknownParameterException& e) {
//    streamlog_out( MESSAGE ) << "  TrackFinderFTF - TPC missing in gear file exit(1)" << std::endl ;
//    exit(1);
//  }
//  
//  try{
//    
//    const gear::ZPlanarParameters& pSETDetMain = marlin::Global::GEAR->getSETParameters();
//    const gear::ZPlanarLayerLayout& pSETLayerLayout = pSETDetMain.getZPlanarLayerLayout();
//    _nlayers_sit = pSETLayerLayout.getNLayers();; 
//    
//  }
//  catch (gear::UnknownParameterException& e) {
//    streamlog_out( MESSAGE ) << "  TrackFinderFTF - SET missing in gear file exit(1) " << std::endl ;
//    exit(1);
//  }

  
  
  _trackFinder = new ftf::TrackFinder;
  _trackFinder->para.setDefaults();
  
  this->setFTFParameters( &(_trackFinder->para) );
  
  int maxHits        = 30000 ;
	int maxTracks      =  1000 ;
  
	_trackFinder->hit        = new Hit[maxHits] ;
	_trackFinder->track      = new ftf::Track[maxTracks] ;
	_trackFinder->maxTracks  = maxTracks ;
  
  _n_run = 0 ;
  _n_evt = 0 ;
 
  
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

  
}

void TrackFinderFTF::processRunHeader( LCRunHeader* ) {
  
  ++_n_run ;
} 

void TrackFinderFTF::processEvent( LCEvent * evt ) { 
  
  
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  
  streamlog_out(DEBUG) << "   processing event: " << _n_evt 
  << std::endl ;
  
  _input_vxd_hits_col = this->GetCollection( evt, _input_vxd_hits_col_name ) ;
  
  _input_sit_hits_col = this->GetCollection( evt, _input_sit_hits_col_name ) ;
  
  _input_ftd_pxl_hits_col = this->GetCollection( evt, _input_ftd_pixel_hits_col_name ) ;

  _input_ftd_spp_hits_col = this->GetCollection( evt, _input_ftd_spacepoint_hits_col_name ) ;
  
//  _input_tpc_hits_col = this->GetCollection( evt, _input_tpc_hits_col_name ) ;
//  
//  _input_set_hits_col = this->GetCollection( evt, _input_set_hits_col_name ) ;
//  
  _trackFinder->reset ( ) ;
  hitmap.clear();
  
  // establish the track collection that will be created 
  LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
  
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackVec->setFlag( trkFlag.getFlag()  ) ;
  
  
  int nhits = 0;
  
  
  nhits = _input_vxd_hits_col ? _input_vxd_hits_col->getNumberOfElements() : 0;
  
  int counter = 0 ;
  
  // loop over the input tacks and refit using KalTest    
  for(int i=0; i< nhits ; ++i) {
    
    TrackerHit* trkhit = dynamic_cast<TrackerHit*>( _input_vxd_hits_col->getElementAt( i ) ) ;
    
    const gear::Vector3D point(trkhit->getPosition()); 
    
    const double phi = this->angular_range_2PI(point.phi());
    
    if ( phi < _trackFinder->para.phiMin ) continue ;
    if ( phi > _trackFinder->para.phiMax ) continue ;
    
    const double eta = -std::log(tan(point.theta()/2.));
    
    if ( eta < _trackFinder->para.etaMin ) continue ;
    if ( eta > _trackFinder->para.etaMax ) continue ;
    
    _trackFinder->hit[counter].id       = counter ;
    _trackFinder->hit[counter].row      = this->getLayerID(trkhit);
    _trackFinder->hit[counter].x        = point.x();
    _trackFinder->hit[counter].y        = point.y();
    _trackFinder->hit[counter].z        = point.z();
    _trackFinder->hit[counter].dx       = 0.01;
    _trackFinder->hit[counter].dy       = 0.01;
    _trackFinder->hit[counter].dz       = 0.01;
    
    hitmap[ _trackFinder->hit[counter].id ] = trkhit;
    
    ++counter;
    
  }    
  
  nhits = _input_sit_hits_col ? _input_sit_hits_col->getNumberOfElements() : 0;

  // loop over the input tacks and refit using KalTest    
  for(int i=0; i< nhits ; ++i) {
    
    TrackerHit* trkhit = dynamic_cast<TrackerHit*>( _input_sit_hits_col->getElementAt( i ) ) ;
    
    const gear::Vector3D point(trkhit->getPosition()); 
    
    const double phi = this->angular_range_2PI(point.phi());
    
    if ( phi < _trackFinder->para.phiMin ) continue ;
    if ( phi > _trackFinder->para.phiMax ) continue ;
    
    const double eta   = -std::log(tan(point.theta()/2.));
    
    if ( eta < _trackFinder->para.etaMin ) continue ;
    if ( eta > _trackFinder->para.etaMax ) continue ;
    
    _trackFinder->hit[counter].id       = counter ;
    _trackFinder->hit[counter].row      = this->getLayerID(trkhit)+_nlayers_vxd;
    _trackFinder->hit[counter].x        = point.x();
    _trackFinder->hit[counter].y        = point.y();
    _trackFinder->hit[counter].z        = point.z();
    _trackFinder->hit[counter].dx       = trkhit->getCovMatrix()[0];
    _trackFinder->hit[counter].dy       = trkhit->getCovMatrix()[1];
    _trackFinder->hit[counter].dz       = trkhit->getCovMatrix()[2];
    
    hitmap[ _trackFinder->hit[counter].id ] = trkhit;
    
    ++counter;
    
  }    

  //  nhits = _input_ftd_hits_col ? _input_ftd_hits_col->getNumberOfElements() : 0;  
  //
//
//  // loop over the input tacks and refit using KalTest    
//  for(int i=0; i< nhits ; ++i) {
//    
//    TrackerHit* trkhit = dynamic_cast<TrackerHit*>( _input_ftd_hits_col->getElementAt( i ) ) ;
//    
//    const gear::Vector3D point(trkhit->getPosition()); 
//    
//    const double phi = this->angular_range_2PI(point.phi());
//    
//    if ( phi < _trackFinder->para.phiMin ) continue ;
//    if ( phi > _trackFinder->para.phiMax ) continue ;
//    
//    const double eta   = -std::log(tan(point.theta()/2.));
//    
//    if ( eta < _trackFinder->para.etaMin ) continue ;
//    if ( eta > _trackFinder->para.etaMax ) continue ;
//    
//    _trackFinder->hit[counter].id       = counter ;
//    _trackFinder->hit[counter].row      = this->getLayerID(trkhit);
//    _trackFinder->hit[counter].x        = point.x();
//    _trackFinder->hit[counter].y        = point.y();
//    _trackFinder->hit[counter].z        = point.z();
//    _trackFinder->hit[counter].dx       = trkhit->getCovMatrix()[0] ;
//    _trackFinder->hit[counter].dy       = trkhit->getCovMatrix()[2] ;
//    _trackFinder->hit[counter].dz       = trkhit->getCovMatrix()[5] ;
//    
//    hitmap[ _trackFinder->hit[counter].id ] = trkhit;
//    
//    ++counter;
//    
//  }    
  
  
//  nhits = _input_tpc_hits_col ? _input_tpc_hits_col->getNumberOfElements() : 0; 
  
//  // loop over the input tacks and refit using KalTest    
//  for(int i=0; i< nhits ; ++i) {
//    
//    TrackerHit* trkhit = dynamic_cast<TrackerHit*>( _input_tpc_hits_col->getElementAt( i ) ) ;
//    
//    const gear::Vector3D point(trkhit->getPosition()); 
//    
//    const double phi = this->angular_range_2PI(point.phi());
//    
//    if ( phi < _trackFinder->para.phiMin ) continue ;
//    if ( phi > _trackFinder->para.phiMax ) continue ;
//    
//    const double eta   = -std::log(tan(point.theta()/2.));
//    
//    if ( eta < _trackFinder->para.etaMin ) continue ;
//    if ( eta > _trackFinder->para.etaMax ) continue ;
//    
//    _trackFinder->hit[counter].id       = counter ;
//    _trackFinder->hit[counter].row      = this->getLayerID(trkhit)+_nlayers_vxd+_nlayers_sit;
//    _trackFinder->hit[counter].x        = point.x();
//    _trackFinder->hit[counter].y        = point.y();
//    _trackFinder->hit[counter].z        = point.z();
//    _trackFinder->hit[counter].dx       = trkhit->getCovMatrix()[0] ;
//    _trackFinder->hit[counter].dy       = trkhit->getCovMatrix()[2] ;
//    _trackFinder->hit[counter].dz       = trkhit->getCovMatrix()[5] ;
//    
//    hitmap[ _trackFinder->hit[counter].id ] = trkhit;
//    
//    ++counter;
//    
//  }    
//  
//  nhits = _input_set_hits_col ? _input_set_hits_col->getNumberOfElements() : 0; 
//  
//  // loop over the input tacks and refit using KalTest    
//  for(int i=0; i< nhits ; ++i) {
//    
//    TrackerHit* trkhit = dynamic_cast<TrackerHit*>( _input_set_hits_col->getElementAt( i ) ) ;
//    
//    const gear::Vector3D point(trkhit->getPosition()); 
//    
//    const double phi = this->angular_range_2PI(point.phi());
//    
//    if ( phi < _trackFinder->para.phiMin ) continue ;
//    if ( phi > _trackFinder->para.phiMax ) continue ;
//    
//    const double eta   = -std::log(tan(point.theta()/2.));
//    
//    if ( eta < _trackFinder->para.etaMin ) continue ;
//    if ( eta > _trackFinder->para.etaMax ) continue ;
//    
//    _trackFinder->hit[counter].id       = counter ;
//    _trackFinder->hit[counter].row      = this->getLayerID(trkhit)+_nlayers_vxd+_nlayers_sit+_nlayers_tpc;
//    _trackFinder->hit[counter].x        = point.x();
//    _trackFinder->hit[counter].y        = point.y();
//    _trackFinder->hit[counter].z        = point.z();
//    _trackFinder->hit[counter].dx       = 0.01 ;
//    _trackFinder->hit[counter].dy       = 0.01 ;
//    _trackFinder->hit[counter].dz       = 0.01 ;
//    
//    hitmap[ _trackFinder->hit[counter].id ] = trkhit;
//    
//    ++counter;
//    
//  }    
  
  
  _trackFinder->nHits = counter ;
  
  for ( int h = 0 ; h < _trackFinder->nHits ; h++ ) {
    _trackFinder->hit[h].track = 0 ;
  }
  
  
  _trackFinder->para.eventReset = 1 ;
  _trackFinder->nTracks         = 0 ;
  
  // double sectorTime =
  _trackFinder->process ( ) ;

  // print input and output for this event
  _trackFinder->para.print();
  //
  std::cout << std::endl << std::endl;
  
  //    printf ( "hit Row    x      y     z\n" ) ;
  //    for(int i=0; i<_trackFinder->nHits; ++i) {
  //      _trackFinder->hit[i].print(1);
  //    }
  
  std::cout << std::endl;
  
  for ( int i = 0 ; i < _trackFinder->nTracks ; i++ ) {
    printf ( " %d pt %f tanl %f nHits %d \n ", i, 
            _trackFinder->track[i].pt, 
            _trackFinder->track[i].tanl,
            _trackFinder->track[i].nHits );
    
    ftf::Track trk = _trackFinder->track[i];
    
    IMPL::TrackImpl* Track = new TrackImpl;

    std::vector<TrackerHit*> hit_list;
    
    for ( trk.startLoop() ; trk.done() ; trk.nextHit()  ) { 
      hit_list.push_back(hitmap[(trk.currentHit)->id]);
    }
    
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
    
    
    
    streamlog_out( DEBUG3 ) << "Create track with " << hit_list.size() << " hits" << std::endl;
    
    // First sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to in. 
    sort( hit_list.begin(), hit_list.end(), TrackFinderFTF::compare_r() );
    
    bool fit_backwards = IMarlinTrack::backward;
    
    MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
    
    try {
      
      int error = 0;
      
        error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, hit_list, Track, fit_backwards, covMatrix, _Bz, _maxChi2PerHit);
        
      
      if( error != IMarlinTrack::success || Track->getNdf() < 0 ) {       
        streamlog_out(DEBUG3) << "TrackFinderFTF:: EVENT: << " << evt->getEventNumber() << " >> Track fit returns error code " << error << " NDF = " << Track->getNdf() <<  ". Number of hits = "<< hit_list.size() << std::endl;       
        continue ;
      }
      
#ifdef MARLINTRK_DIAGNOSTICS_ON
      if (error != IMarlinTrack::success) {        
        void * dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();
      }        
#endif
      
      
      
    } catch (...) {
      
      streamlog_out(ERROR) << "TrackFinderFTF::createTrack: EVENT: << " << evt->getEventNumber() << " >> exception caught and rethown." << std::endl;       
      
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
    
    MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, *_encoder);
    
    marlinTrk->getOutliers(outliers);
    
    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }
    
    MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, *_encoder);
    
    delete marlinTrk;
    
    trackVec->addElement(Track);
    
  }
  
  
  
  
  streamlog_out(DEBUG4) << " \n **************************************************** " ;
  streamlog_out(DEBUG4) << " \n ** Number of reconstructed tracks: " << _trackFinder->nTracks  ;
  streamlog_out(DEBUG4) << " \n **************************************************** "  ;
  streamlog_out(DEBUG4) << " \n "  ;
  
  
  evt->addCollection( trackVec , _output_track_col_name) ;

  ++_n_evt ;
  
}  






void TrackFinderFTF::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrackFinderFTF::end(){ 
  
  //
  //    Destroy objects
  //
  delete[] _trackFinder->track ;
  delete[] _trackFinder->hit   ;
  delete _trackFinder;
  
  streamlog_out(DEBUG) << "TrackFinderFTF::end()  " << name() 
  << " processed " << _n_evt << " events in " << _n_run << " runs "
  << std::endl ;
  
}

LCCollection* TrackFinderFTF::GetCollection( LCEvent * evt, std::string colName ){
  
  LCCollection* col = NULL;
  
  try{
    col = evt->getCollection( colName.c_str() ) ;
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() << " track collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;     
  }
  
  return col; 
  
}

void TrackFinderFTF::setFTFParameters ( ftf::TrackFindingParameters* para ) {
  
  para->segmentRowSearchRange = 5      ;
	para->trackRowSearchRange = 5    ;
	para->infoLevel   = 8;
	para->getErrors        =  0     ;
	//fillTracks
	//ghostFlag
	para->goBackwards      =  1     ;
	para->mergePrimaries   =  1    ;
	para->minHitsPerTrack  = 5      ;
	//modRow
	para->nHitsForSegment  = 4      ;
	para->minHitsForFit    = 4      ;
	para->nEta             = 60     ;
	para->nPhi             = 10     ;
	para->deta              = 0.08F  ;
	para->dphi              = 0.08F  ;
	para->etaMin            = -3.5;
	para->etaMax            = 3.5;
	para->phiMin            = -1.75e-008;
	para->phiMax            = 6.29;
	//phishift
	para->detaMerge        = 0.02F  ;
	//distanceMerge
	para->nEtaTrack        = 70     ;
	para->nPhiTrack        = 40     ;
	para->etaMinTrack      = -3.5F  ;
	para->etaMaxTrack      =  3.5F  ;
	para->phiMinTrack      = -1.75e-008;
	para->phiMaxTrack      = 6.29;
	para->nPrimaryPasses   = 1      ;
	para->nSecondaryPasses = 0      ;
	//vertexconstrainedfit
	//parameterlocation
	para->rowInnerMost =         0;
	para->rowOuterMost =        _nlayers_vxd+_nlayers_sit-1;
	para->rowStart =            _nlayers_vxd+_nlayers_sit-1;
	para->rowEnd =               0;
	//szfitflag
	//maxchi2primary
	para->bField  =              _Bz;
	para->hitChi2Cut        = 10000000.F  ;
	para->goodHitChi2       = 10000000.F  ;
	para->trackChi2Cut      = 50000000.F  ;
	para->goodDistance     = 5000000.F   ;	
	para->ptMinHelixFit     = 100.F  ;	
	para->maxDistanceSegment = 500000.F ;	
	para->xyErrorScale     = 1.0F   ;
	para->szErrorScale     = 1.0F   ;
	para->xVertex          = 0.F    ;
	para->yVertex          = 0.F    ;
	para->dxVertex         = 0.005F ;
	para->dyVertex         = 0.005F ;
	para->zVertex          = 0.F    ;
	//xyWeightVertex
	para->phiVertex        = 0.F    ;
	para->rVertex          = 0.F    ;
  
	//not used?
	para->dphiMerge        = 0.01F  ;
	para->phiClosed        = 0      ;
  
  
}

double TrackFinderFTF::angular_range_2PI( double phi ) const {
  
  //bring phi_point into range 0 < phi < +2PI
  while (phi < 0) {
    phi += 2.0 * M_PI;
  }
  while (phi >= 2.0*M_PI) {
    phi -= 2.0 * M_PI;
  }
  
  return phi;
  
}

