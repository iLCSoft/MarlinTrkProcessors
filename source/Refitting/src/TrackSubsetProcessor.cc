#include "TrackSubsetProcessor.h"


#include <EVENT/LCCollection.h>
#include "EVENT/TrackerHit.h"
#include "EVENT/Track.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/TrackImpl.h"

#include <UTIL/ILDConf.h>

#include <gear/BField.h>

#include "marlin/VerbosityLevels.h"
#include "marlin/Global.h"

#include "KiTrack/SubsetSimple.h"
#include "KiTrack/SubsetHopfieldNN.h"
#include "Tools/Fitter.h"
#include "Tools/KiTrackMarlinTools.h"

#include "MarlinTrk/MarlinTrkUtils.h"

using namespace lcio ;
using namespace marlin ;
using namespace KiTrack;









TrackSubsetProcessor aTrackSubsetProcessor ;


TrackSubsetProcessor::TrackSubsetProcessor() : Processor("TrackSubsetProcessor") {

  // modify processor description
  _description = "TrackSubsetProcessor takes tracks from multiple sources and outputs them (or modified versions, or a subset of them) as one track collection." ;


  // register steering parameters: name, description, class-variable, default value
  
  
  std::vector< std::string > trackInputColNamesDefault;
  trackInputColNamesDefault.push_back( "ForwardTracks" );
  trackInputColNamesDefault.push_back( "SiTracks" );
  
  registerProcessorParameter( "TrackInputCollections" , 
                              "A vector of the input track collections"  ,
                              _trackInputColNames ,
                              trackInputColNamesDefault );
  
  
  registerOutputCollection(LCIO::TRACK,
                           "TrackOutputCollection",
                           "Name of the output track collection",
                           _trackOutputColName,
                           std::string("SubsetTracks"));
  
  
  //For fitting:
  
  registerProcessorParameter("MultipleScatteringOn",
                             "Use MultipleScattering in Fit",
                             _MSOn,
                             bool(true));
  
  registerProcessorParameter("EnergyLossOn",
                             "Use Energy Loss in Fit",
                             _ElossOn,
                             bool(true));
  
  registerProcessorParameter("SmoothOn",
                             "Smooth All Measurement Sites in Fit",
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

  /////////////////////////////////////////////////////////////////////////////////
  
  
  registerProcessorParameter("Omega",
                             "The parameter omega for the HNN. Controls the influence of the quality indicator.",
                             _omega,
                             double( 0.75 ) );
  
}



void TrackSubsetProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
  
  /**********************************************************************************************/
  /*       Initialise the MarlinTrkSystem, needed by the tracks for fitting                     */
  /**********************************************************************************************/
  
  _bField = Global::GEAR->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  // set upt the geometry
  _trkSystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  if( _trkSystem == 0 ) throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
  
  
  // set the options   
  _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;       //multiple scattering
  _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;     //energy loss
  _trkSystem->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;    //smoothing
  
  // initialise the tracking system
  _trkSystem->init() ;
    

}


void TrackSubsetProcessor::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
} 



void TrackSubsetProcessor::processEvent( LCEvent * evt ) { 

   
   
  std::vector< Track* > tracks;
  
  /**********************************************************************************************/
  /*       Read in the collections                                                              */
  /**********************************************************************************************/
  
  streamlog_out( DEBUG4 ) << "Try to load " << _trackInputColNames.size() << " input track collections\n";
  
  unsigned nTrackLoaded = 0;
  
  for( unsigned i=0; i < _trackInputColNames.size(); i++ ){
    
    
    LCCollection* col = NULL;
    
    try {
      
      col = evt->getCollection( _trackInputColNames[i] ) ;
      
    }
    catch(DataNotAvailableException &e) {
      
      streamlog_out( ERROR ) << "Collection " << _trackInputColNames[i] <<  " is not available!\n";     
      continue;
      
    }
    
    
    if( col != NULL ){
      
      int nTracks = col->getNumberOfElements();
      streamlog_out( DEBUG4 ) << "Load track input collection " << _trackInputColNames[i] << " with " << nTracks << " tracks\n";
      
      
      for(int j=0; j< nTracks ; j++){
        
        Track* track = dynamic_cast<Track*>( col->getElementAt( j ) ) ;
        
        if( track == NULL ){
          
          streamlog_out( ERROR ) << "Entry " << j << " from collection " << _trackInputColNames[i] << " couldn't get casted to a Track*\n";
          
        }
        else{
          
          tracks.push_back( track );
          nTrackLoaded++;
          
        }        
        
      } 
      
    }
    else streamlog_out( ERROR ) << "track input collection " << _trackInputColNames[i] << " could not be found!!!\n";
    
  }
  
  
  streamlog_out( DEBUG4 ) << "Loaded all in all " << nTrackLoaded << " tracks, which will now get further processed\n";
  
  

  
  
  
  /**********************************************************************************************/
  /*       Make sure that all tracks are compatible: find the best subset                       */
  /**********************************************************************************************/
  streamlog_out( DEBUG4 ) << "Find the best subset of tracks using the Hopfield Neural Network \n";
  
  TrackQI trackQI( _trkSystem );
  
  streamlog_out( DEBUG3 ) << "The tracks and their qualities (and their hits ): \n";
 
  for( unsigned i=0; i < tracks.size(); i++ ){
    
    double qi = trackQI( tracks[i] );
    streamlog_out( DEBUG3 ) << tracks[i] << "\t" << qi << "( ";
    std::vector< TrackerHit* > hits = tracks[i]->getTrackerHits();
    
    std::sort( hits.begin(), hits.end(), KiTrackMarlin::compare_TrackerHit_z );
    
    for( unsigned j=0; j<hits.size(); j++ ){
      
      streamlog_out( DEBUG3 ) << hits[j] << " ";
      double x = hits[j]->getPosition()[0];
      double y = hits[j]->getPosition()[1];
      double z = hits[j]->getPosition()[2];
      
      streamlog_out(DEBUG2)<< "[" << x << "," << y << "," << z << "]";
      
    }
    
    streamlog_out( DEBUG3 ) << ")\n";
    
  }
  
  TrackCompatibility comp;
  
  
  SubsetHopfieldNN< Track* > subset;
//   SubsetSimple< Track* > subset;
  subset.add( tracks );
  subset.setOmega( _omega );
  subset.calculateBestSet( comp, trackQI );

  std::vector< Track* > accepted = subset.getAccepted();
  std::vector< Track* > rejected = subset.getRejected();
  
  streamlog_out( DEBUG3 ) << "\tThe accepted tracks: \n";
  for( unsigned i=0; i < accepted.size(); i++ ){
    
    streamlog_out( DEBUG3 ) << accepted[i] <<  "\n";
    
  }
  
  streamlog_out( DEBUG3 ) << "\tThe rejected tracks: \n";
  for( unsigned i=0; i < rejected.size(); i++ ){
    
    streamlog_out( DEBUG3 ) << rejected[i] <<  "\n";
    
  }
  
  /**********************************************************************************************/
  /*            Save the tracks to a collection (make new TrackImpls from them)                 */
  /**********************************************************************************************/
  streamlog_out( DEBUG4 ) << "Fitting and saving of the tracks\n";
  
  LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
  
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackVec->setFlag( trkFlag.getFlag()  ) ;
  
  for( unsigned i=0; i < accepted.size(); i++ ){
    
    
    TrackImpl* trackImpl = new TrackImpl();
    
    Track* track = accepted[i];
    std::vector< TrackerHit* > trackerHits = track->getTrackerHits();


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
    
    
    std::vector< std::pair<float, EVENT::TrackerHit*> > r2_values;
    r2_values.reserve(trackerHits.size());
    
    for (TrackerHitVec::iterator it=trackerHits.begin(); it!=trackerHits.end(); ++it) {
      EVENT::TrackerHit* h = *it;
      float r2 = h->getPosition()[0]*h->getPosition()[0]+h->getPosition()[1]*h->getPosition()[1];
      r2_values.push_back(std::make_pair(r2, *it));
    }
    
    sort(r2_values.begin(),r2_values.end());
    
    trackerHits.clear();
    trackerHits.reserve(r2_values.size());
    
    for (std::vector< std::pair<float, EVENT::TrackerHit*> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
      trackerHits.push_back(it->second);
    }

    bool fit_backwards = MarlinTrk::IMarlinTrack::backward;
    
    MarlinTrk::IMarlinTrack* marlinTrk = _trkSystem->createTrack();
    
    
    int error = 0;
    
    try {
      
      error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trackerHits, trackImpl, fit_backwards, covMatrix, _bField, _maxChi2PerHit);
      
    } catch (...) {
      
      //      delete Track;
      //      delete marlinTrk;
      
      throw ;
      
    }
    
    // Add hit numbers 
    
    std::vector<std::pair<EVENT::TrackerHit* , double> > hits_in_fit ;
    std::vector<std::pair<EVENT::TrackerHit* , double> > outliers ;
    std::vector<TrackerHit*> all_hits;
    all_hits.reserve(300);
    
    marlinTrk->getHitsInFit(hits_in_fit);
    
    for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }
    
    UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ;
    
    MarlinTrk::addHitNumbersToTrack(trackImpl, all_hits, true, cellID_encoder);
    
    marlinTrk->getOutliers(outliers);
    
    for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }
    
    MarlinTrk::addHitNumbersToTrack(trackImpl, all_hits, false, cellID_encoder);
    
    delete marlinTrk;
    
    
    
    if( error != MarlinTrk::IMarlinTrack::success ) {
      delete trackImpl;
      streamlog_out(DEBUG3) << "TrackSubsetProcessor:: Track fit failed with error code " << error << " track dropped. Number of hits = "<< trackerHits.size() << std::endl;
      continue ;
    }
    
    if( trackImpl->getNdf() < 0) {
      delete trackImpl;
      streamlog_out(DEBUG3) << "TrackSubsetProcessor:: Track fit returns " << trackImpl->getNdf() << " degress of freedom track dropped. Number of hits = "<< trackerHits.size() << std::endl;
      continue ;
    }
    
    trackVec->addElement(trackImpl);
    
    
    
//    try{
//
//      Fitter fitter( trackImpl , _trkSystem );
//
//      TrackStateImpl* trkStateIP = new TrackStateImpl( fitter.getTrackState( lcio::TrackState::AtIP ) ) ;
//      trackImpl->setChi2( fitter.getChi2( lcio::TrackState::AtIP ) );
//      trackImpl->setNdf( fitter.getNdf( lcio::TrackState::AtIP ) );
//      trkStateIP->setLocation( TrackState::AtIP );
//      trackImpl->addTrackState( trkStateIP );
//      
//      trackVec->addElement(trackImpl);
//      
//    }
//    catch( FitterException e ){
//      
//      delete trackImpl;
//      
//      streamlog_out( ERROR ) << "TrackImpl nr " << i  << " rejected, because fit failed: " <<  e.what() << "\n";
//      continue;
//      
//    }
    
    
  }


streamlog_out( DEBUG4 ) << "Saving " << trackVec->getNumberOfElements() << " tracks\n";

  evt->addCollection( trackVec , _trackOutputColName) ;


  _nEvt ++ ;
  
}



void TrackSubsetProcessor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrackSubsetProcessor::end(){ 

    //   std::cout << "TrackSubsetProcessor::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

}





