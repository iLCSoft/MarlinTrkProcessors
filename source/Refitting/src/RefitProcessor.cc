#include "RefitProcessor.h"

#include <algorithm>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <EVENT/Track.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include "gear/BField.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"

using namespace lcio ;
using namespace marlin ;

using namespace MarlinTrk ;


RefitProcessor aRefitProcessor ;


RefitProcessor::RefitProcessor() : Processor("RefitProcessor") {
  
  // modify processor description
  _description = "RefitProcessor refits an input track collection, producing a new collection of tracks." ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  registerInputCollection( LCIO::TRACK,
                          "InputTrackCollectionName" , 
                          "Name of the input track collection"  ,
                          _input_track_col_name ,
                          std::string("TruthTracks") ) ;
  
  registerInputCollection( LCIO::LCRELATION,
                          "InputTrackRelCollection" , 
                          "Name of the MCParticle-Track Relations collection for input tracks"  ,
                          _input_track_rel_name ,
                          std::string("TruthTracksMCP") ) ;
  
  registerOutputCollection( LCIO::TRACK,
                           "OutputTrackCollectionName" , 
                           "Name of the output track collection"  ,
                           _output_track_col_name ,
                           std::string("RefittedTracks") ) ;
  
  registerOutputCollection( LCIO::LCRELATION,
                           "OutputTrackRelCollection" , 
                           "Name of the MCParticle-Track Relations collection for output tracks"  ,
                           _output_track_rel_name ,
                           std::string("RefittedTracksMCP") ) ;
  
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
  
  
}


void RefitProcessor::init() { 
  
  streamlog_out(DEBUG) << "   init called  " 
  << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  
  _bField = Global::GEAR->getBField().at( gear::Vector3D(0., 0., 0.) ).z();    //The B field in z direction
  
  // set up the geometery needed by KalTest
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  if( _trksystem == 0 ){
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
  
  _n_run = 0 ;
  _n_evt = 0 ;
  
}

void RefitProcessor::processRunHeader( LCRunHeader* run) { 
  
  ++_n_run ;
} 

void RefitProcessor::processEvent( LCEvent * evt ) { 
  
  
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  
  streamlog_out(DEBUG4) << "   processing event: " << _n_evt 
  << std::endl ;
  
  // get input collection and relations 
  LCCollection* input_track_col = this->GetCollection( evt, _input_track_col_name ) ;
  
  LCRelationNavigator* input_track_rels = this->GetRelations( evt, _input_track_rel_name ) ;
  
  if( input_track_col != 0 ){
    
    
    // establish the track collection that will be created 
    LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
    
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    trackVec->setFlag( trkFlag.getFlag()  ) ;
    
    // establish the track relations collection that will be created 
    LCCollectionVec* trackRelVec = new LCCollectionVec( LCIO::LCRELATION )  ;
    
    int nTracks = input_track_col->getNumberOfElements()  ;
    
    streamlog_out(DEBUG4) << "Processing input collection " << _input_track_col_name << " with " << nTracks << " tracks\n";
    
    // loop over the input tacks and refit using KalTest    
    for(int i=0; i< nTracks ; ++i){
      
      
      Track* track = dynamic_cast<Track*>( input_track_col->getElementAt( i ) ) ;
      
      MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
      
      EVENT::TrackerHitVec trkHits = track->getTrackerHits() ;
      
      // sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to in. 
      sort(trkHits.begin(), trkHits.end(), RefitProcessor::compare_r() );
      
      EVENT::TrackerHitVec::iterator it = trkHits.begin();
      
      int number_of_added_hits = 0;
      int ndof_added = 0;
      TrackerHitVec added_hits;
      
      for( it = trkHits.begin() ; it != trkHits.end() ; ++it ){
        
        TrackerHit* trkHit = *it;
        bool isSuccessful = false; 
          
        if( BitSet32( trkHit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ]   ){ //it is a composite spacepoint
          
          //Split it up and add both hits to the MarlinTrk
          const LCObjectVec rawObjects = trkHit->getRawHits();                    
          
          for( unsigned k=0; k< rawObjects.size(); k++ ){
            
            TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
            
            if( marlin_trk->addHit( rawHit ) == IMarlinTrack::success ){
              
              isSuccessful = true; //if at least one hit from the spacepoint gets added
              ++ndof_added;
            }
          }
        }
        else{ // normal non composite hit
         
          if (marlin_trk->addHit( trkHit ) == 0){
            isSuccessful = true;
            ndof_added += 2;
          }
          
        }
        
        if (isSuccessful) {
          added_hits.push_back(trkHit);
          ++number_of_added_hits;
        }
        
        
        else{
          streamlog_out(DEBUG4) << "Hit " << it - trkHits.begin() << " Dropped " << std::endl;          
        }
        
      }
      
      if( ndof_added < 8 ) {
        streamlog_out(DEBUG3) << "RefitProcessor: Cannot fit less with less than 8 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << std::endl;
        delete marlin_trk ;
        continue ;
      }
      
      // initialise with space-points not strips 
      // make a helix from 3 hits to get a trackstate
      const double* x1 = added_hits[0]->getPosition();
      const double* x2 = added_hits[ added_hits.size()/2 ]->getPosition();
      const double* x3 = added_hits.back()->getPosition();
      
      HelixTrack helixTrack( x1, x2, x3, _bField, IMarlinTrack::backward );
      
      helixTrack.moveRefPoint(0.0, 0.0, 0.0);
      
      const float referencePoint[3] = { helixTrack.getRefPointX() , helixTrack.getRefPointY() , helixTrack.getRefPointZ() };
      
      
      EVENT::FloatVec covMatrix;
      
      covMatrix.resize(15);
      
      for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
        covMatrix[icov] = 0;
      }
      
      covMatrix[0]  = ( 1.e4 ); //sigma_d0^2
      covMatrix[2]  = ( 1.e4 ); //sigma_phi0^2
      covMatrix[5]  = ( 1.e4 ); //sigma_omega^2
      covMatrix[9]  = ( 1.e4 ); //sigma_z0^2
      covMatrix[14] = ( 1.e4 ); //sigma_tanl^2
      
      
      TrackStateImpl trackState( TrackState::AtOther, 
                                 helixTrack.getD0(), 
                                 helixTrack.getPhi0(), 
                                 helixTrack.getOmega(), 
                                 helixTrack.getZ0(), 
                                 helixTrack.getTanLambda(), 
                                 covMatrix, 
                                 referencePoint) ;
                                 
      marlin_trk->initialise( trackState, _bField, IMarlinTrack::backward ) ;
//       marlin_trk->initialise( IMarlinTrack::backward ) ;
      
      int fit_status = marlin_trk->fit() ; 
      
      streamlog_out(DEBUG4) << "fit_status = " << fit_status << std::endl ;
      
      if( fit_status != 0 ){ 
        delete marlin_trk ;
        continue;
      }
      
      const gear::Vector3D point(0.,0.,0.); // nominal IP
      int return_code = 0;
      
      double chi2;
      int ndf;
      
      TrackImpl * refittedTrack = new TrackImpl();
      
      TrackStateImpl* trkStateIP = new TrackStateImpl;
      return_code = marlin_trk->propagate(point, *trkStateIP, chi2, ndf ) ;
      
      if (return_code !=MarlinTrk::IMarlinTrack::success ) {
        streamlog_out( ERROR ) << "  >>>>>>>>>>> could not get TrackState at IP: Track Discarded" << std::endl ;
        delete marlin_trk ;
        delete trkStateIP;
        delete refittedTrack;
        continue;
      }
      
      // fitting finished get hit in the fit
      
      std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
      
      // remember the hits are ordered in the order in which they were fitted
      // here we are fitting inwards to the first is the last and vice verse
      
      marlin_trk->getHitsInFit(hits_in_fit);
      
      if( hits_in_fit.size() < 3 ) {
        streamlog_out(DEBUG3) << "RefitProcessor: Less than 3 hits in fit: Track Discarded. Number of hits =  " << trkHits.size() << std::endl;
        delete marlin_trk ;
        delete trkStateIP;
        delete refittedTrack;
        continue ; 
      }
      
      
      EVENT::TrackerHit* last_hit_in_fit = hits_in_fit.front().first;
      if (!last_hit_in_fit) {
        throw EVENT::Exception( std::string("RefitProcessor: TrackerHit pointer == NULL ")  ) ;
      }
      
      EVENT::TrackerHit* first_hit_in_fit = hits_in_fit.back().first;
      if (!first_hit_in_fit) {
        throw EVENT::Exception( std::string("RefitProcessor: TrackerHit pointer == NULL ")  ) ;
      }
      
      
      
      TrackStateImpl* trkStateFirstHit = new TrackStateImpl;
      return_code = marlin_trk->getTrackState(first_hit_in_fit, *trkStateFirstHit, chi2, ndf ) ;
      
      if(return_code !=MarlinTrk::IMarlinTrack::success){
        streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> could not get TrackState at First Hit " << std::endl ;
        //        delete marlin_trk ;
        //        delete trkStateFirstHit;
        //        delete refittedTrack;
        //        continue;
      }
      
      TrackStateImpl* trkStateLastHit = new TrackStateImpl;
      return_code = marlin_trk->getTrackState(last_hit_in_fit, *trkStateLastHit, chi2, ndf ) ;
      
      if (return_code !=MarlinTrk::IMarlinTrack::success ) {
        streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> could not get TrackState at Last Hit " << std::endl ;
        //        delete marlin_trk ;
        //        delete trkStateLastHit;
        //        delete refittedTrack;
        //        continue;
      }
      
      TrackStateImpl* trkStateCalo = new TrackStateImpl;
      
      UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
      encoder.reset() ;  // reset to 0
      
      encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::ECAL ;
      encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::barrel;
      encoder[lcio::ILDCellID0::layer]  = 0 ;
      
      int detElementID = 0;
      return_code = marlin_trk->propagateToLayer(encoder.lowWord(), last_hit_in_fit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
      
      if (return_code == MarlinTrk::IMarlinTrack::no_intersection ) { // try forward or backward
       if (trkStateLastHit->getTanLambda()>0) {
         encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::fwd;
       }
       else{
         encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::bwd;
       }
       return_code = marlin_trk->propagateToLayer(encoder.lowWord(), last_hit_in_fit, *trkStateCalo, chi2, ndf, detElementID, IMarlinTrack::modeForward ) ;
      }
      
      if (return_code !=MarlinTrk::IMarlinTrack::success ) {
        streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> could not get TrackState at Calo Face: return_code = " << return_code << std::endl ;
        //        delete marlin_trk ;
        //        delete trkStateCalo;
        //        delete refittedTrack;
        //        continue;
      }
      
      delete marlin_trk;
      
      trkStateIP->setLocation(  lcio::TrackState::AtIP ) ;
      trkStateFirstHit->setLocation(  lcio::TrackState::AtFirstHit ) ;
      trkStateLastHit->setLocation(  lcio::TrackState::AtLastHit ) ;
      trkStateCalo->setLocation(  lcio::TrackState::AtCalorimeter ) ;
      
      refittedTrack->trackStates().push_back(trkStateIP);
      refittedTrack->trackStates().push_back(trkStateFirstHit);
      refittedTrack->trackStates().push_back(trkStateLastHit);
      refittedTrack->trackStates().push_back(trkStateCalo);
      
      refittedTrack->setChi2(chi2);
      refittedTrack->setNdf(ndf);
      
      const double* pos = trkHits.front()->getPosition();
      
      double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
      refittedTrack->setRadiusOfInnermostHit(r);
      
      //add the hits to the track
      for( unsigned j=0; j < trkHits.size(); j++ ) refittedTrack->addHit( trkHits[j] );
      
      unsigned size_of_vec = track->getSubdetectorHitNumbers().size() ;
      refittedTrack->subdetectorHitNumbers().resize(size_of_vec) ;
      for ( unsigned detIndex = 0 ;  detIndex < size_of_vec ; detIndex++ ) 
      {
        refittedTrack->subdetectorHitNumbers()[detIndex] = track->getSubdetectorHitNumbers()[detIndex] ;
      }
      
      trackVec->addElement(refittedTrack);
      
      // assign the relations previously assigned to the input tracks  
      LCObjectVec objVec = input_track_rels->getRelatedToObjects( track );
      FloatVec weights   = input_track_rels->getRelatedToWeights( track ); 
      
      for( unsigned int irel=0 ; irel < objVec.size() ; ++irel ){
        
        LCRelationImpl* rel = new LCRelationImpl ;
        rel->setFrom (refittedTrack) ;
        rel->setTo ( objVec[irel] ) ;
        rel->setWeight(weights[irel]) ; 
        trackRelVec->addElement( rel );
        
      }
      
      
      
    } 
    
    evt->addCollection( trackVec , _output_track_col_name) ;
    evt->addCollection( trackRelVec , _output_track_rel_name) ;
    
  }  
  ++_n_evt ;
}



void RefitProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void RefitProcessor::end(){ 
  
  streamlog_out(DEBUG) << "RefitProcessor::end()  " << name() 
  << " processed " << _n_evt << " events in " << _n_run << " runs "
  << std::endl ;
  
}

LCCollection* RefitProcessor::GetCollection( LCEvent * evt, std::string colName ){
  
  LCCollection* col = NULL;
  
  int nElements = 0;
  
  try{
    col = evt->getCollection( colName.c_str() ) ;
    nElements = col->getNumberOfElements()  ;
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() << " track collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG4 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;     
  }
  
  return col; 
  
}

LCRelationNavigator* RefitProcessor::GetRelations(LCEvent * evt , std::string RelName ) {
  
  LCRelationNavigator* nav = NULL ;
  
  try{
    nav = new LCRelationNavigator(evt->getCollection( RelName.c_str() ));
    streamlog_out( DEBUG2 ) << "RefitProcessor --> " << RelName << " track relation collection in event = " << nav << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG2 ) << "RefitProcessor --> " << RelName.c_str() << " track relation collection absent in event" << std::endl;     
  }
  
  return nav;
  
}


