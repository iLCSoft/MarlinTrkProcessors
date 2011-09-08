#include "RefitProcessor.h"
#include <iostream>
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


#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"

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

  streamlog_out(DEBUG) << "   processing event: " << _n_evt 
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

    // loop over the input tacks and refit using KalTest    
    for(int i=0; i< nTracks ; ++i)
      {
      
	Track* track = dynamic_cast<Track*>( input_track_col->getElementAt( i ) ) ;
		
	MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

	EVENT::TrackerHitVec trkHits = track->getTrackerHits() ;	

	// sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to in. 
	sort(trkHits.begin(), trkHits.end(), RefitProcessor::compare_r() );

	EVENT::TrackerHitVec::iterator it = trkHits.begin();

	for( it = trkHits.begin() ; it != trkHits.end() ; ++it )
	  {
	    
	    marlin_trk->addHit(*it);
	    
	  }
	
	marlin_trk->initialise( IMarlinTrack::backward ) ;
	int fit_status = marlin_trk->fit() ; 

	if( fit_status == 0 ){ 

	  gear::Vector3D xing_point ; 

	  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
	
	  encoder.reset() ;  // reset to 0
	  
	  encoder[ILDCellID0::subdet] = ILDDetID::TPC ;
	  encoder[ILDCellID0::side] = 0 ;
	  encoder[ILDCellID0::layer]  = 200 ;
	  encoder[ILDCellID0::module] = 0 ;
	  encoder[ILDCellID0::sensor] = 0 ;

	  
	  int layerID = encoder.lowWord() ;  
	  int elementID = 0 ;
	  marlin_trk->intersectionWithLayer( layerID, xing_point, elementID, IMarlinTrack::modeForward ); 

	  encoder[ILDCellID0::subdet] = ILDDetID::VXD ;
	  encoder[ILDCellID0::layer]  = 0 ;
	  layerID = encoder.lowWord() ;  

	  marlin_trk->intersectionWithLayer( layerID, xing_point, elementID ); // first VXD layer	  

	  encoder[ILDCellID0::layer]  = 1 ;
	  layerID = encoder.lowWord() ;  
	  marlin_trk->intersectionWithLayer( layerID, xing_point, elementID, IMarlinTrack::modeForward ); // second VXD layer	  

	  encoder[ILDCellID0::layer]  = 2 ;
	  layerID = encoder.lowWord() ;  
	  marlin_trk->intersectionWithLayer( layerID, xing_point, elementID, IMarlinTrack::modeForward ); // third VXD layer	  

	  encoder[ILDCellID0::subdet] = ILDDetID::SIT ;
	  encoder[ILDCellID0::layer]  = 1 ;
	  layerID = encoder.lowWord() ;  
	  marlin_trk->intersectionWithLayer( layerID, xing_point, elementID, IMarlinTrack::modeBackward ); // first SIT layer	  

	  
	  streamlog_out(DEBUG4) << "elementID = " << elementID << std::endl ;

	  marlin_trk->intersectionWithDetElement( elementID, xing_point, IMarlinTrack::modeBackward ); // just as an example get the intersection again using the elementID returned 




	  double chi2 = 0 ;
	  int ndf = 0 ;
	  // get track state at SIT outer  layer  
	  TrackStateImpl trkState_at_sit;	  
	  marlin_trk->extrapolateToLayer( layerID, trkState_at_sit, chi2, ndf, elementID, IMarlinTrack::modeBackward); 

	  // get track state at the first and last measurement sites
	  TrackStateImpl trkState_at_begin;	  

	  marlin_trk->getTrackState( trkState_at_begin, chi2, ndf ) ;

	  TrackStateImpl trkState_at_end;

	  marlin_trk->getTrackState(  trkHits.back(), trkState_at_end, chi2, ndf ) ;	  

	  const gear::Vector3D point(0.,0.,0.); // nominal IP

	  // use extrapolate, i.e. do not include material during propagation of the cov matrix 
	  TrackStateImpl trkState_extrapolated;
	  int return_code = marlin_trk->extrapolate(point, trkState_extrapolated, chi2, ndf ) ;	  

	  // use propagate, i.e. include material during propagation of the cov matrix 
	  TrackStateImpl* trkState = new TrackStateImpl() ;
	  return_code = marlin_trk->propagate(point, *trkState, chi2, ndf ) ;

	  if ( return_code == 0 ) {
	    IMPL::TrackImpl* refittedTrack = new IMPL::TrackImpl();
	    
	    refittedTrack->addTrackState(trkState);
	    refittedTrack->setChi2(chi2) ;
	    refittedTrack->setNdf(ndf) ;

	    for( it = trkHits.begin() ; it != trkHits.end() ; ++it )
	      {
		
		refittedTrack->addHit(*it);
		
	      }
	    
	    
	    // assign the relations previously assigned to the input tracks  
	    LCObjectVec objVec = input_track_rels->getRelatedToObjects( track );
	    FloatVec weights   = input_track_rels->getRelatedToWeights( track ); 
	    
	    for( unsigned int irel=0 ; irel < objVec.size() ; ++irel )
	      {
		LCRelationImpl* rel = new LCRelationImpl ;
		rel->setFrom (refittedTrack) ;
		rel->setTo ( objVec[irel] ) ;
		rel->setWeight(weights[irel]) ; 
		trackRelVec->addElement( rel );
	      }
	  
	    //	//SJA:FIXME: This has to go away. The use of hardcoded number here is completely error prone ...
	    unsigned int size_of_vec = track->getSubdetectorHitNumbers().size() ;
	    refittedTrack->subdetectorHitNumbers().resize(size_of_vec) ;
	    for ( unsigned int detIndex = 0 ;  detIndex < size_of_vec ; detIndex++ ) 
	      {
		refittedTrack->subdetectorHitNumbers()[detIndex] = track->getSubdetectorHitNumbers()[detIndex] ;
	      }
	    
	    trackVec->addElement( refittedTrack );
	  }
	}
	
	delete marlin_trk;
	
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


