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

  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,    _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,  _ElossOn) ;
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

	sort(trkHits.begin(), trkHits.end(), RefitProcessor::compare_r() );

	EVENT::TrackerHitVec::iterator it = trkHits.begin();

	for( it = trkHits.begin() ; it != trkHits.end() ; ++it )
	  {
	    
	    marlin_trk->addHit(*it);
	    
	  }
	
	bool direction = false ;
	marlin_trk->initialise( direction ) ;
	int fit_status = marlin_trk->fit( direction ) ; // SJA:FIXME: false means from out to in here i.e. backwards. This would be better if had a more meaningful name perhaps fit_fwd and fit_rev

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
	  marlin_trk->intersectionWithLayer( true, layerID, xing_point); 

	  encoder[ILDCellID0::subdet] = ILDDetID::VXD ;
	  encoder[ILDCellID0::layer]  = 0 ;
	  layerID = encoder.lowWord() ;  

	  // note the last hit to be added and filtered will probably be the first VXD hit so that is where the track state will be 
	  // this means it is not clear if the track will be considered to be backwards or forwards from here
	  // instead of a bool perhaps int with values -1,0,1 would be better ... 
	  marlin_trk->intersectionWithLayer( true, layerID, xing_point); // first VXD layer	  
	  marlin_trk->intersectionWithLayer( false, layerID, xing_point); // first VXD layer	  

	  encoder[ILDCellID0::layer]  = 1 ;
	  layerID = encoder.lowWord() ;  
	  marlin_trk->intersectionWithLayer( true, layerID, xing_point); // second VXD layer	  

	  encoder[ILDCellID0::layer]  = 2 ;
	  layerID = encoder.lowWord() ;  
	  marlin_trk->intersectionWithLayer( true, layerID, xing_point); // third VXD layer	  

	  // get track state at VXD layer 3 
	  TrackStateImpl trkState_at_vxd3;	  
	  marlin_trk->extrapolateToLayer( true, layerID, trkState_at_vxd3); 

	  // get track state at the first and last measurement sites
	  TrackStateImpl trkState_at_begin;	  

	  marlin_trk->getTrackState( trkState_at_begin ) ;

	  TrackStateImpl trkState_at_end;

	  marlin_trk->getTrackState(  trkHits.back(), trkState_at_end ) ;	  

	  const gear::Vector3D point(0.,0.,0.); // nominal IP

	  // use extrapolate, i.e. do not include material during propagation of the cov matrix 
	  TrackStateImpl trkState_extrapolated;
	  int return_code = marlin_trk->extrapolate(point, trkState_extrapolated) ;	  

	  // use propagate, i.e. include material during propagation of the cov matrix 
	  TrackStateImpl* trkState = new TrackStateImpl() ;
	  return_code = marlin_trk->propagate(point, *trkState) ;

	  if ( return_code == 0 ) {
	    IMPL::TrackImpl* refittedTrack = new IMPL::TrackImpl();
	    
	    refittedTrack->addTrackState(trkState);

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
	    refittedTrack->subdetectorHitNumbers().resize(12);
	    for ( unsigned int detIndex = 0 ;  detIndex < refittedTrack->subdetectorHitNumbers().size() ; detIndex++ ) 
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


