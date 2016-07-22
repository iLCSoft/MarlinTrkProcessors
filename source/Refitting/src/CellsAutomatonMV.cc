#include "CellsAutomatonMV.h"

//using namespace MarlinTrk ;

#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"


using namespace lcio ;
using namespace marlin ;

using namespace MarlinTrk ;

// Used to fedine the quality of the track output collection
const int CellsAutomatonMV::_output_track_col_quality_GOOD = 1;
const int CellsAutomatonMV::_output_track_col_quality_FAIR = 2;
const int CellsAutomatonMV::_output_track_col_quality_POOR = 3;

const double CellsAutomatonMV::TWOPI = 2*M_PI;

CellsAutomatonMV aCellsAutomatonMV ;

CellsAutomatonMV::CellsAutomatonMV() : Processor("CellsAutomatonMV"){

  registerProcessorParameter("NDivisionsInPhi",
                             "Number of divisions in Phi",
                             _nDivisionsInPhi,
                             int(80));
  
  
  registerProcessorParameter("NDivisionsInTheta",
                             "Number of divisions in Theta",
                             _nDivisionsInTheta,
                             int(80));

  FloatVec resUEx ;
  resUEx.push_back( 0.0040 ) ;
  
  registerProcessorParameter( "PixRes" ,
                              "VXD pixels resolution"  ,
                              _resU ,
                              resUEx) ;

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
                             bool(true));
  /*
  registerProcessorParameter("UseMiddleLayer",
                             "Use middle layer for the cellular automaton",
                             _middleLayer,
                             bool(true));
  */
  registerProcessorParameter("UseSIT",
                             "Use SIT",
                             _useSIT,
                             int(1));


  registerProcessorParameter( "HitsPerTrackMin",
			      "The minimum number of hits to create a track",
			      _hitsPerTrackMin,
			      int( 4 ) );
  
  registerProcessorParameter( "StepMax",
			      "The maximum step between layers",
			      _layerStepMax,
			      int( 2 ) );

  registerProcessorParameter( "LastLayerToIP",
			      "The maximum step between layers",
			      _lastLayerToIP,
			      int( 3 ) );
  
  registerProcessorParameter("HelixFitMax",
			     "The maximum chi2/Ndf that is allowed as result of a helix fit",
			     _helixFitMax,
			     double( 500 ) );

  registerProcessorParameter("kalFitMax",
			     "The maximum chi2/Ndf that is allowed as result of a Kalman fit",
			     _chi2OverNdfCut,
			     double( 500 ) );
  
  registerProcessorParameter("Chi2ProbCut",
			     "Tracks with a chi2 probability below this will get sorted out",
			     _chi2ProbCut,
			     double(0.005));
  
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
			      float(1.e2));
  
  registerProcessorParameter("MaxHitsPerSector",
			     "Maximal number of hits allowed on a sector. More will cause drop of hits in sector",
			     _maxHitsPerSector,
			     int(1000));
  
  registerProcessorParameter( "BestSubsetFinder",
			      "The method used to find the best non overlapping subset of tracks. Available are: SubsetHopfieldNN, SubsetSimple and None",
			      _bestSubsetFinder,
			      std::string( "SubsetHopfieldNN" ) );

  // Security checks to prevent combinatorial disasters
   
  registerProcessorParameter( "MaxConnectionsAutomaton",
			      "If the automaton has more connections than this it will be redone with the next set of cut off parameters",
			      _maxConnectionsAutomaton,
			      int( 100000 ) );

  registerProcessorParameter("MaxDistance",
			     "The maximum distance between two hits in adjacent layers in order to form a minivector",
			     _maxDist,
			     double(20) );

  registerProcessorParameter("MVHitsThetaDifference",
			     "The difference in polar angle (in degrees)  between two hits in adjacent layers in order to form a minivector",
			     _hitPairThDiff,
			     double(0.5) );

  registerProcessorParameter("MVHitsThetaDifferenceInner",
			     "The difference in polar angle (in degrees)  between two hits in the INNER layer in order to form a minivector",
			     _hitPairThDiffInner,
			     double(0.1) );

  registerProcessorParameter("MVHitsThetaDifferenceSIT",
			     "The difference in polar angle (in degrees)  between two hits in SIT layer in order to form a minivector",
			     _hitPairThDiffSIT,
			     double(2.0) );

  registerProcessorParameter("NHitsChi2",
                             "Maximal number of hits for which a track with n hits is better than one with n-1hits. (defaut 5)",
                             _nHitsChi2,
                             int(5));

  registerProcessorParameter( "InitialTrackState",
			      "TrackState to use for initialization of the fit: -1: refit from hits [default], 1: AtIP, 2: AtFirstHit, 3: AtLastHit, 4:AtCalo" ,
			      _initialTrackState,
			      int(-1) );

  registerProcessorParameter( "FitDirection",
			      "Fit direction: -1: backward [default], +1: forward",
			      _fitDirection,
			      int(-1) );
  
  
  // Input Collections
  
  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "VTXHitCollectionName",
                          "VTX Hit Collection Name",
                          _VTXHitCollection,
                          std::string("VTXTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT,
                          "SITHitCollectionName",
                          "SIT Hit Collection Name",
                          _SITHitCollection,
                          std::string("SITSpacePoints"));
  
  // Output Collections
  
  registerOutputCollection(LCIO::TRACK,
			   "CATrackCollection",
			   "Name of the Cellular Automaton Tracking output collection",
			   _CATrackCollection,
			   std::string("CATracks"));
  
  // The Criteria for the Cellular Automaton:
  
  std::vector< std::string > allCriteria = Criteria::getAllCriteriaNamesVec();
  
  
  registerProcessorParameter( "Criteria",
			      "A vector of the criteria that are going to be used. For every criterion a min and max needs to be set!!!",
			      _criteriaNames,
			      allCriteria);
  
  
  // Now set min and max values for all the criteria
  for( unsigned i=0; i < _criteriaNames.size(); i++ ){
    
    std::vector< float > emptyVec;
    
    std::string critMinString = _criteriaNames[i] + "_min";
    
    registerProcessorParameter( critMinString,
				"The minimum of " + _criteriaNames[i],
				_critMinima[ _criteriaNames[i] ],
				emptyVec);
    
    
    std::string critMaxString = _criteriaNames[i] + "_max";
    
    registerProcessorParameter( critMaxString,
				"The maximum of " + _criteriaNames[i],
				_critMaxima[ _criteriaNames[i] ],
				emptyVec);
    
      
  }
  
}

void CellsAutomatonMV::init() {

  this->setupGearGeom(Global::GEAR);
  
  /**********************************************************************************************/
  /*       Make a SectorSystemVXD                                                             */
  /**********************************************************************************************/
  
  // The SectorSystemVXD is the object translating the sectors of the hits into layers, modules etc. and vice versa

  MiniVectors_sectors = 0 ;
  MiniVectors_CutSelection = 0 ;

  int nDivisionsInPhi = _nDivisionsInPhi ;
  int nDivisionsInTheta = _nDivisionsInTheta ;
  if (_useSIT==1){
    _nLayers = _nLayersVTX + _nLayersSIT + 1 ;   // + 1 adding a virtual layer for the IP
  }
  else 
    _nLayers = _nLayersVTX + 1 ;  // + 1 adding a virtual layer for the IP
  _sectorSystemVXD = new SectorSystemVXD( _nLayers, nDivisionsInPhi , nDivisionsInTheta );
  
  
  _dPhi = TWOPI/_nDivisionsInPhi;
  _dTheta = 2.0/_nDivisionsInTheta;

  /**********************************************************************************************/
  /*       Initialise the MarlinTrkSystem, needed by the tracks for fitting                     */
  /**********************************************************************************************/
  
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


void CellsAutomatonMV::processRunHeader( LCRunHeader* run) { 
    
} 




void CellsAutomatonMV::processEvent( LCEvent * evt ) {

  TrackerHitVec HitsTemp; //Hits to be deleted at the end
  std::vector< IHit* > hitsTBD; //Hits to be deleted at the end
  //std::vector< IHit* > MiniVectorsTemp;

  // reset the hit map
  _map_sector_spacepoints.clear();
  _map_sector_hits.clear();

  InitialiseVTX( evt, HitsTemp );

  unsigned round = 0; // the round we are in
  std::vector < RawTrack > rawTracks;
    


  /**********************************************************************************************/
  /*                Check if no sector is overflowing with hits                                 */
  /**********************************************************************************************/
  
  
  std::map< int , EVENT::TrackerHitVec >::iterator it;
  
  for( it=_map_sector_spacepoints.begin(); it != _map_sector_spacepoints.end(); it++ ){
    
    
    int nHits = it->second.size();
    streamlog_out( DEBUG2 ) << "Number of hits in sector " << it->first << " = " << nHits << "\n";
    
    if( nHits > _maxHitsPerSector ){
      
      it->second.clear(); //delete the hits in this sector, it will be dropped
      
      streamlog_out(ERROR)  << " ### EVENT " << evt->getEventNumber() << " :: RUN " << evt->getRunNumber() << " \n ### Number of Hits in VXD Sector " << it->first << ": " << nHits << " > " << _maxHitsPerSector << " (MaxHitsPerSector)\n : This sector will be dropped from track search, and QualityCode set to \"Poor\" " << std::endl;
      
      _output_track_col_quality = _output_track_col_quality_POOR; // We had to drop hits, so the quality of the result is decreased
      
    }
    
  }


  /**********************************************************************************************/
  /*                Add the IP as a virtual hit                                                 */
  /**********************************************************************************************/
  
  IHit* virtualIPHitForward = createVirtualIPHit(_sectorSystemVXD );
  //HitsTemp.push_back( virtualIPHitForward );
  _map_sector_hits[0].push_back( virtualIPHitForward );
  streamlog_out(DEBUG4) << " sector of the IP hit = " << virtualIPHitForward->getSector() << std::endl ;
  
  //*********************************************************************************************************************  


  /**********************************************************************************************/
  /*                        Create Mini-Vectors                                                 */
  /**********************************************************************************************/


  for ( std::map< int , EVENT::TrackerHitVec >::iterator itSecHit = _map_sector_spacepoints.begin(); itSecHit != _map_sector_spacepoints.end(); itSecHit++ ){ //over all sectors
    
    int sector = itSecHit->first;
    CreateMiniVectors( sector ); // Process one VXD sector     
    
  }



  while( setCriteria( round ) ){

    streamlog_out(DEBUG4) << " DO I ENTER IN THE GAME " << std::endl ;
    
    
    round++; // count up the round we are in
    
    SegmentBuilder segBuilder( _map_sector_hits );
    
    segBuilder.addCriteria ( _crit2Vec ); // Add the criteria on when to connect two hits. The vector has been filled by the method setCriteria
    
    unsigned layerStepMax = _layerStepMax ; // how many layers to go at max
    unsigned lastLayerToIP = _lastLayerToIP ;  
    //bool MiddleLayer = _middleLayer ;
    //VXDSectorConnector secCon( _sectorSystemVXD , layerStepMax, lastLayerToIP, MiddleLayer );  
    VXDSectorConnector secCon( _sectorSystemVXD , layerStepMax, lastLayerToIP );
    
    segBuilder.addSectorConnector ( & secCon ); // Add the sector connector (so the SegmentBuilder knows what hits from different sectors it is allowed to look for connections)
    
    // And get out the Cellular Automaton with the 1-segments 
    Automaton automaton = segBuilder.get1SegAutomaton();
    
    // Check if there are not too many connections
    if( automaton.getNumberOfConnections() > unsigned( _maxConnectionsAutomaton ) ){
      
      streamlog_out(DEBUG4) << "Redo the Automaton with different parameters, because there are too many connections:\n"
			      << "\tconnections( " << automaton.getNumberOfConnections() << " ) > MaxConnectionsAutomaton( " << _maxConnectionsAutomaton << " )\n";
      continue;
      
    }


    /**********************************************************************************************/
    /*                Automaton                                                                   */
    /**********************************************************************************************/
    
    if( automaton.getNumberOfConnections() >  0){
    
      streamlog_out( DEBUG4 ) << "\t\t---Automaton---\n" ;
    
      //if( _useCED ) KiTrackMarlin::drawAutomatonSegments( automaton ); // draws the 1-segments (i.e. hits)
      
      
      /*******************************/
      /*      2-hit segments         */
      /*******************************/
      
      streamlog_out( DEBUG4 ) << "\t\t--2-hit-Segments--\n" ;
    
      //          streamlog_out(DEBUG4) << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates\n"; //should be commented out, because it takes time
      
      automaton.clearCriteria();
      automaton.addCriteria( _crit3Vec );  // Add the criteria for 3 hits (i.e. 2 2-hit segments )
      
      
      // Let the automaton lengthen its 1-hit-segments to 2-hit-segments
      automaton.lengthenSegments();
      
      
      // So now we have 2-hit-segments and are ready to perform the Cellular Automaton.
      
      // Perform the automaton
      automaton.doAutomaton();
      
      
      // Clean segments with bad states
      automaton.cleanBadStates();
      
      
      // Reset the states of all segments
      automaton.resetStates();
      
      //          streamlog_out(DEBUG4) << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates\n"; //should be commented out, because it takes time
      
      
      // Check if there are not too many connections
      if( automaton.getNumberOfConnections() > unsigned( _maxConnectionsAutomaton ) ){
	
	streamlog_out( DEBUG4 ) << "Redo the Automaton with different parameters, because there are too many connections:\n"
				<< "\tconnections( " << automaton.getNumberOfConnections() << " ) > MaxConnectionsAutomaton( " << _maxConnectionsAutomaton << " )\n";
	continue;
      
      }
      

      // get the raw tracks (raw track = just a vector of hits, the most rudimentary form of a track)
      rawTracks = automaton.getTracks( _hitsPerTrackMin );
      
      break; // if we reached this place all went well and we don't need another round --> exit the loop
    }
  }


  streamlog_out(DEBUG4) << "Automaton returned " << rawTracks.size() << " raw tracks \n";
 


  //*************************************************************************************************************
  // Track fitting similar to forward tracking
  //*************************************************************************************************************

  std::vector <ITrack*> trackCandidates;

  // for all raw tracks we got from the automaton
  for( unsigned i=0; i < rawTracks.size(); i++){

    RawTrack rawTrack = rawTracks[i];

    if( rawTrack.size() < unsigned( _hitsPerTrackMin ) ){
      
      streamlog_out(DEBUG4) << "Trackversion discarded, too few hits: only " << rawTrack.size() << " < " << _hitsPerTrackMin << "(hitsPerTrackMin)\n";
      continue;
	       
    }
    
    VXDTrack* trackCand = new VXDTrack( _trkSystem );
   
    // add the hits to the track
    for( unsigned k=0; k<rawTrack.size(); k++ ){
               
      IMiniVector* mvHit = dynamic_cast< IMiniVector* >( rawTrack[k] ); // cast to IMiniVectors, as needed for a VXDTrack
      if( mvHit != NULL ) trackCand->addHit( mvHit );
      else streamlog_out( DEBUG4 ) << "Hit " << rawTrack[k] << " could not be casted to IMiniVector\n";

    }

    // Make here a loop on track hits to see what is finally written to VXDTrack candidate


    /*-----------------------------------------------*/
    /*                Helix Fit                      */
    /*-----------------------------------------------*/
           
    streamlog_out( DEBUG2 ) << "Fitting with Helix Fit\n";
    try{
               
      VXDHelixFitter helixFitter( trackCand->getLcioTrack() );

      TrackerHitVec testVec = trackCand->getLcioTrack()->getTrackerHits() ;

      streamlog_out( DEBUG2 ) << " $$$$ fitting track with " << testVec.size() << " hits " << std::endl ;

      float chi2OverNdf = helixFitter.getChi2() / float( helixFitter.getNdf() );
      streamlog_out( DEBUG2 ) << "chi2OverNdf = " << chi2OverNdf << "\n";

      
      if( chi2OverNdf > _helixFitMax ){
                  
	streamlog_out( DEBUG2 ) << "Discarding track because of bad helix fit: chi2/ndf = " << chi2OverNdf << "\n";

	// debug
	//streamlog_out(DEBUG4) << " pre-fitting: deleting track " << trackCand << std::endl ;
	delete trackCand;
	continue;
        
      }
      else streamlog_out( DEBUG2 ) << "Keeping track because of good helix fit: chi2/ndf = " << chi2OverNdf << "\n";
      
    }
    catch( VXDHelixFitterException e ){
      
      
      streamlog_out( DEBUG2 ) << "Track rejected, because fit failed: " <<  e.what() << "\n";
      delete trackCand;
      continue;
      
    }
    
    /*-----------------------------------------------*/
    /*                Kalman Fit                      */
    /*-----------------------------------------------*/
          
    streamlog_out( DEBUG3 ) << "Fitting with Kalman Filter\n";
    try{
      
      trackCand->fit();
      
      streamlog_out( DEBUG3 ) << " Track " << trackCand 
			      << " chi2Prob = " << trackCand->getChi2Prob() 
			      << "( chi2=" << trackCand->getChi2() 
			      <<", Ndf=" << trackCand->getNdf() << " )\n";
      

      double  test = trackCand->getChi2() /  (1.0*trackCand->getNdf()) ;  // FIXME: give a less dull name to this var. YV
      float NoOfHitsTimes2 = 2.0*( trackCand->getHits().size());
      
      //if ( trackCand->getChi2Prob() >= _chi2ProbCut ){
      if ( test < _chi2OverNdfCut ){
	
	streamlog_out( DEBUG2 ) << "Track accepted (chi2prob " << trackCand->getChi2Prob() << " >= " << _chi2ProbCut << " chi2 over ndf " << test << " No of hits (x4) " << NoOfHitsTimes2 << " Overall sorting variable " << NoOfHitsTimes2/test <<   "\n";
	
      }
      else{
	
	streamlog_out( DEBUG2 ) << "Track rejected (chi2prob " << trackCand->getChi2Prob() << " < " << _chi2ProbCut << " chi2 over ndf " << test << "\n";

	// debug
	//streamlog_out(DEBUG4) << " Kalman fitting: deleting track " << trackCand << std::endl ;

	delete trackCand;
        
	continue;
        
      }
      
      
    }
    catch( FitterException e ){
      
               
      streamlog_out( DEBUG4 ) << "Track rejected, because fit failed: " <<  e.what() << "\n";
      delete trackCand;
      continue;
      
    }
    


    // Kalman fitting over
    //____________________________________________________________________________________________________________

    
    trackCandidates.push_back( trackCand );

  }
  
  
    
  // FTD like track fitting over. 
  //_________________________________________________________________________________________________
    

  // establish the track collection that will be created 
  LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
  
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackVec->setFlag( trkFlag.getFlag()  ) ;
  

  // Selection of best consistent track subsample
  //*************************************************************************************************************



 
  std::vector< ITrack* > GoodTracks;
  std::vector< ITrack* > RejectedTracks;
 
  TrackCompatibilityShare1_MV comp;

  // Various ways to define the quality of a track. They are defined in  the header file
  TrackQISpecial_MV  testQty;
  MaxHits TheMoreTheMerrier ;

  streamlog_out(DEBUG4) << " best subset finder = " << _bestSubsetFinder << " no of tracks fed to the nnets " << trackCandidates.size() << std::endl ;

  //FIX ME: implement properly HopfieldNN also
  /*
  if( _bestSubsetFinder == "SubsetHopfieldNN" ){  
   
    streamlog_out(DEBUG4) << "Use SubsetHopfieldNN for getting the best subset\n" << std::endl ;
    
    SubsetHopfieldNN< ITrack* > subset_tracks ;
    subset_tracks.add( trackCandidates );
    subset_tracks.calculateBestSet( comp, trackQI );
    GoodTracks = subset_tracks.getAccepted();
    RejectedTracks = subset_tracks.getRejected();
    
  }
  */


  if( _bestSubsetFinder == "SubsetSimple" ){
    
            
    SubsetSimple< ITrack* > subset_tracks;
    subset_tracks.add( trackCandidates );
    subset_tracks.calculateBestSet( comp, TheMoreTheMerrier );
    GoodTracks = subset_tracks.getAccepted();
    RejectedTracks = subset_tracks.getRejected();
    
  }
  
  if ( _bestSubsetFinder == "NoSelection") { // in any other case take all tracks
    
    streamlog_out( DEBUG3 ) << "Input for subset = \"" << _bestSubsetFinder << "\". All tracks are kept\n" ;
    GoodTracks = trackCandidates ;
    
  }
  
  streamlog_out(DEBUG4) <<  "######### End of Sorting, Good tracks number: " << GoodTracks.size() <<  std::endl;

 
 
  //******************************************************************************************************************
  // best consistent track subsample selection ends here


  // Finalise the tracks

  for (unsigned int i=0; i < GoodTracks.size(); i++){
    
    VXDTrack* myTrack = dynamic_cast< VXDTrack* >( GoodTracks[i] );
    
    
    
    
    if ( myTrack != NULL ){
      
      TrackerHitVec CompSpcpoints ;

      //----------------------------------------------------------------------------
      // turn SIT mini-vectors consisting of 1D hits to spacepoints
      
      std::vector< IMiniVector* > MVTrkHits =   myTrack->getMVs();

      for (int mv = 0 ; mv <  MVTrkHits.size() ; mv++ ){
	
	MiniVector *testMV = MVTrkHits[mv]->getMiniVector() ;
	TrackerHitVec spchits = testMV->getTrackerHitVec() ;
	if ( BitSet32( spchits[0]->getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {

	  TrackerHitPlane * mvhit0 = dynamic_cast<TrackerHitPlane*>(spchits[0]);
	  TrackerHitPlane * mvhit1 = dynamic_cast<TrackerHitPlane*>(spchits[1]);

	  if ( _map_1dhits_spacepoints[mvhit0->id()] == _map_1dhits_spacepoints[mvhit1->id()]){
	    CompSpcpoints.push_back(_map_1dhits_spacepoints[mvhit0->id()]);
	  }
	}
      }

      //--------------------------------------------------------------------------
      
      
      TrackImpl* trackImpl = new TrackImpl( *(myTrack->getLcioTrack()) );   // possible leak

      //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
      // Create a new vector of hits, where we keep the 2D hits coming from 2D hits mini-vectors, but the 1D hits mini-vectors are turned back to the initial spacepoints
      // The reason for doing so is to be compatible with the trackerhit extended classes, that maps the @D hits and composite spacepoints to extended hits

      TrackImpl* newTrackImpl = new TrackImpl();
      // loop on hits of track's implementation , to see what's there
      TrackerHitVec testHitVector = trackImpl->getTrackerHits();
      for( unsigned kk=0; kk<testHitVector.size(); kk++ ){
	if ( BitSet32(testHitVector[kk]->getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
	  streamlog_out(DEBUG2) << " we skip an 1D hit " <<  std::endl;
	}
	else{   // add the 2D hits
	  newTrackImpl->addHit( testHitVector[kk]  ) ;
	}
      }

      // Now add the spacepoints
      for( unsigned kk=0; kk<CompSpcpoints.size(); kk++ ){
	newTrackImpl->addHit( CompSpcpoints[kk]  ) ;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------------------

      try{
	
	finaliseTrack( newTrackImpl, trackVec );
	
	//trackVec->addElement( newTrackImpl );
	
      }
      
      catch( FitterException e ){
	  
	streamlog_out( DEBUG4 ) << "CellsAutomatonMV: track couldn't be finalized due to fitter error: " << e.what() << "\n";
	  delete newTrackImpl;
      }
    }
  }
  // Finalisation ends
  
  evt->addCollection( trackVec , _CATrackCollection) ;



  /**********************************************************************************************/
  /*                Clean up                                                                    */
  /**********************************************************************************************/
  
  // delete all the created IHits
  for ( unsigned i=0; i<HitsTemp.size(); i++ )  delete HitsTemp[i];

  if ( MiniVectorsTemp.size() > 0 ) {
    for ( unsigned i=0; i<MiniVectorsTemp.size(); i++ ){
      delete  MiniVectorsTemp[i] ;
      MiniVectorsTemp[i] = NULL ;
    }
  }

  if ( TestMiniVectorsTemp.size() > 0 ) {
    for ( unsigned i=0; i<TestMiniVectorsTemp.size(); i++ ){
      delete  TestMiniVectorsTemp[i] ;
      TestMiniVectorsTemp[i] = NULL ;
    }
  }

  // cleanup of tracks
  for (unsigned int i=0; i < GoodTracks.size(); i++){ delete GoodTracks[i]; } 
  for ( unsigned i=0; i<RejectedTracks.size(); i++){ delete RejectedTracks[i]; }
  //for ( unsigned i=0; i<trackCandidates.size(); i++){ delete trackCandidates[i]; }
}


void CellsAutomatonMV::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void CellsAutomatonMV::end(){

   for ( unsigned i=0; i< _crit2Vec.size(); i++) delete _crit2Vec[i];
   for ( unsigned i=0; i< _crit3Vec.size(); i++) delete _crit3Vec[i];
   for ( unsigned i=0; i< _crit4Vec.size(); i++) delete _crit4Vec[i];
   _crit2Vec.clear();
   _crit3Vec.clear();
   _crit4Vec.clear();
   
   delete _sectorSystemVXD;
   _sectorSystemVXD = NULL; 

   //streamlog_out(DEBUG4) << " no of candidate minivectors created with sectors apprach " << MiniVectors_sectors << " no of candidate minivectors created applying a cut selection " <<  MiniVectors_CutSelection  << std::endl ;

}


void CellsAutomatonMV::InitialiseVTX( LCEvent * evt, EVENT::TrackerHitVec HitsTemp ) {


  // Reading out VTX Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
  try {
    
    LCCollection * hitCollection = evt->getCollection(_VTXHitCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of VTX hits = " << nelem << std::endl;
    
    for (int ielem=0; ielem<nelem; ++ielem) {
      
      TrackerHitPlane * hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));

      gear::Vector3D U(1.0,hit->getU()[1],hit->getU()[0],gear::Vector3D::spherical);
      gear::Vector3D V(1.0,hit->getV()[1],hit->getV()[0],gear::Vector3D::spherical);
      gear::Vector3D Z(0.0,0.0,1.0);
      
      const float eps = 1.0e-07;
      // V must be the global z axis 
      if( fabs(1.0 - V.dot(Z)) > eps ) {
        streamlog_out(ERROR) << "SiTrkHopfieldNN: VXD Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }
      
      if( fabs(U.dot(Z)) > eps ) {
        streamlog_out(ERROR) << "SiTrkHopfieldNN: VXD Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      double pos[3];
      double radius = 0;

      for (int i=0; i<3; ++i) {
        pos[i] = hit->getPosition()[i];
        radius += pos[i]*pos[i];
      }


      radius = sqrt(radius);

      double cosTheta = pos[2]/radius;
      double Phi = atan2(pos[1],pos[0]);
      

      if (Phi < 0.) Phi = Phi + TWOPI;
      
      // get the layer number

      int celId = hit->getCellID0() ;

      UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
      encoder.setValue(celId) ;
      int layer  = encoder[lcio::ILDCellID0::layer] + 1 ;  // +1 cause IP is considered layer 0

      int iPhi = int(Phi/_dPhi);
      int iTheta = int ((cosTheta + double(1.0))/_dTheta);
      int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;   

      streamlog_out(DEBUG4) << " CA: making the VXD hit "  << hit  << " at layer " << layer << " phi sector " << iPhi << " theta sector " << iTheta << " theta angle " << acos(cosTheta)*(180.0/M_PI) << " sector code " << iCode << " total layers " << _nLayers << " no of phi sectors " << _nDivisionsInPhi << " no of theta sectors " << _nDivisionsInTheta  << std::endl ; 

      
      //Make an VXDHit01 from the TrackerHit 
      //VXDHit01* vxdHit = new VXDHit01 ( hit , _sectorSystemVXD );   // Don't need to create VXDHits, we stick to mini - vectors
      HitsTemp.push_back(hit); //so we can easily delete every created hit afterwards
      
      _map_sector_spacepoints[ iCode ].push_back( hit );         

    }

  }

  catch(DataNotAvailableException &e) {
    streamlog_out( DEBUG4 ) << " collection not found : " << _VTXHitCollection.c_str() << std::endl ;
  }
 
  //________________________________________________________________________________________________________________
  // readout of SIT collection                                                                                     |     
  //_______________________________________________________________________________________________________________|
  
  if ( _useSIT == 1 ){
    
    try {
      LCCollection *hitCollection = evt->getCollection(_SITHitCollection.c_str());
      
      int nelem = hitCollection->getNumberOfElements();
      
      streamlog_out(DEBUG4) << "Number of SIT hits = " << nelem << std::endl;
      
      TrackerHit* trkhit = 0;
      
      for (int ielem=0; ielem<nelem; ++ielem) {
	
	trkhit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
	
	// first check that we have not been given 1D hits by mistake, as they won't work here
	if ( BitSet32( trkhit->getType() )[ UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
	  
	  streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL COMPOSITE SPACEPOINTS must be use instead. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
	  exit(1);
          
	} 
	// most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
	else if ( BitSet32( trkhit->getType() )[ UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ] ) {
	  
	  const LCObjectVec rawObjects = trkhit->getRawHits();
	  
	  for(unsigned k = 0; k < rawObjects.size(); k++){
	    
	    TrackerHit* rawHit = dynamic_cast< TrackerHit* >(rawObjects[k]);
	    TrackerHit* planarhit = dynamic_cast<TrackerHit*>(rawHit);

	    _pairs_1dhits_spacepoints = std::make_pair( planarhit, trkhit ) ;
	    _map_1dhits_spacepoints[planarhit->id()] = trkhit ;

	    int celId_SIT = planarhit->getCellID0() ;
	
	    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
	    encoder.setValue(celId_SIT) ;
	    int layer  = encoder[lcio::ILDCellID0::layer] + 1 ;  // + 1 cause IP is considered layer 0
        
	    // VXD and SIT are treated as one system so SIT layers start from _nLayersVTX
	    layer = layer + _nLayersVTX;
	
	    if (layer < 0 || layer >= _nLayers) {
	      streamlog_out(ERROR) << "SiliconTracking_MarlinTrk => fatal error in SIT : layer is outside allowed range : " << layer << std::endl;
	      exit(1);
	    }

	    double pos[3];
	    double radius = 0;
	

	    for (int i=0; i<3; ++i) {
	      pos[i] = planarhit->getPosition()[i];
	      radius += pos[i]*pos[i];
	    }
	
	    radius = sqrt(radius);
	
	    double cosTheta = pos[2]/radius;
	    double Phi = atan2(pos[1],pos[0]);

	    if (Phi < 0.) Phi = Phi + TWOPI;
	
	    int iPhi = int(Phi/_dPhi);
	    int iTheta = int ((cosTheta + double(1.0))/_dTheta);
	    int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;  

	    streamlog_out(DEBUG2) << " CA: making the SIT hit " << planarhit << " at layer " << layer << " phi sector " << iPhi << " theta sector " << iTheta << " sector code " << iCode << " total layers " << _nLayers << " no of phi sectors " << _nDivisionsInPhi << std::endl ;    

	    _map_sector_spacepoints[ iCode ].push_back( planarhit );
	  }
	} 
      }
    } catch(DataNotAvailableException &e) {
      streamlog_out( DEBUG4 ) << " collection not found : " << _SITHitCollection.c_str() << std::endl ;
    }
  }

}



//________________________________________________________________________________________________________
//                                                                                                       |
// Create all possible mini-vectors of the sector                                                        |   
//                                                                                                       |
//--------------------------------------------------------------------------------------------------------

// Taking the hits on the outer side of each VXD double layer we try to connect them with hits on the inner side.
// The connections considered valid if it satisfies a condition (for example the dtheta difference of the hits or their 
// distance didvided with  the distance of the two sides of the layer. The search is confined in  a number of target sectors
// of the inner side of the VXD layer.

void CellsAutomatonMV::CreateMiniVectors( int sector ) {

  int iTheta = sector/(_nLayers*_nDivisionsInPhi) ;
  int iPhi = ((sector - (iTheta*_nLayers*_nDivisionsInPhi)) / _nLayers) ;
  int layer = sector - (iTheta*_nLayers*_nDivisionsInPhi) - (iPhi*_nLayers) ; 

  streamlog_out(DEBUG4) << " Taking sector " << sector << " of layer " << layer << " Phi sector " << iPhi << " Theta sector " << iTheta << std::endl ;

  int iPhi_Up    = iPhi + 2;
  int iPhi_Low   = iPhi - 2;
  int iTheta_Up  = iTheta + 2; 
  int iTheta_Low = iTheta - 2;
  if (iTheta_Low < 0) iTheta_Low = 0;
  if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;

  TrackerHitVec VXDHits = _map_sector_spacepoints[sector];
  
  for (TrackerHitVec::iterator iter=VXDHits.begin(); iter!=VXDHits.end(); ++iter) {

    TrackerHit *fromHit = *iter ;   // Starting hit

    streamlog_out(DEBUG2) << " hit to initiate a MV: " <<  fromHit << std::endl ;
    
    int celID = fromHit->getCellID0() ;
    int detID = 0 ;
    //int layerID = 0 ;
    
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
    encoder.reset() ;  // reset to 0
    encoder.setValue(celID) ;

    detID = encoder[lcio::ILDCellID0::subdet] ; 
    
    layer = encoder[lcio::ILDCellID0::layer] + 1 ;  // + 1 if we consider the IP hit
    
    if (detID==lcio::ILDDetID::VXD ){

      //if (layer==5 || layer==3 || layer==1) {
      if ( layer==6 || layer==4 || layer==2 ) {     // in case of considering an IP hit
	
	for (int iPhi = iPhi_Low ; iPhi < iPhi_Up ; iPhi++){


	  // construct mini-vectors applying a delta theta cut
	  //***************************************************************************
	  
	  int iTheta_Up_mod  = iTheta + 2; 
	  int iTheta_Low_mod = iTheta - 2;
	  if (iTheta_Low_mod < 0) iTheta_Low_mod = 0;
	  if (iTheta_Up_mod  >= _nDivisionsInTheta) iTheta_Up_mod = _nDivisionsInTheta-1;
	  
	  
	  for (int iTheta = iTheta_Low_mod ; iTheta < iTheta_Up_mod ; iTheta++){
	    
	    int target_sector = ( layer-1) + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta ;
	    TrackerHitVec targetHitsMod = _map_sector_spacepoints[target_sector];

	    streamlog_out(DEBUG2) << " Checking with TARGET sector " << target_sector << " of layer " << layer-1 << " Phi sector " << iPhi << " Theta sector " << iTheta << std::endl ;

	    streamlog_out(DEBUG2) << " No of hits in the sector " << targetHitsMod.size() << std::endl ;

	    for (TrackerHitVec::iterator iterMod=targetHitsMod.begin(); iterMod!=targetHitsMod.end(); ++iterMod) {
	      
	      TrackerHit *toHitMod = *iterMod ;  // Candidate hit to form a mini - vector with the starting hit

	      if (thetaAgreementImproved(toHitMod,fromHit,layer) == true){
		
		MiniVector *testMiniVectorMod = new MiniVector(fromHit,toHitMod) ;
		MiniVectorHit01* miniVectorHitMod = new MiniVectorHit01( testMiniVectorMod , _sectorSystemVXD );  
		MiniVectors_CutSelection++;
		_map_sector_hits[ sector ].push_back( miniVectorHitMod );

		streamlog_out(DEBUG2) << " making the VXD mini-vector hit " << miniVectorHitMod << " at sector " << sector << std::endl ;
		
		MiniVectorsTemp.push_back(miniVectorHitMod);
		TestMiniVectorsTemp.push_back(testMiniVectorMod);
      
	      }
	    }
	    
	    targetHitsMod.clear();
	    
	  }
	}
      }
    }  // end of VXD hits - creation of mini-vectors



    // **** SIT **** 

    // We don't use SIT for time being. FIX ME

    if ( _useSIT == 1 ){
      if (detID==lcio::ILDDetID::SIT){
	
	layer = layer + _nLayersVTX;
      
	if ( layer==10 || layer==8) {
	  
	  int iPhi_Up_Mod    = iPhi + 8;
	  int iPhi_Low_Mod   = iPhi - 8;
	  
	  for (int iPhi = iPhi_Low_Mod ; iPhi < iPhi_Up_Mod ; iPhi++){
	    
	    int iTheta_Up_mod  = iTheta + 8; 
	    int iTheta_Low_mod = iTheta - 8;
	    if (iTheta_Low_mod < 0) iTheta_Low_mod = 0;
	    if (iTheta_Up_mod  >= _nDivisionsInTheta) iTheta_Up_mod = _nDivisionsInTheta-1;
	    
	    for (int iTheta = iTheta_Low ; iTheta < iTheta_Up ; iTheta++){
	      
	      int target_sector = ( layer-1) + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta ;
	      
	      //double ThetaAngle = acos(((2*iTheta)/(_nDivisionsInTheta*1.0)) - 1);
	      
	      TrackerHitVec targetHits = _map_sector_spacepoints[target_sector];
	    
	      streamlog_out(DEBUG2) << " Checking with TARGET sector " << target_sector << " of layer " << layer-1 << " Phi sector " << iPhi << " Theta sector " << iTheta << std::endl ;
	      streamlog_out(DEBUG2) << " No of hits in the sector " << targetHits.size() << std::endl ;	      
	      
	      for (TrackerHitVec::iterator iter2=targetHits.begin(); iter2!=targetHits.end(); ++iter2) {
		
		TrackerHit *toHit = *iter2 ;
		
		//if ( Dist(fromHit,toHit) < _maxDist ){
		if (thetaAgreementImproved(toHit,fromHit,layer) == true){

		  MiniVector *testMiniVectorMod = new MiniVector(fromHit,toHit) ;
		  MiniVectorHit01* miniVectorHitMod = new MiniVectorHit01( testMiniVectorMod , _sectorSystemVXD );  
		  MiniVectors_CutSelection++;
		  _map_sector_hits[ sector ].push_back( miniVectorHitMod );
		  
		  streamlog_out(DEBUG4) << " making the SIT mini-vector hit " << miniVectorHitMod << " at sector " << sector << std::endl ;
		  
		  MiniVectorsTemp.push_back(miniVectorHitMod);
		  TestMiniVectorsTemp.push_back(testMiniVectorMod);
		  
		  
		}
	      }
	      
	      targetHits.clear();
	    }
	  }
	}
      }  
    }
    

   
  } // end of looping on VXD -  SIT trackerhits
  
  VXDHits.clear();
  
}





void CellsAutomatonMV::setupGearGeom( const gear::GearMgr* gearMgr ){
  
  _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  //-- VXD Parameters--
  _nLayersVTX = 0 ;
  const gear::VXDParameters* pVXDDetMain = 0;
  const gear::VXDLayerLayout* pVXDLayerLayout = 0;
  
  try{
    
    streamlog_out( DEBUG9 ) << " filling VXD parameters from gear::SITParameters " << std::endl ;
    
    pVXDDetMain = &Global::GEAR->getVXDParameters();
    pVXDLayerLayout = &(pVXDDetMain->getVXDLayerLayout());
    _nLayersVTX = pVXDLayerLayout->getNLayers();
  }
  catch( gear::UnknownParameterException& e){
    
    streamlog_out( DEBUG9 ) << " ### gear::VXDParameters Not Present in GEAR FILE" << std::endl ;
    
  }
  
  
  
  //-- SIT Parameters--
  _nLayersSIT = 0 ;
  const gear::ZPlanarParameters* pSITDetMain = 0;
  const gear::ZPlanarLayerLayout* pSITLayerLayout = 0;
  
  try{
    
    streamlog_out( DEBUG9 ) << " filling SIT parameters from gear::SITParameters " << std::endl ;
    
    pSITDetMain = &Global::GEAR->getSITParameters();
    pSITLayerLayout = &(pSITDetMain->getZPlanarLayerLayout());
    _nLayersSIT = pSITLayerLayout->getNLayers();
    
  }
  catch( gear::UnknownParameterException& e){
    
    streamlog_out( DEBUG9 ) << " ### gear::SITParameters Not Present in GEAR FILE" << std::endl ;
    
  }
  
}


bool CellsAutomatonMV::setCriteria( unsigned round ){
 
   // delete the old ones
   for ( unsigned i=0; i< _crit2Vec.size(); i++) delete _crit2Vec[i];
   for ( unsigned i=0; i< _crit3Vec.size(); i++) delete _crit3Vec[i];
   for ( unsigned i=0; i< _crit4Vec.size(); i++) delete _crit4Vec[i];
   _crit2Vec.clear();
   _crit3Vec.clear();
   _crit4Vec.clear();
   
   streamlog_out(DEBUG4) << " Calling setCriteria function - criteria vector size =  " << _criteriaNames.size() << std::endl ;
   
   bool newValuesGotUsed = false; // if new values are used
   
   for( unsigned i=0; i<_criteriaNames.size(); i++ ){
      
      std::string critName = _criteriaNames[i];


      float min = _critMinima[critName].back();
      float max = _critMaxima[critName].back();
      
      streamlog_out(DEBUG4) << " iterator " << i << " criterio " << critName << " min. value " << min << " max. value "  << max << std::endl ;
      
      // use the value corresponding to the round, if there are no new ones for this criterion, just do nothing (the previous value stays in place)
      if( round + 1 <= _critMinima[critName].size() ){
         
         min =  _critMinima[critName][round];
         newValuesGotUsed = true;
         
      }
      
      if( round + 1 <= _critMaxima[critName].size() ){
         
         max =  _critMaxima[critName][round];
         newValuesGotUsed = true;
         
      }
      
      ICriterion* crit = Criteria::createCriterion( critName, min , max );

      // Some debug output about the created criterion
      std::string type = crit->getType();
      
      streamlog_out(DEBUG4) <<  "Added: Criterion " << critName << " (type =  " << type 
      << " ). Min = " << min
      << ", Max = " << max
      << ", round " << round << "\n";
      
      
      // Add the new criterion to the corresponding vector
      if( type == "2Hit" ){
         
         _crit2Vec.push_back( crit );
         
      }
      else if( type == "3Hit" ){
         
         _crit3Vec.push_back( crit );
         
      }
      else if( type == "4Hit" ){
         
         _crit4Vec.push_back( crit );
         
      }
      else delete crit;
      
      
   }
   
   return newValuesGotUsed;
   
   
}


void CellsAutomatonMV::finaliseTrack( TrackImpl* trackImpl, LCCollectionVec* trackVec ){
   


  MarlinTrk::IMarlinTrack* marlinTrk = _trkSystem->createTrack();


  //EVENT::FloatVec covMatrix = trackImpl->getCovMatrix();
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


  bool fit_direction =  IMarlinTrack::backward ;

  
  EVENT::TrackerHitVec trkHits = trackImpl->getTrackerHits() ;
  
  std::vector< std::pair<float, EVENT::TrackerHit*> > r2_values;
  r2_values.reserve(trkHits.size());
  
  for (TrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
    EVENT::TrackerHit* h = *it;
    float r2 = h->getPosition()[0]*h->getPosition()[0]+h->getPosition()[1]*h->getPosition()[1];
    r2_values.push_back(std::make_pair(r2, *it));
  }
  
  sort(r2_values.begin(),r2_values.end());
  
  trkHits.clear();
  trkHits.reserve(r2_values.size());
  
  //UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ;
  
  for (std::vector< std::pair<float, EVENT::TrackerHit*> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
    
    trkHits.push_back(it->second);
    
  }

  TrackImpl* refittedTrack = new TrackImpl ; 

  try {
      
    int error = 0;
    
    if( _initialTrackState < 0 ) { // initialize the track from three hits

      // error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_direction, covMatrix, _bField, _maxChi2PerHit);
	  
      // call with empty pre_fit  -> should use default initialisation of the implementation, e.g.
      // use an internal pre fit in aidaTT
      error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_direction, covMatrix , _bField, _maxChi2PerHit);



    } else {  // use the specified track state 

  
      EVENT::TrackState* ts = const_cast<EVENT::TrackState* > ( trackImpl->getTrackState(lcio::TrackState::AtIP) )  ;  
	  
      if( !ts ){

	std::stringstream ess ; ess << "  Could not get track state from track to refit " ;
	throw EVENT::Exception( ess.str() ) ;
      } 

      IMPL::TrackStateImpl pre_fit( *ts ) ;
      pre_fit.setCovMatrix( covMatrix )  ;
	  
      error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_direction, &pre_fit , _bField, _maxChi2PerHit);
    }
        

    if( error != IMarlinTrack::success || refittedTrack->getNdf() < 0 ) {

      streamlog_out(MESSAGE) << "::createTrack: Track fit returns error code " << error << " NDF = " << refittedTrack->getNdf() 
			     <<  ". Number of hits = "<< trkHits.size() << std::endl;

      //fg: to write out also incomplete tracks comment this out 
      delete refittedTrack;
      //continue ;
    }
        
        
  } catch (...) {
        
    streamlog_out(ERROR) << "CellsAutomatonMV::finaliseTrack exception caught and rethown. Track = " 
			 << refittedTrack << std::endl;

    delete refittedTrack;
        
    throw ;
        
  }

  //std::cout << " Track " << refittedTrack << " after the final refit has "  <<   (refittedTrack->getTrackerHits()).size() << " hits " << " and track state " << refittedTrack->getTrackState(lcio::TrackState::AtIP) << std::endl ;


  // fitting finished get hit in the fit
      
  std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
  std::vector<std::pair<EVENT::TrackerHit* , double> > outliers ;
      
  // remember the hits are ordered in the order in which they were fitted
  // here we are fitting inwards to the first is the last and vice verse
      
  marlinTrk->getHitsInFit(hits_in_fit);
      
  if( hits_in_fit.size() < 3 ) {
    streamlog_out(DEBUG3) << "RefitProcessor: Less than 3 hits in fit: Track Discarded. Number of hits =  " << trkHits.size() << std::endl;
    delete marlinTrk ;
    delete refittedTrack;
    //continue ; 
  }
      
    
  std::vector<TrackerHit*> all_hits;
  all_hits.reserve(300);
      
      
  for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
    all_hits.push_back(hits_in_fit[ihit].first);
  }
      
  UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ;
      
  MarlinTrk::addHitNumbersToTrack(refittedTrack, all_hits, true, cellID_encoder);
      
  marlinTrk->getOutliers(outliers);
      
  for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
    all_hits.push_back(outliers[ihit].first);
  }
      
  MarlinTrk::addHitNumbersToTrack(refittedTrack, all_hits, false, cellID_encoder);
      
  delete marlinTrk;
      
      
  int nhits_in_vxd = refittedTrack->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 1 ];
  int nhits_in_ftd = refittedTrack->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 1 ];
  int nhits_in_sit = refittedTrack->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 1 ];
  int nhits_in_tpc = refittedTrack->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ];
  int nhits_in_set = refittedTrack->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 1 ];
      
      
  streamlog_out( DEBUG3 ) << " Hit numbers for Track "<< refittedTrack->id() << ": "
			  << " vxd hits = " << nhits_in_vxd
			  << " ftd hits = " << nhits_in_ftd
			  << " sit hits = " << nhits_in_sit
			  << " tpc hits = " << nhits_in_tpc
			  << " set hits = " << nhits_in_set
			  << std::endl;
      

  if (nhits_in_vxd > 0) refittedTrack->setTypeBit( lcio::ILDDetID::VXD ) ;
  if (nhits_in_ftd > 0) refittedTrack->setTypeBit( lcio::ILDDetID::FTD ) ;
  if (nhits_in_sit > 0) refittedTrack->setTypeBit( lcio::ILDDetID::SIT ) ;
  if (nhits_in_tpc > 0) refittedTrack->setTypeBit( lcio::ILDDetID::TPC ) ;
  if (nhits_in_set > 0) refittedTrack->setTypeBit( lcio::ILDDetID::SET ) ;
 
  trackVec->addElement( refittedTrack );   
   
  return;
   
   
}


double CellsAutomatonMV::Dist( EVENT::TrackerHit *toHit, EVENT::TrackerHit *fromHit ) {

  const double *posOuter = fromHit->getPosition();
  const double *posInner = toHit->getPosition();

  double dist = sqrt( (posOuter[0] - posInner[0])*(posOuter[0] - posInner[0]) + (posOuter[1] - posInner[1])*(posOuter[1] - posInner[1]) +  (posOuter[2] - posInner[2])*(posOuter[2] - posInner[2])  ) ;

  //streamlog_out(DEBUG4) << " distance between minivectors = " << dist << std::endl ;

  return dist;

}


bool CellsAutomatonMV::thetaAgreementImproved( EVENT::TrackerHit *toHit, EVENT::TrackerHit *fromHit, int layer ) {

  // Improved error estimation of polar angle theta, taking into account the uncertainty of the radius as well...
  // also the single point resolution is not hard-coded any more

  bool agreement = false ;

  //double resolution = 0.004 ; // just temporarily here
  double MPS_factor = 0.0 ;  // just for test, and only for vertical muons of 4GeV

  double  pos_outer[3];
  double  pos_inner[3];

  pos_outer[0] = fromHit->getPosition()[0];
  pos_outer[1] = fromHit->getPosition()[1];
  pos_outer[2] = fromHit->getPosition()[2];
  double radout = sqrt(pos_outer[0]*pos_outer[0]+pos_outer[1]*pos_outer[1]);
  //  double theta_out = (180.0 * atan(radout/pos_outer[2])) / M_PI ;
  //BEFORE BAD RANGE - at 90deg possible flip of sign -- fixed track ineff, still bad track
  double theta_out = atan2(radout,pos_outer[2]);
  theta_out = 2*M_PI-theta_out;
  //theta_out = theta_out % (2*M_PI);
  //theta_out = fmod(theta_out,2*M_PI);
  while(theta_out>=2*M_PI){
    theta_out = theta_out - 2*M_PI; 
  }
  theta_out = theta_out*180./M_PI;

  pos_inner[0] = toHit->getPosition()[0];
  pos_inner[1] = toHit->getPosition()[1];
  pos_inner[2] = toHit->getPosition()[2];
  double radinn = sqrt(pos_inner[0]*pos_inner[0]+pos_inner[1]*pos_inner[1]);
  //double theta_inn = (180.0 *atan(radinn/pos_inner[2])) / M_PI;
  //BEFORE BAD RANGE - at 90deg possible flip of sign -- fixed track ineff, still bad track
  double theta_inn = atan2(radinn,pos_inner[2]);
  theta_inn = 2*M_PI-theta_inn;
  while(theta_inn>=2*M_PI){
    theta_inn = theta_inn - 2*M_PI; 
  }
  theta_inn = theta_inn*180./M_PI;

  double  diff_theta = fabs ( theta_out - theta_inn ) ;

  streamlog_out(DEBUG4) << "Theta of outer hit " << theta_out << " Theta of inner hit " << theta_inn << " Theta diff between the adjacent hits " << diff_theta << " maximum value " << _hitPairThDiff << std::endl ;

 

  if ( layer == 6 || layer == 4 ){

    if ( diff_theta < _hitPairThDiff ){
      agreement = true ;
    }
  }

  if ( layer == 2 ){

    if ( diff_theta < _hitPairThDiffInner ){
      agreement = true ;
    }
  }

  if ( layer == 10 || layer == 8 ) {

    if ( diff_theta < _hitPairThDiffSIT ){
      agreement = true ;
    }
  }


  return agreement ;

}


