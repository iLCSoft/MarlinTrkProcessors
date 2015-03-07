#include "ExtrToTracker.h"

#include <IMPL/TrackerHitPlaneImpl.h>




#include "TFile.h"
#include "TTree.h"
#include "LinkDef.h"
#include "MarlinTrk/MarlinTrkUtils.h"



#if defined GEO1
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"
#elif defined GEO2
//GEOMETRY - FOR NIKIFOROS: PUT THE INCLUDE FOR DD4HEP HERE
#else
#error Geometry type not defined
#endif



using namespace lcio ;
using namespace marlin ;
using namespace MarlinTrk ;



// //----------------------------------------------------------------
// struct Distance3D2{
//   gear::Vector3D _pos ;
//   Distance3D2( const gear::Vector3D& pos) : _pos( pos ) {}
//   template <class T>
//   double operator()( const T* t) { 
//     gear::Vector3D p( t->getPosition() ) ;
//     return ( p - _pos ).r2() ; 

//   }
// };
// //----------------------------------------------------------------


//------------------------------------------------------------------------------------------

struct PtSort {  // sort tracks wtr to pt - largest first
  inline bool operator()( const lcio::LCObject* l, const lcio::LCObject* r) {      
    return ( std::abs( ( (const lcio::Track*) l )->getOmega() ) < std::abs( ( (const lcio::Track*) r )->getOmega() )  );  // pt ~ 1./omega  
  }
};

//------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------

struct ZSort {  // sort TPC track segments wtr to Z - smallest first
  inline bool operator()( const lcio::LCObject* l, const lcio::LCObject* r) {      
    return ( std::abs( ( (const lcio::Track*) l )->getZ0() ) < std::abs( ( (const lcio::Track*) r )->getZ0() )  );  // pt ~ 1./omega  
  }
};

//------------------------------------------------------------------------------------------


ExtrToTracker aExtrToTracker ;


ExtrToTracker::ExtrToTracker() : Processor("ExtrToTracker") {
  
  // modify processor description
  _description = "ExtrToTracker refits an input track collection (TPC or VXD), and used IMarlinTrk tools to propagate it to SIT" ;
  
  
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

  registerInputCollection( LCIO::TRACKERHITPLANE,
			   "digitisedVXDHits" , 
			   "Name of the VTXTrackerHit collection"  ,
			   _vxdColName ,
			   std::string("VTXTrackerHits") ) ;

  registerInputCollection( LCIO::TRACKERHITPLANE,
			   "digitisedSITHits" , 
			   "Name of the SITTrackerHit collection"  ,
			   _sitColName ,
			   std::string("SITTrackerHits") ) ;
  
  registerOutputCollection( LCIO::TRACK,
			    "OutputTrackCollectionName" , 
			    "Name of the output track collection"  ,
			    _output_track_col_name ,
			    std::string("RefittedTracks") ) ;

  registerOutputCollection( LCIO::TRACK,
			    "SiliconCollectionName" , 
			    "Name of the output silicon track collection"  ,
			    _siTrkColName ,
			    std::string("SiliconTracks") ) ;
  
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
  

  registerProcessorParameter("Max_Chi2_Incr",
                             "maximum allowable chi2 increment when moving from one site to another",
                             _Max_Chi2_Incr,
                             double(1000));
  
  registerProcessorParameter("PropagateToLayer",
                             "Which layer should the seed be propagated to",
                             _propToLayer,
                             int(4));


  registerProcessorParameter("TPCHitsCut",
                             "minimum acceptable no of hits for the TPC tracks",
                             _tpcHitsCut,
                             int(6));

  registerProcessorParameter("Chi2NDoFCut",
                             "maximum acceptable chi2/ndof for the TPC tracks",
                             _chi2NDoFCut,
                             float(100));

  registerProcessorParameter("DoCut",
                             "maximum acceptable D0 at IP",
                             _DoCut,
                             float(2));

  registerProcessorParameter("ZoCut",
                             "maximum acceptable Z0 at IP",
                             _ZoCut,
                             float(5));

  registerProcessorParameter("SearchSigma",
                             "times d0(Z0) acceptable from track extrapolation point",
                             _searchSigma,
                             double(3));

  registerProcessorParameter("IsSpacePoints",
                             "If we use space points rather than hits (SIT)",
                             _isSpacePoints,
                             bool(false));

  registerProcessorParameter("NHitsChi2",
                             "Maximal number of hits for which a track with n hits is better than one with n-1hits. (defaut 5)",
                             _nHitsChi2,
                             int(5));

  StringVec trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back( "VXDTrackerHitRelations" );
  trackerHitsRelInputColNamesDefault.push_back( "SITTrackerHitRelations" );

  registerInputCollections("LCRelation",
                           "TrackerHitsRelInputCollections",
                           "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits.",
                           _colNamesTrackerHitRelations,
                           trackerHitsRelInputColNamesDefault );

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle input collection"  ,
			   _mcParticleCollectionName ,
			   std::string("MCParticle") ) ;



  registerProcessorParameter( "doNtuple",
			      "flag that say to fill and write the ntuple",
			      _doNtuple,
			      bool(false)
			      );

  registerProcessorParameter( "outFileName",
			      "Name of the output root file",
			      _outFileName,
			      std::string("ExtrToTracker.root")
			      );

  registerProcessorParameter( "treeName",
			      "Name of the tree",
			      _treeName,
			      std::string("extrtree")
			      );






}


void ExtrToTracker::init() { 
  
  streamlog_out(DEBUG) << "   init called  " 
  << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  

#if defined GEO1
  // set up the geometery needed by KalTest
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
#elif defined GEO2
  //GEOMETRY - FOR NIKIFOROS: PUT DD4HEP HERE
#else
#error Geometry type not defined
#endif

  
  if( _trksystem == 0 ){
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
  
  _n_run = 0 ;
  _n_evt = 0 ;
  SITHitsFitted = 0 ;
  SITHitsNonFitted = 0 ;
  TotalSITHits = 0 ;


#if defined GEO1
  //_doNtuple = true;
  _bField = Global::GEAR->getBField().at( gear::Vector3D(0., 0., 0.) ).z();    //The B field in z direction
#elif defined GEO2
  _doNtuple = false;
  //GEOMETRY - FOR NIKIFOROS: PUT DD4HEP HERE
#else
#error Geometry type not defined
#endif


  _maxChi2PerHit = 100;



  if (_doNtuple){  

    //open file and tree
    // //TFile _out(_outFileName.c_str(),"RECREATE");
    _out = new TFile(_outFileName.c_str(),"RECREATE");
    // _out->cd();
    // //_out->ls();
    _tree = new TTree(_treeName.c_str(),_treeName.c_str());

    // init tree variables 
    this->clearEventVar();

    // set tree branches 
    int bufsize = 32000; //default buffer size 32KB
    this->setTreeBranches(bufsize);    
  
  }//end _doNtuple flag
    
}


void ExtrToTracker::processRunHeader( LCRunHeader* run) { 
  
  ++_n_run ;
} 

void ExtrToTracker::processEvent( LCEvent * evt ) { 


  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  
  streamlog_out(DEBUG1) << "   processing event: " << _n_evt 
		       << std::endl ;
  
  if (_doNtuple) clearEventVar();


  // get input collection and relations 
  LCCollection* input_track_col = this->GetCollection( evt, _input_track_col_name ) ;
  LCCollection* sitHitsCol = this->GetCollection( evt, _sitColName ) ;

  if( input_track_col != 0 ){
    
    // establish the track collection that will be created 
    LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    
    //LCCollectionVec* SiTrkCol  = new LCCollectionVec(LCIO::TRACK);

    
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0) ;
    trkFlag.setBit( LCIO::TRBIT_HITS ) ;
    trackVec->setFlag( trkFlag.getFlag()  ) ;
    
    int nTracks = input_track_col->getNumberOfElements()  ;

    streamlog_out(DEBUG4) << " ######### NO OF TRACKS $$$$$$$$$$ " << nTracks << std::endl;

    LCCollectionVec* inputTrackVec = new LCCollectionVec( LCIO::TRACK )  ; 

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    for( int n= 0 ; n < nTracks ; ++n ) {
      Track* testTrack = dynamic_cast<Track*>( input_track_col->getElementAt( n ) ) ;
      inputTrackVec->addElement( testTrack );
    }

    std::sort( inputTrackVec->begin() , inputTrackVec->end() ,  PtSort()  ) ;

    int sitHits = 0 ; 
    if ( sitHitsCol != 0  ) { sitHits = sitHitsCol->getNumberOfElements();   }
    streamlog_out(DEBUG4) << "Number of sitHits in the collection " << sitHits << std::endl ;
	

    //EVENT::TrackerHitVec usedSiHits ;
    std::vector< IMPL::TrackImpl* > trackCandidates ;

    // loop over the input tracks and refit using KalTest    
    for(int i=0; i< nTracks ; ++i) {

      
      if (_doNtuple) this->clearLayerHelperVar();
     

      int SITHitsPerTrk = 0 ;
      
      Track* track = dynamic_cast<Track*>( inputTrackVec->getElementAt( i ) ) ;
      
      MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
     
      EVENT::TrackerHitVec trkHits = track->getTrackerHits() ;
      
      // //apply a chi2/ndf track sample selection
      // float chisquare = track->getChi2() ;
      // int ndof = track->getNdf() ;
      
      // // d0, Z0 values of TPC tracks
	
      // float chi2_ndf = chisquare / ( 1.0 * ndof ) ;
      
      // //accept only tracks with > minimum hits at the TPC & chi2/ndf < maximum cut & Pt cut?
      // if (chi2_ndf < _chi2NDoFCut){

      // 	streamlog_out(DEBUG4) << "%%%% Checking track " << track->id() << " chi2/ndf " << chi2_ndf << " cut " << _chi2NDoFCut << " D0 " << track->getD0() << " Z0 " << track->getZ0() << std::endl ;
	
	
	// sort the hits in R, so here we are assuming that the track came from the IP 
 
	sort(trkHits.begin(), trkHits.end(), ExtrToTracker::compare_r() );
	
	EVENT::TrackerHitVec::iterator it = trkHits.begin();
	
	for( it = trkHits.begin() ; it != trkHits.end() ; ++it ){
	  marlin_trk->addHit(*it);
	}
	
	//int init_status = marlin_trk->initialise( IMarlinTrack::backward ) ;
	//int init_status = FitInit(trkHits, marlin_trk) ;                   // alternative way to initialise the marlin_trk object
	int init_status = FitInit2(track, marlin_trk) ;          
	
	if (init_status==0) {
	  
	  streamlog_out(DEBUG4) << "track initialised " << std::endl ;
	  
	  int fit_status = marlin_trk->fit(); 
	  
	  if ( fit_status == 0 ){
	    
	    int testFlag=0;
	    
	    streamlog_out(DEBUG3) << "###########$$$$$$$$$$##############" << std::endl;
	    
	    
	    // now we have a track constisting only of TPC hits
	    
	    // 1) extrapolate to outer most SIT layer  
	    
	    // 2) select best hit candidate 
	    
	    // 3) marlin_trk->addAndFit(hit,max_chi2_increment);
	    
	    // 3a) select best hit candidate on n+1 layer 
	      
	    // 3b) marlin_trk->addAndFit(hit,max_chi2_increment);
	    
	    // 4) loop and repeat for other layers, plus VXD
	    
	    // 5) all hits have been added, propogate to IP for instance
	    
	    //_________________________________________________________
	    
	      
	    //marlin_trk->smooth(trkHits.back());	  

#if defined GEO1
	    const gear::ZPlanarParameters& gearSIT = Global::GEAR->getSITParameters() ;
	    const gear::ZPlanarLayerLayout& layerSIT = gearSIT.getZPlanarLayerLayout(); 
	    const unsigned int nSITR = layerSIT.getNLayers() ;
#elif defined GEO2
	    //GEOMETRY - FOR NIKIFOROS: PUT DD4HEP HERE
#else
#error Geometry type not defined
#endif
	    
	    
	    double chi2 = 0 ;
	    int ndf = 0 ;
	    TrackStateImpl trkState;	
	    
	    //gear::Vector3D xing_point ; 
	    
	    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
	    
	    encoder.reset() ;  // reset to 0
	    
	    int layerID = encoder.lowWord() ;  
	    int elementID = 0 ;    
	    
	    
	    //________________________________________________________________________________________________________
	    //
	    // starting loop to SIT layers
	    //________________________________________________________________________________________________________
	      

	    //for loop to all SIT layers
	    
	    for (int iL=0;iL<nSITR;iL++){
	      
	      if ( sitHitsCol != 0 ){
		
		streamlog_out(DEBUG4) << " Do I come into the loop " << std::endl;
		
		streamlog_out(DEBUG4) << "LOOP" << iL << " begins "<< std::endl;
		

		encoder[lcio::ILDCellID0::subdet] = ILDDetID::SIT ;
		encoder[lcio::ILDCellID0::layer]  = iL ;   // in case we propagate outwards from VXD
		//encoder[lcio::ILDCellID0::layer]  = 3 - iL ;  //  in case we propagate inwards from TPC
		layerID = encoder.lowWord() ;  
	      
		streamlog_out(DEBUG4) << "-- layerID " << layerID << std::endl;
		

		///////////////////////////////////////////////////////////
   

		if ( marlin_trk->propagateToLayer( layerID, trkState, chi2, ndf, elementID, IMarlinTrack::modeClosest) == MarlinTrk::IMarlinTrack::success) {
		    
		  const FloatVec& covLCIO = trkState.getCovMatrix();
		  const float* pivot = trkState.getReferencePoint();
		  double r = sqrt( pivot[0]*pivot[0]+pivot[1]*pivot[1] ) ;
		  
		  streamlog_out( DEBUG4 ) << " kaltest track parameters: "
					  << " chi2/ndf " << chi2 / ndf  
					  << " chi2 " <<  chi2 << std::endl 
		    
					  << "\t D0 "          <<  trkState.getD0() <<  "[+/-" << sqrt( covLCIO[0] ) << "] " 
					  << "\t Phi :"        <<  trkState.getPhi()<<  "[+/-" << sqrt( covLCIO[2] ) << "] " 
					  << "\t Omega "       <<  trkState.getOmega() <<  "[+/-" << sqrt( covLCIO[5] ) << "] " 
					  << "\t Z0 "          <<  trkState.getZ0() <<  "[+/-" << sqrt( covLCIO[9] ) << "] " 
					  << "\t tan(Lambda) " <<  trkState.getTanLambda() <<  "[+/-" << sqrt( covLCIO[14]) << "] " 
		      
					  << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", "  << pivot[2] 
					  << " - r: " << r << "]" 
					  << std::endl ;
		  
		  
		  //streamlog_out(DEBUG4) << " layer " << 3 - iL << " max search distances Z : Rphi " << _searchSigma*sqrt( covLCIO[9] ) << " : " << _searchSigma*sqrt( covLCIO[0] ) << std::endl ;
		  streamlog_out(DEBUG4) << " layer " << iL << " max search distances Z : Rphi " << _searchSigma*sqrt( covLCIO[9] ) << " : " << _searchSigma*sqrt( covLCIO[0] ) << std::endl ;


		  //_______________________________________________________________________________________
		  //
		  
		  streamlog_out(DEBUG2) << " element ID " << elementID << std::endl;
		  
		  if ( elementID != 0 ){
		    
		    testFlag = 1;
		    
		    //this->findHitsOnElementID(...)
		    ////this->getSinglePointResolution(...)

		    //if (!_pixelSIT){
	
		      float dU_spres = 0.007;
		      float dV_spres = 0.05;
		      //}
		      // ---------------------------------------------------------------------------------------------
		    
		    // if (HitsInLayer.size()!=0){
		      
		      // TrackerHit *BestHit;
		      // bool BestHitFound = false;
		      
		      // int pointer = 0 ;
		      // int PossibleHits = 0 ;	
		      // double DimDist = 0 ;	  
		      
		      // streamlog_out(DEBUG4) << " calling selectbestcandidatelimited with value for possible hits = " << PossibleHits << std::endl ;

		      // SelectBestCandidateLimited(HitsInLayer, pivot, BestHit, covLCIO, r, BestHitFound, _searchSigma, pointer, PossibleHits, dU_spres, dV_spres, DimDist, usedSiHits );		      
		      // if ( BestHitFound ) {      // when selection is restricted to an area maybe we will not find an appropriate hit

		      bool isSuccessfulFit = false; 

		      int nhits=0;
		      //TrackerHit* siHitToBeAdded = getSiHit(sitHitsCol, elementID, marlin_trk, nhits);
		      TrackerHit* BestHit = getSiHit(sitHitsCol, elementID, marlin_trk, nhits);
		      if (BestHit != 0){
		      		
	  
			streamlog_out(DEBUG4) << " Best hit found: call add and fit " << std::endl ; 
						  
			double chi2_increment = 0.;

			
			//smearing on the hit is really needed? has not been done in digi? - turned off for the moment
			bool doSinglePointResolutionSmearing = false; //make it a general parameter

			if (doSinglePointResolutionSmearing){
			  TrackerHitPlaneImpl *TestHitPlane = new TrackerHitPlaneImpl ;   
			  TestHitPlane->setCellID0(BestHit->getCellID0()) ;
			  TestHitPlane->setPosition(BestHit->getPosition());
			  TestHitPlane->setdU(dU_spres);
			  TestHitPlane->setdV(dV_spres);
			  isSuccessfulFit = marlin_trk->addAndFit( TestHitPlane, chi2_increment, _Max_Chi2_Incr ) == IMarlinTrack::success ;
			  delete TestHitPlane ;
			} else {
			  isSuccessfulFit = marlin_trk->addAndFit( BestHit, chi2_increment, _Max_Chi2_Incr ) == IMarlinTrack::success ;
			}

		  
			double hitx = BestHit->getPosition()[0];
			double hity = BestHit->getPosition()[1];
			double hitz = BestHit->getPosition()[2];
			double hitr = sqrt(hitx*hitx+hity*hity);
	

			// isSuccessful2 = marlin_trk->addAndFit( TestHitPlane, chi2_increment, _Max_Chi2_Incr ) == IMarlinTrack::success;
		
			TotalSITHits++;
			
			//if ( marlin_trk->addAndFit( TestHitPlane, chi2_increment, _Max_Chi2_Incr ) == IMarlinTrack::success){
			if ( isSuccessfulFit ){
			  
			  streamlog_out(DEBUG4) << " successful fit " << std::endl ; 
			  streamlog_out(DEBUG4) << " incriment in the chi2 = " << chi2_increment << "  , max chi2 to accept the hit " << _Max_Chi2_Incr  << std::endl;
			  
			  trkHits.push_back(BestHit) ;
			  
			  //usedSiHits.push_back(BestHit) ;
			  
			  streamlog_out(DEBUG4) << " hit added " << BestHit << std::endl ;
			  
			  //DO NTUPLE HERE



			  if (_doNtuple){


			    _fitX_layer_helper.push_back(pivot[0]);
			    _fitY_layer_helper.push_back(pivot[1]);
			    _fitZ_layer_helper.push_back(pivot[2]);
			    _fitR_layer_helper.push_back( sqrt(pivot[0]*pivot[0]+pivot[1]*pivot[1]) );
			    _fitD0_layer_helper.push_back( trkState.getD0() );
			    _fitD0Err_layer_helper.push_back( sqrt(covLCIO[0]) );
			    _fitZ0_layer_helper.push_back( trkState.getZ0() );
			    _fitZ0Err_layer_helper.push_back( sqrt(covLCIO[9]) );

		
#if defined GEO1
			    gear::MeasurementSurface const* surf = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface(elementID);
			    CLHEP::Hep3Vector fit_global(pivot[0],pivot[1],pivot[2]);
			    CLHEP::Hep3Vector fit_local = surf->getCoordinateSystem()->getLocalPoint(fit_global);

#elif defined GEO2
			    //GEOMETRY - FOR NIKIFOROS: do nothing (inside residuals computation) 
#else
#error Geometry type not defined
#endif	    

			    _fitU_layer_helper.push_back(fit_local.x());
			    _fitV_layer_helper.push_back(fit_local.y());
			    _fitW_layer_helper.push_back(fit_local.z());
			    _fitID_layer_helper.push_back(elementID);
		  
			    // double hitx = TestHitPlane->getPosition()[0];
			    // double hity = TestHitPlane->getPosition()[1];
			    // double hitz = TestHitPlane->getPosition()[2];
			    // // double hitx = BestHit->getPosition()[0];
			    // // double hity = BestHit->getPosition()[1];
			    // // double hitz = BestHit->getPosition()[2];
			    // double hitr = sqrt(hitx*hitx+hity*hity);
			    streamlog_out(DEBUG2) << " hit positon - x = " << hitx
						  << " mm  - y = " << hity 
						  << " mm  - z = " << hitz
						  << " mm  - R(mm) = " << hitr
						  << std::endl;

			    _hitX_layer_helper.push_back(hitx);
			    _hitY_layer_helper.push_back(hity);
			    _hitZ_layer_helper.push_back(hitz);
			    _hitR_layer_helper.push_back(hitr);

		
#if defined GEO1
			    //gear::MeasurementSurface const* surf0 = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface(elementID);
			    CLHEP::Hep3Vector hit_global(hitx,hity,hitz);
			    //CLHEP::Hep3Vector hit_local = surf0->getCoordinateSystem()->getLocalPoint(hit_global);
			    CLHEP::Hep3Vector hit_local = surf->getCoordinateSystem()->getLocalPoint(hit_global); //use same surf than extr fit
#elif defined GEO2
			    //GEOMETRY - FOR NIKIFOROS: do nothing (inside residuals computation) 
#else
#error Geometry type not defined
#endif
	    

			    _hitU_layer_helper.push_back(hit_local.x());
			    _hitV_layer_helper.push_back(hit_local.y());
			    _hitW_layer_helper.push_back(hit_local.z());
			    _hitID_layer_helper.push_back(elementID);

			    //_hitN_layer_helper.push_back(nhits);


			  }//end _doNtuple flag




			  SITHitsPerTrk++;
			  SITHitsFitted++;
			  
			  //HitsInLayer.erase( HitsInLayer.begin() + pointer ) ;
			  
			} //end successful fit
      			else{
			  
			  SITHitsNonFitted++;
			  
			}

			// if (doSinglePointResolutionSmearing) delete TestHitPlane ;

		      }   // besthit found
		    //    HitsInLayer.clear();
		    // }  // end of loop on layer i
		    
		      streamlog_out(DEBUG4) << "LOOP" << iL << " ends "<< std::endl;
		      streamlog_out(DEBUG4) << "###########$$$$$$$$$$##############" << std::endl;
		  }  
		  //track - hit association fini
		  
		}  // successful propagation to layer
		
	      } // sitHits collection not empty
	      
	    } // loop to all SIT layers
	    
	    
	    streamlog_out(MESSAGE) << " no of hits in the track (after adding SIT hits) " << trkHits.size() << " SIT hits added " << SITHitsPerTrk  << " event " <<  _n_evt<< std::endl;
	    
	    
	    // refitted track collection creation
	      
	    if (testFlag==1){
	      
	      
	      TrackStateImpl* trkState = new TrackStateImpl() ;
	      double chi2_fin = 0. ;
	      int ndf_fin = 0 ;
	      
	      marlin_trk->getTrackState(*trkState, chi2_fin, ndf_fin);
	      //const FloatVec& covMatrix = trkState->getCovMatrix();
	      
	      IMPL::TrackImpl* lcio_trk = new IMPL::TrackImpl();
	      
	      
	      //this->longSorting(...);
	      

	      sort(trkHits.begin(), trkHits.end(), ExtrToTracker::compare_r() );



	      bool fit_direction = IMarlinTrack::backward;
	      //bool fit_forwards = IMarlinTrack::forward;
	      MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();		
	      //int error = 0;
	      
	      // try {
		
	      // 	//error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, lcio_trk, fit_forwards, covMatrix, _bField, _maxChi2PerHit);  
	      // 	//error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, lcio_trk, fit_backwards, covMatrix, _bField, _maxChi2PerHit);   

	      // int error_fit =  MarlinTrk::createFit(trkHits, marlinTrk, trkState, _bField, fit_backwards, _maxChi2PerHit);
	      // int error = finaliseLCIOTrack(marlinTrk, lcio_trk, trkHits);

	      // } catch (...) {
		
	      // 	delete lcio_trk;
	      // 	delete marlinTrk;
		
	      // 	throw ;
		
	      // }
	      

	      //Kalman filter smoothing - fit track from out to in
	      int error_fit =  MarlinTrk::createFit(trkHits, marlinTrk, trkState, _bField, fit_direction, _maxChi2PerHit);
	      if (error_fit == 0) {
		int error = finaliseLCIOTrack(marlinTrk, lcio_trk, trkHits, fit_direction );
		if (error!=0){	   
		  delete lcio_trk;
		  delete marlinTrk;
		}
	      } else {
		delete lcio_trk;
		delete marlinTrk;
	      }

	      
	      // fitting finished get hit in the fit
	      
	      std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
	      std::vector<std::pair<EVENT::TrackerHit* , double> > outliers;
	      
	      // remember the hits are ordered in the order in which they were fitted
	      // here we are fitting inwards to the first is the last and vice verse
	      
	      marlinTrk->getHitsInFit(hits_in_fit);
	      
	      if( hits_in_fit.size() < 3 ) {
		streamlog_out(DEBUG3) << "RefitProcessor: Less than 3 hits in fit: Track Discarded. Number of hits =  " << trkHits.size() << std::endl;
		delete marlinTrk ;
		delete lcio_trk;
		continue ; 
	      }
	      
	      
	      std::vector<TrackerHit*> all_hits;
	      all_hits.reserve(300);
	      
	      
	      for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
		  all_hits.push_back(hits_in_fit[ihit].first);
	      }
	      
	      UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ;
	      
	      MarlinTrk::addHitNumbersToTrack(lcio_trk, all_hits, true, cellID_encoder);
		
	      marlinTrk->getOutliers(outliers);
	      
	      for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
		all_hits.push_back(outliers[ihit].first);
	      }
	      
	      MarlinTrk::addHitNumbersToTrack(lcio_trk, all_hits, false, cellID_encoder);

	      delete marlinTrk;	      
	      
	      int nhits_in_vxd = lcio_trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ];
	      int nhits_in_ftd = lcio_trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ];
	      int nhits_in_sit = lcio_trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 2 ];
	      int nhits_in_tpc = lcio_trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ];
	      int nhits_in_set = lcio_trk->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 2 ];
	      
	      
	      streamlog_out( DEBUG4 ) << " Hit numbers for Track "<< lcio_trk->id() << ": "
				      << " vxd hits = " << nhits_in_vxd
				      << " ftd hits = " << nhits_in_ftd
				      << " sit hits = " << nhits_in_sit
				      << " tpc hits = " << nhits_in_tpc
				      << " set hits = " << nhits_in_set
				      << std::endl;
	      
	      
	      if (nhits_in_vxd > 0) lcio_trk->setTypeBit( lcio::ILDDetID::VXD ) ;
	      if (nhits_in_ftd > 0) lcio_trk->setTypeBit( lcio::ILDDetID::FTD ) ;
	      if (nhits_in_sit > 0) lcio_trk->setTypeBit( lcio::ILDDetID::SIT ) ;
	      if (nhits_in_tpc > 0) lcio_trk->setTypeBit( lcio::ILDDetID::TPC ) ;
	      if (nhits_in_set > 0) lcio_trk->setTypeBit( lcio::ILDDetID::SET ) ;
	      
	      // trackCandidates.push_back(lcio_trk) ;  // trackCandidates vector stores all the candidate tracks of the event
	      
	      trackVec->addElement( lcio_trk );


	    }// end of the creation of the refitted track collection
	  }  // good fit status
	} // good initialisation status

      // } //minimum acceptable TPC hits
      


	if (_doNtuple) {
	  
	  _hitX.push_back(_hitX_layer_helper);
	  _hitY.push_back(_hitY_layer_helper);
	  _hitZ.push_back(_hitZ_layer_helper);
	  _hitR.push_back(_hitR_layer_helper);
	  _hitU.push_back(_hitU_layer_helper);
	  _hitV.push_back(_hitV_layer_helper);
	  _hitW.push_back(_hitW_layer_helper);
	  _hitID.push_back(_hitID_layer_helper);

	  _hitN.push_back(_hitN_layer_helper);

	  _fitX.push_back(_fitX_layer_helper);
	  _fitY.push_back(_fitY_layer_helper);
	  _fitZ.push_back(_fitZ_layer_helper);
	  _fitR.push_back(_fitR_layer_helper);
	  _fitU.push_back(_fitU_layer_helper);
	  _fitV.push_back(_fitV_layer_helper);
	  _fitW.push_back(_fitW_layer_helper);
	  _fitID.push_back(_fitID_layer_helper);

	  _fitD0Err.push_back(_fitD0Err_layer_helper);
	  _fitD0.push_back(_fitD0_layer_helper);

	  _fitZ0Err.push_back(_fitZ0Err_layer_helper);
	  _fitZ0.push_back(_fitZ0_layer_helper);

	}    


	delete marlin_trk;
      
    }    // for loop to the tracks 
    
    //-------------------------------------------------------------------------------------------------------		
   
    // trackVec already filled with corresponding trackCandidates
    // for (unsigned int i=0; i < trackCandidates.size(); i++){
    //   TrackImpl* SiliconRefTrack = dynamic_cast< TrackImpl* >( trackCandidates[i] );
    //   trackVec->addElement( SiliconRefTrack );
    // }

    
    // for debugging reasons
    /*
    for (int ii=0; ii<usedSiHits.size(); ii++){
      streamlog_out(DEBUG0) << " ii " << ii << " hit added " << usedSiHits.at(ii) << std::endl;
    }
    */

    evt->addCollection( trackVec , _output_track_col_name ) ;

  }// track collection no empty  
  
  ++_n_evt ;
  
  //cout << " event " << _n_evt << std::endl ;


  if (_doNtuple){

    _fitN = _fitX.size();
    _tree->Fill();
  }

  
}



void ExtrToTracker::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ExtrToTracker::end(){ 
  
  streamlog_out(DEBUG) << "ExtrToTracker::end()  " << name() 
  << " processed " << _n_evt << " events in " << _n_run << " runs "
  << std::endl ;

  streamlog_out(DEBUG4) << " SIT hits considered for track-hit association " << TotalSITHits << " how many of them were matched and fitted successfully ? " << SITHitsFitted << " for how many the fit failed ? " << SITHitsNonFitted << std::endl ;

  streamlog_out(MESSAGE) << " _doNtuple " << _doNtuple <<  std::endl ;


  
  if (_doNtuple) {

    //TFile *out = new TFile(_outFileName.c_str(),"RECREATE");
    //_out->ls();
    //_out->cd();
    _tree->Write();
    //out->Close();
    //_out->Close();

    delete _tree;
    //delete out;
    delete _out;

  }  
  
}



void ExtrToTracker::SelectBestCandidateLimited(EVENT::TrackerHitVec &HitsInLayer, const float* &pivot, EVENT::TrackerHit* &BestHit, const FloatVec& covLCIO, double& radius, bool &BestHitFound, double &sigma, int &pointer, int &PossibleHits, float &dU, float &dV, double &DimDist, TrackerHitVec &usedSiHits)
{

  BestHitFound = false ;

  int NoOfHits = HitsInLayer.size();

  double MaxDistZ = sigma*(sqrt( covLCIO[9] + dV*dV)) ;
  double MaxDistRphi = sigma*(sqrt( covLCIO[0] + dU*dU)) ;
  //double 3Ddist = 0;

  streamlog_out(DEBUG4) << " sensor single point resolution Z : Rphi " << dV << " : " << dU << " at radius " << radius << std::endl;

  for (int i=0;i<NoOfHits;i++){

    EVENT::TrackerHit *CandidateHit = HitsInLayer.at(i);

    streamlog_out(DEBUG1) << " Checking candidate hit " <<  CandidateHit << std::endl ;

    // int pointer = 0;
    
    double distZ = 0 ;   double distRphi = 0 ;
    
    double posX = CandidateHit->getPosition()[0];
    double posY = CandidateHit->getPosition()[1];
    double posZ = CandidateHit->getPosition()[2];
    
    distZ = fabs(posZ-pivot[2]);
    distRphi = sqrt( (posX-pivot[0])*(posX-pivot[0]) + (posY-pivot[1])*(posY-pivot[1]) ) ;
    
    streamlog_out(DEBUG2) << " maximum acceptable distances for track-hit matching Z : Rphi " << MaxDistZ << " : " << MaxDistRphi << " at radius " << radius << std::endl;
    
    if ( distZ < MaxDistZ && distRphi < MaxDistRphi ){
      PossibleHits++; 
      BestHit = CandidateHit;
      BestHitFound = true ;
      streamlog_out(DEBUG4) << " found appropriate hit at distance in Z : Rphi " << distZ << " : " << distRphi << std::endl; 
      MaxDistZ = distZ;
      MaxDistRphi = distRphi;
      pointer = i ;
      DimDist = sqrt((distZ*distZ) + (distRphi*distRphi));
    }
  }

  streamlog_out(DEBUG4) << " No of candidate hits to be associated to the track = " << PossibleHits << std::endl ;

}


LCCollection* ExtrToTracker::GetCollection( LCEvent * evt, std::string colName ){
  
  LCCollection* col = NULL;
  
  // int nElements = 0;
  
  try{
    col = evt->getCollection( colName.c_str() ) ;
    //int nElements = col->getNumberOfElements()  ;
    streamlog_out( DEBUG3 ) << " --> " << colName.c_str() << " track collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG3 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;     
  }
  
  return col; 
  
}

LCRelationNavigator* ExtrToTracker::GetRelations(LCEvent * evt , std::string RelName ) {
  
  LCRelationNavigator* nav = NULL ;
  
  try{
    nav = new LCRelationNavigator(evt->getCollection( RelName.c_str() ));
    streamlog_out( DEBUG2 ) << "ExtrToTracker --> " << RelName << " track relation collection in event = " << nav << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG2 ) << "ExtrToTracker --> " << RelName.c_str() << " track relation collection absent in event" << std::endl;     
  }
  
  return nav;
  
}


int ExtrToTracker::FitInit( std::vector < TrackerHit* > trackerHits , MarlinTrk::IMarlinTrack* _marlinTrk ){

  
  if  ( trackerHits.size() < 3 ){
    streamlog_out( ERROR) << "<<<<<< FitInit: Shortage of Hits! nhits = "  << trackerHits.size() << " >>>>>>>" << std::endl;
    return IMarlinTrack::error ;
  }
  
   // initialise with space-points not strips 
   // make a helix from 3 hits to get a trackstate
   const double* x1 = trackerHits[0]->getPosition();
   const double* x2 = trackerHits[ trackerHits.size()/2 ]->getPosition();
   const double* x3 = trackerHits.back()->getPosition();

   
   double r1 = sqrt(x1[0]*x1[0]+x1[1]*x1[1]);
   double r2 = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
   double r3 = sqrt(x3[0]*x3[0]+x3[1]*x3[1]);
   streamlog_out(DEBUG4) << " Radii of hits used for initialisation: " << r1 << ", " << r2 << " and " << r3 << std::endl ;
   

   HelixTrack helixTrack( x1, x2, x3, _bField, HelixTrack::forwards );
   
   helixTrack.moveRefPoint(0.0, 0.0, 0.0);
   
   //const float referencePoint[3] = { helixTrack.getRefPointX() , helixTrack.getRefPointY() , helixTrack.getRefPointZ() };
   const float referencePoint[3] = { float(helixTrack.getRefPointX()) , float(helixTrack.getRefPointY()) , float(helixTrack.getRefPointZ()) };
   
   
   
   EVENT::FloatVec covMatrix;
   
   covMatrix.resize(15);
   
   for (unsigned icov = 0; icov<covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
   }
   
   covMatrix[0]  = ( 1.e6 ); //sigma_d0^2
   covMatrix[2]  = ( 1.e2 ); //sigma_phi0^2
   covMatrix[5]  = ( 1.e-4 ); //sigma_omega^2
   covMatrix[9]  = ( 1.e6 ); //sigma_z0^2
   covMatrix[14] = ( 1.e2 ); //sigma_tanl^2
   
   
   TrackStateImpl trackState( TrackState::AtOther, 
                              helixTrack.getD0(), 
                              helixTrack.getPhi0(), 
                              helixTrack.getOmega(), 
                              helixTrack.getZ0(), 
                              helixTrack.getTanLambda(), 
                              covMatrix, 
                              referencePoint) ;
                              
   _marlinTrk->initialise( trackState, _bField, IMarlinTrack::backward ) ;

   return IMarlinTrack::success ;   

}




int ExtrToTracker::FitInit2( Track* track, MarlinTrk::IMarlinTrack* _marlinTrk ){


  //EVENT::FloatVec covMatrix = track->getCovMatrix();
  
  
  TrackStateImpl trackState( TrackState::AtOther, 
			     track->getD0(), 
			     track->getPhi(), 
			     track->getOmega(), 
			     track->getZ0(), 
			     track->getTanLambda(), 
			     track->getCovMatrix(), 
			     track->getReferencePoint()
			     ) ;
  

  //_marlinTrk->initialise( trackState, _bField, IMarlinTrack::backward ) ;
  _marlinTrk->initialise( trackState, _bField, IMarlinTrack::forward ) ;

  
  return IMarlinTrack::success ;   
  

}












 
TrackerHit* ExtrToTracker::getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk, int& nHitsOnDetEl){

  if ( sitHitsCol != 0  ) { 
    int sitHits = sitHitsCol->getNumberOfElements();   
    std::vector<TrackerHit* > hitsOnDetEl; 
    hitsOnDetEl.clear();

    for(int i=0; i<sitHits; i++){
            
      TrackerHit* hit = dynamic_cast<TrackerHit*>( sitHitsCol->getElementAt( i ) );	

      int cellID0 = hit->getCellID0();
      UTIL::BitField64 encoder0( lcio::ILDCellID0::encoder_string ); 	    
      encoder0.reset();  // reset to 0
      encoder0.setValue(cellID0);
      int hitElID = encoder0.lowWord();  
      
      streamlog_out(DEBUG4) << "-- hit element ID: " << hitElID << ", fit element ID = " << fitElID << std::endl ;

      if (hitElID==fitElID) {
	hitsOnDetEl.push_back(hit);
      }//hit and extr fit on the same detector element

    }//end loop on siHits

    double min = 9999999.;
    double testChi2=0.;
    //size_t nHitsOnDetEl = hitsOnDetEl.size();
    nHitsOnDetEl = hitsOnDetEl.size();
    int index = -1; //index of the selected hits

    streamlog_out(DEBUG2) << "-- number of hits on the same detector element: " << nHitsOnDetEl << std::endl ;

    if (nHitsOnDetEl==0) return 0;

    for(size_t i=0; i<nHitsOnDetEl; i++){
      //if ( marlin_trk->testChi2Increment(hitsOnDetEl.at(i), testChi2) == MarlinTrk::IMarlinTrack::success ) {..} // no do not do this the testChi2 internally call the addandfir setting as max acceptable chi2 a very small number to be sure to never add the the hit to the fit (so it always fails)
      marlin_trk->testChi2Increment(hitsOnDetEl.at(i), testChi2);
      streamlog_out(DEBUG2) << "-- testChi2: " << testChi2 << std::endl ;
      if (min>testChi2) {
	min = testChi2;
	index = i;
      } //end min chi2
    }//end loop on hits on same detector element

    streamlog_out(DEBUG2) << "-- index of the selected hit: " << index << std::endl ;
    
    if (index == -1) return 0;
    else {
      TrackerHit* selectedHit = dynamic_cast<TrackerHit*>( hitsOnDetEl.at(index) ) ;
      return selectedHit;
    }

  }//end hit collection not empty

  else return 0;

}//end getSiHit




void ExtrToTracker::setTreeBranches(int bufsize){
  //event variables
  _tree->Branch("nRun",&_n_run,"nRun/I");
  _tree->Branch("nEvt",&_n_evt,"nEvt/I");
  _tree->Branch("hitVXDN",&_hitVXDN,"hitVXDN/I");
  _tree->Branch("hitN", "std::vector<std::vector<double > >",&_hitN,bufsize,0); 
  _tree->Branch("hitX", "std::vector<std::vector<double > >",&_hitX,bufsize,0); 
  _tree->Branch("hitY", "std::vector<std::vector<double > >",&_hitY,bufsize,0); 
  _tree->Branch("hitZ", "std::vector<std::vector<double > >",&_hitZ,bufsize,0); 
  _tree->Branch("hitR", "std::vector<std::vector<double > >",&_hitR,bufsize,0); 
  _tree->Branch("hitU", "std::vector<std::vector<double > >",&_hitU,bufsize,0); 
  _tree->Branch("hitV", "std::vector<std::vector<double > >",&_hitV,bufsize,0); 
  _tree->Branch("hitW", "std::vector<std::vector<double > >",&_hitW,bufsize,0); 
  _tree->Branch("hitID", "std::vector<std::vector<double > >",&_hitID,bufsize,0); 
  _tree->Branch("fitN",&_fitN,"fitN/I");
  _tree->Branch("fitX", "std::vector<std::vector<double > >",&_fitX,bufsize,0); 
  _tree->Branch("fitY", "std::vector<std::vector<double > >",&_fitY,bufsize,0); 
  _tree->Branch("fitZ", "std::vector<std::vector<double > >",&_fitZ,bufsize,0); 
  _tree->Branch("fitR", "std::vector<std::vector<double > >",&_fitR,bufsize,0); 
  _tree->Branch("fitU", "std::vector<std::vector<double > >",&_fitU,bufsize,0); 
  _tree->Branch("fitV", "std::vector<std::vector<double > >",&_fitV,bufsize,0); 
  _tree->Branch("fitW", "std::vector<std::vector<double > >",&_fitW,bufsize,0); 
  _tree->Branch("fitID", "std::vector<std::vector<double > >",&_fitID,bufsize,0); 
  _tree->Branch("fitD0Err", "std::vector<std::vector<double > >",&_fitD0Err,bufsize,0); 
  _tree->Branch("fitD0", "std::vector<std::vector<double > >",&_fitD0,bufsize,0); 
  _tree->Branch("fitZ0Err", "std::vector<std::vector<double > >",&_fitZ0Err,bufsize,0); 
  _tree->Branch("fitZ0", "std::vector<std::vector<double > >",&_fitZ0,bufsize,0); 
}



void ExtrToTracker::clearEventVar(){

  _hitVXDN = 0;

  _hitN.clear();
  _hitX.clear();
  _hitY.clear();
  _hitZ.clear();
  _hitR.clear();
  _hitU.clear();
  _hitV.clear();
  _hitW.clear();
  _hitID.clear();

  _fitN=0;
  _fitX.clear();
  _fitY.clear();
  _fitZ.clear();
  _fitR.clear();
  _fitU.clear();
  _fitV.clear();
  _fitW.clear();
  _fitID.clear();

  _fitD0Err.clear();
  _fitD0.clear();

  _fitZ0Err.clear();
  _fitZ0.clear();

}

void ExtrToTracker::clearLayerHelperVar(){

  _hitN_layer_helper.clear();
  _hitX_layer_helper.clear();
  _hitY_layer_helper.clear();
  _hitZ_layer_helper.clear();
  _hitR_layer_helper.clear();
  _hitU_layer_helper.clear();
  _hitV_layer_helper.clear();
  _hitW_layer_helper.clear();
  _hitID_layer_helper.clear();

  _fitX_layer_helper.clear();
  _fitY_layer_helper.clear();
  _fitZ_layer_helper.clear();
  _fitR_layer_helper.clear();
  _fitU_layer_helper.clear();
  _fitV_layer_helper.clear();
  _fitW_layer_helper.clear();
  _fitID_layer_helper.clear();

  _fitD0Err_layer_helper.clear();
  _fitD0_layer_helper.clear();

  _fitZ0Err_layer_helper.clear();
  _fitZ0_layer_helper.clear();

}



void ExtrToTracker::fillDummy(){

  _hitN_layer_helper.push_back(-99999.);
  _hitX_layer_helper.push_back(-99999.);
  _hitY_layer_helper.push_back(-99999.);
  _hitZ_layer_helper.push_back(-99999.);
  _hitR_layer_helper.push_back(-99999.);
  _hitU_layer_helper.push_back(-99999.);
  _hitV_layer_helper.push_back(-99999.);
  _hitW_layer_helper.push_back(-99999.);
  _hitID_layer_helper.push_back(-99999.);

  _fitX_layer_helper.push_back(-99999.);
  _fitY_layer_helper.push_back(-99999.);
  _fitZ_layer_helper.push_back(-99999.);
  _fitR_layer_helper.push_back(-99999.);
  _fitU_layer_helper.push_back(-99999.);
  _fitV_layer_helper.push_back(-99999.);
  _fitW_layer_helper.push_back(-99999.);
  _fitID_layer_helper.push_back(-99999.);

  _fitD0Err_layer_helper.push_back(-99999.);
  _fitD0_layer_helper.push_back(-99999.);

  _fitZ0Err_layer_helper.push_back(-99999.);
  _fitZ0_layer_helper.push_back(-99999.);

}





//findHitsOnElementID(...){

		    // EVENT::TrackerHitVec HitsInLayer;
		    
		    // float dU_spres = 0 ;
		    // float dV_spres = 0 ;
		 
		    // Something similar done in getSiHit(...) => commented here

		    // for (int i=0;i<sitHits;i++){
		      
		    //   TrackerHit* hit = dynamic_cast<TrackerHit*>( sitHitsCol->getElementAt( i ) ) ;
		      
		    //   streamlog_out(DEBUG2) << " type = " << hit->getType() << std::endl;
		      
		    //   const int celId = hit->getCellID0() ;
		      
		    //   streamlog_out(DEBUG1) << " hit cellid0 = " << celId << std::endl;
		      
		    //   int layerNumber = 0 ;
		    //   int ladderNumber = 0 ;
		    //   int sideTest = 0 ;
		    //   int sensorNumber = 0 ;
		      
		    //   //UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
		    //   encoder.setValue(celId) ;
		    //   layerNumber  = encoder[lcio::ILDCellID0::layer] ;
		    //   ladderNumber = encoder[lcio::ILDCellID0::module] ;
		    //   sideTest = encoder[lcio::ILDCellID0::side] ;
		    //   sensorNumber = encoder[lcio::ILDCellID0::sensor] ;
		    //   encoder.reset() ;
		      
		    //   // // Just to check the element matching between the hit (coming from the digitiser) and the track extrapolation element (coming from Mokka)
		      
		    //   // int mokkaLayerNumber = 0 ;
		    //   // int mokkaLadderNumber = 0 ;
		    //   // int mokkaSideTest = 0 ;
		    //   // int mokkaSensorNumber = 0 ;
		      
		    //   // encoder.setValue(elementID) ;
		    //   // mokkaLayerNumber  = encoder[lcio::ILDCellID0::layer] ;
		    //   // mokkaLadderNumber = encoder[lcio::ILDCellID0::module] ;
		    //   // mokkaSideTest = encoder[lcio::ILDCellID0::side] ;
		    //   // mokkaSensorNumber = encoder[lcio::ILDCellID0::sensor] ;
		    //   // encoder.reset() ;
		      
		    //   streamlog_out(DEBUG2) << " checking hit : type = " << hit->getType() << " cell ID = " << celId << " side = " << sideTest << " layer = " << layerNumber << " ladder = " << ladderNumber << " sensor = " << sensorNumber << std::endl ;
		      
		    //   streamlog_out(DEBUG2) << " the element id where the hit is found " << celId << " the element id where the track is extrapolated " << elementID << std::endl;
		      
		    //   if (celId==elementID){   // cause elementID does not give sensor-side infos...
		    // 	//if (mokkaLayerNumber==layerNumber && mokkaLadderNumber==ladderNumber){
		    // 	streamlog_out(DEBUG2) << " We found a hit at the right element with type : " << hit->getType() << " side = " << sideTest << " layer = " << layerNumber << " ladder = " << ladderNumber << " sensor = " << sensorNumber << std::endl;
		    // 	streamlog_out(DEBUG3) << " We found a hit at the right element with type : " << hit->getType() << std::endl ;
		    // 	HitsInLayer.push_back(hit);
		    //   }
		    // }   //end loop on hits


//}


//getSinglePointResolution{
		    /*
		    // ------------------ in order to obtain sensors single point resolution ---------------------
		    if (_pixelSIT){
		      for (int jj=0;jj<sitHits;jj++){
			
			TrackerHitPlane* hitplane = dynamic_cast<TrackerHitPlane*>( sitHitsCol->getElementAt( jj ) ) ;
			if (hitplane){
			  const int celId_spres = hitplane->getCellID0() ;
			  
			  if (celId_spres==elementID){
			    dU_spres = hitplane->getdU();
			    dV_spres = hitplane->getdV();
			    //streamlog_out(DEBUG4) << " U resolution = " << dU_spres << " V resolution = " << dV_spres << std::endl ;
			  }
			}
		      }
		    }
		    */
//}



//longSorting(){
//-----------------------------------------------------------------------------------------------------------------------
	      // final refit of the track

	      // std::vector< std::pair<float, EVENT::TrackerHit*> > r2_values;
	      // r2_values.reserve(trkHits.size());
	      
	      // streamlog_out(DEBUG4) << " size of hits vector thats gonna be added to the output track " << trkHits.size() << std::endl ;
	      
	      // for (TrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
	      // 	EVENT::TrackerHit* h = *it;
	      // 	float r2 = h->getPosition()[0]*h->getPosition()[0]+h->getPosition()[1]*h->getPosition()[1];
	      // 	r2_values.push_back(std::make_pair(r2, *it));
	      // }
	      
	      // sort(r2_values.begin(),r2_values.end());
	      
	      // trkHits.clear();
	      // trkHits.reserve(r2_values.size());
	      

	      // for (std::vector< std::pair<float, EVENT::TrackerHit*> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
	      // 	trkHits.push_back(it->second);
	      // }
	      
//}












