#include "ExtrToTracker.h"

#include <IMPL/TrackerHitPlaneImpl.h>


#include "TFile.h"
#include "TTree.h"
#include "LinkDef.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include "UTIL/Operators.h"

#if defined GEO1 //ILD GEAR GEOMETRY
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"
#elif defined GEO2 //CLIC DD4HEP GEOMETRY
#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/DetectorData.h"
#include "DDRec/DDGear.h"

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


  /////////////////////////

  StringVec vecDigiHitsDefault;
  vecDigiHitsDefault.push_back( "ITrackerHits" );
  vecDigiHitsDefault.push_back( "OTrackerHits" );
  vecDigiHitsDefault.push_back( "ITrackerEndcapHits" );
  vecDigiHitsDefault.push_back( "OTrackerEndcapHits" );

  registerInputCollections(LCIO::TRACKERHITPLANE,
                           "vecDigiHits",
                           "vector of name of the digi hits collection - nned to be syncro with vecSubdetName!",
                           _vecDigiHits,
                           vecDigiHitsDefault );


  StringVec vecSubdetNameDefault;
  vecSubdetNameDefault.push_back("InnerTrackerBarrel") ;
  vecSubdetNameDefault.push_back("OuterTrackerBarrel") ;
  vecSubdetNameDefault.push_back("InnerTrackerEndcap") ;
  vecSubdetNameDefault.push_back("OuterTrackerEndcap") ;

  registerProcessorParameter( "vecSubdetName" , 
                              "vector of names of all subdetector to exrapolate to" ,
			      _vecSubdetName ,
                              vecSubdetNameDefault );

  /////////////////////////






  registerInputCollection( LCIO::TRACKERHITPLANE,
			   "digitisedSITHits" , 
			   "Name of the SITTrackerHit collection"  ,
			   _sitColName ,
			   std::string("ITrackerHits") ) ;


  registerInputCollection( LCIO::TRACKERHITPLANE,
			   "digitisedOTHits" , 
			   "Name of the SITTrackerHit collection"  ,
			   _otColName ,
			   std::string("OTrackerHits") ) ;
  

  registerOutputCollection( LCIO::TRACK,
			    "OutputTrackCollectionName" , 
			    "Name of the output track collection"  ,
			    _output_track_col_name ,
			    std::string("ExtrTracks") ) ;

  registerProcessorParameter( "detectorElementName",
			      "Name of the detector element you are extrapolating to",
			      _detElName,
			      std::string("InnerTrackerBarrel")
			      );

  registerProcessorParameter( "detectorID",
			      "ID of the detector element you are extrapolating to",
			      _detID,
			      int(3)
			      );


  registerProcessorParameter( "detectorElementOTName",
			      "Name of the detector element you are extrapolating to",
			      _detElOTName,
			      std::string("")
			      );
  //std::string("OuterTrackerBarrel")


  registerProcessorParameter( "detectorIDOT",
			      "ID of the detector element you are extrapolating to",
			      _detIDOT,
			      int(5)
			      );


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

  registerProcessorParameter("SearchSigma",
                             "times d0(Z0) acceptable from track extrapolation point",
                             _searchSigma,
                             double(3));

  registerProcessorParameter("PerformFinalRefit",
                             "perform a final refit of the extrapolated track",
                             _performFinalRefit,
                             bool(false));

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

  _bField = Global::GEAR->getBField().at( gear::Vector3D(0., 0., 0.) ).z();    //The B field in z direction

  // set up the geometery needed by KalTest
  //FIXME: for now do KalTest only - make this a steering parameter to use other fitters
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
#elif defined GEO2

  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();

  const double pos[3]={0,0,0}; 
  double bFieldVec[3]={0,0,0}; 
  lcdd.field().magneticField(pos,bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2]/dd4hep::tesla; // z component at (0,0,0)

  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , marlin::Global::GEAR , "" ) ;


  // DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();


  // //===========  get the surface map from the SurfaceManager ================                                                                                 

  // DD4hep::DDRec::SurfaceManager& surfMan = *lcdd.extension<DD4hep::DDRec::SurfaceManager>() ;

  // std::string _subDetName = _detElName;
  // DD4hep::Geometry::DetElement det = lcdd.detector( _subDetName ) ;

  // _map = surfMan.map( det.name() ) ;

  // if( ! _map ) {
  //   std::stringstream err  ; err << " Could not find surface map for detector: "
  //                                << _subDetName << " in SurfaceManager " ;
  //   throw EVENT::Exception( err.str() ) ;
  // }

  // streamlog_out( DEBUG3 ) << " DDPlanarDigiProcessor::init(): found " << _map->size()
  //                         << " surfaces for detector:" <<  _subDetName << std::endl ;

  
  // streamlog_out( DEBUG3 ) << " ExtrToTracker::init(): found " << _map->size() << " surfaces for detector with requested ids."<<std::endl ;
  
  // if( _map->empty() ) {   std::stringstream err  ; err << " Could not find any surfaces for requested detector ids in processor " 
  // 						       <<   this->name() ;
  //   throw EVENT::Exception( err.str() ) ;
  // }

  
#else
#error Geometry type not defined
#endif






  //   /////////////////////////////
  // #if defined GEO1
  //   const gear::ZPlanarParameters& gearSIT = Global::GEAR->getSITParameters() ;
  //   const gear::ZPlanarLayerLayout& layerSIT = gearSIT.getZPlanarLayerLayout(); 
  //   const unsigned int nSITR = layerSIT.getNLayers() ;
  // #elif defined GEO2

  //   DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();

  //   _vecSubdetNLayers.clear();
  //   for(size_t i=0; i<_vecSubdetName.size(); i++){
  //     streamlog_out(DEBUG2) << "ExtrToTracker - vecSubdetName.at("<<i<<") = " << _vecSubdetName.at(i) << std::endl;
  //     int nSubdetLayers = 0;
  //     if (_vecSubdetName.at(i)!="" && (_vecSubdetName.at(i).std::string::find("Barrel") != std::string::npos) ){
  //       DD4hep::Geometry::DetElement detEl = lcdd.detector(_vecSubdetName.at(i));
  //       DD4hep::DDRec::ZPlanarData* data = detEl.extension<DD4hep::DDRec::ZPlanarData>(); 
  //       nSubdetLayers=data->layers.size();
  //     } 
  //     else if (_vecSubdetName.at(i).std::string::find("OuterTrackerEndcap") != std::string::npos){
  //       nSubdetLayers=5;
  //     }
  //     else if (_vecSubdetName.at(i).std::string::find("InnerTrackerEndcap") != std::string::npos){
  //       nSubdetLayers=6;
  //     }
  //     _vecSubdetNLayers.push_back(nSubdetLayers);
  //     streamlog_out(DEBUG2) << "ExtrToTracker - vecSubdetName.at("<<i<<") = " << _vecSubdetName.at(i) << std::endl;
  //     streamlog_out(DEBUG2) << "ExtrToTracker - vecSubdetNLayers.at("<<i<<") = " << _vecSubdetNLayers.at(i) << std::endl;

  //   }
   

  // #else
  // #error Geometry type not defined
  // #endif


  //   /////////////////////////////
	    


  
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



  _maxChi2PerHit = 100;



  if (_doNtuple){  

    _out = new TFile(_outFileName.c_str(),"RECREATE");
    _tree = new TTree(_treeName.c_str(),_treeName.c_str());

    // init tree variables 
    this->clearEventVar();

    // set tree branches 
    int bufsize = 32000; //default buffer size 32KB
    //setTreeBranches(bufsize);    
  

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
  LCCollection* otHitsCol = this->GetCollection( evt, _otColName ) ;

  if( input_track_col != 0 ){
    

    // establish the track collection that will be created 
    LCCollectionVec* trackVec = new LCCollectionVec( LCIO::TRACK )  ;    

    
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
	

    // std::vector< IMPL::TrackImpl* > trackCandidates ;

    // loop over the input tracks and refit using KalTest    
    for(int i=0; i< nTracks ; ++i) {


      
      if (_doNtuple) this->clearLayerHelperVar();
     

      int SITHitsPerTrk = 0 ;
      
      Track* track = dynamic_cast<Track*>( inputTrackVec->getElementAt( i ) ) ;
      
      MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
     
      EVENT::TrackerHitVec trkHits = track->getTrackerHits() ;
      
	
      // sort the hits in R, so here we are assuming that the track came from the IP 
 
      sort(trkHits.begin(), trkHits.end(), ExtrToTracker::compare_r() );
	
      EVENT::TrackerHitVec::iterator it = trkHits.begin();
	
      for( it = trkHits.begin() ; it != trkHits.end() ; ++it ){
	marlin_trk->addHit(*it);
      }
	
      int init_status = FitInit2(track, marlin_trk) ;          
	
      if (init_status==0) {
	  

	;

	streamlog_out(DEBUG4) << "track initialised " << std::endl ;
	  
	int fit_status = marlin_trk->fit(); 
	  
	if ( fit_status == 0 ){
	    
	  int testFlag=0;




	    
	  double chi2 = 0 ;
	  int ndf = 0 ;
	  TrackStateImpl trkState;	
	    
	    
	  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
	    
	  encoder.reset() ;  // reset to 0
	    
	  int layerID = encoder.lowWord() ;  
	  int elementID = 0 ;    
	    
	    
	  //________________________________________________________________________________________________________
	  //
	  // starting loop on subdetectors and loop on each subdetector layer
	  //________________________________________________________________________________________________________
	      
	  printParameters();




	  ////////////////////////

	  fillVecSubdet(evt);

	  ////////////////////////


	  for (size_t idet=0; idet<_vecSubdetName.size(); idet++){
	    streamlog_out(DEBUG4) << "LOOP - idet = " << idet << " begins "<< std::endl;
        
	    if( _vecDigiHitsCol.at(idet) != 0 ){ 


	      for (int iL=0;iL<_vecSubdetNLayers.at(idet);iL++){
		streamlog_out(DEBUG4) << "LOOP" << iL << " begins "<< std::endl;
	 	
		encoder[lcio::ILDCellID0::subdet] = _vecSubdetID.at(idet);
		encoder[lcio::ILDCellID0::layer]  = iL;   
		layerID = encoder.lowWord();  
		streamlog_out(DEBUG4) << "layerID = " << layerID << std::endl;
		
		///////////////////////////////////////////////////////////

   

		if ( marlin_trk->propagateToLayer( layerID, trkState, chi2, ndf, elementID, IMarlinTrack::modeClosest) == MarlinTrk::IMarlinTrack::success) {
		    

		  streamlog_out(DEBUG4) << "-- layerID " << layerID << std::endl;


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
		  
		  
		  streamlog_out(DEBUG4) << " layer " << iL << " max search distances Z : Rphi " << _searchSigma*sqrt( covLCIO[9] ) << " : " << _searchSigma*sqrt( covLCIO[0] ) << std::endl ;


		  //_______________________________________________________________________________________
		  //
		  
		  streamlog_out(DEBUG2) << " element ID " << elementID << std::endl;
		  
		  if ( elementID != 0 ){
		    
		    testFlag = 1;
		    
		    float dU_spres = 0.007;
		    float dV_spres = 0.05;

		    bool isSuccessfulFit = false; 

		    int nhits=0;
		    TrackerHitPlane* BestHit;
		    BestHit = getSiHit(_vecDigiHitsCol.at(idet), elementID, marlin_trk, nhits);
		    if (BestHit != 0){
		      			  
		      streamlog_out(DEBUG4) << " --- Best hit found: call add and fit _Max_Chi2_Incr "<< _Max_Chi2_Incr<< std::endl ; 
						  
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
			streamlog_out(DEBUG4) << " --- chi2_increment "<< chi2_increment << std::endl ; 
			streamlog_out(DEBUG4) << " --- isSuccessfulFit "<< isSuccessfulFit << std::endl ; 
		      }

		  
		      double hitx = BestHit->getPosition()[0];
		      double hity = BestHit->getPosition()[1];
		      double hitz = BestHit->getPosition()[2];
		      double hitr = sqrt(hitx*hitx+hity*hity);
	
		
		      TotalSITHits++;
			
		      if ( isSuccessfulFit ){
			  
			streamlog_out(DEBUG4) << " successful fit " << std::endl ; 
			streamlog_out(DEBUG4) << " increment in the chi2 = " << chi2_increment << "  , max chi2 to accept the hit " << _Max_Chi2_Incr  << std::endl;
			  
			trkHits.push_back(BestHit) ;
			  
			  
			streamlog_out(DEBUG4) << " +++ hit added " << BestHit << std::endl ;
	      

			SITHitsPerTrk++;
			SITHitsFitted++;
			  
			//HitsInLayer.erase( HitsInLayer.begin() + pointer ) ;
			  
		      } //end successful fit
		      else{
			  
			SITHitsNonFitted++;
			streamlog_out(DEBUG4) << " +++ HIT NOT ADDED "<< std::endl;
			  
		      }

		      // if (doSinglePointResolutionSmearing) delete TestHitPlane ;

		    }   // besthit found
		    //    HitsInLayer.clear();
		    		  
		    // addHitOnNextElID(elementID+1, marlin_trk, trkHits, sitHitsCol, otHitsCol, iL, nSITR, TotalSITHits, SITHitsPerTrk, SITHitsFitted, SITHitsNonFitted);
		    // addHitOnNextElID(elementID-1, marlin_trk, trkHits, sitHitsCol, otHitsCol, iL, nSITR, TotalSITHits, SITHitsPerTrk, SITHitsFitted, SITHitsNonFitted);


		  }//elementID !=0 
		  
		} // successful propagation to layer
		  	      
	      } // loop to all subdetector layers

	    }//end colDigi not empty 

	  }//end loop on subdetectors
	  streamlog_out(DEBUG4) << " no of hits in the track (after adding SIT hits) " << trkHits.size() << " SIT hits added " << SITHitsPerTrk  << " event " <<  _n_evt<< std::endl;
	    
	    

	  //==============================================================================================================

	  IMPL::TrackImpl* lcio_trk = new IMPL::TrackImpl();

	  IMarlinTrack* marlinTrk = 0 ;

	  if( ! _performFinalRefit ) {

	    //fg: ------ here we just create a final LCIO track from the extrapolation :
	    
	    marlinTrk = marlin_trk ;
	    
	    bool fit_direction = IMarlinTrack::forward ;
	    int return_code =  finaliseLCIOTrack( marlin_trk, lcio_trk, trkHits,  fit_direction ) ;
	    
	    streamlog_out( DEBUG ) << " *** created finalized LCIO track - return code " << return_code  << std::endl 
				   << *lcio_trk << std::endl ;
	    

	  } else { //fg: ------- perform a final refit - does not work right now ...

	    // refitted track collection creation
	    if (  testFlag==1 ){
	      
	      
	      TrackStateImpl* trkState = new TrackStateImpl() ;
	      double chi2_fin = 0. ;
	      int ndf_fin = 0 ;
	      
	      marlin_trk->getTrackState(*trkState, chi2_fin, ndf_fin);
	      //const FloatVec& covMatrix = trkState->getCovMatrix();


	      //////////////////////////////////////////////////////////////////////////////////
	      

	      sort(trkHits.begin(), trkHits.end(), ExtrToTracker::compare_r() );



	      bool fit_backwards = IMarlinTrack::backward;
	      //bool fit_forwards = IMarlinTrack::forward;
	      MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();		


	      std::vector<EVENT::TrackerHit* > vec_hits;
	      vec_hits.clear();
	      for(int ih=0; ih<trkHits.size(); ih++){
		vec_hits.push_back(trkHits.at(ih));
	      }//end loop on hits
	      streamlog_out(DEBUG) << " --- vec_hits.size() = " <<   vec_hits.size()  <<std::endl;	


	      int ndf_test_0;
	      int return_error_0 = marlinTrk->getNDF(ndf_test_0);
	      streamlog_out(DEBUG3) << "++++ 0 - getNDF returns " << return_error_0 << std::endl;
	      streamlog_out(DEBUG3) << "++++ 0 - getNDF returns ndf = " << ndf_test_0 << std::endl;


	      //Kalman filter smoothing - fit track from out to in
	      int error_fit =  createFit(vec_hits, marlinTrk, trkState, _bField, fit_backwards, _maxChi2PerHit);
	      streamlog_out(DEBUG) << "---- createFit - error_fit = " << error_fit << std::endl;

	      bool fit_direction  = fit_backwards ;

	      if (error_fit == 0) {
		int error = finaliseLCIOTrack(marlinTrk, lcio_trk, vec_hits, fit_direction );
		streamlog_out(DEBUG) << "---- finalisedLCIOTrack - error = " << error << std::endl;
		
		int ndf_test;
		int return_error = marlinTrk->getNDF(ndf_test);
		streamlog_out(DEBUG3) << "++++ getNDF returns " << return_error << std::endl;
		streamlog_out(DEBUG3) << "++++ getNDF returns ndf = " << ndf_test << std::endl;


		if (error!=0){
            
		  streamlog_out(DEBUG3) << "Error from finaliseLCIOTrack non zero! deleting tracks. error=" << error <<" noHits: "<<trkHits.size()<<" marlinTrk: "<<marlinTrk<<" lcio_trk: "<<lcio_trk<< std::endl;
            
		  delete lcio_trk;
		  continue ; 
          
		}
	      } else {
		streamlog_out(DEBUG3) << "Error from createFit non zero! deleting tracks. error_fit=" << error_fit << std::endl;
              
		delete lcio_trk;
		continue ; 
        
	      }

	      delete trkState;

	    } // end of the creation of the refitted track collection
	  } // !_perfomFinalRefit 

	      
	  // fit finished - get hits in the fit
	  
	  std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
	  std::vector<std::pair<EVENT::TrackerHit* , double> > outliers;
	  
	  // remember the hits are ordered in the order in which they were fitted
	  
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
	    

	  if( _performFinalRefit ) delete marlinTrk ;

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
   

    evt->addCollection( trackVec , _output_track_col_name ) ;

    //delete trackVec;
    //delete inputTrackVec;
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



  
  if (_doNtuple) {

    _tree->Write();
    _out->Write();
    _out->Close();

    // delete _tree;
    // delete _out;

  }  
  
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







// TrackerHit* ExtrToTracker::getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk, int& nHitsOnDetEl){
TrackerHitPlane* ExtrToTracker::getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk, int& nHitsOnDetEl){

  if ( sitHitsCol != 0  ) { 
    int sitHits = sitHitsCol->getNumberOfElements();   
    //std::vector<TrackerHit* > hitsOnDetEl; 
    std::vector<TrackerHitPlane* > hitsOnDetEl; 
    hitsOnDetEl.clear();

    for(int i=0; i<sitHits; i++){
            
      //TrackerHit* hit = dynamic_cast<TrackerHit*>( sitHitsCol->getElementAt( i ) );	
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( sitHitsCol->getElementAt( i ) );	

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
      //TrackerHit* selectedHit = dynamic_cast<TrackerHit*>( hitsOnDetEl.at(index) ) ;
      TrackerHitPlane* selectedHit = dynamic_cast<TrackerHitPlane*>( hitsOnDetEl.at(index) ) ;
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









void ExtrToTracker::addHitOnNextElID(int elementID, MarlinTrk::IMarlinTrack*& marlin_trk, EVENT::TrackerHitVec& trkHits, LCCollection*& sitHitsCol, LCCollection*& otHitsCol, int& iL, int& nSITR, int& TotalSITHits, int& SITHitsPerTrk, int& SITHitsFitted, int& SITHitsNonFitted){

  if ( elementID != 0 ){
		    
    bool isSuccessfulFit = false; 

    int nhits=0;
    TrackerHitPlane* BestHit;
    if (iL<nSITR) BestHit = getSiHit(sitHitsCol, elementID, marlin_trk, nhits);
    else BestHit = getSiHit(otHitsCol, elementID, marlin_trk, nhits);
    if (BestHit != 0){
		      			  
      streamlog_out(DEBUG4) << " --- Best hit found: call add and fit _Max_Chi2_Incr "<< _Max_Chi2_Incr<< std::endl ; 
						  
      double chi2_increment = 0.;

      isSuccessfulFit = marlin_trk->addAndFit( BestHit, chi2_increment, _Max_Chi2_Incr ) == IMarlinTrack::success ;
      streamlog_out(DEBUG4) << " --- chi2_increment "<< chi2_increment << std::endl ; 
      streamlog_out(DEBUG4) << " --- isSuccessfulFit "<< isSuccessfulFit << std::endl ; 
		  
      double hitx = BestHit->getPosition()[0];
      double hity = BestHit->getPosition()[1];
      double hitz = BestHit->getPosition()[2];
      double hitr = sqrt(hitx*hitx+hity*hity);
	
		
      TotalSITHits++;
			
      if ( isSuccessfulFit ){
			  
	streamlog_out(DEBUG4) << " successful fit " << std::endl ; 
	streamlog_out(DEBUG4) << " incriment in the chi2 = " << chi2_increment << "  , max chi2 to accept the hit " << _Max_Chi2_Incr  << std::endl;
			  
	trkHits.push_back(BestHit) ;
			  			  
	streamlog_out(DEBUG4) << " +++ hit added " << BestHit << std::endl ;
			  

	SITHitsPerTrk++;
	SITHitsFitted++;
			  			  
      } //end successful fit
      else{
			  
	SITHitsNonFitted++;
	streamlog_out(DEBUG4) << " +++ HIT NOT ADDED "<< std::endl;
			  
      }

    }   // besthit found
		    
    streamlog_out(DEBUG4) << "LOOP" << iL << " ends "<< std::endl;
    streamlog_out(DEBUG4) << "###########$$$$$$$$$$##############" << std::endl;
  }  //end elementID != 0


}







void ExtrToTracker::fillVecSubdet(lcio::LCEvent*& evt){


  _vecSubdetNLayers.clear();
  _vecSubdetID.clear();
  _vecDigiHitsCol.clear();

  for (size_t idet=0; idet<_vecSubdetName.size(); idet++){
       
    streamlog_out(DEBUG4) << "idet = " << idet << std::endl;   
 
    int nSubdetLayers = 0;
    int idSubdet = 0;

    if ( _vecDigiHits.at(idet)!="" ){ 
	 	
      LCCollection* colDigi = 0 ;
      try{
	colDigi = evt->getCollection( _vecDigiHits.at(idet) ) ;
      }catch(DataNotAvailableException &e){
	streamlog_out(DEBUG4) << "Collection " << _vecDigiHits.at(idet).c_str() << " is unavailable " << std::endl;
      }	

      if( colDigi != 0 ){ 
	CellIDDecoder<TrackerHitPlane> cellid_decoder(colDigi);
	int last_element = colDigi->getNumberOfElements()-1;
	streamlog_out(DEBUG2) << "last_element = " << last_element << std::endl;
	if (last_element>=0) {
	  _vecDigiHitsCol.push_back(colDigi);
	  TrackerHitPlane* hit_helper2 = dynamic_cast<TrackerHitPlane*>( colDigi->getElementAt(last_element) ) ;
	  nSubdetLayers = cellid_decoder( hit_helper2 )["layer"] + 1; //from 0 to n-1
	  idSubdet = cellid_decoder( hit_helper2 )["subdet"]; 
	}
	else   _vecDigiHitsCol.push_back(0); 
      }//end colDigi!=0 
      else {
	_vecDigiHitsCol.push_back(colDigi);
      }
      streamlog_out(DEBUG2) << "_vecDigiHitsCol.back() = " << _vecDigiHitsCol.back() << std::endl;
      
    }// diginame !=0

    _vecSubdetNLayers.push_back(nSubdetLayers);
    _vecSubdetID.push_back(idSubdet);

    streamlog_out(DEBUG4) << "nSubdetLayers = " << nSubdetLayers << std::endl;
    streamlog_out(DEBUG4) << "idSubdet = " << idSubdet << std::endl;

  }//end loop on subdetectors
	  
}//end 
