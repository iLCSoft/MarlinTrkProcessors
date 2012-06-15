
#include "SiliconTracking_MarlinTrk.h"


#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>

//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_linalg.h>

#include <marlin/Global.h>

//#include "ClusterShapes.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include "gear/FTDLayerLayout.h"
#include "gear/FTDParameters.h"

#include <gear/BField.h>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"

#include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

using namespace lcio ;
using namespace marlin ;
using namespace MarlinTrk ;

using std::min;
using std::max;
using std::abs;

const double SiliconTracking_MarlinTrk::TWOPI = 2*M_PI;

SiliconTracking_MarlinTrk aSiliconTracking_MarlinTrk ;

SiliconTracking_MarlinTrk::SiliconTracking_MarlinTrk() : Processor("SiliconTracking_MarlinTrk") {

 _description = "Pattern recognition in silicon trackers";

 _fastfitter = new MarlinTrk::HelixFit();

 _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);

 std::vector<int> combinations;

 combinations.push_back(6);
 combinations.push_back(4);
 combinations.push_back(3);

 combinations.push_back(6);
 combinations.push_back(4);
 combinations.push_back(2);

 combinations.push_back(6);
 combinations.push_back(3);
 combinations.push_back(2);

 combinations.push_back(5);
 combinations.push_back(4);
 combinations.push_back(3);

 combinations.push_back(5);
 combinations.push_back(4);
 combinations.push_back(2);

 combinations.push_back(5);
 combinations.push_back(3);
 combinations.push_back(2);

 combinations.push_back(5);
 combinations.push_back(3);
 combinations.push_back(1);

 combinations.push_back(5);
 combinations.push_back(2);
 combinations.push_back(1);

 combinations.push_back(4);
 combinations.push_back(3);
 combinations.push_back(2);

 combinations.push_back(4);
 combinations.push_back(3);
 combinations.push_back(1);

 combinations.push_back(4);
 combinations.push_back(2);
 combinations.push_back(1);

 combinations.push_back(3);
 combinations.push_back(2);
 combinations.push_back(1);


 registerProcessorParameter("LayerCombinations",
                            "Combinations of Hits in Layers",
                            _Combinations,
                            combinations);

 std::vector<int> combinationsFTD;

 combinationsFTD.push_back(6);
 combinationsFTD.push_back(5);
 combinationsFTD.push_back(4);

 combinationsFTD.push_back(6);
 combinationsFTD.push_back(5);
 combinationsFTD.push_back(3);

 combinationsFTD.push_back(5);
 combinationsFTD.push_back(4);
 combinationsFTD.push_back(3);

 combinationsFTD.push_back(5);
 combinationsFTD.push_back(4);
 combinationsFTD.push_back(2);

 combinationsFTD.push_back(5);
 combinationsFTD.push_back(4);
 combinationsFTD.push_back(1);

 combinationsFTD.push_back(5);
 combinationsFTD.push_back(3);
 combinationsFTD.push_back(2);

 combinationsFTD.push_back(5);
 combinationsFTD.push_back(3);
 combinationsFTD.push_back(1);

 combinationsFTD.push_back(5);
 combinationsFTD.push_back(2);
 combinationsFTD.push_back(1);

 combinationsFTD.push_back(4);
 combinationsFTD.push_back(3);
 combinationsFTD.push_back(2);

 combinationsFTD.push_back(4);
 combinationsFTD.push_back(3);
 combinationsFTD.push_back(1);

 combinationsFTD.push_back(4);
 combinationsFTD.push_back(3);
 combinationsFTD.push_back(0);

 combinationsFTD.push_back(4);
 combinationsFTD.push_back(2);
 combinationsFTD.push_back(1);

 combinationsFTD.push_back(4);
 combinationsFTD.push_back(2);
 combinationsFTD.push_back(0);

 combinationsFTD.push_back(4);
 combinationsFTD.push_back(1);
 combinationsFTD.push_back(0);

 combinationsFTD.push_back(3);
 combinationsFTD.push_back(2);
 combinationsFTD.push_back(1);

 combinationsFTD.push_back(3);
 combinationsFTD.push_back(2);
 combinationsFTD.push_back(0);

 combinationsFTD.push_back(3);
 combinationsFTD.push_back(1);
 combinationsFTD.push_back(0);

 combinationsFTD.push_back(2);
 combinationsFTD.push_back(1);
 combinationsFTD.push_back(0);


 registerProcessorParameter("LayerCombinationsFTD",
                            "Combinations of Hits in FTD",
                            _CombinationsFTD,
                            combinationsFTD);

 registerProcessorParameter("NDivisionsInPhi",
                            "Number of divisions in Phi",
                            _nDivisionsInPhi,
                            int(80));

 registerProcessorParameter("NDivisionsInPhiFTD",
                            "Number of divisions in Phi for FTD",
                            _nPhiFTD,
                            int(30));

 registerProcessorParameter("NDivisionsInTheta",
                            "Number of divisions in Theta",
                            _nDivisionsInTheta,
                            int(80));

 // Input Collections
 // ^^^^^^^^^^^^^^^^^
 registerInputCollection(LCIO::TRACKERHITPLANE,
                         "VTXHitCollectionName",
                         "VTX Hit Collection Name",
                         _VTXHitCollection,
                         std::string("VTXTrackerHits"));


 registerInputCollection(LCIO::TRACKERHITPLANE,
                         "FTDPixelHitCollectionName",
                         "FTD Pixel Hit Collection Name",
                         _FTDPixelHitCollection,
                         std::string("FTDPixelTrackerHits"));  

  registerInputCollection(LCIO::TRACKERHIT,
                          "FTDSpacePointCollectionName",
                          "FTD FTDSpacePoint Collection Name",
                          _FTDSpacePointCollection,
                          std::string("FTDSpacePoints"));  
  

 registerInputCollection(LCIO::TRACKERHIT,
                         "SITHitCollectionName",
                         "SIT Hit Collection Name",
                         _SITHitCollection,
                         std::string("SITTrackerHits"));  

 // Output Collections
 // ^^^^^^^^^^^^^^^^^^
 registerOutputCollection(LCIO::TRACK,
                          "SiTrackCollectionName",
                          "Silicon track Collection Name",
                          _siTrkCollection,
                          std::string("SiTracks"));



 // Steering parameters
 // ^^^^^^^^^^^^^^^^^^^
 registerProcessorParameter("Chi2WRphiTriplet",
                            "Chi2WRphiTriplet",
                            _chi2WRPhiTriplet,
                            float(1.));

 registerProcessorParameter("Chi2WRphiQuartet",
                            "Chi2WRphiQuartet",
                            _chi2WRPhiQuartet,
                            float(1.));

 registerProcessorParameter("Chi2WRphiSeptet",
                            "Chi2WRphiSeptet",
                            _chi2WRPhiSeptet,
                            float(1.));

 registerProcessorParameter("Chi2WZTriplet",
                            "Chi2WZTriplet",
                            _chi2WZTriplet,
                            float(0.5));

 registerProcessorParameter("Chi2WZQuartet",
                            "Chi2WZQuartet",
                            _chi2WZQuartet,
                            float(0.5));

 registerProcessorParameter("Chi2WZSeptet",
                            "Chi2WZSeptet",
                            _chi2WZSeptet,
                            float(0.5));

 registerProcessorParameter("Chi2FitCut",
                            "Chi2 Fit Cut",
                            _chi2FitCut,
                            float(120.0));

 registerProcessorParameter("Chi2PrefitCut",
                            "Chi2 Prefit Cut",
                            _chi2PrefitCut,
                            float(1.0e+10));

 registerProcessorParameter("AngleCutForMerging",
                            "Angle Cut For Merging",
                            _angleCutForMerging,
                            float(0.1));

 registerProcessorParameter("MinDistCutAttach",
                            "MinDistCutAttach",
                            _minDistCutAttach,
                            float(2.5));

 registerProcessorParameter("MinLayerToAttach",
                            "MinLayerToAttach",
                            _minimalLayerToAttach,
                            int(-1));

 registerProcessorParameter("CutOnD0",
                            "cut on D0 for tracks",
                            _cutOnD0,
                            float(100.0));

 registerProcessorParameter("CutOnZ0",
                            "cut on Z0 for tracks",
                            _cutOnZ0,
                            float(100.0));

 registerProcessorParameter("CutOnPt",
                            "cut on Pt",
                            _cutOnPt,
                            float(0.05));

 registerProcessorParameter("MinimalHits",
                            "minimal hits (default 3)",
                            _minimalHits,
                            int(3));

 registerProcessorParameter("NHitsChi2",
                            "Maximal number of hits for which a track with n hits is better than one with n-1hits. (defaut 5)",
                            _nHitsChi2,
                            int(5));

 registerProcessorParameter("FastAttachment",
                            "Fast attachment",
                            _attachFast,
                            int(0));

 registerProcessorParameter("UseSIT",
                            "Use SIT",
                            _useSIT,
                            int(1));

  
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
  
  
 registerProcessorParameter("Debug",
                            "Print out debugging info?",
                            _debug,
                            int(1));


 registerProcessorParameter("CheckForDelta",
                            "Check for Delta rays hits in hit-to-track assignment",
                            _checkForDelta,
                            int(1));

 registerProcessorParameter("MinDistToDelta",
                            "Minimal distance of track hit to the delta electron hit",
                            _minDistToDelta,
                            float(0.25));

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

 registerProcessorParameter("ReadingLoiData",
                            "Legacy mode for reading loi data",
                            _reading_loi_data,
                            bool(false));

#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  registerOptionalParameter("RunMarlinTrkDiagnostics", "Run MarlinTrk Diagnostics. MarlinTrk must be compiled with MARLINTRK_DIAGNOSTICS_ON defined", _runMarlinTrkDiagnostics, bool(false));
  
  registerOptionalParameter("DiagnosticsName", "Name of the root file and root tree if running Diagnostics", _MarlinTrkDiagnosticsName, std::string("SiliconTrackingDiagnostics"));    
  
#endif

  

}



void SiliconTracking_MarlinTrk::init() { 

 _nRun = -1 ;
 _nEvt = 0 ;
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

#ifdef MARLINTRK_DIAGNOSTICS_ON
  
  void * dcv = _trksystem->getDiagnositicsPointer();
  DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
  dc->init(_MarlinTrkDiagnosticsName,_MarlinTrkDiagnosticsName, _runMarlinTrkDiagnostics);
  
#endif

 this->setupGearGeom(Global::GEAR);

  
  // planar sit is arranged in double layers which creates spacepoints per stereo layer
  // this means that here there are only "half" the number layers
  if ( _reading_loi_data == false ) {
    _nLayersSIT = _nLayersSIT/2 ;
  }
  
 if (_useSIT == 0)
   _nLayers = _nLayersVTX;
 else 
   _nLayers = _nLayersVTX + _nLayersSIT;

 // initialise the container to have separate vectors for up to _nHitsChi2 hits.
 _tracksWithNHitsContainer.resize(_nHitsChi2);

 _dPhi = TWOPI/_nDivisionsInPhi;
 _dTheta = 2.0/_nDivisionsInTheta;
 _dPhiFTD = TWOPI/_nPhiFTD;
 // I leave this for the moment, but 0.3 is c/1e9.
 // For the cut it does not make too much of a difference
 double cutOnR = _cutOnPt/(0.3*_bField);
 cutOnR = 1000.*cutOnR;
 _cutOnOmega = 1/cutOnR;


}


void SiliconTracking_MarlinTrk::processRunHeader( LCRunHeader* run) { 

 _nRun++ ;
 _nEvt = 0;

 streamlog_out(MESSAGE) << "SiliconTracking_MarlinTrk ---> new run : run number = " << _nRun << std::endl;

} 

void SiliconTracking_MarlinTrk::processEvent( LCEvent * evt ) { 

 // Clearing the working containers from the previous event
 // FIXME: partly done at the end of the event, in CleanUp. Make it consistent.
 _tracksWithNHitsContainer.clear();
 _trackImplVec.clear();

 streamlog_out(DEBUG4) << "SiliconTracking_MarlinTrk -> run = " << _nRun 
 << "  event = " << _nEvt << std::endl;

  int successVTX = InitialiseVTX( evt );
  int successFTD = InitialiseFTD( evt );

 if (successVTX == 1) {

   streamlog_out(DEBUG1) << "   phi \t\t theta \t\tlayer \t nh o : m : i : o*m*i " << std::endl; 

   for (int iPhi=0; iPhi<_nDivisionsInPhi; ++iPhi) { 
     for (int iTheta=0; iTheta<_nDivisionsInTheta;++iTheta) {
       ProcessOneSector(iPhi,iTheta); // Process one VXD sector     
     }
   }
   
   streamlog_out(DEBUG4) << "End of Processing VXD and SIT sectors" << std::endl;

 }

 if (successFTD == 1) {
   TrackingInFTD(); // Perform tracking in the FTD
   streamlog_out(DEBUG4) << "End of Processing FTD sectors" << std::endl;
 }

  if (successVTX == 1 || successFTD == 1) {
    //if (successVTX == 1 ) {

    for (int nHits = _nHitsChi2; nHits >= 3 ;// the three is hard coded, sorry.
                                                      // It's the minimal number to form a track
         nHits--) {
      Sorting( _tracksWithNHitsContainer.getTracksWithNHitsVec( nHits ) );
   }
   streamlog_out(DEBUG4) <<  "End of Sorting " << std::endl;

   for (int nHits = _nHitsChi2; nHits >= 3 ;// the three is hard coded, sorry.
                                                      // It's the minimal number to form a track
	 nHits--) {
     TrackExtendedVec &tracksWithNHits = _tracksWithNHitsContainer.getTracksWithNHitsVec( nHits );
     for (TrackExtendedVec::iterator trackIter = tracksWithNHits.begin();
	   trackIter < tracksWithNHits.end(); trackIter++) {
       CreateTrack( *trackIter );
     }
     streamlog_out(DEBUG4) <<  "End of creating "<< nHits << " hits tracks " << std::endl;
   }

   if (_attachFast == 0) {
     AttachRemainingVTXHitsSlow();
     AttachRemainingFTDHitsSlow();
   }
   else {
     AttachRemainingVTXHitsFast();
     AttachRemainingFTDHitsFast();
   }
   streamlog_out(DEBUG4) <<  "End of picking up remaining hits " << std::endl;


   LCCollectionVec * trkCol = new LCCollectionVec(LCIO::TRACK);
   // if we want to point back to the hits we need to set the flag
   LCFlagImpl trkFlag(0) ;
   trkFlag.setBit( LCIO::TRBIT_HITS ) ;
   trkCol->setFlag( trkFlag.getFlag()  ) ;

   LCCollectionVec * relCol = NULL;


   FinalRefit(trkCol, relCol);


   evt->addCollection(trkCol,_siTrkCollection.c_str());     

  }
 CleanUp();
 streamlog_out(DEBUG4) << "Event is done " << std::endl;
 _nEvt++;

}


void SiliconTracking_MarlinTrk::CleanUp() {

 _tracksWithNHitsContainer.clear();

 for (int il=0;il<_nLayers;++il) {
   for (int ip=0;ip<_nDivisionsInPhi;++ip) {
     for (int it=0;it<_nDivisionsInTheta; ++it) {
       unsigned int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      

       if( iCode >= _sectors.size()){          
         std::cerr<< "iCode index out of range: iCode =   " << iCode << " _sectors.size() = " << _sectors.size() << " exit(1) called from file " << __FILE__ << " line " << __LINE__<< std::endl;
         exit(1);
       }


       TrackerHitExtendedVec& hitVec = _sectors[iCode];
       int nH = int(hitVec.size());
       for (int iH=0; iH<nH; ++iH) {
         TrackerHitExtended * hit = hitVec[iH];
         delete hit;
       }
     }
   }
 }

 for (int iS=0;iS<2;++iS) {
   for (unsigned int layer=0;layer<_nlayersFTD;++layer) {
     for (int ip=0;ip<_nPhiFTD;++ip) {
       unsigned int iCode = iS + 2*layer + 2*_nlayersFTD*ip;

       if( iCode >= _sectorsFTD.size()){          
         std::cerr<< "iCode index out of range: iCode =   " << iCode << " _sectorsFTD.size() = " << _sectorsFTD.size() << " exit(1) called from file " << __FILE__ << " line " << __LINE__<< std::endl;
         exit(1);
       }


       TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
       int nH = int(hitVec.size());
       for (int iH=0; iH<nH; ++iH) {
         TrackerHitExtended * hit = hitVec[iH];
         delete hit;
       }
     }
   }
 }

}

int SiliconTracking_MarlinTrk::InitialiseFTD(LCEvent * evt) {

 int success = 1;

 _nTotalFTDHits = 0;
 _sectorsFTD.clear();
 _sectorsFTD.resize(2*_nlayersFTD*_nPhiFTD);

  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    
    LCCollection * hitCollection = evt->getCollection(_FTDPixelHitCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of FTD Pixel Hits = " << nelem << std::endl;
    _nTotalFTDHits = nelem;
    
    for (int ielem=0; ielem<nelem; ++ielem) {
      
      TrackerHitPlane * hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
                  
      double point_res_rphi = sqrt( hit->getdU()*hit->getdU() + hit->getdV()*hit->getdV() );
            
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));
      
      double pos[3];
      
      for (int i=0; i<3; ++i) {
        pos[i] = hit->getPosition()[i];
      }
      
      double Phi = atan2(pos[1],pos[0]);
      if (Phi < 0.) Phi = Phi + TWOPI;
      
      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));
      
      if ( _reading_loi_data == false ) {
        
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals 
        if ( petalIndex % 2 == 0 ) {
          layer = 2*layer;
        }
        else {
          layer = 2*layer + 1;
        }
        
      }
      
      if (layer >= _nlayersFTD) {
        streamlog_out(ERROR) << "SiliconTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nlayersFTD <<  std::endl;
        exit(1);
      }
      
      int iPhi = int(Phi/_dPhiFTD);
      
      int side = getSideID(hit);
      int iSemiSphere = 0;
      
      if (side > 0) 
        iSemiSphere = 1;
      
      int iCode = iSemiSphere + 2*layer + 2*_nlayersFTD*iPhi;
      _sectorsFTD[iCode].push_back( hitExt );
      
    }
  }
  catch(DataNotAvailableException &e ) {
    success = 0;
  }


 // Reading out FTD SpacePoint Collection
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 try {

   LCCollection * hitCollection = evt->getCollection(_FTDSpacePointCollection.c_str());

   int nelem = hitCollection->getNumberOfElements();

   streamlog_out(DEBUG4) << "Number of FTD SpacePoints = " << nelem << std::endl;
   _nTotalFTDHits += nelem;

   for (int ielem=0; ielem<nelem; ++ielem) {
     
     TrackerHit * hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));

     TrackerHitExtended * hitExt = new TrackerHitExtended( hit );

     double point_res_rphi = sqrt( hit->getCovMatrix()[0] );
     
     hitExt->setResolutionRPhi( point_res_rphi );

     // SJA:FIXME why is this needed? 
     hitExt->setResolutionZ(0.1);

     // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
     hitExt->setType(int(INT_MAX));
     // det is no longer used set to INT_MAX to try and catch any missuse
     hitExt->setDet(int(INT_MAX));

     double pos[3];

     for (int i=0; i<3; ++i) {
       pos[i] = hit->getPosition()[i];
     }

     double Phi = atan2(pos[1],pos[0]);
     if (Phi < 0.) Phi = Phi + TWOPI;

     // get the layer number
     unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
     unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));

     if ( _reading_loi_data == false ) {

       // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
       // the layers will follow the even odd numbering of the petals 
       if ( petalIndex % 2 == 0 ) {
         layer = 2*layer;
       }
       else {
         layer = 2*layer + 1;
       }

     }

     if (layer >= _nlayersFTD) {
       streamlog_out(ERROR) << "SiliconTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : " << layer << " number of layers = " << _nlayersFTD <<  std::endl;
       exit(1);
     }

     int iPhi = int(Phi/_dPhiFTD);

     int side = getSideID(hit);
     int iSemiSphere = 0;

     if (side > 0) 
       iSemiSphere = 1;

     int iCode = iSemiSphere + 2*layer + 2*_nlayersFTD*iPhi;
     _sectorsFTD[iCode].push_back( hitExt );

     streamlog_out( DEBUG1 ) << " FTD Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << "  iPhi = " << iPhi <<  " iSemiSphere "  << iSemiSphere << " iCode = " << iCode << "  layer = " << layer << std::endl;  
     
   }
 }
 catch(DataNotAvailableException &e ) {
   success = 0;
 }
  
 return success;
}

int SiliconTracking_MarlinTrk::InitialiseVTX(LCEvent * evt) {

 int success = 1;

 _nTotalVTXHits = 0;
 _nTotalSITHits = 0;
 _sectors.clear();
 _sectors.resize(_nLayers*_nDivisionsInPhi*_nDivisionsInTheta);


 // Reading out VTX Hits Collection
 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
 try {

   LCCollection * hitCollection = evt->getCollection(_VTXHitCollection.c_str());

   int nelem = hitCollection->getNumberOfElements();

   streamlog_out(DEBUG4) << "Number of VTX hits = " << nelem << std::endl;
   _nTotalVTXHits = nelem;

   for (int ielem=0; ielem<nelem; ++ielem) {
     TrackerHitPlane * hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));

     gear::Vector3D U(1.0,hit->getU()[1],hit->getU()[0],gear::Vector3D::spherical);
     gear::Vector3D V(1.0,hit->getV()[1],hit->getV()[0],gear::Vector3D::spherical);
     gear::Vector3D Z(0.0,0.0,1.0);

     const float eps = 1.0e-07;
     // V must be the global z axis 
     if( fabs(1.0 - V.dot(Z)) > eps ) {
       streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: VXD Hit measurment vectors V is not equal to the global Z axis. \n\n  exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
       exit(1);
     }

     if( fabs(U.dot(Z)) > eps ) {
       streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: VXD Hit measurment vectors U is not in the global X-Y plane. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
       exit(1);
     }


     TrackerHitExtended * hitExt = new TrackerHitExtended( hit );


     // SJA:FIXME: just use planar res for now
     hitExt->setResolutionRPhi(hit->getdU());
     hitExt->setResolutionZ(hit->getdV());

     // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
     hitExt->setType(int(INT_MAX));
     // det is no longer used set to INT_MAX to try and catch any missuse
     hitExt->setDet(int(INT_MAX));

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
     int layer = getLayerID(hit);

     if (layer < 0 || layer >= _nLayers) {
       streamlog_out(ERROR) << "SiliconTracking_MarlinTrk => fatal error in VTX : layer is outside allowed range : " << layer << std::endl;
       exit(1);
     }

     int iPhi = int(Phi/_dPhi);
     int iTheta = int ((cosTheta + double(1.0))/_dTheta);
     int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
     _sectors[iCode].push_back( hitExt );

     streamlog_out( DEBUG1 ) << " VXD Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi <<  " iTheta "  << iTheta << " iCode = " << iCode << "  layer = " << layer << std::endl;  
   
   }
 }
 catch(DataNotAvailableException &e) {
   streamlog_out( DEBUG4 ) << " collection not found : " << _VTXHitCollection.c_str() << std::endl ;
   success = 0;
 }

 if (_useSIT > 0 ) {
   // Reading out SIT Hits Collection
   //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   try {
     LCCollection *hitCollection = evt->getCollection(_SITHitCollection.c_str());

     int nelem = hitCollection->getNumberOfElements();

     streamlog_out(DEBUG4) << "Number of SIT hits = " << nelem << std::endl;
     _nTotalSITHits = nelem;

     for (int ielem=0; ielem<nelem; ++ielem) {

       double drphi(NAN);
       double dz(NAN);

       TrackerHit* hit(NULL);

       if (_reading_loi_data) { // uses cylindrical hits

         TrackerHitZCylinder * hit_c = dynamic_cast<TrackerHitZCylinder*>(hitCollection->getElementAt(ielem));

         if (hit_c) {
           hit= hit_c;
         }
         else{
           streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: Dynamic cast to TrackerHitZCylinder failed. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
           exit(1);
         }

         drphi = hit_c->getdRPhi();
         dz    = hit_c->getdZ();

       }


       
       else {
         
         
         TrackerHit * hit_p = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));
         
         if (hit_p) {
           hit= hit_p;
         }
         else{
           streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: Dynamic cast to TrackerHit failed. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
           exit(1);
         }

//         double x = hit->getPosition()[0];
//         double y = hit->getPosition()[1];
         
//         double atan_xy = atan2(y, x);
//         double r = sqrt(x*x+y*y);
         
//         double drphidx = ( x * atan_xy - y ) / r ;
//         double drphidy = ( y * atan_xy + x ) / r ;
//         
//         drphi =  sqrt( drphidx*drphidx*hit->getCovMatrix()[0] + drphidy*drphidy*hit->getCovMatrix()[2] + 2.0 * drphidx * drphidy * hit->getCovMatrix()[1] );

         // SJA:FIXME: fudge for now by a factor of two and ignore covariance
         drphi =  2 * sqrt(hit->getCovMatrix()[0] + hit->getCovMatrix()[2]);

         dz =     sqrt(hit->getCovMatrix()[5]);
                  
       }

       
       TrackerHitExtended * hitExt = new TrackerHitExtended( hit );

       // SJA:FIXME: just use planar res for now
       hitExt->setResolutionRPhi(drphi);
       hitExt->setResolutionZ(dz);

       // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and catch any missuse
       hitExt->setType(int(INT_MAX));
       // det is no longer used set to INT_MAX to try and catch any missuse
       hitExt->setDet(int(INT_MAX));

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

       // VXD and SIT are treated as one system so SIT layers start from _nLayersVTX
       // SIT layer has to be divided by 2 as the space points are from stereo layers
       int layer = getLayerID(hit)/2 + _nLayersVTX;

       if (layer < 0 || layer >= _nLayers) {
         streamlog_out(ERROR) << "SiliconTracking_MarlinTrk => fatal error in SIT : layer is outside allowed range : " << layer << std::endl;
         exit(1);
       }
       int iPhi = int(Phi/_dPhi);
       int iTheta = int ((cosTheta + double(1.0))/_dTheta);
       int iCode = layer + _nLayers*iPhi + _nLayers*_nDivisionsInPhi*iTheta;      
       _sectors[iCode].push_back( hitExt );

       streamlog_out( DEBUG1 ) << " SIT Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi <<  " iTheta "  << iTheta << " iCode = " << iCode << "  layer = " << layer << std::endl;  
     
     }
   }
   catch(DataNotAvailableException &e) {

   }
 }
 return success;

}

void SiliconTracking_MarlinTrk::check(LCEvent * evt) { 

}

void SiliconTracking_MarlinTrk::end() {

  delete _fastfitter ;
  delete _encoder ;
  delete _trksystem ;

}


void SiliconTracking_MarlinTrk::ProcessOneSector(int iPhi, int iTheta) {

 int counter = 0 ;

 int iPhi_Up    = iPhi + 1;
 int iPhi_Low   = iPhi - 1;
 int iTheta_Up  = iTheta + 1; 
 int iTheta_Low = iTheta - 1;
 if (iTheta_Low < 0) iTheta_Low = 0;
 if (iTheta_Up  >= _nDivisionsInTheta) iTheta_Up = _nDivisionsInTheta-1;

 int nComb = int( _Combinations.size() / 3 ); // number of triplet combinations
                                              //  std::cout << iPhi << " " << iTheta << " " << _nEvt << std::endl;
 int iNC = 0;

 for (int iComb=0; iComb < nComb; ++iComb) { // loop over triplets
   
   int nLR[3];

   for (int iS=0; iS<3; ++iS) {
     nLR[iS] = _Combinations[iNC];
     iNC++;
   }    

//   std::cout << iPhi << " " << iTheta << " " << nLR[0] << " " << nLR[1] << " " << nLR[2] << " " << std::endl;

   int iCode = nLR[0] + _nLayers*iPhi +  _nLayers*_nDivisionsInPhi*iTheta;

//   std::cout << "size of vector = " << _sectors.size() << " iCode = " << iCode << std::endl;

   TrackerHitExtendedVec& hitVecOuter =  _sectors[iCode]; 

   int nHitsOuter = int(hitVecOuter.size());
   if (nHitsOuter > 0) {

//     std::cout << " " << iPhi << " " << iTheta << " " << nLR[0] << " " << nLR[1] << " " << nLR[2] << " size of vector = " << hitVecOuter.size() << std::endl;
     
     for (int ipMiddle=iPhi_Low; ipMiddle<iPhi_Up+1;ipMiddle++) { // loop over phi in the Middle
       for (int itMiddle=iTheta_Low; itMiddle<iTheta_Up+1;itMiddle++) { // loop over theta in the Middle 
         int iPhiMiddle = ipMiddle;
         if (ipMiddle < 0) iPhiMiddle = _nDivisionsInPhi-1;
         if (ipMiddle >= _nDivisionsInPhi) iPhiMiddle = ipMiddle - _nDivisionsInPhi;
         iCode = nLR[1] + _nLayers*iPhiMiddle +  _nLayers*_nDivisionsInPhi*itMiddle;
         TrackerHitExtendedVec& hitVecMiddle = _sectors[iCode];
         int nHitsMiddle = int(hitVecMiddle.size());

         int iPhiLowInner = iPhi_Low;
         int iPhiUpInner = iPhi_Up;
         int iThetaLowInner = iTheta_Low;
         int iThetaUpInner = iTheta_Up;        

         if (ipMiddle == iPhi && itMiddle==iTheta) { 
           iPhiLowInner = iPhi_Low;
           iPhiUpInner  = iPhi_Up;
           iThetaLowInner = iTheta_Low;
           iThetaUpInner = iTheta_Up;
         }
         else {
           int difP = abs(ipMiddle-iPhi);
           int difT = abs(itMiddle-iTheta);
           int minP = min(ipMiddle,iPhi);
           int minT = min(itMiddle,iTheta);
           int maxP = max(ipMiddle,iPhi);
           int maxT = max(itMiddle,iTheta);
           if (difP==1 && difT==1) {
             iPhiLowInner = minP;
             iPhiUpInner = maxP;
             iThetaLowInner = minT;
             iThetaUpInner = maxT;
           }
           if (difP==0) {
             iPhiLowInner = iPhi_Low;
             iPhiUpInner  = iPhi_Up;
             iThetaLowInner = minT;
             iThetaUpInner = maxT;
           }
           if (difT==0) {
             iPhiLowInner = minP;
             iPhiUpInner  = maxP;
             iThetaLowInner = iTheta_Low;
             iThetaUpInner = iTheta_Up;            
           }
         }               
         if (nHitsMiddle > 0) {
           for (int ipInner=iPhiLowInner; ipInner<iPhiUpInner+1;ipInner++) { // loop over phi in the Inner
             for (int itInner=iThetaLowInner; itInner<iThetaUpInner+1;itInner++) { // loop over theta in the Inner 
               int iPhiInner = ipInner;
               if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
               if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
               iCode = nLR[2] + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
               TrackerHitExtendedVec& hitVecInner = _sectors[iCode];
               int nHitsInner = int(hitVecInner.size());
               if (nHitsInner > 0) {

                 streamlog_out(DEBUG1) << " " 
                 << iPhi << " " << ipMiddle << " " << ipInner << "\t" 
                 << iTheta << " " << itMiddle << " " << itInner << "\t" 
                 << nLR[0] << " " << nLR[1] << " " << nLR[2] << "\t    " 
                 << nHitsOuter << " : " << nHitsMiddle << " : " << nHitsInner << "  :: " 
                 << nHitsOuter*nHitsMiddle* nHitsInner << std::endl;

                 for (int iOuter=0; iOuter<nHitsOuter; ++iOuter) { // loop over hits in the outer sector
                   TrackerHitExtended * outerHit = hitVecOuter[iOuter];
                   for (int iMiddle=0;iMiddle<nHitsMiddle;iMiddle++) { // loop over hits in the middle sector
                     TrackerHitExtended * middleHit = hitVecMiddle[iMiddle];
                     for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the inner sector
                       TrackerHitExtended * innerHit = hitVecInner[iInner];
                       HelixClass helix;
                       TrackExtended * trackAR = TestTriplet(outerHit,middleHit,innerHit,helix);
                       if ( trackAR != NULL ) {
                         int nHits = BuildTrack(outerHit,middleHit,innerHit,helix,nLR[2],
                                             iPhiLowInner,iPhiUpInner,
                                             iThetaLowInner,iThetaUpInner,trackAR);
			  
                         _tracksWithNHitsContainer.getTracksWithNHitsVec(nHits).push_back(trackAR);

                         counter ++ ;
                       }       
                     } // endloop over hits in the inner sector
                   } // endloop over hits in the middle sector
                 } // endloop over hits in the outer sector
               } // endif nHitsInner > 0
             } // endloop over theta in the Inner
           } // endloop over phi in the Inner      
         } // endif nHitsMiddle > 0
       } // endloop over theta in the Middle
     } // endloop over phi in the Middle
   } // endif nHitsOuter > 0
 } // endloop over triplets


 //  streamlog_out( DEBUG2 ) << " process one sectector theta,phi " << iTheta << ", " << iPhi <<
 //  "  number of loops : " << counter << std::endl  ;
}

TrackExtended * SiliconTracking_MarlinTrk::TestTriplet(TrackerHitExtended * outerHit, 
                                                      TrackerHitExtended * middleHit,
                                                      TrackerHitExtended * innerHit,
                                                      HelixClass & helix) {
 /*
  Methods checks if the triplet of hits satisfies helix hypothesis
  */


  
 TrackExtendedVec& trackOuterVec  = outerHit->getTrackExtendedVec();
 TrackExtendedVec& trackMiddleVec = middleHit->getTrackExtendedVec();
 TrackExtendedVec& trackInnerVec  = innerHit->getTrackExtendedVec();

 if ( (!trackOuterVec.empty())  && (!trackMiddleVec.empty()) && (!trackInnerVec.empty())) {

   TrackExtendedVec::const_iterator middleEndIter = trackMiddleVec.end();
   TrackExtendedVec::const_iterator outerEndIter  = trackOuterVec.end();
   TrackExtendedVec::const_iterator innerEndIter  = trackInnerVec.end();
   TrackExtendedVec::const_iterator outerBeginIter  = trackOuterVec.begin();
   TrackExtendedVec::const_iterator innerBeginIter  = trackInnerVec.begin();
   
   for (TrackExtendedVec::const_iterator middleIter = trackMiddleVec.begin();
        middleIter < middleEndIter;
        ++middleIter) {

     for (TrackExtendedVec::const_iterator outerIter = outerBeginIter;
          outerIter < outerEndIter;
          ++outerIter) {
       
       if ( *outerIter != *middleIter ) continue;

       for (TrackExtendedVec::const_iterator innerIter = innerBeginIter;
            innerIter < innerEndIter;
            ++innerIter) {

         // no need to check against middle, it is idendical to outer here
         if ( *outerIter == *innerIter ) {
           // return a null pointer
//           streamlog_out( DEBUG2 ) << " TestTriplet: Returns NULL " << std::endl ;
           return 0;            
         }

       }// for inner
     }// for outer    
   }// for middle
 }// if vectors not empty
  

 //    float dZ = FastTripletCheck(innerHit, middleHit, outerHit);

 //    if (fabs(dZ) > _minDistCutAttach)
 //      return trackAR;    


  
 double xh[3];
 double yh[3];
 float  zh[3];
 double wrh[3];
 float  wzh[3];
 float  rh[3];
 float  ph[3];

 float par[5];
 float epar[15];

 xh[0] = outerHit->getTrackerHit()->getPosition()[0];
 yh[0] = outerHit->getTrackerHit()->getPosition()[1];
 zh[0] = float(outerHit->getTrackerHit()->getPosition()[2]);
 wrh[0] = double(1.0/(outerHit->getResolutionRPhi()*outerHit->getResolutionRPhi()));
 wzh[0] = 1.0/(outerHit->getResolutionZ()*outerHit->getResolutionZ());

 xh[1] = middleHit->getTrackerHit()->getPosition()[0];
 yh[1] = middleHit->getTrackerHit()->getPosition()[1];
 zh[1] = float(middleHit->getTrackerHit()->getPosition()[2]);
 wrh[1] = double(1.0/(middleHit->getResolutionRPhi()*middleHit->getResolutionRPhi()));
 wzh[1] = 1.0/(middleHit->getResolutionZ()*middleHit->getResolutionZ());


 xh[2] = innerHit->getTrackerHit()->getPosition()[0];
 yh[2] = innerHit->getTrackerHit()->getPosition()[1];
 zh[2] = float(innerHit->getTrackerHit()->getPosition()[2]);
 wrh[2] = double(1.0/(innerHit->getResolutionRPhi()*innerHit->getResolutionRPhi()));
 wzh[2] = 1.0/(innerHit->getResolutionZ()*innerHit->getResolutionZ());

 for (int ih=0; ih<3; ih++) {
   rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
   ph[ih] = atan2(yh[ih],xh[ih]);
   if (ph[ih] < 0.) 
     ph[ih] = TWOPI + ph[ih]; 
 }

 int NPT = 3;
 int iopt = 2;
 float chi2RPhi;
 float chi2Z;

//  streamlog_out( DEBUG2 ) << " TestTriplet: Use fastHelixFit " << std::endl ;  
  
 _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
 par[3] = par[3]*par[0]/fabs(par[0]);

 float omega = par[0];
 float tanlambda = par[1];
 float phi0 = par[2];
 float d0 = par[3];
 float z0 = par[4];

 float Chi2 = chi2RPhi*_chi2WRPhiTriplet+chi2Z*_chi2WZTriplet;
 int ndf = 2*NPT-5;


 // Check if track satisfies all conditions


 //   std::cout << "Chi2/ndf = " << Chi2/float(ndf) << " , cut = " << _chi2FitCut << std::endl;
 //   std::cout << "d0 = " << d0 << " , cut = " << _cutOnD0  << std::endl;
 //   std::cout << "z0 = " << z0 << " , cut = " << _cutOnZ0  << std::endl;
 //   std::cout << "omega = " << omega << " , cut = " << _cutOnOmega << std::endl;

  if ( Chi2/float(ndf) > _chi2FitCut || fabs(d0) > _cutOnD0 || fabs(z0) > _cutOnZ0 || fabs(omega)>_cutOnOmega)
    // return a null pointer
    return 0;
  
//  bool failed = false;
//  
//  if ( Chi2/float(ndf) > _chi2FitCut ) {
//    streamlog_out(DEBUG1) << "Chi2/ndf = " << Chi2/float(ndf) << " , cut = " << _chi2FitCut << std::endl;
//    failed = true;
//  } else if (fabs(d0) > _cutOnD0 ) {
//    streamlog_out(DEBUG1) << "d0 = " << d0 << " , cut = " << _cutOnD0  << std::endl;
//    failed = true;
//  } else if (fabs(z0) > _cutOnZ0 ) {
//    streamlog_out(DEBUG1) << "z0 = " << z0 << " , cut = " << _cutOnZ0  << std::endl;
//    failed = true;
//  } else if ( fabs(omega)>_cutOnOmega)  {
//    streamlog_out(DEBUG1) << "omega = " << omega << " , cut = " << _cutOnOmega << std::endl;
//    failed = true;
//  }
//  
//  if( failed ) {
//    // return a null pointer
//    return 0;
//  }
//  
 helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);

 TrackExtended * trackAR = new TrackExtended();
 trackAR->addTrackerHitExtended(outerHit);
 trackAR->addTrackerHitExtended(middleHit);
 trackAR->addTrackerHitExtended(innerHit);
 outerHit->addTrackExtended(trackAR);
 middleHit->addTrackExtended(trackAR);
 innerHit->addTrackExtended(trackAR);    
 trackAR->setD0(d0);
 trackAR->setZ0(z0);
 trackAR->setPhi(phi0);
 trackAR->setTanLambda(tanlambda);
 trackAR->setOmega(omega);
 trackAR->setChi2( Chi2 );
 trackAR->setNDF( ndf );
 trackAR->setCovMatrix(epar);

//  streamlog_out(DEBUG1) << "Success !!!!!!!" << std::endl;


 return trackAR;

}

int SiliconTracking_MarlinTrk::BuildTrack(TrackerHitExtended * outerHit, 
                                         TrackerHitExtended * middleHit,
                                         TrackerHitExtended * innerHit,
                                         HelixClass & helix,
                                         int innerLayer,
                                         int iPhiLow, int iPhiUp,
                                         int iThetaLow, int iThetaUp, 
                                         TrackExtended * trackAR) {
 /**
  Method for building up track in the VXD. Method starts from the found triplet and performs
  sequential attachment of hits in other layers 
  */


 for (int layer = innerLayer-1; layer>=0; layer--) { // loop over remaining layers
   float distMin = 1.0e+20;
   TrackerHitExtended * assignedhit = NULL;

   for (int ipInner=iPhiLow; ipInner<iPhiUp+1;ipInner++) { // loop over phi in the Inner region

     for (int itInner=iThetaLow; itInner<iThetaUp+1;itInner++) { // loop over theta in the Inner region 

       int iPhiInner = ipInner;
       if (ipInner < 0) iPhiInner = _nDivisionsInPhi-1;
       if (ipInner >= _nDivisionsInPhi) iPhiInner = ipInner - _nDivisionsInPhi;
       int iCode = layer + _nLayers*iPhiInner +  _nLayers*_nDivisionsInPhi*itInner;
       TrackerHitExtendedVec& hitVecInner = _sectors[iCode];

       int nHitsInner = int(hitVecInner.size());

       for (int iInner=0;iInner<nHitsInner;iInner++) { // loop over hits in the Inner sector
         TrackerHitExtended * currentHit = hitVecInner[iInner];
         float pos[3];
         float distance[3];
         //      float ref[3];
         //      float point[6];
         for (int i=0; i<3; ++i) {
           pos[i] = float(currentHit->getTrackerHit()->getPosition()[i]);
           //      ref[i] = helix.getReferencePoint()[i];
         }
         //      float radius = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
         float time = helix.getDistanceToPoint(pos,distance);    
         //      float time = helix.getPointOnCircle(radius,ref,point);
         if (time < 1.0e+10) {
           //      float distRPhi1 = sqrt((point[0]-pos[0])*(point[0]-pos[0])+
           //                             (point[1]-pos[1])*(point[1]-pos[1]));
           //      float distRPhi2 = sqrt((point[3]-pos[0])*(point[3]-pos[0])+
           //                             (point[4]-pos[1])*(point[4]-pos[1]));
           //      float distZ1 = fabs(point[2]-pos[2]);
           //      float distZ2 = fabs(point[5]-pos[2]);           
           //      float dist1 = distRPhi1/_resolutionRPhi + distZ1/_resolutionZ;
           //      float dist2 = distRPhi2/_resolutionRPhi + distZ2/_resolutionZ;
           //      float dist = fmin(dist1,dist2);
           //      _distRPhi = fmin(distRPhi1,distRPhi2);
           //      _distZ = fmin(distZ1,distZ2);
           if (distance[2] < distMin) {
             distMin = distance[2];
             assignedhit = currentHit;
           }
         }
       } // endloop over hits in the Inner sector
     } // endloop over theta in the Inner region 
   } // endloop over phi in the Inner region
   if (distMin < _minDistCutAttach) {
     TrackerHitExtendedVec& hvec = trackAR->getTrackerHitExtendedVec();
     int  nHits = int(hvec.size());
     double * xh = new double[nHits+1];
     double * yh = new double[nHits+1];
     float * zh = new float[nHits+1];
     double * wrh = new double[nHits+1];
     float * wzh = new float[nHits+1];
     float * rh = new float[nHits+1];
     float * ph = new float[nHits+1];
     float par[5];
     float epar[15];

     for (int ih=0;ih<nHits;++ih) {
       TrackerHit * trkHit = hvec[ih]->getTrackerHit();
       xh[ih] = trkHit->getPosition()[0];
       yh[ih] = trkHit->getPosition()[1];
       zh[ih] = float(trkHit->getPosition()[2]);
       wrh[ih] = double(1.0/(hvec[ih]->getResolutionRPhi()*hvec[ih]->getResolutionRPhi()));
       wzh[ih] = 1.0/(hvec[ih]->getResolutionZ()*hvec[ih]->getResolutionZ());
       rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
       ph[ih] = float(atan2(yh[ih],xh[ih]));
       if (ph[ih] < 0.) 
         ph[ih] = TWOPI + ph[ih]; 
     }      
     TrackerHit * assignedTrkHit = assignedhit->getTrackerHit();
     xh[nHits] = assignedTrkHit->getPosition()[0];
     yh[nHits] = assignedTrkHit->getPosition()[1];
     zh[nHits] = float(assignedTrkHit->getPosition()[2]);
     rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));
     ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
     if (ph[nHits] < 0.) 
       ph[nHits] = TWOPI + ph[nHits]; 
     wrh[nHits] = double(1.0/(assignedhit->getResolutionRPhi()*assignedhit->getResolutionRPhi()));
     wzh[nHits] = 1.0/(assignedhit->getResolutionZ()*assignedhit->getResolutionZ());

     int NPT = nHits + 1;
     int iopt = 2;
     float chi2RPhi;
     float chi2Z;

     //std::cout << "######## number of hits to fit with _fastfitter = " << NPT << std::endl; 

     _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
     par[3] = par[3]*par[0]/fabs(par[0]);


     delete[] xh;
     delete[] yh;
     delete[] zh;
     delete[] wrh;
     delete[] wzh;
     delete[] rh;
     delete[] ph;

     bool validCombination = 0;
     float Chi2 = NAN;

     if ((nHits+1) == 4) {
       Chi2 = chi2RPhi*_chi2WRPhiQuartet+chi2Z*_chi2WZQuartet;
     }         
     if ((nHits+1) >= 5) {
       Chi2 = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
     }
     int ndf = 2*NPT-5;


     validCombination = Chi2/float(ndf) < _chi2FitCut;

     if ( validCombination ) {
       trackAR->addTrackerHitExtended(assignedhit);
       assignedhit->addTrackExtended(trackAR);
       float omega = par[0];
       float tanlambda = par[1];
       float phi0 = par[2];
       float d0 = par[3];
       float z0 = par[4];
       helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
       trackAR->setD0(d0);
       trackAR->setZ0(z0);
       trackAR->setPhi(phi0);
       trackAR->setTanLambda(tanlambda);
       trackAR->setOmega(omega);
       trackAR->setChi2( Chi2 );
       trackAR->setCovMatrix(epar);
       trackAR->setNDF( ndf );
     }

   }
 } // endloop over remaining layers

 TrackerHitExtendedVec& hvec = trackAR->getTrackerHitExtendedVec();  
 int nTotalHits = int(hvec.size());
 return nTotalHits;

}


void SiliconTracking_MarlinTrk::Sorting(TrackExtendedVec & trackVec) {
 /**
  Sorting of Track Vector in ascending order of chi2
  */

 std::sort(trackVec.begin(), trackVec.end(), compare_TrackExtended() );

 // also clean up? what does this do here?
 for (size_t i=0, sizeOfVector=trackVec.size(); i<sizeOfVector; ++i) {
   TrackerHitExtendedVec& hitVec = trackVec[i]->getTrackerHitExtendedVec();
   int nHits = int(hitVec.size());
   for (int ih=0;ih<nHits;ih++) {
     hitVec[ih]->clearTrackVec();
   }
 }

}

void SiliconTracking_MarlinTrk::CreateTrack(TrackExtended * trackAR ) {

 /**
  Method which creates Track out of TrackExtended objects. Checks for possible
  track splitting (separate track segments in VXD and FTD).
  */


 TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
 int nHits = int(hitVec.size());

 for (int i=0; i<nHits; ++i) {
   TrackExtendedVec& trackVec = hitVec[i]->getTrackExtendedVec();
   if (trackVec.size() != 0) 
     return ;
 }

 // First check if the current track is piece of the split one
 // look for matching track segment

 int found = 0;
 int nTrk = int(_trackImplVec.size());
 for (int i=0; i<nTrk; ++i) {
   TrackExtended * trackOld = _trackImplVec[i];
   TrackerHitExtendedVec& hitVecOld = trackOld->getTrackerHitExtendedVec();

   float phiNew = trackAR->getPhi();
   float phiOld = trackOld->getPhi();
   float thetaNew = M_PI_2 - atan(trackAR->getTanLambda());
   float thetaOld = M_PI_2 - atan(trackOld->getTanLambda());

   float angle = (cos(phiNew)*cos(phiOld)+sin(phiNew)*sin(phiOld))*sin(thetaNew)*sin(thetaOld)+cos(thetaNew)*cos(thetaOld);
   angle = acos(angle);

   if (angle < _angleCutForMerging) {
     int nHitsOld = int(hitVecOld.size());
     int nTotHits = nHits + nHitsOld;
     double * xh = new double[nTotHits];
     double * yh = new double[nTotHits];
     float * zh = new float[nTotHits];
     double * wrh = new double[nTotHits];
     float * wzh = new float[nTotHits];
     float * rh = new float[nTotHits];
     float * ph = new float[nTotHits];
     float par[5];
     float epar[15];
     float refPoint[3] = {0.,0.,0.};
     for (int ih=0;ih<nHits;++ih) {
       TrackerHit * trkHit = hitVec[ih]->getTrackerHit();
       float rR = hitVec[ih]->getResolutionRPhi();
       float rZ = hitVec[ih]->getResolutionZ();
       if (int(hitVec[ih]->getTrackExtendedVec().size()) != 0)
         streamlog_out(DEBUG2) << "WARNING : HIT POINTS TO TRACK " << std::endl;
       xh[ih] = trkHit->getPosition()[0];
       yh[ih] = trkHit->getPosition()[1];
       zh[ih] = float(trkHit->getPosition()[2]);
       wrh[ih] = double(1.0/(rR*rR));
       wzh[ih] = 1.0/(rZ*rZ);
       rh[ih] = float(sqrt(xh[ih]*xh[ih]+yh[ih]*yh[ih]));
       ph[ih] = float(atan2(yh[ih],xh[ih]));
     }      
     for (int ih=0;ih<nHitsOld;++ih) {
       TrackerHit * trkHit = hitVecOld[ih]->getTrackerHit();
       xh[ih+nHits] = trkHit->getPosition()[0];
       yh[ih+nHits] = trkHit->getPosition()[1];
       zh[ih+nHits] = float(trkHit->getPosition()[2]);
       float rR = hitVecOld[ih]->getResolutionRPhi();
       float rZ = hitVecOld[ih]->getResolutionZ();     
       wrh[ih+nHits] = double(1.0/(rR*rR));
       wzh[ih+nHits] = 1.0/(rZ*rZ);
       rh[ih+nHits] = float(sqrt(xh[ih+nHits]*xh[ih+nHits]+yh[ih+nHits]*yh[ih+nHits]));
       ph[ih+nHits] = float(atan2(yh[ih+nHits],xh[ih+nHits]));

     }
     int NPT = nTotHits;
     int iopt = 3;
     float chi2RPhi;
     float chi2Z;
     int ndf = 2*NPT - 5;

     _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
     par[3] = par[3]*par[0]/fabs(par[0]);

     float omega = par[0];
     float tanlambda = par[1];
     float phi0 = par[2];
     float d0 = par[3];
     float z0 = par[4];

     float eparmin[15];
     for (int iparam=0;iparam<15;++iparam)
       eparmin[iparam] = epar[iparam];      

     float refPointMin[3];
     for (int ipp=0;ipp<3;++ipp)
       refPointMin[ipp] = refPoint[ipp];

     float chi2Min = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
     chi2Min = chi2Min/float(ndf);

     float chi2MinRPhi = chi2RPhi;
     float chi2MinZ = chi2Z;

     int iBad = -1;
     if (chi2Min < _chi2FitCut) {
       found = 1;
     }
     else { // SJA:FIXME: UH What is going on here? setting weights to 0 and refitting?
       float * wzhOld = new float[nTotHits];
       double * wrhOld = new double[nTotHits];
       for (int i=0;i<nTotHits;++i) {
         wzhOld[i] = wzh[i];
         wrhOld[i] = wrh[i];
       }
       for (int i=0; i<nTotHits; ++i) {
         for (int j=0;j<nTotHits;++j) {
           if (i == j) {
             wrh[j] = 0.0;
             wzh[j] = 0.0;
           } 
           else {
             wrh[j] = wrhOld[j];
             wzh[j] = wzhOld[j];
           }
         }

         _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
         par[3] = par[3]*par[0]/fabs(par[0]);

         float chi2Cur = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
         chi2Cur = chi2Cur/float(ndf);

         if (chi2Cur < chi2Min) {
           chi2Min = chi2Cur;
           chi2MinRPhi = chi2RPhi;
           chi2MinZ = chi2Z;
           omega = par[0];
           tanlambda = par[1];
           phi0 = par[2];
           d0 = par[3];
           z0 = par[4];
           for (int iparam=0;iparam<15;++iparam)
             eparmin[iparam] = epar[iparam];
           for (int ipp=0;ipp<3;++ipp)
             refPointMin[ipp] = refPoint[ipp];
           iBad = i;
         }
       }
       if (chi2Min < _chi2FitCut) {
         found = 1;
       }
       delete[] wzhOld;
       delete[] wrhOld;
     }

     // Split track is found.
     // Attach hits belonging to the current track segment to  
     // the track already created
     if (found == 1) {
       trackOld->ClearTrackerHitExtendedVec();
       for (int i=0;i<nHits;++i) {
         TrackerHitExtended * trkHit = hitVec[i];
         trkHit->clearTrackVec();
         if (i == iBad) {          
         }
         else {
           trackOld->addTrackerHitExtended(trkHit);
           trkHit->addTrackExtended( trackOld );
         }
       }  
       for (int i=0;i<nHitsOld;++i) {
         int icur = i+nHits;
         TrackerHitExtended * trkHit = hitVecOld[i];
         trkHit->clearTrackVec();
         if (icur == iBad) {
         }
         else {
           trackOld->addTrackerHitExtended(trkHit);
           trkHit->addTrackExtended( trackOld );
         }
       }
       trackOld->setOmega(omega);
       trackOld->setTanLambda(tanlambda);
       trackOld->setPhi(phi0);
       trackOld->setD0(d0);
       trackOld->setZ0(z0);

       //      std::cout << "Split track found " << d0 << " " << z0 << std::endl;
	
       // killeb:  In the original SiliconTracking this was in the NOT simple helix branch.
	// The rest of the code uses the simple helix branch, where ndf_D is never set.
	// In fact it has never been initialised or used anywhere. I think this line should not be executed.
	// ndf = ndf_D;

       trackOld->setChi2(chi2Min*float(ndf));  
       trackOld->setNDF(ndf);
       trackOld->setCovMatrix(eparmin);
       //      trackOld->setReferencePoint(refPointMin);
     }

     delete[] xh;
     delete[] yh;
     delete[] zh;
     delete[] wrh;
     delete[] wzh;
     delete[] rh;
     delete[] ph;

   }
   if (found == 1)
     break;
 }

 // Candidate is a unique track
 // No other segments are found
 if (found == 0 ) {
   _trackImplVec.push_back(trackAR);
   for (int i=0;i<nHits;++i) {
     TrackerHitExtended * hit = hitVec[i];
     hit->addTrackExtended( trackAR );
   }
 }


}

void SiliconTracking_MarlinTrk::AttachRemainingVTXHitsFast() {

 std::vector<TrackerHitExtendedVec> nonAttachedHits;
 nonAttachedHits.resize(_nDivisionsInPhi*_nDivisionsInTheta);
 std::vector<TrackExtendedVec> trackVector;
 trackVector.resize(_nDivisionsInPhi*_nDivisionsInTheta);
 int nTrk = int(_trackImplVec.size());

 for (int iTrk=0;iTrk<nTrk;++iTrk) {
   TrackExtended * track = _trackImplVec[iTrk];
   double Phi = double(track->getPhi());
   if (Phi < 0)
     Phi = Phi + TWOPI;
   float tanlambda = track->getTanLambda();
   double cosTheta = double(tanlambda/sqrt(1+tanlambda*tanlambda));
   int iPhi = int(Phi/_dPhi);
   int iTheta = int ((cosTheta + double(1.0))/_dTheta);
   int iCode = iPhi + _nDivisionsInPhi*iTheta; 
   trackVector[iCode].push_back( track );
 }

 for (int il=0;il<_nLayers;++il) {
   for (int ip=0;ip<_nDivisionsInPhi;++ip) {
     for (int it=0;it<_nDivisionsInTheta; ++it) {
       int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
       TrackerHitExtendedVec& hitVec = _sectors[iCode];
       int nH = int(hitVec.size());
       for (int iH=0; iH<nH; ++iH) {
         TrackerHitExtended * hitExt = hitVec[iH];
         TrackExtendedVec& trackVec = hitExt->getTrackExtendedVec();
         if (trackVec.size()==0) {
           TrackerHit * hit = hitExt->getTrackerHit();
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
           int iPhi = int(Phi/_dPhi);
           int iTheta = int ((cosTheta + double(1.0))/_dTheta);
           int iCode = iPhi + _nDivisionsInPhi*iTheta;      
           nonAttachedHits[iCode].push_back( hitExt );
         }
       }
     }
   }
 }

 for (int iT=0; iT<_nDivisionsInTheta; ++iT) {
   for (int iP=0; iP<_nDivisionsInPhi; ++iP) {
     int iCode = iP + _nDivisionsInPhi*iT; 
     int nHits = int(nonAttachedHits[iCode].size());
     int iT1 = iT - 1;
     int iT2 = iT + 1; 
     if (iT == 0) {
       iT1 = iT;
       iT2 = iT1 + 1;
     }
     if (iT == _nDivisionsInTheta - 1) {
       iT2 = iT;
       iT1 = iT2 - 1;
     }
     int iPHI[3];
     iPHI[0] = iP - 1;
     iPHI[1] = iP;
     iPHI[2] = iP + 1;
     if (iP == 0) 
       iPHI[0] = _nDivisionsInPhi - 1;
     if (iP == _nDivisionsInPhi - 1 )
       iPHI[2] = 0;

     for (int ihit = 0; ihit<nHits; ++ihit) {

       TrackerHitExtended * hit = nonAttachedHits[iCode][ihit];
       TrackExtended * trackToAttach = NULL;
       float minDist = 1.0e+6;

       for (int iTheta = iT1; iTheta <iT2+1; ++iTheta) {
         for (int indexP=0;indexP<3;++indexP) {
           int iPhi = iPHI[indexP];        
           int iCodeForTrack = iPhi + _nDivisionsInPhi*iTheta;
           int nTrk = int(trackVector[iCodeForTrack].size());
           for (int iTrk=0; iTrk<nTrk; ++iTrk) {         
             TrackExtended * trackAR = trackVector[iCodeForTrack][iTrk];
             bool consider = true;
             if (_checkForDelta) {
               TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
               int NHITS = int(hitVector.size());
               for (int IHIT=0;IHIT<NHITS;++IHIT) {

                 // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
                 TrackerHit* trkhit1 = hit->getTrackerHit();
                 TrackerHit* trkhit2 = hitVector[IHIT]->getTrackerHit();                  

                 if ( trkhit1->getCellID0() == trkhit2->getCellID0() ){ // i.e. they are in the same sensor
                   float distance = 0.;
                   for (int iC=0;iC<3;++iC) {
                     float posFirst = float(hit->getTrackerHit()->getPosition()[iC]);
                     float posSecond = float(hitVector[IHIT]->getTrackerHit()->getPosition()[iC]);
                     float deltaPos = posFirst - posSecond;
                     distance += deltaPos*deltaPos;
                   }
                   distance = sqrt(distance);
                   if (distance<_minDistToDelta) {
                     consider = false;
                     break;
                   }
                 }
               }
             }
             if (consider) {   
               float phi0 = trackAR->getPhi();
               float d0 = trackAR->getD0();
               float z0 = trackAR->getZ0();
               float omega = trackAR->getOmega();
               float tanlambda = trackAR->getTanLambda();
               HelixClass helix;
               helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
               int layer = getLayerID(hit->getTrackerHit());
               if (layer > _minimalLayerToAttach) {
                 float pos[3];
                 for (int i=0; i<3; ++i) 
                   pos[i] = hit->getTrackerHit()->getPosition()[i];      
                 float distance[3];
                 float time = helix.getDistanceToPoint(pos,distance);
                 if (time < 1.0e+10) {
                   if (distance[2] < minDist) {
                     minDist = distance[2];
                     trackToAttach = trackAR;
                   }                      
                 }    
               }
             }
           }
         }
       }
       if (minDist < _minDistCutAttach && trackToAttach != NULL) {
         int iopt = 3;
         AttachHitToTrack(trackToAttach,hit,iopt);
       }      
     }
   }
 }
}

void SiliconTracking_MarlinTrk::AttachRemainingVTXHitsSlow() {
 TrackerHitExtendedVec nonAttachedHits;
 nonAttachedHits.clear();

 for (int il=0;il<_nLayers;++il) {
   for (int ip=0;ip<_nDivisionsInPhi;++ip) {
     for (int it=0;it<_nDivisionsInTheta; ++it) {
       int iCode = il + _nLayers*ip + _nLayers*_nDivisionsInPhi*it;      
       TrackerHitExtendedVec& hitVec = _sectors[iCode];
       int nH = int(hitVec.size());
       for (int iH=0; iH<nH; ++iH) {
         TrackerHitExtended * hit = hitVec[iH];
         TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
         // if (trackVec.size()==0)
         // nonAttachedHits.push_back( hit );
         //-- allow hits that are only used in triplets to be re-attached 
         unsigned int maxTrackSize = 0;
         for(unsigned int itrack = 0;itrack < trackVec.size();itrack++){
           TrackerHitExtendedVec hitVec = trackVec[itrack]->getTrackerHitExtendedVec();
           unsigned int isize = hitVec.size();
           if(isize>maxTrackSize)maxTrackSize = isize;
         }     
         if (maxTrackSize<=3)nonAttachedHits.push_back( hit );


       }
     }
   }
 }

 int nNotAttached = int(nonAttachedHits.size());

 int nTrk = int(_trackImplVec.size()); 
 for (int iHit=0; iHit<nNotAttached; ++iHit) {
   TrackerHitExtended * hit = nonAttachedHits[iHit];
   int layer = getLayerID( hit->getTrackerHit() );
   if (layer > _minimalLayerToAttach) {
     float pos[3];
     for (int i=0; i<3; ++i) 
       pos[i] = hit->getTrackerHit()->getPosition()[i];      
     float minDist = 1e+10;
     TrackExtended * trackToAttach = NULL;
     for (int iTrk=0; iTrk<nTrk; ++iTrk) {
       TrackExtended * trackAR = _trackImplVec[iTrk];
       bool consider = true;
       if (_checkForDelta) {
         TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
         int NHITS = int(hitVector.size());
         for (int IHIT=0;IHIT<NHITS;++IHIT) {

           // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
           TrackerHit* trkhit1 = hit->getTrackerHit();
           TrackerHit* trkhit2 = hitVector[IHIT]->getTrackerHit();                  

           if ( trkhit1->getCellID0() == trkhit2->getCellID0() ){ // i.e. they are in the same sensor

             float distance = 0.;
             for (int iC=0;iC<3;++iC) {
               float posFirst = float(hit->getTrackerHit()->getPosition()[iC]);
               float posSecond = float(hitVector[IHIT]->getTrackerHit()->getPosition()[iC]);
               float deltaPos = posFirst - posSecond;
               distance += deltaPos*deltaPos;
             }
             distance = sqrt(distance);
             if (distance<_minDistToDelta) {
               consider = false;
               break;
             }
           }       
         }
       }
       if (consider) {
         HelixClass helix;
         float phi0 = trackAR->getPhi();
         float d0 = trackAR->getD0();
         float z0 = trackAR->getZ0();
         float omega = trackAR->getOmega();
         float tanlambda = trackAR->getTanLambda();
         helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
         float distance[3];
         float time = helix.getDistanceToPoint(pos,distance);
         if (time < 1.0e+10) {
           if (distance[2] < minDist) {
             minDist = distance[2];
             trackToAttach = trackAR;
           }
         }
       }
     }
     if (minDist < _minDistCutAttach && trackToAttach != NULL) {
       int iopt = 3;
       AttachHitToTrack(trackToAttach,hit,iopt);
     }      
   }
 }  
}

void SiliconTracking_MarlinTrk::AttachRemainingFTDHitsSlow() {
 TrackerHitExtendedVec nonAttachedHits;
 nonAttachedHits.clear();

 for (int iS=0;iS<2;++iS) {
   for (unsigned int layer=0;layer<_nlayersFTD;++layer) {
     for (int ip=0;ip<_nPhiFTD;++ip) {
       int iCode = iS + 2*layer + 2*_nlayersFTD*ip;      
       TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
       int nH = int(hitVec.size());
       for (int iH=0; iH<nH; ++iH) {
         TrackerHitExtended * hit = hitVec[iH];
         TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
         if (trackVec.size()==0)
           nonAttachedHits.push_back( hit );
       }
     }
   }
 }

 int nNotAttached = int(nonAttachedHits.size());

 int nTrk = int(_trackImplVec.size()); 
 for (int iHit=0; iHit<nNotAttached; ++iHit) {
   TrackerHitExtended * hit = nonAttachedHits[iHit];
   float pos[3];
   for (int i=0; i<3; ++i) 
     pos[i] = hit->getTrackerHit()->getPosition()[i];      
   float minDist = 1e+10;
   TrackExtended * trackToAttach = NULL;
   for (int iTrk=0; iTrk<nTrk; ++iTrk) {
     TrackExtended * trackAR = _trackImplVec[iTrk];
     bool consider = true;
     TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
     int NHITS = int(hitVector.size());

     for (int IHIT=0;IHIT<NHITS;++IHIT) {

       // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works 
       //        if (hit->getTrackerHit()->getType() == hitVector[IHIT]->getTrackerHit()->getType()) {
       if (hit->getTrackerHit()->getCellID0() == hitVector[IHIT]->getTrackerHit()->getCellID0()) {

         consider = false;
         break;
       }
     }


     if (consider) {
       HelixClass helix;
       float phi0 = trackAR->getPhi();
       float d0 = trackAR->getD0();
       float z0 = trackAR->getZ0();
       float omega = trackAR->getOmega();
       float tanlambda = trackAR->getTanLambda();
       if (tanlambda*float(getSideID(hit->getTrackerHit())) > 0) {
         helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
         float distance[3];
         float time = helix.getDistanceToPoint(pos,distance);
         if (time < 1.0e+10) {
           if (distance[2] < minDist) {
             minDist = distance[2];
             trackToAttach = trackAR;
           }
         }
       }
     }
   }
   if (minDist < _minDistCutAttach && trackToAttach != NULL) {
     int iopt = 2;
     AttachHitToTrack(trackToAttach,hit,iopt);
   }      
 }  
}


void SiliconTracking_MarlinTrk::AttachRemainingFTDHitsFast() {
 int nTrk = _trackImplVec.size();

 for (int iTrk=0; iTrk<nTrk; ++iTrk) {
   TrackExtended * trackAR = _trackImplVec[iTrk];
   HelixClass helix;
   float phi0 = trackAR->getPhi();
   float d0 = trackAR->getD0();
   float z0 = trackAR->getZ0();
   float omega = trackAR->getOmega();
   float tanlambda = trackAR->getTanLambda();
   helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
   int iSemiSphere = 0;
   if (tanlambda > 0) 
     iSemiSphere = 1;
   float ref[3];
   for (int i=0;i<3;++i) 
     ref[i] = helix.getReferencePoint()[i];
   // Start loop over FTD layers
   for (unsigned int layer=0;layer<_nlayersFTD;layer++) {
     float ZL = _zLayerFTD[layer];
     if (iSemiSphere == 0)
       ZL = - ZL;
     float point[3];
     helix.getPointInZ(ZL,ref,point);
     float Phi = atan2(point[1],point[0]);
     if (Phi < 0) 
       Phi = Phi + TWOPI;
     int iPhi = int(Phi/_dPhiFTD);
     float distMin = 1e+10;
     TrackerHitExtended * attachedHit = NULL;     
     for (int iP=iPhi-1;iP<=iPhi+1;++iP) {
       int iPP = iP;
       if (iP < 0) 
         iPP = iP + _nPhiFTD;
       if (iP >= _nPhiFTD)
         iPP = iP - _nPhiFTD;  
       int iCode = iSemiSphere + 2*layer + 2*_nlayersFTD*iPP;
       int nHits = int(_sectorsFTD[iCode].size());
       for (int iHit=0;iHit<nHits;++iHit) {
         TrackerHitExtended * hit = _sectorsFTD[iCode][iHit];
         bool consider = true;
         TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
         int NHITS = int(hitVector.size());

         // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works 
         for (int IHIT=0;IHIT<NHITS;++IHIT) {
           //            if (hit->getTrackerHit()->getType() == hitVector[IHIT]->getTrackerHit()->getType()) {
           if (hit->getTrackerHit()->getCellID0() == hitVector[IHIT]->getTrackerHit()->getCellID0()) {
             consider = false;
             break;
           }
         }


         if (consider) {
           float pos[3];
           for (int i=0;i<3;++i) {
             pos[i] = hit->getTrackerHit()->getPosition()[i];
           }
           float distance[3];
           float time = helix.getDistanceToPoint(pos,distance);
           if (time < 1.0e+10) {
             if (distance[2] < distMin) {
               distMin = distance[2];
               attachedHit = hit;
             }
           }
         }
       }
     }
     if (distMin < _minDistCutAttach && attachedHit != NULL) {
       int iopt = 2;
       AttachHitToTrack(trackAR,attachedHit, iopt);
     }
   }
 }
}

void SiliconTracking_MarlinTrk::TrackingInFTD() {
 int nComb = int(_CombinationsFTD.size()) / 3;
 for (int iComb=0;iComb<nComb;++iComb) {
   int nLS[3];
   nLS[0] = _CombinationsFTD[3*iComb];
   nLS[1] = _CombinationsFTD[3*iComb+1];
   nLS[2] = _CombinationsFTD[3*iComb+2];
   for (int iS=0;iS<2;++iS) {
     //      std::cout << "Combinations : " << iS << " " << nLS[0] << " " << nLS[1] << " " << nLS[2] << std::endl;
     //      int iC = iS + 2*nLS[0];
     //      TrackerHitExtendedVec& hitVec = _sectorsFTD[iC];
     //      int nO = int(hitVec.size());
     //      iC = iS + 2*nLS[1];
     //      hitVec = _sectorsFTD[iC];
     //      int nM = int(hitVec.size());
     //      iC = iS + 2*nLS[2];
     //      hitVec = _sectorsFTD[iC];
     //      int nI = int(hitVec.size());
     //      std::cout << nO << " " << nM << " " << nI << std::endl;
     for (int ipOuter=0;ipOuter<_nPhiFTD;++ipOuter) {
       int ipMiddleLow = ipOuter - 1;
       int ipMiddleUp  = ipOuter + 1;

       unsigned int iCodeOuter = iS + 2*nLS[0] + 2*_nlayersFTD*ipOuter;

       if( iCodeOuter >= _sectorsFTD.size()){          
         std::cout << "iCodeOuter index out of range: iCodeOuter =   " << iCodeOuter << " _sectorsFTD.size() = " << _sectorsFTD.size() << " exit(1) called from file " << __FILE__ << " line " << __LINE__<< std::endl;
         exit(1);
       }

       TrackerHitExtendedVec& hitVecOuter = _sectorsFTD[iCodeOuter];

       int nOuter = int(hitVecOuter.size());
       for (int iOuter=0;iOuter<nOuter;++iOuter) {

         TrackerHitExtended * hitOuter = hitVecOuter[iOuter];
         for (int ipMiddle=ipMiddleLow;ipMiddle<=ipMiddleUp;++ipMiddle) {
           //for(int ipMiddle=0;ipMiddle<_nPhiFTD;++ipMiddle) {
           int ipM = ipMiddle;
           if (ipM < 0) 
             ipM = ipMiddle + _nPhiFTD;
           if (ipM >= _nPhiFTD) 
             ipM = ipMiddle - _nPhiFTD;
           int iCodeMiddle = iS + 2*nLS[1] + 2*_nlayersFTD*ipM;
           TrackerHitExtendedVec& hitVecMiddle = _sectorsFTD[iCodeMiddle];
           int ipInnerLow,ipInnerUp;       
           ipInnerLow = ipMiddle - 1;
           ipInnerUp =  ipMiddle + 1;
           int nMiddle = int(hitVecMiddle.size());
           for (int iMiddle=0;iMiddle<nMiddle;++iMiddle) {
             TrackerHitExtended * hitMiddle = hitVecMiddle[iMiddle];
             for (int ipInner=ipInnerLow;ipInner<=ipInnerUp;++ipInner) {
               //for (int ipInner=0;ipInner<_nPhiFTD;++ipInner) {
               int ipI = ipInner;
               if (ipI < 0)
                 ipI = ipInner + _nPhiFTD;
               if (ipI >= _nPhiFTD) 
                 ipI = ipInner - _nPhiFTD;
               int iCodeInner = iS + 2*nLS[2] + 2*_nlayersFTD*ipI;
               TrackerHitExtendedVec& hitVecInner = _sectorsFTD[iCodeInner];
               int nInner = int(hitVecInner.size());
               for (int iInner=0;iInner<nInner;++iInner) {
                 TrackerHitExtended * hitInner = hitVecInner[iInner];
                 HelixClass helix;
                 //                  std::cout << std::endl;
                 //                  std::cout << "Outer Hit Type " << hitOuter->getTrackerHit()->getType() << " z = " << hitOuter->getTrackerHit()->getPosition()[2] 
                 //                  << "\nMiddle Hit Type "<< hitMiddle->getTrackerHit()->getType() << " z = " << hitMiddle->getTrackerHit()->getPosition()[2]  
                 //                  << "\nInner Hit Type "<< hitInner->getTrackerHit()->getType() << " z = " << hitInner->getTrackerHit()->getPosition()[2]  << std::endl;
                 TrackExtended * trackAR = TestTriplet(hitOuter,hitMiddle,hitInner,helix);
                 if (trackAR != NULL) {
                   //                    std::cout << "FTD triplet found" << std::endl;
                   int nHits = BuildTrackFTD(trackAR,nLS,iS);

		    _tracksWithNHitsContainer.getTracksWithNHitsVec( nHits ).push_back( trackAR );
                 }
               }
             }
           }
         }       
       }
     }
   }
 }
}


int SiliconTracking_MarlinTrk::BuildTrackFTD(TrackExtended * trackAR, int * nLR, int iS) {
 //  std::cout << "BuildTrackFTD: Layers = " << nLR[0] << " " << nLR[1] << " " << nLR[2] << std::endl;

 // initialise a helix from the track
 HelixClass helix;
 const float d0 = trackAR->getD0();
 const float z0 = trackAR->getZ0();
 const float phi0 = trackAR->getPhi();
 const float tanlambda = trackAR->getTanLambda();
 const float omega = trackAR->getOmega();
 helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
 float ref[3] = {helix.getReferencePoint()[0],
		  helix.getReferencePoint()[1],
		  helix.getReferencePoint()[2]};

 for (int iL=0; iL < static_cast<int>(_nlayersFTD); ++iL) {
   if (iL != nLR[0] && iL != nLR[1] && iL != nLR[2]) {
     float point[3];
     float ZL = _zLayerFTD[iL];
     if (iS == 0) 
       ZL = - ZL;
     helix.getPointInZ(ZL,ref,point);
     //      float Phi = atan2(point[1],point[0]);
     //      int iPhi = int(Phi/_dPhiFTD);
     float distMin = 1e+6;
     TrackerHitExtended * attachedHit = NULL;
     for (int ip=0;ip<=_nPhiFTD;++ip) {
       int iP = ip;
       if (iP < 0)
         iP = ip + _nPhiFTD;
       if (iP >= _nPhiFTD)
         iP = ip - _nPhiFTD;   
       int iCode = iS + 2*iL + 2*_nlayersFTD*iP;
       TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
       int nH = int(hitVec.size());
       for (int iH=0; iH<nH; ++iH) {
         TrackerHitExtended * hit = hitVec[iH];
         TrackerHit * trkHit = hit->getTrackerHit();
         float pos[3];
         for (int i=0;i<3;++i)
           pos[i] = float(trkHit->getPosition()[i]);
         float distance[3];
         float time = helix.getDistanceToPoint(pos,distance);
         if (time < 1.0e+10) {
           if (distance[2] < distMin) {
             distMin = distance[2];
             attachedHit = hit;
           }
         }
       }
     }
     //      std::cout << "Layer = " << iL << "  distMin = " << distMin << std::endl;
     if (distMin < _minDistCutAttach && attachedHit != NULL) {
       int iopt = 2;
       AttachHitToTrack( trackAR, attachedHit, iopt);
     }
   }
 }
 TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
 int nH = int (hitVec.size());
 return nH;
}

int SiliconTracking_MarlinTrk::AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit, int iopt) {

 int attached = 0;
 TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
 int nHits = int(hitVec.size());

 double * xh = new double[nHits+1];
 double * yh = new double[nHits+1];
 float  * zh = new float[nHits+1];
 double * wrh = new double[nHits+1];
 float * wzh = new float[nHits+1];
 float * rh = new float[nHits+1];
 float * ph = new float[nHits+1];
 float par[5];
 float epar[15];

 for (int i=0; i<nHits; ++i) {
   TrackerHit * trkHit = hitVec[i]->getTrackerHit();
   xh[i] = double(trkHit->getPosition()[0]);
   yh[i] = double(trkHit->getPosition()[1]);
   zh[i] = float(trkHit->getPosition()[2]);
   ph[i] = float(atan2(yh[i],xh[i]));
   rh[i] = float(sqrt(xh[i]*xh[i]+yh[i]*yh[i]));
   float rR = hitVec[i]->getResolutionRPhi();
   float rZ = hitVec[i]->getResolutionZ();
   wrh[i] = double(1.0/(rR*rR));
   wzh[i] = 1.0/(rZ*rZ);
 }

 TrackerHit * trkHit = hit->getTrackerHit();
 xh[nHits] = double(trkHit->getPosition()[0]);
 yh[nHits] = double(trkHit->getPosition()[1]);
 zh[nHits] = float(trkHit->getPosition()[2]);
 ph[nHits] = float(atan2(yh[nHits],xh[nHits]));
 rh[nHits] = float(sqrt(xh[nHits]*xh[nHits]+yh[nHits]*yh[nHits]));

 float rR = hit->getResolutionRPhi();
 float rZ = hit->getResolutionZ();
 wrh[nHits] = double(1.0/(rR*rR));
 wzh[nHits] = 1.0/(rZ*rZ);


 int NPT = nHits + 1;

 // SJA:FIXME the newtonian part is giving crazy results for FTD so just use iopt 2 for simply attaching hits 
 // using SIT and VXD doesn't seem to give any problems, so make it a function parameter and let the caller decide
 //  int iopt = 3;

 float chi2RPhi;
 float chi2Z;


 _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh,iopt, par, epar, chi2RPhi, chi2Z);
 par[3] = par[3]*par[0]/fabs(par[0]);


 float omega = par[0];
 float tanlambda = par[1];
 float phi0 = par[2];
 float d0 = par[3];
 float z0 = par[4];
 float chi2 = NAN;
 int ndf = INT_MAX;

 if (NPT == 3) {
   chi2 = chi2RPhi*_chi2WRPhiTriplet+chi2Z*_chi2WZTriplet;
 }
 if (NPT == 4) {
   chi2 = chi2RPhi*_chi2WRPhiQuartet+chi2Z*_chi2WZQuartet;
 }
 if (NPT > 4) {
   chi2 = chi2RPhi*_chi2WRPhiSeptet+chi2Z*_chi2WZSeptet;
 }
 ndf = 2*NPT-5;

 if (chi2/float(ndf) < _chi2FitCut) {
   trackAR->addTrackerHitExtended(hit);
   hit->addTrackExtended( trackAR );
   trackAR->setChi2( chi2 );
   trackAR->setOmega( omega );
   trackAR->setTanLambda( tanlambda );
   trackAR->setD0( d0 );
   trackAR->setZ0( z0 );
   trackAR->setPhi( phi0 );
   trackAR->setNDF( ndf );
   trackAR->setCovMatrix( epar );
   attached = 1;
 }     

 delete[] xh;
 delete[] yh;
 delete[] zh;
 delete[] wrh;
 delete[] wzh;
 delete[] rh;
 delete[] ph;

 return attached;


}

void SiliconTracking_MarlinTrk::FinalRefit(LCCollectionVec* trk_col, LCCollectionVec* rel_col) {

 int nTracks = int(_trackImplVec.size());

 int nSiSegments = 0;        
 float eTot = 0.;
 float pxTot = 0.;
 float pyTot = 0.;
 float pzTot = 0.;

 for (int iTrk=0;iTrk<nTracks;++iTrk) {

   TrackExtended * trackAR = _trackImplVec[iTrk];    
   TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();

   int nHits = int(hitVec.size());

   if( nHits >= _minimalHits) {
     //    int * lh = new int[nHits];
     std::vector<int> lh;
     lh.resize(nHits);

     for (int i=0; i<nHits; ++i) {
       lh[i]=0;
     }

     float d0 = trackAR->getD0();
     float z0 = trackAR->getZ0();
     float omega = trackAR->getOmega();
     float tanlambda = trackAR->getTanLambda();
     float phi0 = trackAR->getPhi();

     HelixClass * helix = new HelixClass();
     helix->Initialize_Canonical(phi0, d0, z0, omega, 
                                 tanlambda, _bField);


     // get the point of closest approach to the reference point
     // here it is implicitly assumed that the reference point is the origin 
     float Pos[3];
     Pos[0] = -d0*sin(phi0);
     Pos[1] = d0*cos(phi0);
     Pos[2] = z0;


     // at this point is is possible to have hits from the same layer ...
     // so a check is made to ensure that the hit with the smallest distance to the 
     // current helix hypothosis is used, the other hit has lh set to 0 

     // start loop over the hits to
     for (int ihit=0;ihit<nHits;++ihit) {

       lh[ihit] = 1; // only hits which have lh=1 will be used for the fit

       // get the pointer to the lcio trackerhit for this hit
       TrackerHit * trkHit = hitVec[ihit]->getTrackerHit();

       int det = getDetectorID(trkHit);

       if (det == lcio::ILDDetID::VXD || det == lcio::ILDDetID::FTD || det == lcio::ILDDetID::SIT) { // only accept VXD, FTD or SIT


         //        int layer = getLayerID(trkHit);
         //        int moduleIndex = getModuleID(trkHit);

         // start a double loop over the hits which have already been checked 
         for (int lhit=0;lhit<ihit;++lhit) {

           // get the pointer to the lcio trackerhit for the previously checked hit
           TrackerHit * trkHitS = hitVec[lhit]->getTrackerHit();


           //          int layerS = getLayerID(trkHitS);
           //          int moduleIndexS = getModuleID(trkHitS);

           // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works 
           // if they are on the same layer and the previously checked hits has been declared good for fitting
           //          if ((trkHitS->getType() == trkHit->getType()) && (lh[lhit] == 1)) {
           // check if the hits have the same layer and petal number
           //          hitVec[ihit]->
           //          if ((layer == layerS) && (moduleIndex==moduleIndexS) && (lh[lhit] == 1)) {
           if ( (trkHit->getCellID0() == trkHitS->getCellID0()) && (lh[lhit] == 1)) {

             // get the position of the hits 
             float xP[3];
             float xPS[3];
             for (int iC=0;iC<3;++iC) {
               xP[iC] = float(trkHit->getPosition()[iC]);
               xPS[iC] = float(trkHitS->getPosition()[iC]);
             }

             // get the intersection of the helix with the either the cylinder or plane containing the hit
             float Point[3];
             float PointS[3];
             if (det == lcio::ILDDetID::FTD) {
               float time = helix->getPointInZ(xP[2],Pos,Point);
               time = helix->getPointInZ(xPS[2],Pos,PointS);
             }
             else {
               float RAD = sqrt(xP[0]*xP[0]+xP[1]*xP[1]);
               float RADS = sqrt(xPS[0]*xPS[0]+xPS[1]*xPS[1]);
               float time = helix->getPointOnCircle(RAD,Pos,Point);
               time = helix->getPointOnCircle(RADS,Pos,PointS);
             }

             float DIST = 0;
             float DISTS = 0;

             // get the euclidean distance between the hit and the point of intersection
             for (int iC=0;iC<3;++iC) {
               DIST += (Point[iC]-xP[iC])*(Point[iC]-xP[iC]);
               DISTS += (PointS[iC]-xPS[iC])*(PointS[iC]-xPS[iC]);
             }
             if (DIST < DISTS) {
               lh[lhit] = 0;
             }
             else {
               lh[ihit] = 0;
             }
             break;
           }
         }
       }
     }

     delete helix;

     EVENT::TrackerHitVec trkHits;
     EVENT::TrackerHitVec trkHits_used_inFit;

     int nFit = 0;
     for (int i=0; i<nHits; ++i) {
       // check if the hit has been rejected as being on the same layer and further from the helix lh==0
       if (lh[i] == 1) {
         TrackerHit * trkHit = hitVec[i]->getTrackerHit();
         nFit++;
         if(trkHit) { 
           trkHits.push_back(trkHit);   
         }
         else{
           throw EVENT::Exception( std::string("SiliconTracking_MarlinTrk::FinalRefit: TrackerHit pointer == NULL ")  ) ;
         }
       }
       else { // reject hit 
              // SJA:FIXME missuse of type find a better way to signal rejected hits
         hitVec[i]->setType(int(0));
       }
     }



     if( trkHits.size() < 3 ) {
       streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Cannot fit less than 3 hits. Number of hits =  " << trkHits.size() << std::endl;
       continue ; 
     }

     TrackImpl* Track = new TrackImpl ;
     
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
     r2_values.reserve(trkHits.size());
     
     for (TrackerHitVec::iterator it=trkHits.begin(); it!=trkHits.end(); ++it) {
       EVENT::TrackerHit* h = *it;
       float r2 = h->getPosition()[0]*h->getPosition()[0]+h->getPosition()[1]*h->getPosition()[1];
       r2_values.push_back(std::make_pair(r2, *it));
     }
     
     sort(r2_values.begin(),r2_values.end());
     
     trkHits.clear();
     trkHits.reserve(r2_values.size());
     
     for (std::vector< std::pair<float, EVENT::TrackerHit*> >::iterator it=r2_values.begin(); it!=r2_values.end(); ++it) {
       trkHits.push_back(it->second);
     }
           
     
     bool fit_backwards = IMarlinTrack::backward;
          
     MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
     
     
     int error = 0;
     
     try {
              
       error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, Track, fit_backwards, covMatrix, _bField, _maxChi2PerHit);                              
       
     } catch (...) {
       
       //      delete Track;
       //      delete marlinTrk;
       
       throw ;
       
     }
     
     
#ifdef MARLINTRK_DIAGNOSTICS_ON
     if ( error != IMarlinTrack::success && _runMarlinTrkDiagnostics ) {        
       void * dcv = _trksystem->getDiagnositicsPointer();
       DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
       dc->skip_current_track();
     }        
#endif

     
     std::vector<std::pair<EVENT::TrackerHit* , double> > hits_in_fit ;  
     std::vector<std::pair<EVENT::TrackerHit* , double> > outliers ;
     std::vector<TrackerHit*> all_hits;    
     all_hits.reserve(300);
     
     marlinTrk->getHitsInFit(hits_in_fit);
     
     for ( unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
       all_hits.push_back(hits_in_fit[ihit].first);
     }
     
     UTIL::BitField64 cellID_encoder( lcio::ILDCellID0::encoder_string ) ; 
     
     MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);
     
     marlinTrk->getOutliers(outliers);
     
     for ( unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
       all_hits.push_back(outliers[ihit].first);
     }
     
     MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);
     
     delete marlinTrk;

     
     
     if( error != IMarlinTrack::success ) {       
            
       delete Track;
       streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Track fit failed with error code " << error << " track dropped. Number of hits = "<< trkHits.size() << std::endl;       
       continue ;
     }
     
     if( Track->getNdf() < 0) {       
       delete Track;
       streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Track fit returns " << Track->getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << std::endl;       
       continue ;
     }

     trk_col->addElement(Track);     
     
     const TrackState* trkStateIP = Track->getTrackState(lcio::TrackState::AtIP);
     
     if (trkStateIP == 0) {
       streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Track fit returns " << Track->getNdf() << " degress of freedom track dropped. Number of hits = "<< trkHits.size() << std::endl;       
       throw EVENT::Exception( std::string("SiliconTracking_MarlinTrk::FinalRefit: trkStateIP pointer == NULL ")  ) ;
     }
     
     // note trackAR which is of type TrackExtended, only takes fits set for ref point = 0,0,0 
     trackAR->setOmega(trkStateIP->getOmega());
     trackAR->setTanLambda(trkStateIP->getTanLambda());
     trackAR->setPhi(trkStateIP->getPhi());
     trackAR->setD0(trkStateIP->getD0());
     trackAR->setZ0(trkStateIP->getZ0());

     float cov[15];

     for (int i = 0 ; i<15 ; ++i) {
       cov[i] = trkStateIP->getCovMatrix().operator[](i);
     }
               
     trackAR->setCovMatrix(cov);
     trackAR->setChi2(Track->getChi2());
     trackAR->setNDF(Track->getNdf());

     nSiSegments++;

     HelixClass helix_final;

     helix_final.Initialize_Canonical(trkStateIP->getPhi(),trkStateIP->getD0(),trkStateIP->getZ0(),trkStateIP->getOmega(),trkStateIP->getTanLambda(),_bField);

     float trkPx = helix_final.getMomentum()[0];
     float trkPy = helix_final.getMomentum()[1];
     float trkPz = helix_final.getMomentum()[2];
     float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
     eTot += trkP;
     pxTot += trkPx;
     pyTot += trkPy;
     pzTot += trkPz;

   }
 }

 streamlog_out(DEBUG4) << "SiliconTracking_MarlinTrk -> run " << _nRun
 << " event " << _nEvt << std::endl;
 streamlog_out(DEBUG4) << "Number of Si tracks = " << nSiSegments << std::endl;
 streamlog_out(DEBUG4) << "Total 4-momentum of Si tracks : E = " << eTot
 << " Px = " << pxTot
 << " Py = " << pyTot
 << " Pz = " << pzTot << std::endl;


}


void SiliconTracking_MarlinTrk::setupGearGeom( const gear::GearMgr* gearMgr ){

 _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

 if( _reading_loi_data ){

   streamlog_out( MESSAGE ) << "  MarlinKalTest - Simple Disc Based FTD using parameters defined by SFtd05 Mokka driver " << std::endl ;

   // VXD

   const gear::ZPlanarParameters& pVXDDetMain = gearMgr->getVXDParameters();
   const gear::ZPlanarLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();

   _nLayersVTX = pVXDLayerLayout.getNLayers(); 




   // SIT

   const gear::GearParameters& pSIT = gearMgr->getGearParameters("SIT");

   const EVENT::DoubleVec& SIT_r   =  pSIT.getDoubleVals("SITLayerRadius" )  ;
   const EVENT::DoubleVec& SIT_hl  =  pSIT.getDoubleVals("SITSupportLayerHalfLength" )  ;

   _nLayersSIT = SIT_r.size() ; 

   if (_nLayersSIT != SIT_r.size() || _nLayersSIT != SIT_hl.size()) {

     streamlog_out( ERROR ) << "ILDSITCylinderKalDetector miss-match between DoubleVec and nlayers exit(1) called from file " << __FILE__ << " line " << __LINE__  << std::endl ;
     exit(1);

   }


   // FTD

   const gear::GearParameters& pFTD = gearMgr->getGearParameters("FTD");

   const EVENT::DoubleVec* pFTD_si  = NULL;
   const EVENT::DoubleVec* pFTD_sp  = NULL;
   const EVENT::DoubleVec* pFTD_ri  = NULL;
   const EVENT::DoubleVec* pFTD_ro  = NULL;
   const EVENT::DoubleVec* pFTD_z   = NULL;

   streamlog_out( MESSAGE ) << " For FTD using parameters defined by SFtd05 Mokka driver " << std::endl ;

   pFTD_si  =  &pFTD.getDoubleVals("FTDDiskSiThickness" )  ;
   pFTD_sp  =  &pFTD.getDoubleVals("FTDDiskSupportThickness" )  ;
   pFTD_ri  =  &pFTD.getDoubleVals("FTDInnerRadius" )  ;
   pFTD_ro  =  &pFTD.getDoubleVals("FTDOuterRadius" )  ;
   pFTD_z   =  &pFTD.getDoubleVals("FTDZCoordinate" )  ;

   _nlayersFTD = pFTD_z->size();

   for (unsigned int i = 0; i<_nlayersFTD; ++i) {
     _zLayerFTD.push_back((*pFTD_z)[i]);
   }

 }

 else {

   try {

     const gear::ZPlanarParameters& pVXDDetMain = gearMgr->getVXDParameters();
     const gear::ZPlanarLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();

     _nLayersVTX = pVXDLayerLayout.getNLayers(); 

     _VXDgeo.resize(_nLayersVTX);

     //SJA:FIXME: for now the support is taken as the same size the sensitive
     //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
     //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
     //           for a significant distance 

     for( unsigned int layer=0; layer < _nLayersVTX; ++layer){
       _VXDgeo[layer].nLadders = pVXDLayerLayout.getNLadders(layer); 
       _VXDgeo[layer].phi0 = pVXDLayerLayout.getPhi0(layer); 
       _VXDgeo[layer].dphi = 2*M_PI / _VXDgeo[layer].nLadders; 
       _VXDgeo[layer].senRMin = pVXDLayerLayout.getSensitiveDistance(layer); 
       _VXDgeo[layer].supRMin = pVXDLayerLayout.getLadderDistance(layer); 
       _VXDgeo[layer].length = pVXDLayerLayout.getSensitiveLength(layer); 
       _VXDgeo[layer].width = pVXDLayerLayout.getSensitiveWidth(layer); 
       _VXDgeo[layer].offset = pVXDLayerLayout.getSensitiveOffset(layer); 
       _VXDgeo[layer].senThickness = pVXDLayerLayout.getSensitiveThickness(layer); 
       _VXDgeo[layer].supThickness = pVXDLayerLayout.getLadderThickness(layer); 
     }

   } catch (gear::UnknownParameterException& e) {
     streamlog_out( MESSAGE ) << "  MarlinKalTest - VXD missing in gear file: VXD Not Built " << std::endl ;
   }

   try {

     const gear::ZPlanarParameters& pSITDetMain = gearMgr->getSITParameters();
     const gear::ZPlanarLayerLayout& pSITLayerLayout = pSITDetMain.getZPlanarLayerLayout();

     _nLayersSIT = pSITLayerLayout.getNLayers(); 
     _SITgeo.resize(_nLayersSIT);

     //SJA:FIXME: for now the support is taken as the same size the sensitive
     //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
     //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
     //           for a significant distance 

     for( unsigned int layer=0; layer < _nLayersSIT; ++layer){
       _SITgeo[layer].nLadders = pSITLayerLayout.getNLadders(layer); 
       _SITgeo[layer].phi0 = pSITLayerLayout.getPhi0(layer); 
       _SITgeo[layer].dphi = 2*M_PI / _SITgeo[layer].nLadders; 
       _SITgeo[layer].senRMin = pSITLayerLayout.getSensitiveDistance(layer); 
       _SITgeo[layer].supRMin = pSITLayerLayout.getLadderDistance(layer); 
       _SITgeo[layer].length = pSITLayerLayout.getSensitiveLength(layer); 
       _SITgeo[layer].width = pSITLayerLayout.getSensitiveWidth(layer); 
       _SITgeo[layer].offset = pSITLayerLayout.getSensitiveOffset(layer); 
       _SITgeo[layer].senThickness = pSITLayerLayout.getSensitiveThickness(layer); 
       _SITgeo[layer].supThickness = pSITLayerLayout.getLadderThickness(layer); 
     }

   } 
   catch (gear::UnknownParameterException& e) {
     streamlog_out( MESSAGE ) << "  MarlinKalTest - SIT missing in gear file: SIT Not Built " << std::endl ;
   }

   try {


     const gear::FTDParameters& ftdParams = gearMgr->getFTDParameters() ;
     const gear::FTDLayerLayout& ftdlayers = ftdParams.getFTDLayerLayout() ;


     _nlayersFTD = ftdlayers.getNLayers() ; 
     _FTDgeo.resize(_nlayersFTD);

     //SJA:FIXME: for now the support is taken as the same size the sensitive
     //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
     //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
     //           for a significant distance 

     for(unsigned int disk=0; disk< _nlayersFTD; ++disk){

       // numbers taken from the ILD_01 gear file for the sensitive part 
       _FTDgeo[disk].nPetals = ftdlayers.getNPetals(disk) ;    
       _FTDgeo[disk].dphi = 2*M_PI /  _FTDgeo[disk].nPetals ;
       _FTDgeo[disk].phi0 = ftdlayers.getPhi0(disk) ;
       _FTDgeo[disk].alpha = ftdlayers.getAlpha(disk) ;
       _FTDgeo[disk].rInner = ftdlayers.getSensitiveRinner(disk) ;
       _FTDgeo[disk].height = ftdlayers.getSensitiveWidth(disk) ;
       _FTDgeo[disk].innerBaseLength =  ftdlayers.getSensitiveLengthMin(disk) ;
       _FTDgeo[disk].outerBaseLength =  ftdlayers.getSensitiveLengthMax(disk) ;
       _FTDgeo[disk].senThickness =  ftdlayers.getSensitiveThickness(disk) ;
       _FTDgeo[disk].supThickness =  ftdlayers.getSupportThickness(disk) ;

//         _FTDgeo[disk].senZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
//         _FTDgeo[disk].senZPos_even_petal2 = ftdlayers.getSensitiveZposition(disk, 0, 2) ; 
//         _FTDgeo[disk].senZPos_even_petal3 = ftdlayers.getSensitiveZposition(disk, 0, 3) ; 
//         _FTDgeo[disk].senZPos_even_petal4 = ftdlayers.getSensitiveZposition(disk, 0, 4) ; 
//         
//         // currently the design assumes that the petal on the same side are at the same z
//         assert(_FTDgeo[disk].senZPos_even_petal1==_FTDgeo[disk].senZPos_even_petal2);
//         assert(_FTDgeo[disk].senZPos_even_petal3==_FTDgeo[disk].senZPos_even_petal4);
//         
//         _FTDgeo[disk].senZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
//         _FTDgeo[disk].senZPos_odd_petal2 = ftdlayers.getSensitiveZposition(disk, 1, 2) ; 
//         _FTDgeo[disk].senZPos_odd_petal3 = ftdlayers.getSensitiveZposition(disk, 1, 3) ; 
//         _FTDgeo[disk].senZPos_odd_petal4 = ftdlayers.getSensitiveZposition(disk, 1, 4) ; 
//         
//         // currently the design assumes that the petal on the same side are at the same z
//         assert(_FTDgeo[disk].senZPos_odd_petal1==_FTDgeo[disk].senZPos_odd_petal2);
//         assert(_FTDgeo[disk].senZPos_odd_petal3==_FTDgeo[disk].senZPos_odd_petal4);
//         
//         _FTDgeo[disk].supZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
//         _FTDgeo[disk].supZPos_even_petal2 = ftdlayers.getSensitiveZposition(disk, 0, 2) ; 
//         _FTDgeo[disk].supZPos_even_petal3 = ftdlayers.getSensitiveZposition(disk, 0, 3) ; 
//         _FTDgeo[disk].supZPos_even_petal4 = ftdlayers.getSensitiveZposition(disk, 0, 4) ; 
//         
//         assert(_FTDgeo[disk].supZPos_even_petal1==_FTDgeo[disk].supZPos_even_petal2);
//         assert(_FTDgeo[disk].supZPos_even_petal3==_FTDgeo[disk].supZPos_even_petal4);
//         
//         _FTDgeo[disk].supZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
//         _FTDgeo[disk].supZPos_odd_petal2 = ftdlayers.getSensitiveZposition(disk, 1, 2) ; 
//         _FTDgeo[disk].supZPos_odd_petal3 = ftdlayers.getSensitiveZposition(disk, 1, 3) ; 
//         _FTDgeo[disk].supZPos_odd_petal4 = ftdlayers.getSensitiveZposition(disk, 1, 4) ; 
//         
//         assert(_FTDgeo[disk].supZPos_odd_petal1==_FTDgeo[disk].supZPos_odd_petal2);
//         assert(_FTDgeo[disk].supZPos_odd_petal3==_FTDgeo[disk].supZPos_odd_petal4);

       ///////////////////////////////////////////////////////////////////////////////
       double thickness= _FTDgeo[disk].senThickness + _FTDgeo[disk].supThickness;

       _FTDgeo[disk].senZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
       _FTDgeo[disk].senZPos_even_petal2 = _FTDgeo[disk].senZPos_even_petal1;
       _FTDgeo[disk].senZPos_even_petal3 = _FTDgeo[disk].senZPos_even_petal1 + thickness;
       _FTDgeo[disk].senZPos_even_petal4 = _FTDgeo[disk].senZPos_even_petal3;

       _FTDgeo[disk].senZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
       _FTDgeo[disk].senZPos_odd_petal2 = _FTDgeo[disk].senZPos_odd_petal1;
       _FTDgeo[disk].senZPos_odd_petal3 = _FTDgeo[disk].senZPos_odd_petal1 + thickness;
       _FTDgeo[disk].senZPos_odd_petal4 = _FTDgeo[disk].senZPos_odd_petal3;


       _FTDgeo[disk].supZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
       _FTDgeo[disk].supZPos_even_petal2 = _FTDgeo[disk].senZPos_even_petal1;
       _FTDgeo[disk].supZPos_even_petal3 = _FTDgeo[disk].senZPos_even_petal1 + thickness;
       _FTDgeo[disk].supZPos_even_petal4 = _FTDgeo[disk].senZPos_even_petal3;

       _FTDgeo[disk].supZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
       _FTDgeo[disk].supZPos_odd_petal2 = _FTDgeo[disk].senZPos_odd_petal1;
       _FTDgeo[disk].supZPos_odd_petal3 = _FTDgeo[disk].senZPos_odd_petal1 + thickness;
       _FTDgeo[disk].supZPos_odd_petal4 = _FTDgeo[disk].senZPos_odd_petal3;


       ///////////////////////////////////////////////////////////////////////////////




       // rough check to see if the petal is rotated
       if( fabs(_FTDgeo[disk].alpha) > fabs(FLT_MIN)  ) { 
         streamlog_out( ERROR ) << "  SiliconTracking_MarlinTrk - tilted design not supported exit(1) " << std::endl ;
         exit(1);
       }

     }

     for (unsigned int disk=0; disk < _nlayersFTD; ++disk) {
       _zLayerFTD.push_back(_FTDgeo[disk].senZPos_even_petal1); // front petal even numbered
       _zLayerFTD.push_back(_FTDgeo[disk].senZPos_odd_petal1);  // front petal odd numbered
     }

     // SJA as disks are staggered lets treat them internally as 2*ndisksFTD
     _nlayersFTD =_zLayerFTD.size() ;

   } 

   catch (gear::UnknownParameterException& e) {
     streamlog_out( MESSAGE ) << "  SiliconTracking_MarlinTrk - FTD missing in gear file: FTD Not Built " << std::endl ;
   }

 }



}

void SiliconTracking_MarlinTrk::TracksWithNHitsContainer::clear()
{
 for (std::vector< TrackExtendedVec >::iterator trackVecIter = _tracksNHits.begin();
      trackVecIter < _tracksNHits.end(); trackVecIter++)
 {
   for (TrackExtendedVec::iterator trackIter = trackVecIter->begin();
	 trackIter < trackVecIter->end(); trackIter++)
   {
     delete *trackIter;
   }

   trackVecIter->clear();
 }
}
