#include "FullLDCTracking_MarlinTrk.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <map>
#include <marlin/Global.h>
#include "ClusterShapes.h"

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include "gear/FTDLayerLayout.h"
#include "gear/FTDParameters.h"
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>

#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include <climits>
#include <cmath>

using namespace lcio ;
using namespace marlin ;
using namespace MarlinTrk ;

FullLDCTracking_MarlinTrk aFullLDCTracking_MarlinTrk ;

FullLDCTracking_MarlinTrk::FullLDCTracking_MarlinTrk() : Processor("FullLDCTracking_MarlinTrk") {  
  _description = "Performs full tracking in ILD detector" ;  
  
  _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
  
  // Input tracker hit collections
  
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
  
  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "VTXHitCollection",
                          "VTX Hit Collection Name",
                          _VTXTrackerHitCollection,
                          std::string("VTXTrackerHits"));  
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "SITHitCollection",
                          "SIT Hit Collection Name",
                          _SITTrackerHitCollection,
                          std::string("SITTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "SETHitCollection",
                          "SET Hit Collection Name",
                          _SETTrackerHitCollection,
                          std::string("SETTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "ETDHitCollection",
                          "ETD Hit Collection Name",
                          _ETDTrackerHitCollection,
                          std::string("ETDTrackerHits"));
  
  registerInputCollection(LCIO::TRACKERHIT,
                          "TPCHitCollection",
                          "TPC Hit Collection Name",
                          _TPCTrackerHitCollection,
                          std::string("TPCTrackerHits"));
  
  
  // Input track collections
  
  registerInputCollection(LCIO::TRACK,
                          "TPCTracks",
                          "TPC Track Collection",
                          _TPCTrackCollection,
                          std::string("TPCTracks"));
  
  registerInputCollection(LCIO::TRACK,
                          "SiTracks",
                          "Si Track Collection",
                          _SiTrackCollection,
                          std::string("SiTracks"));
  
  // Input relation collections
  
  registerInputCollection(LCIO::LCRELATION,
                          "TPCTracksMCPRelColl",
                          "TPC Track to MCP Relation Collection Name",
                          _TPCTrackMCPCollName,
                          std::string("TPCTracksMCP"));
  
  registerInputCollection(LCIO::LCRELATION,
                          "SiTracksMCPRelColl",
                          "Si Track to Collection",
                          _SiTrackMCPCollName,
                          std::string("SiTracksMCP"));
  
  // Output track collection
  registerOutputCollection(LCIO::TRACK,
                           "LDCTrackCollection",
                           "LDC track collection name",
                           _LDCTrackCollection,
                           std::string("LDCTracks"));
  
  
  // Output relation collection
  registerOutputCollection(LCIO::LCRELATION,
                           "LDCTrackMCPRelCollection",
                           "Collection name for the LDC track to MCParticle relations",
                           _LDCTrackMCPCollection,
                           std::string("LDCTracksMCP"));
  
  
  
  // steering parameters
  
  registerProcessorParameter("D0CutForMerging",
                             "Cut on D0 difference for merging of Si and TPC segments",
                             _d0CutForMerging,
                             float(500.0));
  
  registerProcessorParameter("Z0CutForMerging",
                             "Cut on Z0 difference for merging of Si and TPC segments",
                             _z0CutForMerging,
                             float(1000.0));
  
  registerProcessorParameter("OmegaCutForMerging",
                             "Cut on Omega difference for merging Si and TPC segments",
                             _dOmegaForMerging,
                             float(0.25));
  
  registerProcessorParameter("AngleCutForMerging",
                             "Cut on Opening Angle for merging Si and TPC segments",
                             _angleForMerging,
                             float(0.10));
  
  registerProcessorParameter("Chi2FitCut",
                             "Cut on fit Chi2",
                             _chi2FitCut,
                             float(100.0));
  
  registerProcessorParameter("Chi2PrefitCut",
                             "Cut on fit Chi2",
                             _chi2PrefitCut,
                             float(1.0e+5));
  
  registerProcessorParameter("CreateMap",
                             "Create Track to MCP Relations",
                             _createMap,
                             int(1));
  
  registerProcessorParameter("Debug",
                             "Activate debugging?",
                             _debug,
                             int(0));
  
  
  registerProcessorParameter("ForceSiTPCMerging",
                             "Force merging of Si and TPC segments?",
                             _forceMerging,
                             int(0));
  
  registerProcessorParameter("D0CutForForcedMerging",
                             "Cut on D0 difference for forced merging of Si and TPC segments",
                             _d0CutForForcedMerging,
                             float(50.));
  
  registerProcessorParameter("Z0CutForForcedMerging",
                             "Cut on Z0 difference for forced merging of Si and TPC segments",
                             _z0CutForForcedMerging,
                             float(200.));
  
  registerProcessorParameter("OmegaCutForForcedMerging",
                             "Cut on Omega difference for forced merging of Si and TPC segments",
                             _dOmegaForForcedMerging,
                             float(0.15));
  
  registerProcessorParameter("AngleCutForForcedMerging",
                             "Cut on Opening Angle for forced merging of Si and TPC segments",
                             _angleForForcedMerging,
                             float(0.05));
  
  registerProcessorParameter("ForceTPCSegmentsMerging",
                             "Force merging of TPC Segments?",
                             _mergeTPCSegments,
                             int(1));
  
  registerProcessorParameter("D0CutToMergeTPCSegments",
                             "Cut on D0 difference for merging TPC segments",
                             _d0CutToMergeTPC,
                             float(100.));
  
  registerProcessorParameter("Z0CutToMergeTPCSegments",
                             "Cut on Z0 difference for merging TPC segments",
                             _z0CutToMergeTPC,
                             float(5000.0));
  
  registerProcessorParameter("DeltaPCutToMergeTPCSegments",
                             "Cut on dP/P difference for merging TPC segments",
                             _dPCutToMergeTPC,
                             float(0.1));
  
  registerProcessorParameter("PtCutToMergeTPCSegments",
                             "Cut on Pt of tracks for merging TPC segments",
                             _PtCutToMergeTPC,
                             float(1.2));
  
  
  
  registerProcessorParameter("cosThetaCutHighPtMerge",
                             "Cut on cos theta between the two momentum vectors when considering merger of high Pt tracks",
                             _cosThetaCutHighPtMerge,
                             float(0.99));
  
  registerProcessorParameter("cosThetaCutSoftHighPtMerge",
                             "cut on cos theta between the two momentum vectors when considering merger of high Pt tracks for softer dp/p cut",
                             _cosThetaCutSoftHighPtMerge,
                             float(0.998));
  
  registerProcessorParameter("momDiffCutHighPtMerge",
                             "cut on dp/p when considering merger of high Pt tracks",
                             _momDiffCutHighPtMerge,
                             float(0.01));
  
  registerProcessorParameter("momDiffCutSoftHighPtMerge",
                             "softer cut on dp/p when considering merger of high Pt tracks",
                             _momDiffCutSoftHighPtMerge,
                             float(0.25));
  
  registerProcessorParameter("hitDistanceCutHighPtMerge",
                             "cut on 3D distance between hit and helix extrapolation when considering merger of high Pt tracks",
                             _hitDistanceCutHighPtMerge,
                             float(25.0));
  
  registerProcessorParameter("maxHitDistanceCutHighPtMerge",
                             "cut for max 3D distance between any hit and helix extrapolation when considering merger of high Pt tracks",
                             _maxHitDistanceCutHighPtMerge,
                             float(50.0));
  
  registerProcessorParameter("maxFractionOfOutliersCutHighPtMerge",
                             "cut on maximum fraction of outliers when considering merger of high Pt tracks",
                             _maxFractionOfOutliersCutHighPtMerge,
                             float(0.95));
  
  
  
  
  
  registerProcessorParameter("CutOnTPCHits",
                             "Cut on the number of the TPC hits for tracks with no Si hits",
                             _cutOnTPCHits,
                             int(35));
  registerProcessorParameter("CutOnSiHits",
                             "Cut on the number of the Si hits for tracks with no TPC hits",
                             _cutOnSiHits,
                             int(4));
  
  
  registerProcessorParameter("AssignVTXHits",
                             "Assign left over VTX hits",
                             _assignVTXHits,
                             int(1));
  
  registerProcessorParameter("AssignFTDHits",
                             "Assign left over FTD hits",
                             _assignFTDHits,
                             int(1));
  
  registerProcessorParameter("AssignSITHits",
                             "Assign left over SIT hits",
                             _assignSITHits,
                             int(1));
  
  registerProcessorParameter("AssignTPCHits",
                             "Assign left over TPC hits",
                             _assignTPCHits,
                             int(1));
  
  registerProcessorParameter("AssignSETHits",
                             "Assign SET Hits",
                             _assignSETHits,
                             int(1));
  
  
  registerProcessorParameter("AssignETDHits",
                             "Assign ETD Hits",
                             _assignETDHits,
                             int(1));
  
  registerProcessorParameter("NHitsExtrapolation",
                             "number of hits for outward extrapolation",
                             _nHitsExtrapolation,
                             int(35));
  
  registerProcessorParameter("StoreHitsInFit",
                             "Store only hits used in fit?",
                             _storeHitsInFit,
                             int(0));
  
  
  registerProcessorParameter("VTXHitToTrackDistance",
                             "Cut on distance between track and VTX hits",
                             _distCutForVTXHits,
                             float(1.5));
  
  
  registerProcessorParameter("FTDHitToTrackDistance",
                             "Cut on distance between track and FTD hits",
                             _distCutForFTDHits,
                             float(2.0));
  
  
  registerProcessorParameter("SITHitToTrackDistance",
                             "Cut on distance between track and SIT hits",
                             _distCutForSITHits,
                             float(2.0));
  
  registerProcessorParameter("SETHitToTrackDistance",
                             "Cut on distance between track and SET hits",
                             _distCutForSETHits,
                             float(2.0));
  
  
  registerProcessorParameter("ETDHitToTrackDistance",
                             "Cut on distance between track and ETD hits",
                             _distCutForETDHits,
                             float(10.0));
  
  
  registerProcessorParameter("TPCHitToTrackDistance",
                             "Cut on distance between track and TPC hits",
                             _distCutForTPCHits,
                             float(15.0));
  
  
  registerProcessorParameter("CutOnTrackD0",
                             "Cut on the track parameter D0",
                             _d0TrkCut,
                             float(500.));
  
  
  registerProcessorParameter("CutOnTrackZ0",
                             "Cut on the track parameter Z0",
                             _z0TrkCut,
                             float(500.));
  
  
  registerProcessorParameter("ForbidOverlapInZTPC",
                             "Forbid overlap in Z for the merged TPC segments",
                             _forbidOverlapInZTPC,
                             int(0));
  
  registerProcessorParameter("ForbidOverlapInZComb",
                             "Forbid overlap in Z for combining TPC segments with tracks having Si hits",
                             _forbidOverlapInZComb,
                             int(0));
  
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

  
}



void FullLDCTracking_MarlinTrk::init() { 
  
  printParameters();  
  _nRun = -1 ;
  _nEvt = 0 ;
  PI = acos(-1.);
  PIOVER2 = 0.5*PI;
  TWOPI = 2*PI;
  
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
  
  this->setupGearGeom(Global::GEAR);
  
}

void FullLDCTracking_MarlinTrk::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
  _nEvt = 0;
  streamlog_out(DEBUG4) << std::endl;
  streamlog_out(DEBUG4) << "FullLDCTracking_MarlinTrk ---> new run : run number = " << _nRun << std::endl;
  
} 

void FullLDCTracking_MarlinTrk::processEvent( LCEvent * evt ) { 
  
  _evt = evt;
  
  
  streamlog_out(DEBUG4) << std::endl;
  streamlog_out(DEBUG4) << "FullLDCTracking_MarlinTrk -> run = " << _nRun 
  << "  event = " << _nEvt << std::endl;
  streamlog_out(DEBUG4) << std::endl;
  
  
  prepareVectors( evt );
  streamlog_out(DEBUG3) << "prepareVectors done..." << std::endl;
  streamlog_out(DEBUG3) << "************************************Merge TPC/Si" << std::endl;
  MergeTPCandSiTracks();
  streamlog_out(DEBUG3) << "************************************Merging done..." << std::endl;
  MergeTPCandSiTracksII();
  streamlog_out(DEBUG3) << "************************************Merging II done..." << std::endl;
  Sorting(_allCombinedTracks);
  streamlog_out(DEBUG3) << "************************************Sorting done..." << std::endl;
  SelectCombinedTracks();
  streamlog_out(DEBUG3) << "************************************Selection of combined tracks done..." << std::endl;
  AddNotCombinedTracks( );
  streamlog_out(DEBUG3) << "************************************Not combined tracks added..." << std::endl;
  //CheckTracks( );
  
  AddNotAssignedHits();
  streamlog_out(DEBUG3) << "***********************************Not assigned hits added..." << std::endl;
  AddTrackColToEvt(evt,_trkImplVec,
                   _LDCTrackCollection,_LDCTrackMCPCollection);
  streamlog_out(DEBUG3) << "Collections added to event..." << std::endl;
  CleanUp();
  streamlog_out(DEBUG3) << "Cleanup is done..." << std::endl;
  _nEvt++;
  //  getchar();
  streamlog_out(DEBUG3) << std::endl;
  streamlog_out(DEBUG3) << std::endl;
  
}

void FullLDCTracking_MarlinTrk::AddTrackColToEvt(LCEvent * evt, TrackExtendedVec & trkVec, 
                                                 std::string TrkColName, std::string RelColName) {
  
  LCCollectionVec * colTRK = new LCCollectionVec(LCIO::TRACK);
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  colTRK->setFlag( trkFlag.getFlag()  ) ;  
  
  streamlog_out(DEBUG4)<< "Collection " << TrkColName << " is being added to event " << std::endl;
  
  LCCollectionVec * colRel = NULL;
  
  if (_createMap) {
    colRel = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    colRel->setFlag( lcFlag.getFlag()  ) ;
  }
  
  int nTrkCand = int(trkVec.size());
  
  int nTotTracks = 0;
  float eTot = 0.0;
  float pxTot = 0.0;
  float pyTot = 0.0;
  float pzTot = 0.0;
  
  //SJA:FIXME: So here we are going to do one final refit. This can certainly be optimised, but rather than worry about the mememory management right now lets make it work, and optimise it later ...
  
  
  for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
    
    TrackExtended * trkCand = trkVec[iTRK];
    TrackerHitExtendedVec& hitVec = trkCand->getTrackerHitExtendedVec();
    
    EVENT::TrackerHitVec trkHits;
    
    int nHits = int(hitVec.size());
    for (int ihit=0;ihit<nHits;++ihit) {

      EVENT::TrackerHit* trkHit = hitVec[ihit]->getTrackerHit();
      if(trkHit) { 
        trkHits.push_back(trkHit);   
      }
      else{
        throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: TrackerHit pointer == NULL ")  ) ;
      }
      
    }
    
    
    if( trkHits.size() < 3 ) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Cannot fit less than 3 hits. Number of hits =  " << trkHits.size() << std::endl;
      continue ; 
    }
    
    MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
    
    // hits are in reverse order 
    
    sort(trkHits.begin(), trkHits.end(), FullLDCTracking_MarlinTrk::compare_r() );
    
    EVENT::TrackerHitVec::iterator it = trkHits.begin();
    
    streamlog_out(DEBUG2) << "Start Fitting: AddHits: number of hits to fit " << trkHits.size() << std::endl;
    
//    int number_of_added_hits = 0;
//    for( it = trkHits.begin() ; it != trkHits.end() ; ++it )
//        {
//      
//      if (marlin_trk->addHit(*it) == 0){
//        ++number_of_added_hits;
//      }
//      else{
//        streamlog_out(DEBUG4) << "Hit " << it - trkHits.begin() << " Dropped " << std::endl;
//      }
//      
//        }
//    
//    if( number_of_added_hits < 3 ) {
//      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Cannot fit less than 3 hits. Number of hits =  " << number_of_added_hits << std::endl;
//      delete marlin_trk ;
//      continue ;
//    }
//    
//    marlin_trk->initialise( IMarlinTrack::backward ) ;

    
    int number_of_added_hits = 0;
    int ndof_added = 0;
    TrackerHitVec added_hits;
    
    for( it = trkHits.begin() ; it != trkHits.end() ; ++it ) {
      
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
      else { // normal non composite hit
        
        if (marlin_trk->addHit( trkHit ) == 0) {
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
      streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Cannot fit less with less than 8 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << std::endl;
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
    
    for (int icov = 0; icov<covMatrix.size(); ++icov) {
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

    
    int fit_status = marlin_trk->fit() ; 
    
    if( fit_status != 0 ){ 
      delete marlin_trk ;
      continue;
    }
    
    const gear::Vector3D point(0.,0.,0.); // nominal IP
    int return_code = 0;
    
    double chi2;
    int ndf;
    
    TrackImpl * track_lcio = new TrackImpl();
    
    TrackStateImpl* trkStateIP = new TrackStateImpl;
    return_code = marlin_trk->propagate(point, *trkStateIP, chi2, ndf ) ;
    
    if (return_code !=MarlinTrk::IMarlinTrack::success ) {
      streamlog_out( ERROR ) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at IP: Track Discarded" << std::endl ;
      delete marlin_trk ;
      delete trkStateIP;
      delete track_lcio;
      continue;
    }
    
    // fitting finished get hit in the fit
    
    std::vector<std::pair<EVENT::TrackerHit*, double> > hits_in_fit;
    
    // remember the hits are ordered in the order in which they were fitted
    // here we are fitting inwards to the first is the last and vice verse
    
    marlin_trk->getHitsInFit(hits_in_fit);
    
    if( hits_in_fit.size() < 3 ) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Less than 3 hits in fit: Track Discarded. Number of hits =  " << trkHits.size() << std::endl;
      delete marlin_trk ;
      delete trkStateIP;
      delete track_lcio;
      continue ; 
    }

    
    EVENT::TrackerHit* last_hit_in_fit = hits_in_fit.front().first;
    if (!last_hit_in_fit) {
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::FinalRefit: TrackerHit pointer == NULL ")  ) ;
    }
    
    EVENT::TrackerHit* first_hit_in_fit = hits_in_fit.back().first;
    if (!first_hit_in_fit) {
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::FinalRefit: TrackerHit pointer == NULL ")  ) ;
    }
    
    
    // SJA:FIXME: As we are not smoothing back then the fits at the last hit will probably still be rubbish 
    // must find a method to control smoothing ...

    
    
    TrackStateImpl* trkStateFirstHit = new TrackStateImpl;
    return_code = marlin_trk->getTrackState(first_hit_in_fit, *trkStateFirstHit, chi2, ndf ) ;
    
    if (return_code !=MarlinTrk::IMarlinTrack::success ) {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at First Hit return_code = " << return_code << std::endl ;
//      delete marlin_trk ;
//      delete trkStateFirstHit;
//      delete track_lcio;
//      continue;
    }
    
    TrackStateImpl* trkStateLastHit = new TrackStateImpl;
    return_code = marlin_trk->getTrackState(last_hit_in_fit, *trkStateLastHit, chi2, ndf ) ;
    
    if (return_code !=MarlinTrk::IMarlinTrack::success ) {
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at Last Hit : return_code = " << return_code << std::endl ;
//      delete marlin_trk ;
//      delete trkStateLastHit;
//      delete track_lcio;
//      continue;
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
      streamlog_out( DEBUG5 ) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at Calo Face : return_code = " << 
      return_code << std::endl ;
      
      //      delete marlin_trk ;
      //      delete trkStateCalo;
      //      delete track_lcio;
      //      continue;
    }
    
    delete marlin_trk;
    
    trkStateIP->setLocation(  lcio::TrackState::AtIP ) ;
    trkStateFirstHit->setLocation(  lcio::TrackState::AtFirstHit ) ;
    trkStateLastHit->setLocation(  lcio::TrackState::AtLastHit) ;
    trkStateCalo->setLocation(  lcio::TrackState::AtCalorimeter ) ;
    
    track_lcio->trackStates().push_back(trkStateIP);
    track_lcio->trackStates().push_back(trkStateFirstHit);
    track_lcio->trackStates().push_back(trkStateLastHit);
    track_lcio->trackStates().push_back(trkStateCalo);
    
    track_lcio->setChi2(chi2);
    track_lcio->setNdf(ndf);
    
    const double* pos = trkHits.front()->getPosition();
    
    double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    track_lcio->setRadiusOfInnermostHit(r);
    
    std::map<int, int> hitNumbers; 
    
    hitNumbers[lcio::ILDDetID::VXD] = 0;
    hitNumbers[lcio::ILDDetID::SIT] = 0;
    hitNumbers[lcio::ILDDetID::FTD] = 0;
    hitNumbers[lcio::ILDDetID::TPC] = 0;
    hitNumbers[lcio::ILDDetID::SET] = 0;
    hitNumbers[lcio::ILDDetID::ETD] = 0;
    
    
    std::vector<MCParticle*> mcPointers ;
    std::vector<int> mcHits ;
    mcPointers.clear();
    mcHits.clear();
    
    for(int j=trkHits.size()-1; j>=0; --j) {
      track_lcio->addHit(trkHits.at(j)) ;
      ++hitNumbers[ getDetectorID(trkHits.at(j)) ];
      
      if (_createMap > 0) {
        int nSH = int(trkHits.at(j)->getRawHits().size());
        for (int ish=0;ish<nSH;++ish) {
          SimTrackerHit * simHit = dynamic_cast<SimTrackerHit*>(trkHits.at(j)->getRawHits()[ish]);
          MCParticle * mcp = simHit->getMCParticle();
          bool found = false;
          int nMCP = int(mcPointers.size());
          for (int iMCP=0;iMCP<nMCP;++iMCP) {
            if (mcp == mcPointers[iMCP]) {
              found = true;
              mcHits[iMCP]++;
              break;
            }
          }
          if (!found) {
            mcPointers.push_back(mcp);
            mcHits.push_back(1);
          }
        }
      }
    }
    
    //SJA:FIXME no distiction made for hits in fit or not
    track_lcio->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ] = hitNumbers[lcio::ILDDetID::VXD];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 2 ] = hitNumbers[lcio::ILDDetID::FTD];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 2 ] = hitNumbers[lcio::ILDDetID::SIT];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 2 ] = hitNumbers[lcio::ILDDetID::TPC];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 2 ] = hitNumbers[lcio::ILDDetID::SET];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 2 ] = hitNumbers[lcio::ILDDetID::ETD];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 1 ] = hitNumbers[lcio::ILDDetID::VXD];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::FTD - 1 ] = hitNumbers[lcio::ILDDetID::FTD];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SIT - 1 ] = hitNumbers[lcio::ILDDetID::SIT];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::TPC - 1 ] = hitNumbers[lcio::ILDDetID::TPC];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::SET - 1 ] = hitNumbers[lcio::ILDDetID::SET];
    track_lcio->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::ETD - 1 ] = hitNumbers[lcio::ILDDetID::ETD];
    
    
    GroupTracks * group = trkCand->getGroupTracks();
    
    if (group != NULL) {
      TrackExtendedVec trkVecGrp = group->getTrackExtendedVec();
      int nGrTRK = int(trkVecGrp.size());
      for (int iGr=0;iGr<nGrTRK;++iGr) {
        TrackExtended * subTrack = trkVecGrp[iGr];
        track_lcio->addTrack(subTrack->getTrack());
      }
    }
    
    float d0TrkCand = trkCand->getD0();
    float z0TrkCand = trkCand->getZ0();
    //    float phi0TrkCand = trkCand->getPhi();
    
    
    int nHitsSi = hitNumbers[lcio::ILDDetID::VXD]+hitNumbers[lcio::ILDDetID::FTD]+hitNumbers[lcio::ILDDetID::SIT];
    
    bool rejectTrack = (hitNumbers[lcio::ILDDetID::TPC]<_cutOnTPCHits) && (nHitsSi<=0);
    
    rejectTrack = rejectTrack || ( (hitNumbers[lcio::ILDDetID::TPC]<=0) && (nHitsSi<_cutOnSiHits) );
    rejectTrack = rejectTrack || ( fabs(d0TrkCand) > _d0TrkCut ) || ( fabs(z0TrkCand) > _z0TrkCut );
    
    if ( rejectTrack ) {
      delete track_lcio;
    }
    
    else {
      
      float omega = trkStateIP->getOmega();
      float tanLambda = trkStateIP->getTanLambda();
      float phi0 = trkStateIP->getPhi();
      float d0 = trkStateIP->getD0();
      float z0 = trkStateIP->getZ0();
      
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      
      float trkPx = helix.getMomentum()[0];
      float trkPy = helix.getMomentum()[1];
      float trkPz = helix.getMomentum()[2];
      float trkP = sqrt(trkPx*trkPx+trkPy*trkPy+trkPz*trkPz);
      
      eTot += trkP;
      pxTot += trkPx;
      pyTot += trkPy;
      pzTot += trkPz;   
      nTotTracks++;
      
      colTRK->addElement(track_lcio);
      
      if (_createMap > 0) {
        int nRel = int(mcPointers.size());
        for (int k=0;k<nRel;++k) {
          LCRelationImpl* rel = new LCRelationImpl;
          MCParticle * mcp = mcPointers[k];
          rel->setFrom (track_lcio);
          rel->setTo (mcp);
          float weight = (float)(mcHits[k])/(float)(track_lcio->getTrackerHits().size());
          rel->setWeight(weight);
          colRel->addElement(rel);
        }
      }
    }
  }
  
  streamlog_out(DEBUG4) << std::endl;
  streamlog_out(DEBUG4) << "Number of accepted " << TrkColName << " = " 
  << nTotTracks << std::endl;
  streamlog_out(DEBUG4) << "Total 4-momentum of " << TrkColName << " : E = " << eTot
  << " Px = " << pxTot
  << " Py = " << pyTot
  << " Pz = " << pzTot << std::endl;
  streamlog_out(DEBUG4) << std::endl;
  
  evt->addCollection(colTRK,TrkColName.c_str());
  if (_createMap)
    evt->addCollection(colRel,RelColName.c_str());
  
}


void FullLDCTracking_MarlinTrk::prepareVectors(LCEvent * event ) {
  
  
  
  
  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allSETHits.clear();
  _allETDHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _allNonCombinedTPCTracks.clear();
  _allNonCombinedSiTracks.clear();
  _trkImplVec.clear();
  _candidateCombinedTracks.clear();
  
  
  std::map <TrackerHit*,TrackerHitExtended*> mapTrackerHits;
  
  // Reading TPC hits
  try {
    
    LCCollection * col = event->getCollection(_TPCTrackerHitCollection.c_str());
    
    int nelem = col->getNumberOfElements();
    
    for (int ielem=0;ielem<nelem;++ielem) {
      
      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      
      // Covariance Matrix in LCIO is defined in XYZ convert to R-Phi-Z
      // For no error in r
      
      double tpcRPhiRes = sqrt(hit->getCovMatrix()[0] + hit->getCovMatrix()[2]);
      double tpcZRes = sqrt(hit->getCovMatrix()[5]);
      
      hitExt->setResolutionRPhi(float(tpcRPhiRes));
      hitExt->setResolutionZ(float(tpcZRes));
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));
      _allTPCHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _TPCTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  };
  

  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    
    LCCollection * hitCollection = event->getCollection(_FTDPixelHitCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of FTD Pixel hits = " << nelem << std::endl;
    
    for (int ielem=0; ielem<nelem; ++ielem) {
      TrackerHitPlane * hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));
      
      TrackerHitExtended * hitExt = new TrackerHitExtended( hit );
      
      double point_res_rphi = sqrt( hit->getdU()*hit->getdU() + hit->getdV()*hit->getdV() );
      hitExt->setResolutionRPhi( point_res_rphi );
      
      // SJA:FIXME why is this needed? 
      hitExt->setResolutionZ(0.1);
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));            
      hitExt->setDet(int(INT_MAX));
      
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
      
    }
  }
  catch(DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _FTDPixelHitCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  // Reading in FTD SpacePoint Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    
    LCCollection * hitCollection = event->getCollection(_FTDSpacePointCollection.c_str());
    
    int nelem = hitCollection->getNumberOfElements();
    
    streamlog_out(DEBUG4) << "Number of FTD SpacePoints hits = " << nelem << std::endl;
    
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
            
      _allFTDHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
      
    }
  }
  catch(DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _FTDSpacePointCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  
  //  // Reading ETD Hits
  //  try {
  //    LCCollection * col = event->getCollection(_ETDTrackerHitCollection.c_str());
  //    int nelem = col->getNumberOfElements();
  //    for (int ielem=0;ielem<nelem;++ielem) {
  //      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
  //      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
  //      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[0])));
  //      hitExt->setResolutionZ(0.1);
  //      //      hitExt->setResolutionRPhi(_resolutionRPhi_FTD);
  //      //      hitExt->setResolutionZ(_resolutionZ_FTD);
  //      // type and det are no longer used, set to INT_MAX to try and catch any missuse
  //      hitExt->setType(int(INT_MAX));
  //      hitExt->setDet(int(INT_MAX));
  //      _allETDHits.push_back( hitExt );
  //      mapTrackerHits[hit] = hitExt;
  //    }
  //  }
  //  catch( DataNotAvailableException &e ) {
  //      streamlog_out(DEBUG4) << _ETDTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  //  }
  
  // Reading SIT Hits
  try {

    LCCollection * col = event->getCollection(_SITTrackerHitCollection.c_str());
    
    int nelem = col->getNumberOfElements();
    
    for (int ielem=0;ielem<nelem;++ielem) {
    
      double drphi(NAN);
      double dz(NAN);
      
      TrackerHit* hit(NULL);

      if (_reading_loi_data) { // uses cylindrical hits
        
        TrackerHitZCylinder * hit_c = dynamic_cast<TrackerHitZCylinder*>(col->getElementAt(ielem));
        
        if (hit_c) {
          hit= hit_c;
        }
        else{
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: Dynamic cast to TrackerHitZCylinder failed. \n\n exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }
        
        drphi = hit_c->getdRPhi();
        dz    = hit_c->getdZ();
        
      }
      
      else {
        
        
        TrackerHit * hit_p = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
        
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


      
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(drphi);
      hitExt->setResolutionZ(dz);
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));
      _allSITHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _SITTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  // Reading SET Hits
  //    sprintf(strg,"SET hits ----->\n");
  //    streamlog_out(DEBUG4) << strg;
  //    sprintf(strg," id    r_hit   phi_hit    z_hit  e(r-phi)  e(z)  \n");
  //    streamlog_out(DEBUG4) << strg;
  
  //     "  0   1807.64   -1.97     26.45   0.010   0.010
  //     "  1   1812.64   -1.97     26.51   0.010   0.010
  //  try {
  //    LCCollection * col = event->getCollection(_SETTrackerHitCollection.c_str());
  //    int nelem = col->getNumberOfElements();
  //    for (int ielem=0;ielem<nelem;++ielem) {
  //      TrackerHit * hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
  //      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
  //      hitExt->setResolutionRPhi(float(sqrt(hit->getCovMatrix()[2])));
  //      hitExt->setResolutionZ(float(sqrt(hit->getCovMatrix()[5])));
  //      //       float x_hit = float(hit->getPosition()[0]);
  //      //       float y_hit = float(hit->getPosition()[1]);
  //      //       float z_hit = float(hit->getPosition()[2]);
  //      //       float phi_hit = atan2( y_hit, x_hit );
  //      //       float r_hit = sqrt( x_hit*x_hit + y_hit*y_hit );
  //      //       sprintf(strg,"%3i  %8.2f  %6.3f  %8.2f  %6.3f  %6.3f   ", 
  //      //             ielem, r_hit, phi_hit, z_hit, 
  //      //             hitExt->getResolutionRPhi(), 
  //      //             hitExt->getResolutionZ() );
  //      //       streamlog_out(DEBUG4) << strg;
  //
  //      //       std::cout << hit << std::endl;
  //      // type and det are no longer used, set to INT_MAX to try and catch any missuse
  //      hitExt->setType(int(INT_MAX));
  //      hitExt->setDet(int(INT_MAX));
  //      _allSETHits.push_back( hitExt );
  //      mapTrackerHits[hit] = hitExt;
  //    }
  //  }
  //  catch( DataNotAvailableException &e ) {
  //      streamlog_out(DEBUG4) << _SETTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  //  }
  
  // Reading VTX Hits
  try {
    LCCollection * col = event->getCollection(_VTXTrackerHitCollection.c_str());
    int nelem = col->getNumberOfElements();
    for (int ielem=0;ielem<nelem;++ielem) {
      TrackerHitPlane * hit = dynamic_cast<TrackerHitPlane*>(col->getElementAt(ielem));
      TrackerHitExtended * hitExt = new TrackerHitExtended(hit);
      
      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(hit->getdU());
      hitExt->setResolutionZ(hit->getdV());
      
      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));      
      hitExt->setDet(int(INT_MAX));
      _allVTXHits.push_back( hitExt );
      mapTrackerHits[hit] = hitExt;
    }
  }
  catch( DataNotAvailableException &e ) {
    streamlog_out(DEBUG4) << _VTXTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  
  // Reading TPC Tracks
  try {
    LCCollection * col = event->getCollection(_TPCTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG3) << std::endl;
    streamlog_out(DEBUG3) << "Number of TPC Tracks = " << nelem << std::endl;
    streamlog_out(DEBUG3) << " Trk       p          D0         Z0       Px       Py       Pz    ntpc ndf Chi2/ndf" << std::endl;
    //           "  0  1.111   0.059      0.022    -0.54     0.61    -0.45    0.185
    
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * tpcTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( tpcTrack );
      TrackerHitVec hitVec = tpcTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(tpcTrack->getOmega());
      trackExt->setTanLambda(tpcTrack->getTanLambda());
      trackExt->setPhi(tpcTrack->getPhi());
      trackExt->setD0(tpcTrack->getD0());
      trackExt->setZ0(tpcTrack->getZ0());
      float cov[15];
      float param[5];
      //      float reso[4];
      param[0] = tpcTrack->getOmega();
      param[1] = tpcTrack->getTanLambda();
      param[2] = tpcTrack->getPhi();
      param[3] = tpcTrack->getD0();
      param[4] = tpcTrack->getZ0();
      
      
      const FloatVec Cov = tpcTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
        cov[ic] =  Cov[ic];
      }
      
      
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(tpcTrack->getNdf());
      trackExt->setChi2(tpcTrack->getChi2());            
      char strg[200];
      HelixClass helixTPC;
      
      for (int iHit=0;iHit<nHits;++iHit) {
        TrackerHit * hit = hitVec[iHit];
        TrackerHitExtended * hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended( trackExt );
        trackExt->addTrackerHitExtended( hitExt );      
      }      
      
      
      float d0TPC = trackExt->getD0();
      float z0TPC = trackExt->getZ0();
      float omegaTPC = trackExt->getOmega();
      float phi0TPC = trackExt->getPhi();
      float tanLTPC = trackExt->getTanLambda();
      float Chi2TPC = trackExt->getChi2()/float(trackExt->getNDF());
      const int ndfTPC = trackExt->getNDF();
      
      helixTPC.Initialize_Canonical(phi0TPC,d0TPC,z0TPC,omegaTPC,tanLTPC,_bField);
      
      float pxTPC = helixTPC.getMomentum()[0];
      float pyTPC = helixTPC.getMomentum()[1];
      float pzTPC = helixTPC.getMomentum()[2];
      const float ptot = sqrt(pxTPC*pxTPC+pyTPC*pyTPC+pzTPC*pzTPC);
      sprintf(strg,"%3i  %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk,
              ptot, d0TPC,z0TPC,pxTPC,pyTPC,pzTPC,nHits,ndfTPC,Chi2TPC);
      streamlog_out(DEBUG3) << strg << std::endl;
      
      _allTPCTracks.push_back( trackExt );                
    }      
  }
  catch ( DataNotAvailableException &e) {
    streamlog_out(DEBUG4) << _TPCTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  // Reading Si Tracks
  try {
    LCCollection * col = event->getCollection(_SiTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG3) << std::endl;
    streamlog_out(DEBUG3) << "Number of Si Tracks = " << nelem << std::endl;
    streamlog_out(DEBUG3) << " Trk       p          D0         Z0       Px       Py       Pz   hitsSi ndf Chi2/ndf" << std::endl;
    
    for (int iTrk=0; iTrk<nelem; ++iTrk) {
      Track * siTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));
      TrackExtended * trackExt = new TrackExtended( siTrack );
      TrackerHitVec hitVec = siTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(siTrack->getOmega());
      trackExt->setTanLambda(siTrack->getTanLambda());
      trackExt->setPhi(siTrack->getPhi());
      trackExt->setD0(siTrack->getD0());
      trackExt->setZ0(siTrack->getZ0());
      float cov[15];
      float param[5];
      
      param[0] = siTrack->getOmega();
      param[1] = siTrack->getTanLambda();
      param[2] = siTrack->getPhi();
      param[3] = siTrack->getD0();
      param[4] = siTrack->getZ0();      
      
      
      
      const FloatVec Cov = siTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic=0;ic<NC;ic++) {
        cov[ic] =  Cov[ic];
      } 
      //      }
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(siTrack->getNdf());
      trackExt->setChi2(siTrack->getChi2());      
      char strg[200];
      HelixClass helixSi;
      for (int iHit=0;iHit<nHits;++iHit) {
        TrackerHit * hit = hitVec[iHit];
        TrackerHitExtended * hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended( trackExt );
        
        trackExt->addTrackerHitExtended( hitExt );      
      }
      
      
      float d0Si = trackExt->getD0();
      float z0Si = trackExt->getZ0();
      float omegaSi = trackExt->getOmega();
      float phi0Si = trackExt->getPhi();
      float tanLSi = trackExt->getTanLambda();
      helixSi.Initialize_Canonical(phi0Si,d0Si,z0Si,omegaSi,tanLSi,_bField);
      float pxSi = helixSi.getMomentum()[0];
      float pySi = helixSi.getMomentum()[1];
      float pzSi = helixSi.getMomentum()[2];
      const float pTot = sqrt(pxSi*pxSi+pySi*pySi+pzSi*pzSi);
      const int ndfSi = trackExt->getNDF();
      float Chi2Si = trackExt->getChi2()/float(trackExt->getNDF());
      sprintf(strg,"%3i  %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f",iTrk,
              pTot, d0Si,z0Si,pxSi,pySi,pzSi,nHits, ndfSi, Chi2Si);
      streamlog_out(DEBUG3) << strg << std::endl;
      
      if(nHits>0){
        _allSiTracks.push_back( trackExt );
      }else{
        delete trackExt;
      }
    }
    
    streamlog_out(DEBUG3) << std::endl;
  }
  catch ( DataNotAvailableException &e) {
    streamlog_out(DEBUG4) << _SiTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }
  
  
  
}

void FullLDCTracking_MarlinTrk::CleanUp(){
  
  int nNonCombTpc = int(_allNonCombinedTPCTracks.size());  
  for (int i=0;i<nNonCombTpc;++i) {
    TrackExtended * trkExt = _allNonCombinedTPCTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedTPCTracks.clear();
  
  int nNonCombSi = int(_allNonCombinedSiTracks.size());  
  for (int i=0;i<nNonCombSi;++i) {
    TrackExtended * trkExt = _allNonCombinedSiTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedSiTracks.clear();
  
  int nSITHits = int(_allSITHits.size());
  for (int i=0;i<nSITHits;++i) {
    TrackerHitExtended * hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();
  
  int nSETHits = int(_allSETHits.size());
  for (int i=0;i<nSETHits;++i) {
    TrackerHitExtended * hitExt = _allSETHits[i];
    delete hitExt;
  }
  _allSETHits.clear();
  
  int nTPCHits = int(_allTPCHits.size());
  for (int i=0;i<nTPCHits;++i) {
    TrackerHitExtended * hitExt = _allTPCHits[i];
    delete hitExt;
  }
  _allTPCHits.clear();
  
  int nFTDHits = int(_allFTDHits.size());
  for (int i=0;i<nFTDHits;++i) {
    TrackerHitExtended * hitExt = _allFTDHits[i];
    delete hitExt;
  }
  _allFTDHits.clear();
  
  int nETDHits = int(_allETDHits.size());
  for (int i=0;i<nETDHits;++i) {
    TrackerHitExtended * hitExt = _allETDHits[i];
    delete hitExt;
  }
  _allETDHits.clear();
  
  int nVTXHits = int(_allVTXHits.size());
  for (int i=0;i<nVTXHits;++i) {
    TrackerHitExtended * hitExt = _allVTXHits[i];
    delete hitExt;
  }
  _allVTXHits.clear();
  
  int nSiTrk = int(_allSiTracks.size());
  for (int i=0;i<nSiTrk;++i) {
    TrackExtended * trkExt = _allSiTracks[i];
    delete trkExt;
  }
  _allSiTracks.clear();
  
  int nTPCTrk = int(_allTPCTracks.size());
  for (int i=0;i<nTPCTrk;++i) {
    TrackExtended * trkExt = _allTPCTracks[i];
    delete trkExt;
  }
  _allTPCTracks.clear();
  
  int nCombTrk = int(_allCombinedTracks.size());
  for (int i=0;i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    delete group;
    delete trkExt;    
  }
  _allCombinedTracks.clear();
  
  //   int nImplTrk = int(_trkImplVec.size());
  //   for (int i=0;i<nImplTrk;++i) {
  //     TrackExtended * trkImpl = _trkImplVec[i];
  //     delete trkImpl;
  //   }
  _trkImplVec.clear();
  
  //AS: Dont delete the individual entries, some of them are cleared elsewhere, I think
  _candidateCombinedTracks.clear();
}

void FullLDCTracking_MarlinTrk::MergeTPCandSiTracks() {
  
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());
  
  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      int iComp = 0;
      float angle = 0;
      float dOmega = CompareTrkII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angle);
      if ( (dOmega<_dOmegaForMerging) && (angle<_angleForMerging) && !VetoMerge(siTrackExt,tpcTrackExt)) {
        TrackExtended *combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);       
        if (combinedTrack != NULL) {
          _allCombinedTracks.push_back( combinedTrack );
          _candidateCombinedTracks.insert(tpcTrackExt);
          _candidateCombinedTracks.insert(siTrackExt);
          if (_debug >= 3 ) {
            int iopt = 1;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }else{
          if (_debug >= 3 ) {
            int iopt = 6;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
      }
      else {
        if (_debug >= 3 ) {
          int iopt = 6;
          PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
        }
      }
    }
  }
  
  
}



void FullLDCTracking_MarlinTrk::MergeTPCandSiTracksII() {
  
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks  = int(_allSiTracks.size());
  
  for (int iTPC=0;iTPC<nTPCTracks;++iTPC) {
    TrackExtended * tpcTrackExt = _allTPCTracks[iTPC];
    if(_candidateCombinedTracks.find(tpcTrackExt) != _candidateCombinedTracks.end() )continue;
    
    for (int iSi=0;iSi<nSiTracks;++iSi) {
      
      TrackExtended * siTrackExt = _allSiTracks[iSi];
      if(_candidateCombinedTracks.find(siTrackExt)!= _candidateCombinedTracks.end() )continue;
      int iComp = 0;
      float angleSignificance = 0;
      float significance = CompareTrkIII(siTrackExt,tpcTrackExt,_d0CutForMerging,_z0CutForMerging,iComp,angleSignificance);
      if ( (significance<10) && (angleSignificance<5) && !VetoMerge(tpcTrackExt,siTrackExt) ) {
        TrackExtended * combinedTrack = CombineTracks(tpcTrackExt,siTrackExt);
        
        if (combinedTrack != NULL) {
          
          _allCombinedTracks.push_back( combinedTrack );
          if (_debug >= 3 ) {
            int iopt = 1;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }else{
          if (_debug >= 3 ) {
            int iopt = 6;
            PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
          }
        }
      }
      else {
        if (_debug >= 3 ) {
          int iopt = 6;
          PrintOutMerging(tpcTrackExt,siTrackExt,iopt);
        }
      }
    }
  }
}




TrackExtended * FullLDCTracking_MarlinTrk::CombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack) {
  
  TrackExtended * OutputTrack = NULL;
  
  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();
  
  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());
  int nHits = nTPCHits + nSiHits;
  
  //std::cout << "FullLDCTracking_MarlinTrk::CombineTracks nSiHits = " << nSiHits << std::endl;
  //std::cout << "FullLDCTracking_MarlinTrk::CombineTracks nTPCHits = " << nTPCHits << std::endl;
  
  EVENT::TrackerHitVec trkHits;
  trkHits.reserve(nHits);
  
  for (int ih=0;ih<nSiHits;++ih) {
    TrackerHit * trkHit = siHitVec[ih]->getTrackerHit();
    if(trkHit) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: TrackerHit pointer == NULL ")  ) ;
    }
  }  
  
  for (int ih=0;ih<nTPCHits;++ih) {

    TrackerHit * trkHit = tpcHitVec[ih]->getTrackerHit();
    if(trkHit) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: TrackerHit pointer == NULL ")  ) ;
    }
  }      
  
  double chi2_D;
  int ndf_D;
  
  if( trkHits.size() < 3 ) { 
    
    return NULL ;
    
  }
  
  
  MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
  
  sort(trkHits.begin(), trkHits.end(), FullLDCTracking_MarlinTrk::compare_r() );
  
  EVENT::TrackerHitVec::iterator it = trkHits.begin();
  
  streamlog_out(DEBUG2) << "Start Fitting: AddHits: number of hits to fit " << trkHits.size() << std::endl;
  
//  int number_of_added_hits = 0;
//  for( it = trkHits.begin() ; it != trkHits.end() ; ++it )
//      {
//    if (marlin_trk->addHit(*it) == 0){
//      ++number_of_added_hits;
//    }
//    else{
//      streamlog_out(DEBUG2) << "Hit " << it - trkHits.begin() << " Dropped " << std::endl;
//    }
//      }
//  
//  if( number_of_added_hits < 3 ) {
//    
//    delete marlin_trk ;
//    
//    return NULL;
//    
//  }
//  
//  streamlog_out(DEBUG2) << "Start Fitting: number_of_added_hits  = " << number_of_added_hits << std::endl;
//  
//  // SJA:FIXME: Here we could initialise the fit using the previous fits 
//  marlin_trk->initialise( IMarlinTrack::backward ) ;

  
  int number_of_added_hits = 0;
  int ndof_added = 0;
  TrackerHitVec added_hits;
  
  for( it = trkHits.begin() ; it != trkHits.end() ; ++it ) {
    
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
    else { // normal non composite hit
      
      if (marlin_trk->addHit( trkHit ) == 0) {
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
    streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Cannot fit less with less than 8 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << std::endl;
    delete marlin_trk ;
    return NULL;
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
  
  for (int icov = 0; icov<covMatrix.size(); ++icov) {
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
  
  
  int fit_status = marlin_trk->fit() ; 
  
  if( fit_status != 0 ) {
    
    delete marlin_trk ;
    
    return NULL;
    
  }
  
  
  const gear::Vector3D point(0.,0.,0.); // nominal IP
  int return_code = 0;
  
  TrackStateImpl trkState ;
  return_code = marlin_trk->propagate(point, trkState, chi2_D, ndf_D ) ;
  
  if( return_code != 0 ) {
    
    delete marlin_trk ;
    
    return NULL;
    
  }
  
  float chiQ = chi2_D/float(ndf_D);
  
  // SJA:FIXME remove hardcoded chiQ cut
  if (fit_status == 0 && chiQ > 0.001) {
    
    float chi2Fit = chi2_D/float(ndf_D);
    
    if (chi2Fit > _chi2FitCut) {
      
      delete marlin_trk ;
      
      return NULL;
      
    }    
    
    float omega = trkState.getOmega();
    float tanlambda = trkState.getTanLambda();
    float phi0 = trkState.getPhi();
    float d0 = trkState.getD0();
    float z0 = trkState.getZ0();    
    
    OutputTrack = new TrackExtended();
    GroupTracks * group = new GroupTracks();
    group->addTrackExtended(siTrack);
    group->addTrackExtended(tpcTrack);
    
    // note OutputTrack which is of type TrackExtended, only takes fits set for ref point = 0,0,0 
    OutputTrack->setGroupTracks(group);
    OutputTrack->setOmega(omega);
    OutputTrack->setTanLambda(tanlambda);
    OutputTrack->setPhi(phi0);
    OutputTrack->setZ0(z0);
    OutputTrack->setD0(d0);
    OutputTrack->setChi2(chi2_D);
    OutputTrack->setNDF(ndf_D);
    
    float cov[15];
    
    for (int i = 0 ; i<15 ; ++i) {
      cov[i] = trkState.getCovMatrix().operator[](i);
    }
    
    OutputTrack->setCovMatrix(cov);
    
    
    // SJA:FIXME: for now just assume all hits were used in the fit
    
    for (int i=0;i<nSiHits;++i) {
      TrackerHitExtended * hitExt = siHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(true);
    }
    for (int i=0;i<nTPCHits;++i) {
      TrackerHitExtended * hitExt = tpcHitVec[i];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(true);
    }
  }
  
  
  delete marlin_trk ;
  
  return OutputTrack;
}


TrackExtended * FullLDCTracking_MarlinTrk::TrialCombineTracks(TrackExtended * tpcTrack, TrackExtended * siTrack) {
  
  TrackExtended * OutputTrack = NULL;
  
  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();
  
  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());
  int nHits = nTPCHits + nSiHits;
  
  EVENT::TrackerHitVec trkHits;
  trkHits.reserve(nHits);
  
  for (int ih=0;ih<nSiHits;++ih) {
    TrackerHit * trkHit = siHitVec[ih]->getTrackerHit();
    if(trkHit) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::TrialCombineTracks: TrackerHit pointer == NULL ")  ) ;
    }
  }      
  for (int ih=0;ih<nTPCHits;++ih) {
    TrackerHit * trkHit = tpcHitVec[ih]->getTrackerHit();
    if(trkHit) { 
      trkHits.push_back(trkHit);   
    }
    else{
      throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::TrialCombineTracks: TrackerHit pointer == NULL ")  ) ;
    }
  }      
  
  double chi2_D;
  int ndf_D;
  
  if( trkHits.size() < 3 ) { 
    
    return NULL ;
    
  }
  
  MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
  
  sort(trkHits.begin(), trkHits.end(), FullLDCTracking_MarlinTrk::compare_r() );
  
  EVENT::TrackerHitVec::iterator it = trkHits.begin();
  
  streamlog_out(DEBUG2) << "Start Fitting: AddHits: number of hits to fit " << trkHits.size() << std::endl;
  
//  int number_of_added_hits = 0;
//  for( it = trkHits.begin() ; it != trkHits.end() ; ++it )
//      {
//    if (marlin_trk->addHit(*it) == 0){
//      ++number_of_added_hits;
//    }
//    else{
//      streamlog_out(DEBUG2) << "Hit " << it - trkHits.begin() << " Dropped " << std::endl;
//    }
//      }
//  
//  if( number_of_added_hits < 3 ) {
//    
//    delete marlin_trk ;
//    
//    return NULL;
//    
//  }
//  
//  // SJA:FIXME: Here we could initialise the fit using the previous fits 
//  marlin_trk->initialise( IMarlinTrack::backward ) ;
  
  
  int number_of_added_hits = 0;
  int ndof_added = 0;
  TrackerHitVec added_hits;
  
  for( it = trkHits.begin() ; it != trkHits.end() ; ++it ) {
    
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
    else { // normal non composite hit
      
      if (marlin_trk->addHit( trkHit ) == 0) {
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
    streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Cannot fit less with less than 8 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << std::endl;
    delete marlin_trk ;
    return NULL;
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
  
  for (int icov = 0; icov<covMatrix.size(); ++icov) {
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

  
  int fit_status = marlin_trk->fit() ; 
  
  if( fit_status != 0 ) {
    
    delete marlin_trk ;
    
    return NULL;
    
  }
  
  
  const gear::Vector3D point(0.,0.,0.); // nominal IP
  int return_code = 0;
  
  TrackStateImpl trkState ;
  return_code = marlin_trk->propagate(point, trkState, chi2_D, ndf_D ) ;
  
  if( return_code != 0 ) {
    
    delete marlin_trk ;
    
    return NULL;
    
  }
  
  float chiQ = chi2_D/float(ndf_D);
  
  
  // SJA:FIXME remove hard coded chiQ cut 
  if (fit_status == 0 && chiQ > 0.001) {
    
    float chi2Fit = chi2_D/float(ndf_D);
    
    if (chi2Fit > _chi2FitCut) {
      
      delete marlin_trk ;
      
      return NULL;
      
    }    
    
    float omega = trkState.getOmega();
    float tanlambda = trkState.getTanLambda();
    float phi0 = trkState.getPhi();
    float d0 = trkState.getD0();
    float z0 = trkState.getZ0();    
    
    OutputTrack = new TrackExtended();
    GroupTracks * group = new GroupTracks();
    group->addTrackExtended(siTrack);
    group->addTrackExtended(tpcTrack);
    
    // note OutputTrack which is of type TrackExtended, only takes fits set for ref point = 0,0,0    
    OutputTrack->setGroupTracks(group);
    OutputTrack->setOmega(omega);
    OutputTrack->setTanLambda(tanlambda);
    OutputTrack->setPhi(phi0);
    OutputTrack->setZ0(z0);
    OutputTrack->setD0(d0);
    OutputTrack->setChi2(chi2_D);
    OutputTrack->setNDF(ndf_D);
    
    float cov[15];
    
    for (int i = 0 ; i<15 ; ++i) {
      cov[i] = trkState.getCovMatrix().operator[](i);
    }
    
    OutputTrack->setCovMatrix(cov);
    
  }
  
  delete marlin_trk ;
  
  return OutputTrack;
}


void FullLDCTracking_MarlinTrk::SortingTrackHitPairs(TrackHitPairVec & trackHitPairVec) {
  
  int sizeOfVector = int(trackHitPairVec.size());
  TrackHitPair *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
      one = trackHitPairVec[j];
      two = trackHitPairVec[j+1];
      float oneQ = one->getDistance();
      float twoQ = two->getDistance();
      if( oneQ > twoQ ) {
        Temp = trackHitPairVec[j];
        trackHitPairVec[j] = trackHitPairVec[j+1];
        trackHitPairVec[j+1] = Temp;
      }
    }  
  
  
}

void FullLDCTracking_MarlinTrk::Sorting(TrackExtendedVec & trackVec) {
  
  int sizeOfVector = int(trackVec.size());
  TrackExtended *one,*two,*Temp;
  
  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
      one = trackVec[j];
      two = trackVec[j+1];
      float oneQ = one->getChi2()/float(one->getNDF());
      float twoQ = two->getChi2()/float(two->getNDF());
      if( oneQ > twoQ ) {
        Temp = trackVec[j];
        trackVec[j] = trackVec[j+1];
        trackVec[j+1] = Temp;
      }
    }  
}

void FullLDCTracking_MarlinTrk::SelectCombinedTracks() {
  
  int nCombTrk = int(_allCombinedTracks.size());
  
  for (int i=0; i<nCombTrk;++i) {
    TrackExtended * trkExt = _allCombinedTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    TrackExtendedVec tracks = group->getTrackExtendedVec();
    int nTracks = int(tracks.size());
    if (nTracks == 2) {
      TrackExtended * firstTrack = tracks[0];
      TrackExtended * secondTrack = tracks[1];
      if ((firstTrack->getGroupTracks() == NULL) &&
          (secondTrack->getGroupTracks() == NULL) ) {
        firstTrack->setGroupTracks(group);
        secondTrack->setGroupTracks(group);     
        TrackerHitExtendedVec firstVec = firstTrack->getTrackerHitExtendedVec();
        TrackerHitExtendedVec secondVec = secondTrack->getTrackerHitExtendedVec();
        int nFirst = int(firstVec.size());
        int nSecond = int(secondVec.size());
        float edges[2];
        edges[0] = 1.0e+20;
        edges[1] = -1.0e+20;
        for (int iF=0;iF<nFirst;++iF) {
          TrackerHitExtended * trkHitExt = firstVec[iF];
          TrackerHit * trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit->getPosition()[2]);
          if (zpos>edges[1])
            edges[1] = zpos;
          if (zpos<edges[0])
            edges[0] = zpos;      
        }
        for (int iS=0;iS<nSecond;++iS) {
          TrackerHitExtended * trkHitExt = secondVec[iS];
          TrackerHit * trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit->getPosition()[2]);
          if (zpos>edges[1])
            edges[1] = zpos;
          if (zpos<edges[0])
            edges[0] = zpos;
        }
        group->setEdges(edges);
        _trkImplVec.push_back(trkExt);  
        
        if (_debug >= 3) {
          int iopt = 1;
          PrintOutMerging(secondTrack,firstTrack,iopt);
        }       
      }
    }else{
      if(nTracks>2) streamlog_out(DEBUG3) << " MORE THAN TWO TRACKS " << std::endl;
    }
  }
  
  
}

void FullLDCTracking_MarlinTrk::AddNotCombinedTracks() {  
  
  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());
  
  // we need some buffer vector
  TrackExtendedVec allMergedTracks;
  allMergedTracks.clear();
  
  // forcing merging of Si and TPC track segments
  if (_forceMerging==1) { 
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExtTPC = _allTPCTracks[i];
      GroupTracks * groupTPC = trkExtTPC->getGroupTracks();
      if (groupTPC == NULL) {
        float diffMin = 1.0e+20;  
        TrackExtended * siTrkToAttach = NULL;
        for (int j=0;j<nSiTrk;++j) {
          TrackExtended * trkExtSi = _allSiTracks[j];
          GroupTracks * groupSi = trkExtSi->getGroupTracks();
          if (groupSi == NULL) {
            int iComp = 0;
            //      float deltaP = CompareTrk(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp);
            float angle(0.);
            float angleSignificance(0.);
            
            float dOmega = CompareTrkII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angle);
            float significance = CompareTrkIII(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp,angleSignificance);
            //      if (deltaP < _dPCutForForcedMerging) {
            if ( ((dOmega<_dOmegaForForcedMerging) && (angle<_angleForForcedMerging)) ||
                ((significance<5)                 && (angleSignificance<5))
                ) {
              float chi2O = dOmega/_dOmegaForForcedMerging;
              float chi2A = angle/_angleForForcedMerging;
              float deltaP = chi2O*chi2O + chi2A*chi2A; 
              if (deltaP<diffMin) {
                diffMin = deltaP;
                siTrkToAttach = trkExtSi;
              }
            }else{
              if (_debug==3) {
                int  iopt = 7;
                streamlog_out(DEBUG2) << significance << " " << angleSignificance << std::endl;
                PrintOutMerging(trkExtTPC,trkExtSi,iopt);
              }
            }
          }
        }
        
        if (siTrkToAttach!=NULL) {
          TrackExtended * trkExtSi = siTrkToAttach; 
          TrackExtended * OutputTrack = new TrackExtended();
          GroupTracks * group = new GroupTracks();
          group->addTrackExtended(trkExtSi);
          group->addTrackExtended(trkExtTPC);
          OutputTrack->setGroupTracks(group);
          //        trkExtSi->setGroupTracks(group);
          //        trkExtTPC->setGroupTracks(group);         
          OutputTrack->setOmega(trkExtTPC->getOmega());
          OutputTrack->setTanLambda(trkExtSi->getTanLambda());
          OutputTrack->setPhi(trkExtSi->getPhi());
          OutputTrack->setZ0(trkExtSi->getZ0());
          OutputTrack->setD0(trkExtSi->getD0());
          float covMatTPC[15];
          float covMatSi[15];
          float covMat[15];
          for (int iCov=0;iCov<15;++iCov) {
            covMatTPC[iCov] = trkExtTPC->getCovMatrix()[iCov];
            covMatSi[iCov] = trkExtSi->getCovMatrix()[iCov];                
            covMat[iCov] = covMatSi[iCov];
          }
          float scaling = sqrt(covMatTPC[5]/covMatSi[5]);
          covMat[5] = covMatTPC[5];
          covMat[3] = scaling*covMatSi[3];
          covMat[4] = scaling*covMatSi[4];
          covMat[8] = scaling*covMatSi[8];
          covMat[12] = scaling*covMatSi[12];          
          OutputTrack->setCovMatrix(covMat);
          TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
          TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();              
          int nTPCHits = int( tpcHitVec.size());
          int nSiHits = int( siHitVec.size());        
          float edges[2];
          edges[0] = 1.0e+20;
          edges[1] = -1.0e+20;
          for (int iH=0;iH<nSiHits;++iH) {
            TrackerHitExtended * hitExt = siHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true);
            TrackerHit * hit = hitExt->getTrackerHit();
            float zpos = float(hit->getPosition()[2]);
            if (zpos<edges[0])
              edges[0] = zpos;
            if (zpos>edges[1])
              edges[1] = zpos;
          }       
          for (int iH=0;iH<nTPCHits;++iH) {
            TrackerHitExtended * hitExt = tpcHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true); 
            TrackerHit * hit = hitExt->getTrackerHit();
            float zpos = float(hit->getPosition()[2]);
            if (zpos<edges[0])
              edges[0] = zpos;
            if (zpos>edges[1])
              edges[1] = zpos;
          }
          group->setEdges(edges);
          OutputTrack->setChi2(diffMin); // will be replaced if necessary
          OutputTrack->setNDF(int(1));   // will be replaced if necessary
          _allCombinedTracks.push_back( OutputTrack );
          allMergedTracks.push_back( OutputTrack );
        }           
      }
    }
    
    
    int nMerged = int(allMergedTracks.size());
    if (nMerged>0) {
      Sorting(allMergedTracks); 
      for (int iM=0;iM<nMerged;++iM) {
        TrackExtended * mergedTrack = allMergedTracks[iM];
        GroupTracks * grpTrk = mergedTrack->getGroupTracks();
        TrackExtendedVec trkVec = grpTrk->getTrackExtendedVec();
        TrackExtended * trkTPC = NULL;
        TrackExtended * trkSi = NULL;
        int nT = int(trkVec.size());
        if (nT==2) {
          trkTPC = trkVec[0];
          trkSi = trkVec[1];
          GroupTracks * groupTPC = trkTPC->getGroupTracks();
          GroupTracks * groupSi  = trkSi->getGroupTracks();
          if (groupTPC == NULL && groupSi == NULL) {
            trkTPC->setGroupTracks( grpTrk );
            trkSi->setGroupTracks( grpTrk );
            TrackerHitExtendedVec hitVec = mergedTrack->getTrackerHitExtendedVec();
            int nhits = int(hitVec.size());
            int totNdf = 2*nhits - 5;
            float totChi2 = trkTPC->getChi2() + trkSi->getChi2();
            mergedTrack->setNDF( totNdf );
            mergedTrack->setChi2( totChi2 );
            if (_debug >= 3) {
              int iopt = 2;
              PrintOutMerging(trkTPC,trkSi,iopt);
            }
            _trkImplVec.push_back( mergedTrack );
          }
        }
      }
    }
  }
  
  
  
  // clear buffer vector
  allMergedTracks.clear();
  
  // merging splitted TPC segments
  if (_mergeTPCSegments) {
    std::vector<GroupTracks*> TPCSegments;
    TPCSegments.clear();
    int nNonAssignedTPCSeg = 0;
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
        TrackerHitExtendedVec currentVec = trkExt->getTrackerHitExtendedVec();
        int nCur = int(currentVec.size());
        float zmin = 1e+20;
        float zmax = -1e+20;
        for (int iCur=0;iCur<nCur;++iCur) {
          TrackerHitExtended * curTrkHitExt = currentVec[iCur];
          TrackerHit * curTrkHit = curTrkHitExt->getTrackerHit();
          float zpos = float(curTrkHit->getPosition()[2]);
          if (zpos < zmin)
            zmin = zpos;
          if (zpos > zmax)
            zmax = zpos;
        }
        nNonAssignedTPCSeg++;
        int nGroups = int(TPCSegments.size());
        float dPtMin = 1.0e+10;
        GroupTracks * groupToAttach = NULL;
        TrackExtended * trkToAttach = NULL;
        for (int iG=0;iG<nGroups;++iG) {
          GroupTracks * segments = TPCSegments[iG];
          TrackExtendedVec segVec = segments->getTrackExtendedVec();
          int nTrk = int(segVec.size());
          bool consider = true;
          if (_forbidOverlapInZTPC==1) { // if overlap in Z of the two segments is forbidden
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trkInGroup = segVec[iTrk];
              TrackerHitExtendedVec hitInGroupVec = trkInGroup->getTrackerHitExtendedVec();
              int nHitsInGrp = int(hitInGroupVec.size());
              for (int iHitInGrp=0;iHitInGrp<nHitsInGrp;iHitInGrp++) {
                TrackerHitExtended * xTrkExt = hitInGroupVec[iHitInGrp];
                TrackerHit * xTrk = xTrkExt->getTrackerHit();
                float xZ = float(xTrk->getPosition()[2]);
                if (xZ>zmin&&xZ<zmax) {
                  consider = false;
                  break;
                }
              }
              if (!consider)
                break;
            }
          }
          if (consider) {
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trkInGroup = segVec[iTrk];
              int iComp = 1;
              float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);
              if (dPt < dPtMin && !VetoMerge(trkExt,trkInGroup)) {
                dPtMin = dPt;
                groupToAttach = segments;
                trkToAttach = trkInGroup;
                if (_debug==3) {
                  int iopt = 6;
                  PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
              else {
                if (_debug==3) {
                  int iopt = 9;
                  PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
            }
          }
          else {
            if (_debug >= 3) {
              int iopt = 9;
              for (int iTrk=0;iTrk<nTrk;++iTrk) {
                TrackExtended * trkInGroup = segVec[iTrk];
                int iComp = 1;
                float dPt = CompareTrk(trkExt,trkInGroup,_d0CutToMergeTPC,_z0CutToMergeTPC,iComp);              
                if (dPt >= dPtMin) {          
                  PrintOutMerging(trkExt,trkInGroup,iopt);
                }
              }
            }
          }
        }
        if (dPtMin < _dPCutToMergeTPC && groupToAttach != NULL) {
          
          
          
          groupToAttach->addTrackExtended(trkExt);
          trkExt->setGroupTracks(groupToAttach);
          float zminGroup = groupToAttach->getEdges()[0]; 
          float zmaxGroup = groupToAttach->getEdges()[1];
          float edges[2];
          edges[0] = zmin;
          if (zminGroup<zmin)
            edges[0] = zminGroup;
          edges[1] = zmax;
          if (zmaxGroup>zmax)
            edges[1] = zmaxGroup;
          groupToAttach->setEdges(edges);
          if (_debug==3) {
            int iopt = 3;
            PrintOutMerging(trkExt,trkToAttach,iopt);
          }
        }
        else {
          GroupTracks * newSegment = new GroupTracks(trkExt);
          trkExt->setGroupTracks(newSegment);
          TPCSegments.push_back(newSegment);
          float edges[2];
          edges[0] = zmin;
          edges[1] = zmax;
          newSegment->setEdges(edges);
        }
      }
    }
    
    // combining splitted TPC segments with the reconstructed tracks having Si hits
    int nCombTrk = int(_trkImplVec.size());
    int nSegments = int(TPCSegments.size());
    //    std::cout << "Combined tracks = " << nCombTrk << std::endl;
    //    std::cout << "nSegments = " << nSegments << std::endl;
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      TrackExtendedVec segVec = segments->getTrackExtendedVec();
      float zminTPCSeg = segments->getEdges()[0];
      float zmaxTPCSeg = segments->getEdges()[1];
      int nTrk = int(segVec.size());
      TrackExtended * CombTrkToAttach = NULL;
      TrackExtended * keyTrack = NULL;
      float deltaPtMin = _dPCutToMergeTPC;
      for (int iCTrk=0;iCTrk<nCombTrk;++iCTrk) {
        TrackExtended * combTrk = _trkImplVec[iCTrk];
        GroupTracks * groupComb = combTrk->getGroupTracks();
        bool consider = true;
        if (_forbidOverlapInZComb==1) { // if overlap in Z of the two segments is forbidden
          float zminComb = groupComb->getEdges()[0];
          float zmaxComb = groupComb->getEdges()[1];
          consider = (zminTPCSeg>zmaxComb) || (zmaxTPCSeg<zminComb);
        }
        if (consider) {
          for (int iTrk=0;iTrk<nTrk;++iTrk) {
            TrackExtended * trk = segVec[iTrk];
            int iopt = 0;
            float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
            float angleSignificance(0.);
            float significance = CompareTrkIII(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt,angleSignificance);
            if ( (dPt<deltaPtMin || significance <5 ) ) {
              if(VetoMerge(trk,combTrk)==false){
                CombTrkToAttach = combTrk;
                keyTrack = trk;
                deltaPtMin = dPt;
              }
            }
            else {
              if (_debug==3) {
                GroupTracks * groupCur = combTrk->getGroupTracks();
                TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt = 8;
                PrintOutMerging(trk,dummySi,iopt);
              }
            }
          }
        }
        else {
          if (_debug==3) {
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trk = segVec[iTrk];
              int iopt = 0;
              float dPt = CompareTrk(trk,combTrk,_d0CutToMergeTPC,_z0CutToMergeTPC,iopt);
              if (dPt>deltaPtMin) {
                GroupTracks * groupCur = combTrk->getGroupTracks();
                TrackExtended * dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt = 8;
                PrintOutMerging(trk,dummySi,iopt);
              }
            }
          }
        }
      }
      
      if (CombTrkToAttach != NULL) { // attach TPC segment to existing Comb Track
        GroupTracks * groupToAttach = CombTrkToAttach->getGroupTracks();          
        TrackExtended * SiCombTrk = groupToAttach->getTrackExtendedVec()[0];
        TrackExtended * TpcCombTrk = groupToAttach->getTrackExtendedVec()[1];
        if (_debug==3) {
          int iopt = 4;
          PrintOutMerging(keyTrack,SiCombTrk,iopt);
          iopt = 5;
          PrintOutMerging(keyTrack,TpcCombTrk,iopt);      
        }
        for (int iTrk=0;iTrk<nTrk;iTrk++) {
          TrackExtended * trk = segVec[iTrk];
          groupToAttach->addTrackExtended( trk );
          trk->setGroupTracks( groupToAttach );
          TrackerHitExtendedVec hitVec = trk->getTrackerHitExtendedVec();
          int nHitS = int(hitVec.size());              
          for (int iHS=0;iHS<nHitS;++iHS) {
            TrackerHitExtended * hitExt = hitVec[iHS];
            hitExt->setUsedInFit(false);
            CombTrkToAttach->addTrackerHitExtended( hitExt );
          }
        }
      }
      else {
        if (nTrk==1) { //
          GroupTracks * newGrp = new GroupTracks();
          segVec[0]->setGroupTracks(newGrp);
          newGrp->addTrackExtended(segVec[0]);
          TrackerHitExtendedVec TpcHitVec = segVec[0]->getTrackerHitExtendedVec();
          int nTpcH = int(TpcHitVec.size());
          for (int iTpcH=0;iTpcH<nTpcH;++iTpcH) {
            TpcHitVec[iTpcH]->setUsedInFit( true );
          }
          _trkImplVec.push_back(segVec[0]);
          _allNonCombinedTPCTracks.push_back(segVec[0]);
        }
        else {
          float zMin = 1.0e+20;
          TrackExtended * chosenTrack = NULL;
          for (int iTrk=0;iTrk<nTrk;++iTrk) {
            TrackExtended * trk = segVec[iTrk];
            Track * track = trk->getTrack();
            TrackerHitVec hitVec = track->getTrackerHits();
            int nHits = int(hitVec.size());
            for (int iH=0;iH<nHits;++iH) {
              float zPosi = fabs(hitVec[iH]->getPosition()[2]);
              if (zPosi<zMin) {
                chosenTrack = trk;
                zMin = zPosi;
                break;
              }
            }
          }
          if (chosenTrack!=NULL) {
            GroupTracks * newGroup = new GroupTracks();
            for (int iTrk=0;iTrk<nTrk;++iTrk) {
              TrackExtended * trk = segVec[iTrk];
              trk->setGroupTracks( newGroup );
              newGroup->addTrackExtended( trk );                        
              TrackerHitExtendedVec hitVecS = trk->getTrackerHitExtendedVec();
              int nHitS = int(hitVecS.size());                  
              for (int iH=0;iH<nHitS;++iH) {
                TrackerHitExtended * trkHitExt = hitVecS[iH];
                if (trk!=chosenTrack) {
                  trkHitExt->setUsedInFit( false );
                  chosenTrack->addTrackerHitExtended( trkHitExt );                              
                }
                else {
                  trkHitExt->setUsedInFit( true );
                }
              }
            }
            _allNonCombinedTPCTracks.push_back(chosenTrack);
            _trkImplVec.push_back(chosenTrack);
          }
        }
      }
    }
    for (int iS=0;iS<nSegments;++iS) {
      GroupTracks * segments = TPCSegments[iS];
      delete segments;
    }
    TPCSegments.clear();
  }
  else { // adding all TPC segments to the list of tracks (track splitting is allowed)
    for (int i=0;i<nTPCTrk;++i) {
      TrackExtended * trkExt = _allTPCTracks[i];
      Track * track = trkExt->getTrack();
      GroupTracks * group = trkExt->getGroupTracks();
      if (group == NULL) {
        TrackerHitVec hitVec = track->getTrackerHits();       
        _trkImplVec.push_back(trkExt);
        _allNonCombinedTPCTracks.push_back( trkExt );
        GroupTracks * newGrp = new GroupTracks();
        newGrp->addTrackExtended( trkExt );
        trkExt->setGroupTracks( newGrp );
      }
    }    
  }
  
  for (int i=0;i<nSiTrk;++i) { // adding left-over Si segments to the list of tracks
    TrackExtended * trkExt = _allSiTracks[i];
    GroupTracks * group = trkExt->getGroupTracks();
    if (group == NULL) {
      TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
      int nHSi = int(hitVec.size());
      for (int iHSi=0;iHSi<nHSi;++iHSi) {
        hitVec[iHSi]->setUsedInFit(true);
      }
      _trkImplVec.push_back(trkExt);
      GroupTracks * newGrp = new GroupTracks();
      newGrp->addTrackExtended( trkExt );
      trkExt->setGroupTracks( newGrp );   
      _allNonCombinedSiTracks.push_back( trkExt );
    }
  }
  
}

void FullLDCTracking_MarlinTrk::CheckTracks() {  
  
  for(unsigned int i = 0; i< _trkImplVec.size();i++){
    TrackExtended *first = _trkImplVec[i];
    if(first==NULL)continue;
    float d0First = first->getD0();
    float z0First = first->getZ0();
    float omegaFirst = first->getOmega();
    float tanLFirst = first->getTanLambda();
    float phiFirst = first->getPhi();
    HelixClass helixFirst;
    helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float momFirst[3];
    momFirst[0]= helixFirst.getMomentum()[0];
    momFirst[1]= helixFirst.getMomentum()[1];
    momFirst[2]= helixFirst.getMomentum()[2];
    float pFirst    = sqrt(momFirst[0]*momFirst[0]+momFirst[1]*momFirst[1]+momFirst[2]*momFirst[2]);
    if(std::isnan(pFirst))continue;
    TrackerHitExtendedVec firstHitVec  = first->getTrackerHitExtendedVec();
    if(firstHitVec.size()<1)continue;
    
    for(unsigned int j = i+1; j<_trkImplVec.size();j++){
      TrackExtended *second = _trkImplVec[j];
      if(second==NULL)continue;
      float d0Second = second->getD0();
      float z0Second = second->getZ0();
      float omegaSecond = second->getOmega();
      float tanLSecond = second->getTanLambda();
      float phiSecond = second->getPhi();
      HelixClass helixSecond;
      helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
      float momSecond[3];
      momSecond[0] = helixSecond.getMomentum()[0];
      momSecond[1] = helixSecond.getMomentum()[1];
      momSecond[2] = helixSecond.getMomentum()[2];
      float pSecond    = sqrt(momSecond[0]*momSecond[0]+momSecond[1]*momSecond[1]+momSecond[2]*momSecond[2]);
      if(std::isnan(pSecond))continue;
      TrackerHitExtendedVec secondHitVec  = second->getTrackerHitExtendedVec();
      if(secondHitVec.size()<1)continue;
      if(firstHitVec.size()+secondHitVec.size()<10)continue;
      
      
      float pdot = (momFirst[0]*momSecond[0]+momFirst[1]*momSecond[1]+momFirst[2]*momSecond[2])/pFirst/pSecond;
      if(pdot<0.999)continue;
      // const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      // const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      // const float deltaP = fabs(pFirst-pSecond);
      // const float sigmaDeltaP = sqrt(pFirst*sigmaPOverPFirst*pFirst*sigmaPOverPFirst+pSecond*sigmaPOverPSecond*pSecond*sigmaPOverPSecond);
      //      const float significance = deltaP/sigmaDeltaP;
      
      TrackExtended * combinedTrack = TrialCombineTracks(first,second);
      if(combinedTrack != NULL){
        const int minHits = std::min(firstHitVec.size(),secondHitVec.size());
        const int maxHits = std::max(firstHitVec.size(),secondHitVec.size());
        
        if( combinedTrack->getNDF() <= 2*maxHits+minHits-5){
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
          continue;
        }
        
        float d0    = combinedTrack->getD0();
        float z0    = combinedTrack->getZ0();
        float omega = combinedTrack->getOmega();
        float tanL  = combinedTrack->getTanLambda();
        float phi   = combinedTrack->getPhi();
        
        HelixClass helix;
        helix.Initialize_Canonical(phi,d0,z0,omega,tanL,_bField);
        float mom[3];
        mom[0]  = helix.getMomentum()[0];
        mom[1]  = helix.getMomentum()[1];
        mom[2]  = helix.getMomentum()[2];
        // float p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
        // float chi2Sig =  (combinedTrack->getChi2() - combinedTrack->getNDF());
        // chi2Sig = chi2Sig/sqrt(combinedTrack->getNDF()*2);
        
        
        
        int nTpcFirst(0);
        int nUsedFirst(0);
        for(unsigned int ihit = 0;ihit<firstHitVec.size();ihit++){
          
          if( getDetectorID(firstHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC) nTpcFirst++;
          
          if(firstHitVec[ihit]->getUsedInFit()==true)nUsedFirst++;
        }
        
        
        int nTpcSecond(0);
        int nUsedSecond(0);
        for(unsigned int ihit = 0;ihit<secondHitVec.size();ihit++){
          if( getDetectorID(secondHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC) ++nTpcSecond;
          if( secondHitVec[ihit]->getUsedInFit()==true ) ++nUsedSecond;
        }
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }
    }
  }
  
  
}





float FullLDCTracking_MarlinTrk::CompareTrkII(TrackExtended * first, TrackExtended * second, 
                                              float d0Cut, float z0Cut,int iopt,float & Angle) {
  
  
  float result = 1.0e+20;
  Angle  = 1.0e+20; 
  float omegaFirst = first->getOmega();
  float omegaSecond = second->getOmega();
  float deltaOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
  if(deltaOmega> 2*_dOmegaForMerging)return result;
  
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  if(!isCloseInIP)return result;
  
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  float qSecond = PIOVER2 - atan(tanLSecond);
  
  Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
  sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  Angle = acos(Angle);
  
  result = deltaOmega;
  
  return result;
  
}




float FullLDCTracking_MarlinTrk::CompareTrkIII(TrackExtended * first, TrackExtended * second, 
                                               float d0Cut, float z0Cut,int iopt, float & AngleSignificance) {
  
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qSecond = PIOVER2 - atan(tanLSecond);
  
  
  //MB 2010 03
  float d0ErrFirst = sqrt(first->getCovMatrix()[0]);
  float z0ErrFirst = sqrt(first->getCovMatrix()[9]);
  //  float omegaErrFirst = sqrt(first->getCovMatrix()[5]);
  float phiErrFirst = sqrt(first->getCovMatrix()[2]);
  float qErrFirst = sqrt(cos(qFirst)*cos(qFirst)*first->getCovMatrix()[14]);
  //MB END
  //MB 2010 03
  float d0ErrSecond = sqrt(second->getCovMatrix()[0]);
  float z0ErrSecond = sqrt(second->getCovMatrix()[9]);
  //  float omegaErrSecond = sqrt(second->getCovMatrix()[5]);
  float phiErrSecond = sqrt(second->getCovMatrix()[2]);
  float qErrSecond = sqrt(cos(qSecond)*cos(qSecond)*second->getCovMatrix()[14]);
  //MB END
  
  
  //  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  //isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  //MB 2010 03
  bool isCloseInIP = (fabs(d0First-d0Second)/sqrt(d0ErrFirst*d0ErrFirst+d0ErrSecond*d0ErrSecond)<d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)/sqrt(z0ErrFirst*z0ErrFirst+z0ErrSecond*z0ErrSecond)<z0Cut);
  
  if (!isCloseInIP)return result;
  
  float Angle = (cos(phiFirst)*cos(phiSecond)+sin(phiFirst)*sin(phiSecond))*
  sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float momFirst = 0;
  float momSecond = 0;
  
  for (int iC=0;iC<3;++iC) {
    pFirst[iC] = helixFirst.getMomentum()[iC];
    pSecond[iC] = helixSecond.getMomentum()[iC];
    momFirst += pFirst[iC]* pFirst[iC];
    momSecond += pSecond[iC]*pSecond[iC];
  }
  momFirst = sqrt(momFirst);
  momSecond = sqrt(momSecond);
  
  
  float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
  
  
  const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
  const float deltaP = fabs(momFirst-momSecond);
  const float sigmaPFirst = momFirst*sigmaPOverPFirst;
  const float sigmaPSecond = momSecond*sigmaPOverPSecond;
  const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
  const float significance = deltaP/sigmaDeltaP;
  
  //MB 2010 03
  float errorAngle =sin(phiFirst)*sin(phiFirst)*phiErrFirst*phiErrFirst*cos(phiSecond)*cos(phiSecond)+
  sin(phiSecond)*sin(phiSecond)*phiErrSecond*phiErrSecond*cos(phiFirst)*cos(phiFirst)+
  sin(qFirst)*sin(qFirst)*qErrFirst*qErrFirst*cos(qSecond)*cos(qSecond)+
  sin(qSecond)*sin(qSecond)*qErrSecond*qErrSecond*cos(qFirst)*cos(qFirst)+
  cos(phiFirst)*cos(phiFirst)*phiErrFirst*phiErrFirst*(sin(phiSecond)*sin(qFirst)*sin(qSecond))*(sin(phiSecond)*sin(qFirst)*sin(qSecond))+
  cos(phiSecond)*cos(phiSecond)*phiErrSecond*phiErrSecond*(sin(phiFirst)*sin(qFirst)*sin(qSecond))*(sin(phiFirst)*sin(qFirst)*sin(qSecond))+
  cos(qFirst)*cos(qFirst)*qErrFirst*qErrFirst*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))*(sin(phiFirst)*sin(phiSecond)*sin(qSecond))+
  cos(qSecond)*cos(qSecond)*qErrSecond*qErrSecond*(sin(phiFirst)*sin(phiSecond)*sin(qFirst))*(sin(phiFirst)*sin(phiSecond)*sin(qFirst));
  
  if(Angle<1.){
    errorAngle = sqrt(1./(1.-Angle*Angle)*errorAngle);
  }else{
    errorAngle = sqrt(errorAngle);
  }
  
  if(errorAngle<1.e-6)errorAngle=1.e-6;
  
  AngleSignificance = fabs(acos(Angle)/errorAngle);
  if(pdot<0.999)return result;
  
  return significance;
  
}


float FullLDCTracking_MarlinTrk::CompareTrk(TrackExtended * first, TrackExtended * second, 
                                            float d0Cut, float z0Cut,int iopt) {
  
  float result = 1.0e+20;
  
  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  
  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  
  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  
  if (iopt>0) isCloseInIP = isCloseInIP || (fabs(d0First+d0Second)<d0Cut);
  
  isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);
  
  
  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst,d0First,z0First,omegaFirst,tanLFirst,_bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  
  float pFirst[3];
  float pSecond[3];
  float dPminus[3];
  float dPplus[3];
  float momFirst = 0;
  float momSecond = 0;
  float momMinus = 0;
  float momPlus = 0;
  
  if ( isCloseInIP ) {
    
    for (int iC=0;iC<3;++iC) {
      pFirst[iC] = helixFirst.getMomentum()[iC];
      pSecond[iC] = helixSecond.getMomentum()[iC];
      momFirst += pFirst[iC]* pFirst[iC];
      momSecond += pSecond[iC]*pSecond[iC];
      dPminus[iC] = pFirst[iC] - pSecond[iC];
      dPplus[iC] = pFirst[iC] + pSecond[iC];
      momMinus += dPminus[iC]*dPminus[iC];
      momPlus += dPplus[iC]*dPplus[iC];
    }
    momFirst = sqrt(momFirst);
    momSecond = sqrt(momSecond);
    
    float ptFirst = sqrt(pFirst[0]*pFirst[0]+pFirst[1]*pFirst[1]);
    float ptSecond = sqrt(pSecond[0]*pSecond[0]+pSecond[1]*pSecond[1]);
    
    
    if ( (ptFirst<_PtCutToMergeTPC) && (ptSecond<_PtCutToMergeTPC) ) {
      
      momMinus = sqrt(momMinus);
      momPlus = sqrt(momPlus);
      float nom = momMinus;
      if (momPlus<nom && iopt>0)
        nom = momPlus;
      float den = momFirst;
      if (momSecond<momFirst)
        den = momSecond;
      
      result = nom/den;     
      
    }
    
    else {
      
      
      // check for cases where PatRec splits non-looping TPC tracks 
      // look for two tracks where total tpc hits are not more than total number
      // of pad rows and that the hits on one track are close to the helix of the
      // other track
      
      float dpOverP = 2.0*fabs(momFirst-momSecond)/(momFirst+momSecond);
      const float pdot = (pFirst[0]*pSecond[0]+pFirst[1]*pSecond[1]+pFirst[2]*pSecond[2])/momFirst/momSecond;
      const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      const float deltaP = fabs(momFirst-momSecond);
      const float sigmaPFirst = momFirst*sigmaPOverPFirst;
      const float sigmaPSecond = momSecond*sigmaPOverPSecond;
      const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
      const float significance = deltaP/sigmaDeltaP;
      
      
      //compare angle between the two vectors (cos theta) and their momentum
      //    if( ( pdot>0.99 && dpOverP<0.01 ) || ( pdot>0.998 && dpOverP<0.25 ) ){
      if( (pdot>_cosThetaCutHighPtMerge && dpOverP<_momDiffCutHighPtMerge) 
         || 
         (pdot>_cosThetaCutSoftHighPtMerge && dpOverP<_momDiffCutSoftHighPtMerge) 
         || (pdot > 0.9999 && significance <10) 
         ){
        
        
        int nTrkGrpFirst = 0;
        int nTrkGrpSecond = 0;
        TrackerHitVec hitvecFirst;
        TrackerHitVec hitvecSecond;
        GroupTracks * groupFirst = first->getGroupTracks();
        GroupTracks * groupSecond = second->getGroupTracks();
        
        if(groupFirst!=NULL){
          
          TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
          nTrkGrpFirst = int(tracksInGroupFirst.size());
          
          for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
            
            TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
            TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
            
            for(unsigned int i =0; i<hitVec.size(); ++i){
              hitvecFirst.push_back(hitVec[i]->getTrackerHit());          
            }
          }
        }
        
        if(groupSecond!=NULL){
          
          TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
          nTrkGrpSecond = int(tracksInGroupSecond.size());
          
          for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
            TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
            TrackerHitExtendedVec hitVec = 
            trkGrp->getTrackerHitExtendedVec();
            
            for(unsigned int i=0;i<hitVec.size();++i){
              hitvecSecond.push_back(hitVec[i]->getTrackerHit());
            }
          }
        }
        
        
        
        // for non-looping tracks 
        int nhitsFirst  = (int)hitvecFirst.size();
        int nhitsSecond = (int)hitvecSecond.size();
        int ntpcFirst   = 0;
        int ntpcSecond  = 0;
        float hitxyz[3];
        float dist[3];
        float maxdistFirst=0;
        float maxdistSecond=0;
        int ncloseFirst = 0;
        int ncloseSecond = 0;
        float zminFirst = 99999;
        float zminSecond = 99999;
        float zmaxFirst = -99999;
        float zmaxSecond = -99999;
        
        for(int ih =0;ih<nhitsFirst;++ih){
          
          float x = (float) hitvecFirst[ih]->getPosition()[0];
          float y = (float) hitvecFirst[ih]->getPosition()[1];
          float z = (float) hitvecFirst[ih]->getPosition()[2];
          
          if(fabs(z)<zminFirst) zminFirst=fabs(z);
          if(fabs(z)>zmaxFirst) zmaxFirst=fabs(z);
          
          float r = sqrt(x*x+y*y);
          
          if(r>_tpc_inner_r) ntpcFirst++;
          
          hitxyz[0]=x;
          hitxyz[1]=y;
          hitxyz[2]=z;
          helixSecond.getDistanceToPoint(hitxyz, dist);
          
          // compare 3D distance between hit and extrapolation
          if(dist[2]>maxdistFirst) maxdistFirst=dist[2];
          if(dist[2]<_hitDistanceCutHighPtMerge) ncloseFirst++;
        }
        
        for(int ih =0;ih<nhitsSecond;++ih){
          
          float x = (float)hitvecSecond[ih]->getPosition()[0];
          float y = (float)hitvecSecond[ih]->getPosition()[1];
          float z = (float)hitvecSecond[ih]->getPosition()[2];
          
          if(fabs(z)<zminSecond) zminSecond=fabs(z);
          if(fabs(z)>zmaxSecond) zmaxSecond=fabs(z);
          
          float r = sqrt(x*x+y*y);
          
          if(r>_tpc_inner_r) ntpcSecond++;
          
          hitxyz[0]=x;
          hitxyz[1]=y;
          hitxyz[2]=z;
          helixFirst.getDistanceToPoint(hitxyz, dist);
          
          // compare 3D distance between hit and extrapolation
          if(dist[2]>maxdistSecond) maxdistSecond=dist[2];
          if(dist[2]<_hitDistanceCutHighPtMerge) ncloseSecond++;
        }
        
        float fcloseFirst  = (float)ncloseFirst/(float)nhitsFirst;
        float fcloseSecond = (float)ncloseSecond/(float)nhitsSecond;
        
        
        
        bool split = false;
        //std::cout << "Momenta = " << momFirst << " " << momSecond << std::endl;
        //std::cout << "MaxDist = " << maxdistSecond << " " << maxdistFirst << " " << _maxHitDistanceCutHighPtMerge << std::endl;
        //std::cout << "close   = " << fcloseSecond << " " << fcloseFirst << " " << _maxFractionOfOutliersCutHighPtMerge << std::endl;
        //std::cout << "ntpc    = " << ntpcFirst << " " << ntpcSecond << " " << _tpc_pad_height+10 << std::endl;
        //std::cout << "pdot    = " << pdot << " significance " << significance << std::endl;
        
        TrackExtended * combinedTrack = TrialCombineTracks(first,second);
        
        if(combinedTrack != NULL){
          //std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << first->getNDF()+second->getNDF()+5 << std::endl;
          if(combinedTrack->getNDF()+10>first->getNDF()+second->getNDF()+5){
            split = true;
            dpOverP = 0;
            //std::cout << " Forcing MERGE " << std::endl;
          }
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
        } 
        else {
          //std::cout << "Could not combine track " << std::endl;
          if(significance<5 && fcloseFirst>_maxFractionOfOutliersCutHighPtMerge){
            split = true;
            dpOverP = 0;
            //      int overlap = SegmentRadialOverlap(first,second);
            //std::cout << " Forcing MERGE " << overlap << std::endl;
          }
        }
        
        // criteria for split track
        // old criterion
        if( maxdistSecond < _maxHitDistanceCutHighPtMerge && maxdistFirst < _maxHitDistanceCutHighPtMerge 
           && 
           (fcloseSecond > _maxFractionOfOutliersCutHighPtMerge || fcloseFirst > _maxFractionOfOutliersCutHighPtMerge) 
           && 
           ntpcFirst+ntpcSecond < _tpc_pad_height+10.) {
          
          split = true;        
          
        }
        
        if(split){
          result = dpOverP;
        }
        
      }
    }
  }
  
  return result;
  
}

void FullLDCTracking_MarlinTrk::AddNotAssignedHits() {
  
  
  // Creating helix extrapolations
  //  if (_assignSETHits>0||_assignETDHits>0)
  //    CreateExtrapolations();
  
  //  if (_assignSETHits>0) { // Assignment of SET Hits
  //    
  //    const gear::GearParameters& pSETDet = Global::GEAR->getGearParameters("SET");  
  //    int nLayersSET = int(pSETDet.getDoubleVals("SETLayerRadius").size());
  //    
  //    int nSETHits = _allSETHits.size();
  //    std::vector<TrackerHitExtendedVec> SETHits;
  //    SETHits.resize(nLayersSET);
  //    
  //    for (int iSET=0;iSET<nSETHits;++iSET) {
  //      TrackerHitExtended * trkHit = _allSETHits[iSET];
  //      TrackerHit * hit = trkHit->getTrackerHit();
  //      int layer = getLayerID(trkHit);
  //      if (layer>=0&&layer<nLayersSET) 
  //        SETHits[layer].push_back(trkHit);   
  //    }
  //    for (int iL=0; iL<nLayersSET; ++iL) { // loop over SET layers
  //      TrackerHitExtendedVec hitVec = SETHits[iL];
  //      int refit = 1;
  //      AssignOuterHitsToTracks(hitVec,_distCutForSETHits,refit);
  //    }
  //  }
  
  //  if (_assignETDHits>0) { // Assignment of ETD Hits
  //    
  //    const gear::GearParameters& pETDDet = Global::GEAR->getGearParameters("ETD");  
  //    int nLayersETD = int(pETDDet.getDoubleVals("ETDLayerZ").size());
  //    
  //    int nETDHits = _allETDHits.size();
  //    std::vector<TrackerHitExtendedVec> ETDHits;
  //    ETDHits.resize(nLayersETD);
  //    
  //    for (int iETD=0;iETD<nETDHits;++iETD) {
  //      TrackerHitExtended * trkHit = _allETDHits[iETD];
  //      TrackerHit * hit = trkHit->getTrackerHit();
  //      int layer = getLayerID(trkHit);
  //      if (layer>=0 && layer < nLayersETD) 
  //        ETDHits[layer].push_back(trkHit);
  //    }
  //    for (int iL=0; iL<nLayersETD; ++iL) { // loop over ETD layers
  //      TrackerHitExtendedVec hitVec = ETDHits[iL];
  //      int refit = 0;
  //      AssignOuterHitsToTracks( hitVec, _distCutForETDHits, refit );
  //    }
  //    
  //  }
  
  //  // Cleaning up helix extrapolations
  //  if (_assignSETHits>0||_assignETDHits>0)
  //    CleanUpExtrapolations();
  
  if (_assignSITHits>0) { // Treatment of left-over SIT hits 
    
    std::vector<TrackerHitExtendedVec> nonAssignedSITHits;    
    nonAssignedSITHits.resize(_nLayersSIT);
    
    int nSITHits = int(_allSITHits.size());
    
    // loop over all SIT hits ...
    for (int iH=0;iH<nSITHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allSITHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit * trkHit = trkHitExt->getTrackerHit();
      
        int layer = getLayerID(trkHit);

        if (  _reading_loi_data == false ) { // divide by two as we are treating the SIT as TWO stereo layers 
          layer = layer / 2 ;
        }
        
        if (layer >=0 && layer < _nLayersSIT) {
          nonAssignedSITHits[layer].push_back(trkHitExt);
        }
      }
    }       
    
    for (int iL=_nLayersSIT-1;iL>=0;--iL) { // reverse loop over layers in Si
      TrackerHitExtendedVec hitVec = nonAssignedSITHits[iL];
      AssignSiHitsToTracks(hitVec,
                           _distCutForSITHits);
    }
  }
  
  if (_assignFTDHits>0) { // Treatment of left-over FTD hits
    
    std::vector<TrackerHitExtendedVec> nonAssignedFTDHits;
    nonAssignedFTDHits.resize(_nLayersFTD);
    
    int nFTDHits = int(_allFTDHits.size());
    
    // loop over all FTD hits ...
    for (int iH=0;iH<nFTDHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allFTDHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit * trkHit = trkHitExt->getTrackerHit();
        
        
        // get the layer number
        int layer = getLayerID(trkHit);
        int petalIndex = getModuleID(trkHit);
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals 
        if ( petalIndex % 2 == 0 ) {
          layer = 2*layer;
        }
        else {
          layer = 2*layer + 1;
        }
        
        
        if (layer >=0 && layer < _nLayersFTD)
          nonAssignedFTDHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL=_nLayersFTD-1;iL>=0;--iL) {
      if ( nonAssignedFTDHits[iL].size()!=0 ) {
        
        TrackerHitExtendedVec hitVec = nonAssignedFTDHits[iL];
        AssignSiHitsToTracks(hitVec,
                             _distCutForFTDHits);     
        
      }
    }
  }
  
  
  
  if (_assignVTXHits>0) { // Treatment of left-over VTX hits
    
     streamlog_out(DEBUG4) << "Assign VXD hits *********************************" << std::endl;
    
    std::vector<TrackerHitExtendedVec> nonAssignedVTXHits;
    nonAssignedVTXHits.resize(_nLayersVTX);
    
    int nVTXHits = int(_allVTXHits.size());
    
    // loop over all VXD hits ...
    for (int iH=0;iH<nVTXHits;++iH) {
      
      TrackerHitExtended * trkHitExt = _allVTXHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      
      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit * trkHit = trkHitExt->getTrackerHit();
        
        int layer = getLayerID(trkHit);
        
        if (layer >=0 && layer < _nLayersVTX)
          nonAssignedVTXHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL=_nLayersVTX-1;iL>=0;--iL) {
      TrackerHitExtendedVec hitVec = nonAssignedVTXHits[iL];
      AssignSiHitsToTracks(hitVec,
                           _distCutForVTXHits);     
    }
  }
  
  streamlog_out(DEBUG4) << "Assign TPC hits *********************************" << std::endl;
  
  if (_assignTPCHits) {// Treatment of left-over TPC hits
    TrackerHitExtendedVec nonAssignedTPCHits;
    int nTPCHits = int(_allTPCHits.size());
    for (int iH=0;iH<nTPCHits;++iH) {
      TrackerHitExtended * trkHitExt = _allTPCHits[iH];
      TrackExtended * trkExt = trkHitExt->getTrackExtended();
      if (trkExt == NULL) {
        nonAssignedTPCHits.push_back(trkHitExt);
      }
    }
    AssignTPCHitsToTracks(nonAssignedTPCHits,
                          _distCutForTPCHits);
  }
  
  
}


void FullLDCTracking_MarlinTrk::CreateExtrapolations() {
  
  _trackExtrapolatedHelix.clear();
  
  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix = GetExtrapolationHelix( trk );
    _trackExtrapolatedHelix[trk] = helix;
  }
  
}

void FullLDCTracking_MarlinTrk::CleanUpExtrapolations() {
  
  int nTrk = int(_trkImplVec.size());
  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    TrackExtended * trk = _trkImplVec[iTrk];
    HelixClass * helix =  _trackExtrapolatedHelix[trk];
    delete helix;
  }  
  _trackExtrapolatedHelix.clear();
  
}


//void FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit) {
//  
//  int nHits = int(hitVec.size());
//  int nTrk = int(_trkImplVec.size());
//  
//  std::map <TrackExtended*,bool> flagTrack;
//  std::map <TrackerHitExtended*,bool> flagHit;
//  TrackHitPairVec pairs;
//  flagTrack.clear();
//  flagHit.clear();
//  pairs.clear();
//  
//  for (int iH=0;iH<nHits;++iH) {
//    float pos[3];
//    TrackerHitExtended * trkHitExt = hitVec[iH];
//    TrackerHit * hit = trkHitExt->getTrackerHit();
//    for (int ip=0;ip<3;++ip)
//      pos[ip] = float(hit->getPosition()[ip]);
//    //      float r_hit = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
//    //      float phi_hit = atan2(pos[1],pos[0]);
//    //      float z_hit = pos[2];
//    //      std::cout << r_hit << " " << phi_hit << " "
//    //                << z_hit << " " 
//    //                << trkHitExt->getResolutionRPhi() << " "
//    //                << trkHitExt->getResolutionZ() << std::endl;
//    
//    for (int iT=0;iT<nTrk;++iT) {
//      TrackExtended * trkExt = _trkImplVec[iT];
//      float tanLambda = trkExt->getTanLambda();           
//      float product = pos[2]*tanLambda;
//      if (product>0) {
//        HelixClass * helix = _trackExtrapolatedHelix[trkExt];
//        float distance = helix->getDistanceToPoint(pos,dcut);
//        //      std::cout << "Dist = " << dist[2] << std::endl;
//        if (distance<dcut) {
//          TrackHitPair * trkHitPair = 
//          new TrackHitPair(trkExt,trkHitExt,distance);
//          pairs.push_back(trkHitPair);
//          flagTrack[trkExt] = true;
//          flagHit[trkHitExt] = true;
//        }
//      }
//    }
//  }
//  
//  int nPairs = int(pairs.size());
//  if (nPairs>0) {
//    SortingTrackHitPairs(pairs);
//    for (int iP=0;iP<nPairs;++iP) {
//      TrackHitPair * trkHitPair = pairs[iP];
//      TrackExtended * trkExt = trkHitPair->getTrackExtended();
//      TrackerHitExtended * trkHitExt = 
//      trkHitPair->getTrackerHitExtended();       
//      if (flagTrack[trkExt] && flagHit[trkHitExt]) {
//        if (refit==0) {
//          trkExt->addTrackerHitExtended( trkHitExt );
//          trkHitExt->setUsedInFit( false );
//          trkHitExt->setTrackExtended( trkExt );
//        }
//        else {
//          TrackerHitExtendedVec hitsInTrack = 
//          trkExt->getTrackerHitExtendedVec();
//          int nTotH = int(hitsInTrack.size());
//          int nHitsInFit = 0;
//          for (int iTH=0;iTH<nTotH;++iTH) {
//            TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
//            if (hitInTrack->getUsedInFit())
//              nHitsInFit++;
//          }
//          float * x_h = new float[nHitsInFit+1];
//          float * y_h = new float[nHitsInFit+1];
//          float * z_h = new float[nHitsInFit+1];
//          int * idet_h = new int[nHitsInFit+1];
//          int * ityp_h = new int[nHitsInFit+1];
//          int * lhits = new int[nHitsInFit+1];
//          float * rR_h = new float[nHitsInFit+1];
//          float * rZ_h = new float[nHitsInFit+1];
//          int iHitInFit = 0;
//          for (int iHit=0;iHit<nTotH;++iHit) {
//            TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
//            if (hitInTrack->getUsedInFit()) {
//              TrackerHit * hit = hitInTrack->getTrackerHit();
//              x_h[iHitInFit] = float(hit->getPosition()[0]);
//              y_h[iHitInFit] = float(hit->getPosition()[1]);
//              z_h[iHitInFit] = float(hit->getPosition()[2]);
//              idet_h[iHitInFit] = hitInTrack->getDet();
//              ityp_h[iHitInFit] = hitInTrack->getType();
//              rR_h[iHitInFit] = hitInTrack->getResolutionRPhi();
//              rZ_h[iHitInFit] = hitInTrack->getResolutionZ();
//              iHitInFit++;
//            }
//          }
//          TrackerHit * remainHit = trkHitExt->getTrackerHit();
//          x_h[iHitInFit] = float(remainHit->getPosition()[0]);
//          y_h[iHitInFit] = float(remainHit->getPosition()[1]);
//          z_h[iHitInFit] = float(remainHit->getPosition()[2]);
//          idet_h[iHitInFit] = trkHitExt->getDet();
//          ityp_h[iHitInFit] = trkHitExt->getType();
//          rR_h[iHitInFit] = trkHitExt->getResolutionRPhi();
//          rZ_h[iHitInFit] = trkHitExt->getResolutionZ();          
//          iHitInFit++;
//          
//          int NPT = iHitInFit;
//          float chi2_D;
//          int ndf_D;
//          float chi2rphi,chi2z;
//          float par[5];
//          float epar[15];
//          float refPoint[3];
//          
//          int ierr = _trackFit.DoFitting(_useExtraPoint,_optFit,NPT,
//                                         _bField,idet_h,ityp_h,
//                                         _chi2PrefitCut,
//                                         x_h,y_h,z_h,rR_h,rZ_h,
//                                         par,epar,refPoint,chi2_D,ndf_D,
//                                         chi2rphi,chi2z,lhits);              
//          
//          float chi2ndf = chi2_D/float(ndf_D);
//          if (ierr==0 && chi2ndf<_chi2FitCut) {
//            trkExt->addTrackerHitExtended( trkHitExt );
//            trkExt->setCovMatrix(epar);
//            trkExt->setOmega(par[0]);   
//            trkExt->setTanLambda(par[1]);
//            trkExt->setPhi(par[2]);
//            trkExt->setD0(par[3]);
//            trkExt->setZ0(par[4]);
//            trkExt->setChi2(chi2_D);
//            trkExt->setNDF(ndf_D);
//            trkHitExt->setTrackExtended( trkExt );
//            trkHitExt->setUsedInFit( true );
//            flagTrack[trkExt] = false;
//            flagHit[trkHitExt] = false;
//            if (_optFit == 4) {
//              GroupTracks * groupTracks = trkExt->getGroupTracks();
//              if (groupTracks!=NULL) {
//                TrackExtendedVec group = groupTracks->getTrackExtendedVec();
//                int nTracks = int(group.size());
//                if (nTracks>0) {
//                  TrackExtended * siTrack = group[0];
//                  // preserving positiveness of matrix
//                  float eD0 = siTrack->getCovMatrix()[0];
//                  float eZ0 = siTrack->getCovMatrix()[9];
//                  float sD0 = sqrt(eD0/epar[0]);
//                  float sZ0 = sqrt(eZ0/epar[9]);
//                  epar[0]  = eD0;
//                  epar[9]  = eZ0;
//                  epar[1]  = sD0*epar[1];
//                  epar[3]  = sD0*epar[3];
//                  epar[10] = sD0*epar[10];
//                  epar[7]  = sZ0*epar[7];
//                  epar[8]  = sZ0*epar[8];
//                  epar[13] = sZ0*epar[13];
//                  epar[6]  = sD0*sZ0*epar[6];
//                  trkExt->setZ0(siTrack->getZ0());
//                  trkExt->setD0(siTrack->getD0());
//                  trkExt->setCovMatrix(epar);
//                }
//              }
//            }
//          }
//          delete[] x_h;
//          delete[] y_h;
//          delete[] z_h;
//          delete[] rR_h;
//          delete[] rZ_h;
//          delete[] idet_h;
//          delete[] ityp_h;
//          delete[] lhits;
//        }
//      }
//    }
//    
//    for (int iP=0;iP<nPairs;++iP) {
//      TrackHitPair * trkHitPair = pairs[iP];
//      delete trkHitPair;
//    }
//    
//    pairs.clear();
//    
//  }
//  
//  
//}

HelixClass * FullLDCTracking_MarlinTrk::GetExtrapolationHelix( TrackExtended * track) {
  
  TrackerHitExtendedVec hitVec = track->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());
  
  HelixClass * helix = new HelixClass();
  
  int nHitsFit = nHits;
  if (nHits>_nHitsExtrapolation)
    nHitsFit = _nHitsExtrapolation;
  
  float * ampl = new float[nHitsFit];
  float * xhit = new float[nHitsFit];
  float * yhit = new float[nHitsFit];
  float * zhit = new float[nHitsFit];
  
  
  if (nHits<=_nHitsExtrapolation) {
    for (int iH=0;iH<nHitsFit;++iH) {
      TrackerHit * trkHit = hitVec[iH]->getTrackerHit();
      ampl[iH] = 1.0;
      xhit[iH] = (trkHit->getPosition()[0]);
      yhit[iH] = (trkHit->getPosition()[1]);
      zhit[iH] = (trkHit->getPosition()[2]);
    }    
    
  }
  else {
    int * index = new int[nHits];
    float * zcoor = new float[nHits];        
    for (int iH=0;iH<nHits;++iH) {
      TrackerHit * trkHit = hitVec[iH]->getTrackerHit();
      zcoor[iH] = fabs(float(trkHit->getPosition()[2]));
    }
    int order = 1;
    GeneralSorting(index,zcoor,order,nHits);
    for (int iH=0;iH<nHitsFit;++iH) {
      int idx = index[iH];
      TrackerHit * trkHit = hitVec[idx]->getTrackerHit();
      ampl[iH] = 1.0;
      xhit[iH] = (trkHit->getPosition()[0]);
      yhit[iH] = (trkHit->getPosition()[1]);
      zhit[iH] = (trkHit->getPosition()[2]);
    }
    delete[] index;
    delete[] zcoor;
  }
  
  ClusterShapes * shapes = new ClusterShapes(nHitsFit,ampl,xhit,yhit,zhit);
  
  float parSh[5];
  float dparSh[5];
  float distmax = 0;
  float chi2Sh = 0;
  // do fitting
  shapes->FitHelix(500, 0, 1, parSh, dparSh, chi2Sh, distmax);  
  
  float x0Sh = parSh[0];
  float y0Sh = parSh[1];
  float r0Sh = parSh[2];
  float bzSh = parSh[3];
  float phi0Sh = parSh[4];
  float signPz = 1;
  float zBegin = zhit[0];
  float zEnd   = zhit[0];
  float zMax   = fabs(zBegin);
  float zMin   = fabs(zEnd);
  for (int iH=1;iH<nHitsFit;++iH) {
    float zCurrent = fabs(zhit[iH]);
    if (zCurrent>zMax) {
      zMax = zCurrent;
      zEnd = zhit[iH];
    }
    if (zCurrent<zMin) {
      zMin = zCurrent;
      zBegin = zhit[iH];
    }      
  }
  
  delete shapes;
  
  if (zEnd<zBegin)
    signPz = -1;
  helix->Initialize_BZ(x0Sh, y0Sh, r0Sh, 
                       bzSh, phi0Sh, _bField,signPz,
                       zBegin);
  
  delete[] ampl;
  delete[] xhit;
  delete[] yhit;
  delete[] zhit;
  
  return helix;
  
}


void FullLDCTracking_MarlinTrk::AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec,
                                                      float dcut) {
  
  int nHits = int(hitVec.size());
  int nTrk = int(_trkImplVec.size());
  
  for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
    TrackExtended * foundTrack = _trkImplVec[iT];
    GroupTracks * group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
      TrackerHitExtendedVec hitVecGrp = trkGrp->getTrackerHitExtendedVec();
      int nHits = int(hitVecGrp.size());
      float zMin = 1.0e+20;
      float zMax = -1.0e+20;
      float startPoint[3] = {0.,0.,0.};
      float endPoint[3]   = {0.,0.,0.};
      for (int iH=0;iH<nHits;++iH) {
        TrackerHitExtended * trkHitExt = hitVecGrp[iH];
        float pos[3] = {0.,0.,0.};
        for (int iC=0;iC<3;++iC) 
          pos[iC] = float(trkHitExt->getTrackerHit()->getPosition()[iC]);         
        if (pos[2]>zMax) {
          zMax = pos[2];
          for (int iC=0;iC<3;++iC)
            endPoint[iC] = pos[iC];           
        }
        if (pos[2]<zMin) {
          zMin = pos[2];
          for (int iC=0;iC<3;++iC)
            startPoint[iC] = pos[iC];
        }
      }
      trkGrp->setStart(startPoint);
      trkGrp->setEnd(endPoint);
    }
  }
  
  
  // replace previous version with faster loop ordering
  
  std::vector<float> minDistances(nHits, dcut);
  std::vector<TrackExtended*> tracksToAttach(nHits);
  std::vector< std::vector<float> > HitPositions(nHits);
  std::vector<int> HitSign(nHits);//Positive or Negative side
  for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
    tracksToAttach[iH]=NULL;
    //Get all TrackerHit positions, so we only have to get them once
    TrackerHit* temphit = hitVec[iH]->getTrackerHit();
    const double *temppos = temphit->getPosition();
    HitPositions[iH].push_back(float(temppos[0]));
    HitPositions[iH].push_back(float(temppos[1]));
    HitPositions[iH].push_back(float(temppos[2]));
    HitSign[iH]=std::signbit(temppos[2]);
  }    
  
  streamlog_out(DEBUG4) << " Starting loop " << nTrk << " tracks   and  " << nHits << " hits" << std::endl;
  
  for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
    TrackExtended * foundTrack = _trkImplVec[iT];
    int tanlambdaSign = std::signbit(foundTrack->getTanLambda());//we only care about positive or negative
    GroupTracks * group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
      float tanLambda = trkGrp->getTanLambda();
      float omega = trkGrp->getOmega();
      float d0 = trkGrp->getD0();
      float z0 = trkGrp->getZ0();
      float phi0 = trkGrp->getPhi();
      float startPointZ = trkGrp->getStart()[2];
      float endPointZ   = trkGrp->getEnd()[2];
      HelixClass helix;
      helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
      float OnePFivehalfPeriodZ = 1.5*fabs(acos(-1.)*tanLambda/omega);
      
      for (int iH=0;iH<nHits;++iH) { // loop over leftover TPC hits
        
        //check if the hit and the track or on the same side
        //xor return 1, if hits are different
        if ( tanlambdaSign^HitSign[iH] ) continue;
        
        float DeltaStart = fabs(HitPositions[iH][2]-startPointZ);
        float DeltaEnd = fabs(HitPositions[iH][2]-endPointZ);
        bool consider = DeltaStart <= OnePFivehalfPeriodZ;
        consider = consider || (DeltaEnd <= OnePFivehalfPeriodZ);
        consider = consider || ( (HitPositions[iH][2]>=startPointZ) && (HitPositions[iH][2]<=endPointZ) );
        
        if(consider){
          float distance = helix.getDistanceToPoint(HitPositions[iH], minDistances[iH]);
          if (distance < minDistances[iH]) {
            minDistances[iH] = distance;
            tracksToAttach[iH] = foundTrack;
          }
        }
      } // loop over leftover TPC hits
    } //groups in tracks
  } // loop over all tracks
  
  for (int iH=0;iH<nHits;++iH) {
    TrackerHitExtended * trkHitExt = hitVec[iH];
    if (tracksToAttach[iH]!=NULL) {
      tracksToAttach[iH]->addTrackerHitExtended(trkHitExt);
      trkHitExt->setTrackExtended( tracksToAttach[iH] );
      trkHitExt->setUsedInFit( false );
    }
  }
  
  streamlog_out(DEBUG4) << " Fast loop done " << std::endl;
  
  
  //     for (int iH=0;iH<nHits;iH++) { // loop over leftover TPC hits
  //    TrackerHitExtended * hitExt = hitVec[iH];
  //    float pos[3];
  //    for (int ip=0;ip<3;++ip) 
  //        pos[ip] = float(hitExt->getTrackerHit()->getPosition()[ip]);
  //    float minDist = 1.0e+20;
  //    TrackExtended * trackToAttach = NULL;
  //    for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
  //        TrackExtended * foundTrack = _trkImplVec[iT];
  //        float tanLambdaFound = foundTrack->getTanLambda();
  //        float product = tanLambdaFound*pos[2];
  //        if (product>0) {
  //          GroupTracks * group = foundTrack->getGroupTracks();
  //          TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
  //          int nTrkGrp = int(tracksInGroup.size());
  //          for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
  //            TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
  //            float tanLambda = trkGrp->getTanLambda();
  //            float omega = trkGrp->getOmega();
  //            float d0 = trkGrp->getD0();
  //            float z0 = trkGrp->getZ0();
  //            float phi0 = trkGrp->getPhi();
  //            float dist[3];
  //            float startPoint[3];
  //            float endPoint[3];
  //            for (int iC=0;iC<3;++iC) {
  //              startPoint[iC] = trkGrp->getStart()[iC];
  //              endPoint[iC] = trkGrp->getEnd()[iC];
  //            }
  //            HelixClass helix;
  //            helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
  //            float halfPeriodZ = fabs(acos(-1.)*tanLambda/omega);
  //            helix.getDistanceToPoint(pos,dist);
  //            float DeltaStart = fabs(pos[2]-startPoint[2]);
  //            float DeltaEnd = fabs(pos[2]-endPoint[2]);
  //            bool consider = DeltaStart <= 1.5*halfPeriodZ;
  //            consider = consider || (DeltaEnd <= 1.5*halfPeriodZ);
  //            consider = consider || ( (pos[2]>=startPoint[2]) && (pos[2]<=endPoint[2]) );
  // //                 float ZMin = DeltaStart;
  // //                 if (DeltaEnd<ZMin)
  // //                   ZMin = DeltaEnd;
  //            if (dist[2]<dcut && consider && dist[2]<minDist) {
  //              minDist = dist[2];
  //              trackToAttach = foundTrack;
  //            }
  //          }
  //        }
  //    }
  
  //    if (trackToAttach!=NULL) {
  //   trackToAttach->addTrackerHitExtended(hitExt);
  //  hitExt->setTrackExtended( trackToAttach );
  //  hitExt->setUsedInFit( false );
  //  if(trackToAttach!=tracksToAttach[iH])std::cout << " Check Failed" << trackToAttach << "  " << tracksToAttach[iH] << std::endl;
  //
  //}
  //else {
  ///     std::cout << iH << " hit is not assigned : distance to closest track = " << minDist << std::endl;
  ///}
  //}
  //std::cout << " Slow loop done " << std::endl;
  
}

void FullLDCTracking_MarlinTrk::AssignSiHitsToTracks(TrackerHitExtendedVec hitVec,
                                                     float dcut) {
  
  int nHits = int(hitVec.size());
  int nTrk = int(_allNonCombinedTPCTracks.size());
  
  std::map <TrackExtended*,bool> flagTrack;
  std::map <TrackerHitExtended*,bool> flagHit;
  TrackHitPairVec pairs;
  flagTrack.clear();
  flagHit.clear();
  pairs.clear();
  
  for (int iH=0;iH<nHits;++iH) {
    
    float pos[3];
    TrackerHitExtended * trkHitExt = hitVec[iH];
    TrackerHit * hit = trkHitExt->getTrackerHit();
    
    for (int ip=0;ip<3;++ip) {
      pos[ip] = float(hit->getPosition()[ip]);
    }
    
    for (int iT=0;iT<nTrk;++iT) {
      
      TrackExtended * trkExt = _allNonCombinedTPCTracks[iT];
      
      float tanLambda = trkExt->getTanLambda();       
      float product = pos[2]*tanLambda;
      
      if (product>0) {
        
        float d0 = trkExt->getD0();
        float z0 = trkExt->getZ0();
        float tanLambda = trkExt->getTanLambda();
        float phi0 = trkExt->getPhi();
        float omega = trkExt->getOmega();
        
        HelixClass helix;
        helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
        float distance = helix.getDistanceToPoint(pos,dcut);
        
        if (distance<dcut) {
          TrackHitPair * trkHitPair = 
          new TrackHitPair(trkExt,trkHitExt,distance);
          pairs.push_back(trkHitPair);
          flagTrack[trkExt] = true;
          flagHit[trkHitExt] = true;
        }
      }
    }
  }
  
  int nPairs = int(pairs.size());
  
  if (nPairs>0) {
    
    SortingTrackHitPairs(pairs);
    
    for (int iP=0;iP<nPairs;++iP) {
      
      TrackHitPair * trkHitPair = pairs[iP];
      TrackExtended * trkExt = trkHitPair->getTrackExtended();
      TrackerHitExtended * trkHitExt = 
      
      trkHitPair->getTrackerHitExtended();
      
      if (flagTrack[trkExt] && flagHit[trkHitExt]) {              
        
        // get the hits already assigned to the track
        TrackerHitExtendedVec hitsInTrack = trkExt->getTrackerHitExtendedVec();
        
        int nTotH = int(hitsInTrack.size());
        int nHitsInFit = 0;
        
        for (int iTH=0;iTH<nTotH;++iTH) {
          
          TrackerHitExtended * hitInTrack = hitsInTrack[iTH];
          
          // count the number of hits used in the fit
          if (hitInTrack->getUsedInFit()) {
            nHitsInFit++; 
          }
        }
        
        int iHitInFit = 0;
        
        // add the previously used hits from the track to the vectors 
        
        EVENT::TrackerHitVec trkHits;
        
        for (int iHit=0;iHit<nTotH;++iHit) {
          
          TrackerHitExtended * hitInTrack = hitsInTrack[iHit];
          if (hitInTrack->getUsedInFit()) {
            TrackerHit * hit = hitInTrack->getTrackerHit();
            iHitInFit++;
            if(hit) { 
              trkHits.push_back(hit);   
            }
            else{
              throw EVENT::Exception( std::string("FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: TrackerHit pointer == NULL ")  ) ;
            }
          }
        }
        
        // add the hit to be attached to the vectors 
        TrackerHit * remainHit = trkHitExt->getTrackerHit();
        iHitInFit++;
        trkHits.push_back(remainHit);
        
        
        double chi2;
        int ndf;
        
        if( trkHits.size() < 3 ) return ;
        
        MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
        sort(trkHits.begin(), trkHits.end(), FullLDCTracking_MarlinTrk::compare_r() );
        
        EVENT::TrackerHitVec::iterator it = trkHits.begin();
        
        streamlog_out(DEBUG2) << "Start Fitting: AddHits: number of hits to fit " << trkHits.size() << std::endl;
      
//        
//        int number_of_added_hits = 0;
//        for( it = trkHits.begin() ; it != trkHits.end() ; ++it ) {
//          
//          if (marlin_trk->addHit(*it) == 0){
//            ++number_of_added_hits;
//          }
//          else{
//            streamlog_out(DEBUG4) << "Hit " << it - trkHits.begin() << " Dropped " << std::endl;
//          }
//          
//        }
//        
//        if( number_of_added_hits < 3 ) {
//          delete marlin_trk ;
//          return ;
//        }
                            
        int number_of_added_hits = 0;
        int ndof_added = 0;
        TrackerHitVec added_hits;
        
        for( it = trkHits.begin() ; it != trkHits.end() ; ++it ) {
          
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
          else { // normal non composite hit
            
            if (marlin_trk->addHit( trkHit ) == 0) {
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
          streamlog_out(DEBUG3) << "SiliconTracking_MarlinTrk::FinalRefit: Cannot fit less with less than 8 degrees of freedom. Number of hits =  " << number_of_added_hits << " ndof = " << ndof_added << std::endl;
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
        
        for (int icov = 0; icov<covMatrix.size(); ++icov) {
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
        
        
        int fit_status = marlin_trk->fit() ; 
        
        if( fit_status != 0 ){ 
          delete marlin_trk ;
          return;
        }
        
        const gear::Vector3D point(0.,0.,0.); // nominal IP
        int return_code = 0;
        
        TrackStateImpl trkState ;
        return_code = marlin_trk->propagate(point, trkState, chi2, ndf ) ;
        
        delete marlin_trk ;
        
        if (return_code !=0 || (chi2/(double)ndf) > _chi2FitCut) {
          return;
        }
        
        
        // note trackAR which is of type TrackExtended, only takes fits set for ref point = 0,0,0 
        trkExt->setOmega(trkState.getOmega());
        trkExt->setTanLambda(trkState.getTanLambda());
        trkExt->setPhi(trkState.getPhi());
        trkExt->setD0(trkState.getD0());
        trkExt->setZ0(trkState.getZ0());
        
        
        float cov[15];
        
        for (int i = 0 ; i<15 ; ++i) {
          cov[i] = trkState.getCovMatrix().operator[](i);
        }
        
        trkExt->setCovMatrix(cov);
        trkExt->setChi2(chi2);
        trkExt->setNDF(ndf);
        
        trkExt->addTrackerHitExtended( trkHitExt );
        trkHitExt->setTrackExtended( trkExt );
        trkHitExt->setUsedInFit( true );
        flagTrack[trkExt] = false;
        flagHit[trkHitExt] = false;
        
        
      }
    }
    
    for (int iP=0;iP<nPairs;++iP) {
      TrackHitPair * trkHitPair = pairs[iP];
      delete trkHitPair;
    }
    
    pairs.clear();
    
  }
}

void FullLDCTracking_MarlinTrk::PrintOutMerging(TrackExtended * firstTrackExt, TrackExtended * secondTrackExt, int iopt) {
  // iopt = 1 false Si and TPC merging
  // iopt = 2 false Si and TPC forced merging
  // iopt = 3 false TPC segments merging
  // iopt = 4 false Comb Si and TPC merging
  // iopt = 5 false Comb TPC and TPC merging
  // iopt = 6 unmerged TPC and Si segments (soft merging)
  // iopt = 7 unmerged TPC and Si segments (forced merging)
  // iopt = 8 unmerged Comb and TPC
  // iopt = 9 unmerged TPC segments
  
  char strg[200];
  
  try {
    
    Track * firstTrack = firstTrackExt->getTrack();
    Track * secondTrack = secondTrackExt->getTrack();
    
    std::string firstColName = _TPCTrackMCPCollName;
    std::string secondColName = _TPCTrackMCPCollName;
    
    if (iopt==1) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==2) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==3) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==4) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==5) {
      secondColName = _TPCTrackMCPCollName;
    }
    else if (iopt==6) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==7) {
      secondColName = _SiTrackMCPCollName;
    }
    else if (iopt==8) {
      secondColName = _SiTrackMCPCollName;
    }    
    else {
      secondColName = _TPCTrackMCPCollName;
    }
    
    
    LCCollection * firstCol = _evt->getCollection(firstColName.c_str());
    LCCollection * secondCol = _evt->getCollection(secondColName.c_str());
    
    
    LCRelationNavigator firstNav(firstCol);
    LCRelationNavigator secondNav(secondCol);
    LCObjectVec firstVec = firstNav.getRelatedToObjects(firstTrack);
    FloatVec firstWeights = firstNav.getRelatedToWeights(firstTrack);
    LCObject * firstMCP = NULL;
    float firstWght = 0;
    int nObj = firstVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (firstWeights[iObj]>firstWght) {
        firstWght = firstWeights[iObj];
        firstMCP = firstVec[iObj];
      }
    }
    
    LCObjectVec secondVec = secondNav.getRelatedToObjects(secondTrack);
    FloatVec secondWeights = secondNav.getRelatedToWeights(secondTrack);
    LCObject * secondMCP = NULL;
    float secondWght = 0;
    nObj = secondVec.size();
    for (int iObj=0;iObj<nObj;++iObj) {
      if (secondWeights[iObj]>secondWght) {
        secondWght = secondWeights[iObj];
        secondMCP = secondVec[iObj];
      }
    }
    
    
    float d0First = firstTrackExt->getD0();
    float z0First = firstTrackExt->getZ0();
    float omegaFirst = firstTrackExt->getOmega();
    float tanLFirst = firstTrackExt->getTanLambda();
    float phi0First = firstTrackExt->getPhi();
    
    float d0Second = secondTrackExt->getD0();
    float z0Second = secondTrackExt->getZ0();
    float omegaSecond = secondTrackExt->getOmega();
    float tanLSecond = secondTrackExt->getTanLambda();
    float phi0Second = secondTrackExt->getPhi();            
    
    HelixClass firstHelix;
    firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
    float pxFirst = firstHelix.getMomentum()[0];
    float pyFirst = firstHelix.getMomentum()[1];
    float pzFirst = firstHelix.getMomentum()[2];            
    
    HelixClass secondHelix;
    secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
    float pxSecond = secondHelix.getMomentum()[0];
    float pySecond = secondHelix.getMomentum()[1];
    float pzSecond = secondHelix.getMomentum()[2];          
    
    float dPx = pxFirst - pxSecond;
    float dPy = pyFirst - pySecond;
    float dPz = pzFirst - pzSecond;
    
    float dPplus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    dPx = pxFirst + pxSecond;
    dPy = pyFirst + pySecond;
    dPz = pzFirst + pzSecond;
    
    float dPminus = sqrt(dPx*dPx+dPy*dPy+dPz*dPz);
    
    float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
    float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
    
    //SJA:FIXME Hardcoded cut here should be removed 
    if(pFirst < 1.5 && pSecond < 1.5 )return; 
    
    const float sigmaPOverPFirst  = sqrt(firstTrackExt->getCovMatrix()[5])/fabs(omegaFirst);
    const float sigmaPOverPSecond = sqrt(secondTrackExt->getCovMatrix()[5])/fabs(omegaSecond);
    //    const float deltaP = fabs(pFirst-pSecond);
    const float sigmaPFirst = pFirst*sigmaPOverPFirst;
    const float sigmaPSecond = pSecond*sigmaPOverPSecond;
    //const float sigmaDeltaP = sqrt(sigmaPFirst*sigmaPFirst+sigmaPSecond*sigmaPSecond);
    //    const float significance = deltaP/sigmaDeltaP;
    
    
    float den = pFirst;
    if (pSecond<pFirst)
      den = pSecond;
    
    dPplus  = dPplus/den;
    dPminus = dPminus/den; 
    
    if (firstMCP!=secondMCP && iopt < 6) {
      
      if (iopt==1) {
        streamlog_out(DEBUG4) << "Erroneous combining Si and TPC segments (iopt=1) --->" << std::endl;
      }
      else if (iopt==2) {
        streamlog_out(DEBUG4) << "Erroneous merging of Si and TPC segments (iopt=2) --->" << std::endl;
      }
      else if (iopt==3) {
        streamlog_out(DEBUG4) << "Erroneous merging of TPC segments (iopt=3) ---> " << std::endl; 
      }
      else if (iopt==4) {
        streamlog_out(DEBUG4) << "Erroneous merging of combSi segment with uncombTPC segment (iopt=4) ---> " << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Erroneous merging of combTPC segment with uncombTPC segment (iopt=5) --->" << std::endl;
      }
      
      
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      
      
      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      streamlog_out(DEBUG4) << std::endl;
      
    }
    else if (firstMCP==secondMCP && ( (iopt==8) || (iopt==9) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      if (iopt==8) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Comb segments (iopt=8) --->" << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Unmerged TPC segments (iopt=9) --->" << std::endl;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
      sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega
      << " AngleCut = " << deltaAngle
      << " dOmega = " << dOmega
      << " angle = " << angle << std::endl;
      
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      
      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      streamlog_out(DEBUG4) << std::endl;
      
    }    
    else if (firstMCP==secondMCP && ( (iopt == 6) || (iopt == 7) ) ) {
      
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      
      if (iopt ==6) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=6) --->" << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=7) --->" << std::endl;
        deltaOmega = _dOmegaForForcedMerging;
        deltaAngle = _angleForForcedMerging;
      }
      
      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      
      float dOmega = fabs((omegaFirst-omegaSecond)/omegaSecond);
      float angle = (cos(phi0First)*cos(phi0Second)+sin(phi0First)*sin(phi0Second))*
      sin(qFirst)*sin(qSecond)+cos(qFirst)*cos(qSecond);
      angle = acos(angle);
      
      
      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega
      << " AngleCut = " << deltaAngle
      << " dOmega = " << dOmega
      << " angle = " << angle << std::endl;
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;      
      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      streamlog_out(DEBUG4) << std::endl;      
      
    }else if (firstMCP==secondMCP && iopt < 6) {
      //      return;
      if (iopt==1) {
        streamlog_out(DEBUG4) << "Correct combining Si and TPC segments (iopt=1) --->" << std::endl;
      }
      else if (iopt==2) {
        streamlog_out(DEBUG4) << "Correct merging of Si and TPC segments (iopt=2) --->" << std::endl;
      }
      else if (iopt==3) {
        streamlog_out(DEBUG4) << "Correct merging of TPC segments (iopt=3) ---> " << std::endl; 
      }
      else if (iopt==4) {
        streamlog_out(DEBUG4) << "Correct merging of combSi segment with uncombTPC segment (iopt=4) ---> " << std::endl;
      }
      else {
        streamlog_out(DEBUG4) << "Correct merging of combTPC segment with uncombTPC segment (iopt=5) --->" << std::endl;
      }
      
      streamlog_out(DEBUG3) << "    p         error       D0      Z0     Px      Py      Pz      wieght" << std::endl;
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pFirst, sigmaPFirst, d0First,z0First,pxFirst,pyFirst,pzFirst);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",firstWght);
      streamlog_out(DEBUG4) << strg;
      
      sprintf(strg,"%7.2f +- %7.2f    %7.1f %7.1f %7.2f %7.2f %7.2f  ",
              pSecond, sigmaPSecond, d0Second,z0Second,pxSecond,pySecond,pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg,"  %5.3f\n",secondWght);
      streamlog_out(DEBUG4) << strg;
      
      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;
      TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
      
      if(combinedTrack != NULL){
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << firstTrackExt->getNDF()+secondTrackExt->getNDF()+5 << std::endl;
      }else{
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      delete combinedTrack->getGroupTracks();
      delete combinedTrack;
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt) << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      
      streamlog_out(DEBUG4) << std::endl;
      
    }
  }
  
  catch(DataNotAvailableException &e){};
  
}    

void FullLDCTracking_MarlinTrk::GeneralSorting(int * index, float * val, int direct, int nVal) {
  /**
   Sorting of index vector in ascending (0) /descending (!=0) order of val
   */
  
  float valOne, valTwo, valTemp;
  int   indTemp;
  for (int i=0; i<nVal; ++i) {
    index[i] = i;
  }
  
  for (int i = 0 ; i < nVal-1; i++) {
    for (int j = 0; j < nVal-i-1; j++) {      
      valOne = val[j];
      valTwo = val[j+1];
      bool order = valOne > valTwo;
      if (direct>0) 
        order = valOne <= valTwo;
      if( order )
          {
        valTemp = val[j];
        val[j] = val[j+1];
        val[j+1] = valTemp;
        indTemp = index[j];
        index[j] = index[j+1];
        index[j+1] = indTemp;
          }
    }  
  }
  
  
}

int FullLDCTracking_MarlinTrk::SegmentRadialOverlap(TrackExtended* first, TrackExtended* second){
  
  
  int nTrkGrpFirst = 0;
  int nTrkGrpSecond = 0;
  TrackerHitVec hitvecFirst;
  TrackerHitVec hitvecSecond;
  GroupTracks * groupFirst = first->getGroupTracks();
  GroupTracks * groupSecond = second->getGroupTracks();
  
  if(groupFirst!=NULL){
    
    TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
    nTrkGrpFirst = int(tracksInGroupFirst.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpFirst;++iTrkGrp) {
      
      TrackExtended * trkGrp = tracksInGroupFirst[iTrkGrp];
      TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i =0; i<hitVec.size(); ++i){
        hitvecFirst.push_back(hitVec[i]->getTrackerHit());        
      }
    }
  }
  if(groupSecond!=NULL){
    
    TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
    nTrkGrpSecond = int(tracksInGroupSecond.size());
    
    for (int iTrkGrp=0;iTrkGrp<nTrkGrpSecond;++iTrkGrp) {
      TrackExtended * trkGrp = tracksInGroupSecond[iTrkGrp];
      TrackerHitExtendedVec hitVec = 
      trkGrp->getTrackerHitExtendedVec();
      
      for(unsigned int i=0;i<hitVec.size();++i){
        hitvecSecond.push_back(hitVec[i]->getTrackerHit());
      }
    }
  }
  
  
  int nhitsFirst = (int)hitvecFirst.size();
  int nhitsSecond = (int)hitvecSecond.size();
  int count = 0;
  for(int i =0;i<nhitsFirst;++i){
    float xi = (float)hitvecFirst[i]->getPosition()[0];
    float yi = (float)hitvecFirst[i]->getPosition()[1];
    float ri = sqrt(xi*xi+yi*yi);
    if(ri < _tpc_inner_r || ri > _tpc_pad_height)continue;
    for(int j =0;j<nhitsSecond;++j){
      float xj = (float)hitvecSecond[j]->getPosition()[0];
      float yj = (float)hitvecSecond[j]->getPosition()[1];
      float rj = sqrt(xj*xj+yj*yj);
      if(fabs(ri-rj)<_tpc_pad_height/2.0)count++;
    }
  }  
  return count;
}

bool FullLDCTracking_MarlinTrk::VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt){
  
  
  const float d0First = firstTrackExt->getD0();
  const float z0First = firstTrackExt->getZ0();
  const float omegaFirst = firstTrackExt->getOmega();
  const float tanLFirst = firstTrackExt->getTanLambda();
  const float phi0First = firstTrackExt->getPhi();
  
  const float d0Second = secondTrackExt->getD0();
  const float z0Second = secondTrackExt->getZ0();
  const float omegaSecond = secondTrackExt->getOmega();
  const float tanLSecond = secondTrackExt->getTanLambda();
  const float phi0Second = secondTrackExt->getPhi();        
  
  HelixClass firstHelix;
  firstHelix.Initialize_Canonical(phi0First,d0First,z0First,omegaFirst,tanLFirst,_bField);
  const float pxFirst = firstHelix.getMomentum()[0];
  const float pyFirst = firstHelix.getMomentum()[1];
  const float pzFirst = firstHelix.getMomentum()[2];        
  
  HelixClass secondHelix;
  secondHelix.Initialize_Canonical(phi0Second,d0Second,z0Second,omegaSecond,tanLSecond,_bField);
  const float pxSecond = secondHelix.getMomentum()[0];
  const float pySecond = secondHelix.getMomentum()[1];
  const float pzSecond = secondHelix.getMomentum()[2];      
  const float pFirst = sqrt(pxFirst*pxFirst+pyFirst*pyFirst+pzFirst*pzFirst);
  const float pSecond = sqrt(pxSecond*pxSecond+pySecond*pySecond+pzSecond*pzSecond);
  
  //SJA:FIXME hardcoded cut 
  if(pFirst<2.5 || pSecond<2.5)return false;
  
  TrackExtended * combinedTrack = TrialCombineTracks(firstTrackExt,secondTrackExt);
  bool veto = false;
  if(combinedTrack!=NULL){
    //SJA:FIXME hardcoded cut 
    if(combinedTrack->getNDF()+15<firstTrackExt->getNDF()+secondTrackExt->getNDF()+5)veto=true;
    delete combinedTrack->getGroupTracks();
    delete combinedTrack;
  }else{
    veto = true;
  }
  if(SegmentRadialOverlap(firstTrackExt,secondTrackExt)>10)veto=true;
  return veto;
  
}


void FullLDCTracking_MarlinTrk::check(LCEvent * evt) { }

void FullLDCTracking_MarlinTrk::end() { 
  
  delete _encoder ;
  
}


void FullLDCTracking_MarlinTrk::setupGearGeom( const gear::GearMgr* gearMgr ){
  
  _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  if( _reading_loi_data ){
    
    
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
    
    _nLayersFTD = pFTD_z->size();
    
    for (int i = 0; i<_nLayersFTD; ++i) {
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
      
      // divide by two as we are treating the SIT as TWO stereo layers 
      _nLayersSIT = pSITLayerLayout.getNLayers() / 2.0 ; 
//      _SITgeo.resize(_nLayersSIT);
//      
//      //SJA:FIXME: for now the support is taken as the same size the sensitive
//      //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
//      //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
//      //           for a significant distance 
//      
//      for( unsigned int layer=0; layer < _nLayersSIT; ++layer){
//        _SITgeo[layer].nLadders = pSITLayerLayout.getNLadders(layer); 
//        _SITgeo[layer].phi0 = pSITLayerLayout.getPhi0(layer); 
//        _SITgeo[layer].dphi = 2*M_PI / _SITgeo[layer].nLadders; 
//        _SITgeo[layer].senRMin = pSITLayerLayout.getSensitiveDistance(layer); 
//        _SITgeo[layer].supRMin = pSITLayerLayout.getLadderDistance(layer); 
//        _SITgeo[layer].length = pSITLayerLayout.getSensitiveLength(layer); 
//        _SITgeo[layer].width = pSITLayerLayout.getSensitiveWidth(layer); 
//        _SITgeo[layer].offset = pSITLayerLayout.getSensitiveOffset(layer); 
//        _SITgeo[layer].senThickness = pSITLayerLayout.getSensitiveThickness(layer); 
//        _SITgeo[layer].supThickness = pSITLayerLayout.getLadderThickness(layer); 
//      }
      
    } 
    catch (gear::UnknownParameterException& e) {
      streamlog_out( MESSAGE ) << "  MarlinKalTest - SIT missing in gear file: SIT Not Built " << std::endl ;
    }
    
    try {
      
      
      const gear::FTDParameters& ftdParams = gearMgr->getFTDParameters() ;
      const gear::FTDLayerLayout& ftdlayers = ftdParams.getFTDLayerLayout() ;
      
      
      _nLayersFTD = ftdlayers.getNLayers() ; 
      _FTDgeo.resize(_nLayersFTD);
      
      //SJA:FIXME: for now the support is taken as the same size the sensitive
      //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
      //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
      //           for a significant distance 
      
      for(int disk=0; disk< _nLayersFTD; ++disk){
        
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
          streamlog_out( ERROR ) << "  FullLDCTracking_MarlinTrk - tilted design not supported exit(1) " << std::endl ;
          exit(1);
        }
        
      }
      
      for (int disk=0; disk < _nLayersFTD; ++disk) {
        _zLayerFTD.push_back(_FTDgeo[disk].senZPos_even_petal1); // front petal even numbered
        _zLayerFTD.push_back(_FTDgeo[disk].senZPos_odd_petal1);  // front petal odd numbered
      }
      
      // SJA as disks are staggered lets treat them internally as 2*ndisksFTD
      _nLayersFTD =_zLayerFTD.size() ;
      
    } 
    
    catch (gear::UnknownParameterException& e) {
      streamlog_out( MESSAGE ) << "  FullLDCTracking_MarlinTrk - FTD missing in gear file: FTD Not Built " << std::endl ;
    }
    
  }
  
  
  
}


//void FullLDCTracking_MarlinTrk::setupGearGeom( const gear::GearMgr* gearMgr ){
//  
//  _bField = gearMgr->getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
//  
//  
//  try {
//    
//    const gear::ZPlanarParameters& pVXDDetMain = gearMgr->getVXDParameters();
//    const gear::ZPlanarLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();
//    
//    _nLayersVTX = pVXDLayerLayout.getNLayers(); 
//    _VXDgeo.resize(_nLayersVTX);
//    
//    //SJA:FIXME: for now the support is taken as the same size the sensitive
//    //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
//    //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
//    //           for a significant distance 
//    
//    for( unsigned int layer=0; layer < _nLayersVTX; ++layer){
//      _VXDgeo[layer].nLadders = pVXDLayerLayout.getNLadders(layer); 
//      _VXDgeo[layer].phi0 = pVXDLayerLayout.getPhi0(layer); 
//      _VXDgeo[layer].dphi = 2*M_PI / _VXDgeo[layer].nLadders; 
//      _VXDgeo[layer].senRMin = pVXDLayerLayout.getSensitiveDistance(layer); 
//      _VXDgeo[layer].supRMin = pVXDLayerLayout.getLadderDistance(layer); 
//      _VXDgeo[layer].length = pVXDLayerLayout.getSensitiveLength(layer); 
//      _VXDgeo[layer].width = pVXDLayerLayout.getSensitiveWidth(layer); 
//      _VXDgeo[layer].offset = pVXDLayerLayout.getSensitiveOffset(layer); 
//      _VXDgeo[layer].senThickness = pVXDLayerLayout.getSensitiveThickness(layer); 
//      _VXDgeo[layer].supThickness = pVXDLayerLayout.getLadderThickness(layer); 
//    }
//    
//  } catch (gear::UnknownParameterException& e) {
//    streamlog_out( MESSAGE ) << "  MarlinKalTest - VXD missing in gear file: VXD Not Built " << std::endl ;
//  }
//  
//  try {
//    
//    const gear::ZPlanarParameters& pSITDetMain = gearMgr->getSITParameters();
//    const gear::ZPlanarLayerLayout& pSITLayerLayout = pSITDetMain.getZPlanarLayerLayout();
//    
//    _nLayersSIT = pSITLayerLayout.getNLayers(); 
//    _SITgeo.resize(_nLayersSIT);
//    
//    //SJA:FIXME: for now the support is taken as the same size the sensitive
//    //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
//    //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
//    //           for a significant distance 
//    
//    for( unsigned int layer=0; layer < _nLayersSIT; ++layer){
//      _SITgeo[layer].nLadders = pSITLayerLayout.getNLadders(layer); 
//      _SITgeo[layer].phi0 = pSITLayerLayout.getPhi0(layer); 
//      _SITgeo[layer].dphi = 2*M_PI / _SITgeo[layer].nLadders; 
//      _SITgeo[layer].senRMin = pSITLayerLayout.getSensitiveDistance(layer); 
//      _SITgeo[layer].supRMin = pSITLayerLayout.getLadderDistance(layer); 
//      _SITgeo[layer].length = pSITLayerLayout.getSensitiveLength(layer); 
//      _SITgeo[layer].width = pSITLayerLayout.getSensitiveWidth(layer); 
//      _SITgeo[layer].offset = pSITLayerLayout.getSensitiveOffset(layer); 
//      _SITgeo[layer].senThickness = pSITLayerLayout.getSensitiveThickness(layer); 
//      _SITgeo[layer].supThickness = pSITLayerLayout.getLadderThickness(layer); 
//    }
//    
//  } 
//  catch (gear::UnknownParameterException& e) {
//    streamlog_out( MESSAGE ) << "  MarlinKalTest - SIT missing in gear file: SIT Not Built " << std::endl ;
//  }
//  
//  try {
//    
//    
//    const gear::FTDParameters& ftdParams = gearMgr->getFTDParameters() ;
//    const gear::FTDLayerLayout& ftdlayers = ftdParams.getFTDLayerLayout() ;
//    
//    
//    _nLayersFTD = ftdlayers.getNLayers() ; 
//    _FTDgeo.resize(_nLayersFTD);
//    
//    //SJA:FIXME: for now the support is taken as the same size the sensitive
//    //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
//    //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
//    //           for a significant distance 
//    
//    for(int disk=0; disk< _nLayersFTD; ++disk){
//      
//      // numbers taken from the ILD_01 gear file for the sensitive part 
//      _FTDgeo[disk].nPetals = ftdlayers.getNPetals(disk) ;    
//      _FTDgeo[disk].dphi = 2*M_PI /  _FTDgeo[disk].nPetals ;
//      _FTDgeo[disk].phi0 = ftdlayers.getPhi0(disk) ;
//      _FTDgeo[disk].alpha = ftdlayers.getAlpha(disk) ;
//      _FTDgeo[disk].rInner = ftdlayers.getSensitiveRinner(disk) ;
//      _FTDgeo[disk].height = ftdlayers.getSensitiveWidth(disk) ;
//      _FTDgeo[disk].innerBaseLength =  ftdlayers.getSensitiveLengthMin(disk) ;
//      _FTDgeo[disk].outerBaseLength =  ftdlayers.getSensitiveLengthMax(disk) ;
//      _FTDgeo[disk].senThickness =  ftdlayers.getSensitiveThickness(disk) ;
//      _FTDgeo[disk].supThickness =  ftdlayers.getSupportThickness(disk) ;
//      
//      _FTDgeo[disk].senZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
//      _FTDgeo[disk].senZPos_even_petal2 = ftdlayers.getSensitiveZposition(disk, 0, 2) ; 
//      _FTDgeo[disk].senZPos_even_petal3 = ftdlayers.getSensitiveZposition(disk, 0, 3) ; 
//      _FTDgeo[disk].senZPos_even_petal4 = ftdlayers.getSensitiveZposition(disk, 0, 4) ; 
//      
//      // currently the design assumes that the petal on the same side are at the same z
//      assert(_FTDgeo[disk].senZPos_even_petal1==_FTDgeo[disk].senZPos_even_petal2);
//      assert(_FTDgeo[disk].senZPos_even_petal3==_FTDgeo[disk].senZPos_even_petal4);
//      
//      _FTDgeo[disk].senZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
//      _FTDgeo[disk].senZPos_odd_petal2 = ftdlayers.getSensitiveZposition(disk, 1, 2) ; 
//      _FTDgeo[disk].senZPos_odd_petal3 = ftdlayers.getSensitiveZposition(disk, 1, 3) ; 
//      _FTDgeo[disk].senZPos_odd_petal4 = ftdlayers.getSensitiveZposition(disk, 1, 4) ; 
//      
//      // currently the design assumes that the petal on the same side are at the same z
//      assert(_FTDgeo[disk].senZPos_odd_petal1==_FTDgeo[disk].senZPos_odd_petal2);
//      assert(_FTDgeo[disk].senZPos_odd_petal3==_FTDgeo[disk].senZPos_odd_petal4);
//      
//      _FTDgeo[disk].supZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
//      _FTDgeo[disk].supZPos_even_petal2 = ftdlayers.getSensitiveZposition(disk, 0, 2) ; 
//      _FTDgeo[disk].supZPos_even_petal3 = ftdlayers.getSensitiveZposition(disk, 0, 3) ; 
//      _FTDgeo[disk].supZPos_even_petal4 = ftdlayers.getSensitiveZposition(disk, 0, 4) ; 
//      
//      assert(_FTDgeo[disk].supZPos_even_petal1==_FTDgeo[disk].supZPos_even_petal2);
//      assert(_FTDgeo[disk].supZPos_even_petal3==_FTDgeo[disk].supZPos_even_petal4);
//      
//      _FTDgeo[disk].supZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
//      _FTDgeo[disk].supZPos_odd_petal2 = ftdlayers.getSensitiveZposition(disk, 1, 2) ; 
//      _FTDgeo[disk].supZPos_odd_petal3 = ftdlayers.getSensitiveZposition(disk, 1, 3) ; 
//      _FTDgeo[disk].supZPos_odd_petal4 = ftdlayers.getSensitiveZposition(disk, 1, 4) ; 
//      
//      assert(_FTDgeo[disk].supZPos_odd_petal1==_FTDgeo[disk].supZPos_odd_petal2);
//      assert(_FTDgeo[disk].supZPos_odd_petal3==_FTDgeo[disk].supZPos_odd_petal4);
//      
//      
//      
//      
//      
//      // rough check to see if the petal is rotated
//      if( fabs(_FTDgeo[disk].alpha) > fabs(FLT_MIN)  ) { 
//        streamlog_out( ERROR ) << "  FullLDCTracking_MarlinTrk - tilted design not supported exit(1) " << std::endl ;
//        exit(1);
//      }
//      
//    }
//    
//    for (int disk=0; disk < _nLayersFTD; ++disk) {
//      _zLayerFTD.push_back(_FTDgeo[disk].senZPos_even_petal1); // front petal even numbered
//      _zLayerFTD.push_back(_FTDgeo[disk].senZPos_odd_petal1);  // front petal odd numbered
//    }
//    
//    // SJA as disks are staggered lets treat them internally as 2*ndisksFTD
//    _nLayersFTD =_zLayerFTD.size() ;
//    
//  } 
//  
//  catch (gear::UnknownParameterException& e) {
//    streamlog_out( MESSAGE ) << "  FullLDCTracking_MarlinTrk - FTD missing in gear file: FTD Not Built " << std::endl ;
//  }
//  
//  try {
//    const gear::TPCParameters& pTPC = Global::GEAR->getTPCParameters();
//    const gear::PadRowLayout2D& pTPCpads = pTPC.getPadLayout();
//    _tpc_inner_r = pTPCpads.getPlaneExtent()[0];
//    _tpc_outer_r = pTPCpads.getPlaneExtent()[1];
//    _tpc_nrows = pTPCpads.getNRows();
//    _tpc_pad_height = (_tpc_outer_r-_tpc_inner_r)/(double)_tpc_nrows;
//    
//  } catch (gear::UnknownParameterException& e) {
//    streamlog_out( MESSAGE ) << "  MarlinKalTest - TPC missing in gear file: TPC Not Built " << std::endl ;
//  }
//  
//  
//  
//  
//}
//
//
