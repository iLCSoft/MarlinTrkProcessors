#ifndef ExtrToTracker_h
#define ExtrToTracker_h 1

#include <marlin/Processor.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <time.h>
#include <math.h>
#include <sstream>
#include <assert.h>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/Track.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackImpl.h>

#include "KiTrack/SubsetHopfieldNN.h"
#include "KiTrack/SubsetSimple.h"

//---- GEAR ----
#include "marlin/Global.h"
#include "gear/GEAR.h"
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>
#include "gear/BField.h"

#if defined GEO2
#include "DDRec/Surface.h"

#endif

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"




class TTree;
class TFile;



using namespace KiTrack;
namespace MarlinTrk{
  class IMarlinTrkSystem ;
}


class ExtrToTracker : public marlin::Processor {

#if defined GEO2
    typedef std::map< unsigned long, const DD4hep::DDRec::Surface* > SurfaceMap ;
    
#endif

public:


  virtual marlin::Processor*  newProcessor() { return new ExtrToTracker ; }
  
  ExtrToTracker() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( lcio::LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( lcio::LCEvent * evt ) ; 
  
  
  virtual void check( lcio::LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  void  SelectBestCandidate(EVENT::TrackerHitVec &HitsInLayer, const float* &pivot, EVENT::TrackerHit* &BestHit, bool &BestHitFound, int &pointer) ; 

  void SelectBestCandidateLimited(EVENT::TrackerHitVec &HitsInLayer, const float* &pivot, EVENT::TrackerHit* &BestHit, const FloatVec& covLCIO, double& radius, bool &BestHitFound, double &sigma, int &pointer, int &PossibleHits, float &dU, float &dV, double &DimDist, TrackerHitVec &usedSiHits) ;
  
  int FitInit( std::vector < TrackerHit* > trackerHits , MarlinTrk::IMarlinTrack* _marlinTrk ) ;
  int FitInit2( Track* track , MarlinTrk::IMarlinTrack* _marlinTrk ) ;

  struct compare_r {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      return ( r_a_sqd < r_b_sqd ) ; 
    }
  } ;
  
  
  void setTreeBranches(int bufsize=32000);
  void clearEventVar();
  void clearLayerHelperVar();
  void fillDummy();
  TrackerHit* getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk, int& nHitsOnDetEl);
  
  //bool getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk, TrackerHit*& selectedHit);
  
  
protected:
  
  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection( lcio::LCEvent * evt, std::string colName ) ;
  
  /* helper function to get relations using try catch block */
  lcio::LCRelationNavigator* GetRelations(lcio::LCEvent * evt, std::string RelName ) ;
  
  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name ;
  
  /** Input track relations name for refitting.
   */
  std::string _input_track_rel_name ;

  /** Input SIT tracker summer hit collection.
   */
  std::string _sitColName ;

  /** Input VXD tracker summer hit collection.
   */
  std::string _vxdColName ;
  
  /** refitted track collection name.
   */
  std::string _output_track_col_name ;
  
  /** Output track relations name for refitting.
   */
  std::string _output_track_rel_name ;

  /** Output silicon track collection.
   */
  std::string _siTrkColName ;
  
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;
  
  std::string _mcParticleCollectionName ;

  bool _MSOn ;
  bool _ElossOn ;
  bool _SmoothOn ;
  double _Max_Chi2_Incr ;
  int _tpcHitsCut ;
  float _chi2NDoFCut ;
  float _DoCut ;
  float _ZoCut ;
  double _searchSigma ;
  bool _isSpacePoints ;
  int _propToLayer ;
  
  int _n_run ;
  int _n_evt ;
  int SITHitsFitted ;
  int SITHitsNonFitted ;
  int TotalSITHits ;
  int _nHitsChi2 ;

  float _maxChi2PerHit;
  float _bField;

  StringVec  _colNamesTrackerHitRelations ;
  bool _performFinalRefit ;
 

#if defined GEO2
  SurfaceMap _map ;
  
#endif

  

  //processor parameters
  bool _doNtuple;
  std::string _outFileName;
  std::string _treeName;
  
 //for tree

  TFile* _out; 
  TTree* _tree;

  int _hitVXDN;

  std::vector<std::vector<double > > _hitN;
  std::vector<std::vector<double > > _hitX;
  std::vector<std::vector<double > > _hitY;
  std::vector<std::vector<double > > _hitZ;
  std::vector<std::vector<double > > _hitR;
  std::vector<std::vector<double > > _hitU;
  std::vector<std::vector<double > > _hitV;
  std::vector<std::vector<double > > _hitW;
  std::vector<std::vector<double > > _hitID;

  int _fitN;
  std::vector<std::vector<double > > _fitX;
  std::vector<std::vector<double > > _fitY;
  std::vector<std::vector<double > > _fitZ;
  std::vector<std::vector<double > > _fitR;
  std::vector<std::vector<double > > _fitU;
  std::vector<std::vector<double > > _fitV;
  std::vector<std::vector<double > > _fitW;
  std::vector<std::vector<double > > _fitID;

  std::vector<std::vector<double > > _fitD0Err;
  std::vector<std::vector<double > > _fitD0;


  std::vector<std::vector<double > > _fitZ0Err;
  std::vector<std::vector<double > > _fitZ0;



  std::vector<double > _hitN_layer_helper;
  std::vector<double > _hitX_layer_helper;
  std::vector<double > _hitY_layer_helper;
  std::vector<double > _hitZ_layer_helper;
  std::vector<double > _hitR_layer_helper;
  std::vector<double > _hitU_layer_helper;
  std::vector<double > _hitV_layer_helper;
  std::vector<double > _hitW_layer_helper;
  std::vector<double > _hitID_layer_helper;

  std::vector<double > _fitX_layer_helper;
  std::vector<double > _fitY_layer_helper;
  std::vector<double > _fitZ_layer_helper;
  std::vector<double > _fitR_layer_helper;
  std::vector<double > _fitU_layer_helper;
  std::vector<double > _fitV_layer_helper;
  std::vector<double > _fitW_layer_helper;
  std::vector<double > _fitID_layer_helper;

  std::vector<double > _fitD0Err_layer_helper;
  std::vector<double > _fitD0_layer_helper;

  std::vector<double > _fitZ0Err_layer_helper;
  std::vector<double > _fitZ0_layer_helper;

  
  
} ;





#endif



