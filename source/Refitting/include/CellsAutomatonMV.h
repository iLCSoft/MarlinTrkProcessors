#ifndef CellsAutomatonMV_h
#define CellsAutomatonMV_h 1

#include <algorithm>
#include "Math/ProbFunc.h"

#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "marlin/VerbosityLevels.h"
#include <marlin/Exceptions.h>

#include <EVENT/LCCollection.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

// KiTrack tools
#include "KiTrack/SubsetHopfieldNN.h"
#include "KiTrack/SubsetSimple.h"
#include "KiTrack/SegmentBuilder.h"
#include "KiTrack/Automaton.h"
#include "KiTrack/Segment.h"

// KiTrackMarlin toools
#include "ILDImpl/VXDTrack.h"
#include "ILDImpl/VXDHit01.h"
#include "ILDImpl/VXDSectorConnector.h"
#include "Tools/KiTrackMarlinTools.h"
#include "Tools/KiTrackMarlinCEDTools.h"
#include "KiTrack/ITrack.h"
#include "Criteria/Criteria.h"
#include "ILDImpl/SectorSystemVXD.h"
#include "Tools/VXDHelixFitter.h"
#include "ILDImpl/MiniVectorHit01.h"

// IMarlin tools
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/MarlinTrkUtils.h"

// GEAR tools
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include <gear/BField.h>
#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

//#include "SpacePointBuilder.h"
// CLHEP tools
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"


using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace gear ;
using namespace KiTrack;
using namespace KiTrackMarlin;
namespace MarlinTrk{
  class IMarlinTrkSystem ;
}

namespace gear{
class GearMgr ;
}
/** a simple typedef, making writing shorter. And it makes sense: a track consists of hits. But as a real track
 * has more information, a vector of hits can be considered as a "raw track". */
typedef std::vector< IHit* > RawTrack;

class CellsAutomatonMV : public Processor {
  
 public:
 
  virtual Processor*  newProcessor() { return new CellsAutomatonMV ; }
  
  
  CellsAutomatonMV() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  /*
    private:  

    unsigned _nOutOfBoundary;
    unsigned _nStripsTooParallel;
    unsigned _nPlanesNotParallel;
    
    // functions copied from spacepointbuilder
    static int calculatePointBetweenTwoLines_UsingVertex( 
    const CLHEP::Hep3Vector& PA, 
    const CLHEP::Hep3Vector& PB, 
    const CLHEP::Hep3Vector& PC, 
    const CLHEP::Hep3Vector& PD,
    const CLHEP::Hep3Vector& Vertex,
    CLHEP::Hep3Vector& point);
    
    
    // @return a spacepoint (in the form of a TrackerHitImpl* ) created from two TrackerHitPlane* which stand for si-strips 
    TrackerHitImpl* createSpacePoint( TrackerHitPlane* a , TrackerHitPlane* b, double stripLength );
*/
  
 protected:

  int nEvt;

  int _nDivisionsInPhi;
  int _nDivisionsInTheta;
  int _nDivisionsInPhiMV;
  int _nDivisionsInThetaMV;
  int _nLayers;

  float _bField;

  // two pi is not a constant in cmath. Calculate it, once!
  static const double TWOPI;
  double _dPhi;
  double _dTheta;

  UTIL::BitField64* _encoder;
  int getDetectorID(TrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::subdet]; }
  int getSideID(TrackerHit* hit)     { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::side]; };
  int getLayerID(TrackerHit* hit)    { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::layer]; };
  int getModuleID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::module]; };
  int getSensorID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::sensor]; };
  

  void InitialiseVTX(LCEvent * evt, EVENT::TrackerHitVec HitsTemp);
  void setupGearGeom( const gear::GearMgr* gearMgr ) ;
  bool setCriteria( unsigned round );
  void RawTrackFit( std::vector < MarlinTrk::IMarlinTrack* > candMarlinTracks, std::vector< IMPL::TrackImpl* > &finalTracks ) ;
  void FitFunc2( std::vector < RawTrack > rawTracks, std::vector < MarlinTrk::IMarlinTrack* > &candMarlinTracks ) ;
  void finaliseTrack( TrackImpl* newTrackImpl, LCCollectionVec* trackVec ) ;   
  //void finaliseTrack( TrackImpl* trackImpl, LCCollectionVec* trackVec, const EVENT::TrackState * init_ts,  const EVENT::FloatVec covMatrix ) ;
  void CreateMiniVectors( int sector ) ;
  bool thetaAgreement( EVENT::TrackerHit *toHit, EVENT::TrackerHit *fromHit ) ;
  bool thetaAgreementImproved( EVENT::TrackerHit *toHit, EVENT::TrackerHit *fromHit, int layer ) ;
  double Dist( EVENT::TrackerHit *toHit, EVENT::TrackerHit *fromHit ) ;
  //IMPL::TrackImpl* Refit( TrackImpl* track_to_refit ) ;


  unsigned int _nLayersVTX;
  unsigned int _nLayersSIT;
  
  /** A map to store the hits according to their sectors */
  std::map< int , EVENT::TrackerHitVec > _map_sector_spacepoints;
  std::map< int , std::vector< IHit* > > _map_sector_hits;
  std::map< int , EVENT::TrackerHit* > _map_1dhits_spacepoints ;

  /** A pair to keep information for composite spacepoints */
  std::pair <EVENT::TrackerHit*, EVENT::TrackerHit*  > _pairs_1dhits_spacepoints ;
  
  /** Names of the used criteria */
  std::vector< std::string > _criteriaNames;
   
  /** Map containing the name of a criterion and a vector of the minimum cut offs for it */
  std::map< std::string , std::vector<float> > _critMinima;
  
  /** Map containing the name of a criterion and a vector of the maximum cut offs for it */
  std::map< std::string , std::vector<float> > _critMaxima;
  
  /** Minimum number of hits a track has to have in order to be stored */
  int _hitsPerTrackMin;
  
  /** A vector of criteria for 2 hits (2 1-hit segments) */
  std::vector <ICriterion*> _crit2Vec;
  
  /** A vector of criteria for 3 hits (2 2-hit segments) */
  std::vector <ICriterion*> _crit3Vec;
  
  /** A vector of criteria for 4 hits (2 3-hit segments) */
  std::vector <ICriterion*> _crit4Vec;

  std::vector< IHit* > MiniVectorsTemp;
  std::vector< IHit* > TestMiniVectorsTemp;

  /** Cut for the Kalman Fit (the chi squared probability) */
  double _chi2ProbCut;  

  double _helixFitMax ;

  double _chi2OverNdfCut ;

  const SectorSystemVXD * _sectorSystemVXD;
   
  /** the maximum number of connections that are allowed in the automaton, if this value is surpassed, rerun
   * the automaton with tighter cuts or stop it entirely. */
  int _maxConnectionsAutomaton;

  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trkSystem;
  std::string _trkSystemName ;
  
  bool _MSOn, _ElossOn, _SmoothOn, _middleLayer ;

  int _useSIT ;
  int _ipHit ;

  int _layerStepMax ;

  int _lastLayerToIP ;

  int _nHitsChi2 ;

  int MiniVectors_sectors ;
  int MiniVectors_CutSelection ;

  FloatVec _resU ;

  double _maxDist ;
  double _hitPairThDiff ;
  double _hitPairThDiffInner ;
  double _hitPairThDiffSIT ;
  
  //std::vector< MarlinTrk::IMarlinTrack* > GoodTracks;
  //std::vector< MarlinTrk::IMarlinTrack* > RejectedTracks;


  
  float _initialTrackError_d0;
  float _initialTrackError_phi0;
  float _initialTrackError_omega;
  float _initialTrackError_z0;
  float _initialTrackError_tanL;
  float _maxChi2PerHit;

  int _initialTrackState;
  int _fitDirection ; 

  int _maxHitsPerSector ;
  
  /** Input collection name.
   */
  std::string _VTXHitCollection;
  std::string _SITHitCollection;
  std::string _CATrackCollection;
  
  std::map< LCCollection*, std::string > _colNamesTrackerHits;
  
  std::string _bestSubsetFinder;

   /** The quality of the output track collection */
   int _output_track_col_quality ; 
  
   static const int _output_track_col_quality_GOOD;
   static const int _output_track_col_quality_FAIR;
   static const int _output_track_col_quality_POOR;


} ;



//****************************************************************************************************
// Quality - compatibility for ITrack version


/** A functor to return whether two tracks are compatible: The criterion is if they share a MiniVector or more */
class TrackCompatibilityShare1_MV{
  
public:
   
   inline bool operator()( ITrack* trackA, ITrack* trackB ){
      
      
      std::vector< IHit* > hitsA = trackA->getHits();
      std::vector< IHit* > hitsB = trackB->getHits();
     
      for( unsigned i=0; i < hitsA.size(); i++){
         
         for( unsigned j=0; j < hitsB.size(); j++){
            
            if ( hitsA[i] == hitsB[j] ) return false;      // a hit is shared -> incompatible
            
         }
         
      }
      
      return true;      
      
   }
   
};


/** A functor to return the quality of a track, which is currently the chi2 probability. */
class TrackQIChi2Prob_MV{
   
public:
   
   inline double operator()( ITrack* track ){ return track->getChi2Prob(); }
  
   
};

/** A functor to return the quality of a track, which is the ratio chi2 over degrees of freedom,
    weighted with the number of associated hits. 
*/
class TrackQI{
   
public:
   
  inline double operator()( ITrack* track ){ 

    float NoOfHits = 4.0*( track->getHits().size());
    
    return (1.0*NoOfHits)/(track->getChi2()/track->getNdf()); }
  
   
};


/** A functor to return the quality of a track, which is the number of associated hits. 
*/
class MaxHits{
   
public:
   
  inline double operator()( ITrack* track ){
    
    return track->getHits().size();
  }
   
};


/** A functor to return the quality of a track.
 For tracks with 4 hits or more the chi2prob is mapped to* 0.5-1, with x = prob/2 + 0.5.
 Tracks with 3 hits get the chi2 mapped to 0-0.5 by 1/(ln( e^2 + chi2 );
 That gives 0 for an infinite chi2 and 0.5 for a chi2 of 0.
 
 Reason: now 3-hit-tracks can be compared as well
 */
class TrackQISpecial_MV{
   
public:
   
   inline double operator()( ITrack* track ){ 
      
      if( track->getHits().size() > 2 ){
         
         return track->getChi2Prob()/2. +0.5; 
         
      }
      else{
         
         return 1/( log( 7.3890561 + track->getChi2() ) ); //e^2 = 7.3890561
         
      }
      
   }
   
   
};



#endif
