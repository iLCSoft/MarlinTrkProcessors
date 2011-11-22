#ifndef TruthTracker_h
#define TruthTracker_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>

using namespace lcio ;
using namespace marlin ;

namespace EVENT{
  class MCParticle ;
  class Track ;
}

namespace IMPL {
  class TrackImpl ;
}

namespace UTIL{
  class LCRelationNavigator ;
  class BitField64 ;
}


namespace MarlinTrk{
  class IMarlinTrkSystem ;
}

/**  Track creation based on MC truth. 
 * 
 
 *  <h4>Input - Prerequisites</h4>
 *  Needs a collections of LCIO TrackerHits. 
 *
 *  <h4>Output</h4> 
 *  LCIO Track Collection
 * 
 * @param InputTrackerHitCollectionName Name of the TrackerHit collections 
 * @param OutputTrackCollectionName Name of the output Track collection
 * 
 * @author S. J. Aplin, DESY
 */

class TruthTracker : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new TruthTracker ; }
  
  
  TruthTracker() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
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
  
  struct compare_time {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      SimTrackerHit* hita = dynamic_cast<SimTrackerHit*>(a->getRawHits().at(0));
      SimTrackerHit* hitb = dynamic_cast<SimTrackerHit*>(b->getRawHits().at(0));
      return ( hita->getTime() < hitb->getTime() ) ; 
    }
  };
  
  
  struct compare_r {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      return ( r_a_sqd < r_b_sqd ) ; 
    }
  } ;

  
protected:
  
  /* helper function to get collection using try catch block */
  LCCollection* GetCollection(  LCEvent * evt, std::string colName ) ;
  
  /* helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations( LCEvent * evt, std::string RelName ) ;
  
  /* helper function to get collection using try catch block */
  void SetupInputCollections( LCEvent * evt ) ;
  
  TrackImpl* createTrack( MCParticle* mcp, UTIL::BitField64& cellID_encoder );
  
  /** input MCParticle collection name.
   */
  std::string _colNameMC ;
  
  LCCollection* _colMCP;
  
  /** input SimTrackerHit collections
   */
  std::string _colNameSimTrkHitsVTX;
  std::string _colNameSimTrkHitsSIT;
  std::string _colNameSimTrkHitsFTD;
  std::string _colNameSimTrkHitsTPC;
  
  LCCollection* _VTXSimHits;
  LCCollection* _SITSimHits;
  LCCollection* _FTDSimHits;
  LCCollection* _TPCSimHits;
  
  /** input TrackerHit collections
   */
  std::string _colNameTrkHitsVTX ;
  std::string _colNameTrkHitsFTD ;
  std::string _colNameTrkHitsSIT ;
  std::string _colNameTrkHitsTPC ;
  
  LCCollection* _VTXTrkHits;
  LCCollection* _SITTrkHits;
  LCCollection* _FTDTrkHits;
  LCCollection* _TPCTrkHits;
  
  
  
  /** output track collection 
   */
  std::string _output_track_col_name ;
  
  /** Output track relations
   */
  std::string _output_track_rel_name ;
  
  std::map< MCParticle*, std::vector<TrackerHit*> > _MCParticleTrkHitMap;
  
  std::vector<TrackerHit*> _hit_list;
  
  int _nMCP;
  
  int _nVTXSimHits ;
  int _nSITSimHits ;
  int _nFTDSimHits ;
  int _nTPCSimHits ;
  
  int _nVTXTrkHits ;
  int _nSITTrkHits ;
  int _nFTDTrkHits ;
  int _nTPCTrkHits ;
  
  int _nEventPrintout ;
  int _n_run ;
  int _n_evt ;
  
  float _MCpThreshold ;
  bool  _FitTracksWithMarlinTrk;
  
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;
  
  bool _MSOn ;
  bool _ElossOn ;
  bool _SmoothOn ;


  
} ;

#endif



