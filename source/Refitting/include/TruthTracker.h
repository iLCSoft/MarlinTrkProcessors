#ifndef TruthTracker_h
#define TruthTracker_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

namespace EVENT{
  class MCParticle ;
  class Track ;
}

namespace IMPL {
  class TrackImpl ;
}

namespace UTIL{
  class LCRelationNavigator ;
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

class TruthTracker : public marlin::Processor {
  
public:
  
  virtual marlin::Processor*  newProcessor() { return new TruthTracker ; }
  
  
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
  
  struct SimTrackerHitSortPredicate {
    bool operator()(std::pair<SimTrackerHit*, TrackerHit* > p1, std::pair<SimTrackerHit*, TrackerHit* > p2) {
      
      SimTrackerHit* simHit1 = p1.first;
      SimTrackerHit* simHit2 = p1.first;
      
      
      if( simHit1->getMCParticle() == simHit2->getMCParticle() ) {
        return simHit1->getTime() < simHit2->getTime() ;
      }
      else { 
        return simHit1->getMCParticle() < simHit2->getMCParticle() ;
      }
    }
  };

  
//  struct compare_time {
//    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
//      SimTrackerHit* hita = dynamic_cast<SimTrackerHit*>(a->getRawHits().at(0));
//      SimTrackerHit* hitb = dynamic_cast<SimTrackerHit*>(b->getRawHits().at(0));
//      return ( hita->getTime() < hitb->getTime() ) ; 
//    }
//  };
  
  
  struct compare_r {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      return ( r_a_sqd < r_b_sqd ) ; 
    }
  } ;
  
  struct compare_time_reverse {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      return ( r_a_sqd > r_b_sqd ) ; 
    }
  } ;

  


  
protected:
  
  
  const LCObjectVec* getSimHits( TrackerHit* trkhit, const FloatVec* weights = NULL);
  
  UTIL::BitField64* _encoder;
  int getDetectorID(TrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::subdet]; }
  int getSideID(TrackerHit* hit)     { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::side]; };
  int getLayerID(TrackerHit* hit)    { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::layer]; };
  int getModuleID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::module]; };
  int getSensorID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::ILDCellID0::sensor]; };

  
  /* helper function to get collection using try catch block */
  LCCollection* GetCollection(  LCEvent * evt, std::string colName ) ;
  
  /* helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations( LCEvent * evt, std::string RelName ) ;
  
  /* helper function to get collection using try catch block */
  void SetupInputCollections( LCEvent * evt ) ;
  
  void createTrack( MCParticle* mcp, UTIL::BitField64& cellID_encoder );
  


  
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
  
 
  std::string _vxdTrackerHitRelInputColName;
  std::string _ftdTrackerHitRelInputColName;
  std::string _sitTrackerHitRelInputColName;
  std::string _tpcTrackerHitRelInputColName;
  std::string _setTrackerHitRelInputColName;
  std::string _etdTrackerHitRelInputColName;

  LCRelationNavigator* _navVXDTrackerHitRel;
  LCRelationNavigator* _navSITTrackerHitRel;
  LCRelationNavigator* _navFTDTrackerHitRel;
  LCRelationNavigator* _navTPCTrackerHitRel;
  LCRelationNavigator* _navSETTrackerHitRel;
  LCRelationNavigator* _navETDTrackerHitRel;

  
  /** output track collection 
   */
  std::string _output_track_col_name ;
  LCCollectionVec* _trackVec;
  
  /** Output track relations
   */
  std::string _output_track_rel_name ;
  LCCollectionVec* _trackRelVec;
  
  std::map< MCParticle*, std::vector<TrackerHit*> > _MCParticleTrkHitMap;
  
  std::vector<TrackerHit*> _hit_list;
  
  std::map<int, LCRelationNavigator*> _hit_rels_map;
  
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



