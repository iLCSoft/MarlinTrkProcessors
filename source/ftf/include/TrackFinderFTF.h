#ifndef TrackFinderFTF_h
#define TrackFinderFTF_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <UTIL/LCRelationNavigator.h>

#include <EVENT/TrackerHit.h>

#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"


namespace MarlinTrk{
  class IMarlinTrkSystem ;
}

namespace ftf {
  class TrackFinder ;
  class TrackFindingParameters ;
}

/**  Track Finder using FTF processor for marlin. 
 * 
 
 *  <h4>Input - Prerequisites</h4>
 *  Needs a collection of LCIO TrackerHits. 
 *
 *  <h4>Output</h4> 
 *  LCIO Track Collection
 * 
 * @param InputTrackerHitCollectionName Name of the Track collection to be refitted
 * @param OutputTrackCollectionName Name of the Track collection found
 * 
 * @author S. J. Aplin, DESY
 */

class TrackFinderFTF : public marlin::Processor {
  
public:
  
  virtual marlin::Processor*  newProcessor() { return new TrackFinderFTF ; }
  
  TrackFinderFTF() ;
  
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
  
  
  struct compare_r {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      //       double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      //       double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      //       return ( r_a_sqd < r_b_sqd ) ; 
      //fg: try to work around a compiler issue with optimization that leads to seg faults with above code...
      //    needs further investigation
      return ( a == b  ?  false : 
              ( a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ) < 
              ( b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ) ) ; 
    }
  } ;
  
  struct compare_r_reverse {
    bool operator()( EVENT::TrackerHit* a, EVENT::TrackerHit* b)  const { 
      //       double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ; 
      //       double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ; 
      //       return ( r_a_sqd < r_b_sqd ) ; 
      //fg: try to work around a compiler issue with optimization that leads to seg faults with above code...
      //    needs further investigation
      return ( a == b  ?  false : 
              ( a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1] ) > 
              ( b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1] ) ) ; 
    }
  } ;

  
  
  
protected:

  UTIL::BitField64* _encoder;
  int getDetectorID(TrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::subdet()]; }
  int getSideID(TrackerHit* hit)     { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::side()]; };
  int getLayerID(TrackerHit* hit)    { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::layer()]; };
  int getModuleID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::module()]; };
  int getSensorID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::sensor()]; };

  
  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection( lcio::LCEvent * evt, std::string colName ) ;

  double angular_range_2PI( double phi ) const ;
  
  void setFTFParameters ( ftf::TrackFindingParameters* para ) ;
  

  /** Input tpc tracker hit collection name
   */
  std::string _input_vxd_hits_col_name ;

  /** Input sit tracker hit collection name
   */
  std::string _input_sit_hits_col_name ;

  /** Input ftd pixel tracker hit collection name
   */
  std::string _input_ftd_pixel_hits_col_name ;
  
  /** Input ftd spacepoint tracker hit collection name
   */
  std::string _input_ftd_spacepoint_hits_col_name ;

//  /** Input set tracker hit collection name
//   */
//  std::string _input_set_hits_col_name ;
//
//  /** Input tpc tracker hit collection name
//   */
//  std::string _input_tpc_hits_col_name ;

  
  /** output track collection name.
   */
  std::string _output_track_col_name ;

  
  lcio::LCCollection* _input_vxd_hits_col;
  lcio::LCCollection* _input_sit_hits_col;
  lcio::LCCollection* _input_ftd_pxl_hits_col;
  lcio::LCCollection* _input_ftd_spp_hits_col;

//  lcio::LCCollection* _input_tpc_hits_col;
//  lcio::LCCollection* _input_set_hits_col;

  int _n_run ;
  int _n_evt ;
  
  bool _parameter1;
  
  ftf::TrackFinder* _trackFinder;

  std::map<int, TrackerHit*> hitmap;

  int _nlayers_vxd;
  int _nlayers_sit;
  int _nlayers_ftd;
  
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;
  
  bool _MSOn ;
  bool _ElossOn ;
  bool _SmoothOn ;

  float _initialTrackError_d0;
  float _initialTrackError_phi0;
  float _initialTrackError_omega;
  float _initialTrackError_z0;
  float _initialTrackError_tanL;
  
  double _maxChi2PerHit;
  
  bool _runMarlinTrkDiagnostics;
  std::string _MarlinTrkDiagnosticsName;
  
  double _Bz;
  
  //  int _nlayers_tpc;
//  int _nlayers_set;
  
} ;

#endif



