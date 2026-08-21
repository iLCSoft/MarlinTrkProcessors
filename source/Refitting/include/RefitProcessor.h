#ifndef RefitProcessor_h
#define RefitProcessor_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

#include <UTIL/LCRelationNavigator.h>

#include <EVENT/TrackerHit.h>

#include <memory>

namespace MarlinTrk {
class IMarlinTrkSystem;
}

namespace dd4hep {
class Detector;
}

/**  Track Refitter processor for marlin. Refits an input track collection, producing a new collection of tracks
 *

 *  <h4>Input - Prerequisites</h4>
 *  Needs a collection of LCIO Tracks.
 *
 *  <h4>Output</h4>
 *  Refitted LCIO Track Collection
 *
 * @param InputTrackCollectionName Name of the Track collection to be refitted
 * @param OutputTrackCollectionName Name of the refitted Track collection to be refitted
 * @param TrackSystemName name of the track fitting system to be used (KalTest, DDKalTest, aidaTT, ... )
 *
 * @author S. J. Aplin, DESY
 * @history
 *   Nov 2014 F.Gaede CERN/DESY added processor parameter to instantiate other implementations of IMarlinTrk
 */

class RefitProcessor : public marlin::Processor {
public:
  virtual marlin::Processor* newProcessor() { return new RefitProcessor; }

  RefitProcessor();
  RefitProcessor(const RefitProcessor&) = delete;
  RefitProcessor& operator=(const RefitProcessor&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(lcio::LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(lcio::LCEvent* evt);

  virtual void check(lcio::LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  struct compare_r {
    bool operator()(EVENT::TrackerHit* a, EVENT::TrackerHit* b) const {
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1];
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1];
      return (r_a_sqd < r_b_sqd);
    }
  };

protected:
  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection(lcio::LCEvent* evt, std::string colName);

  /* helper function to get relations using try catch block */
  std::unique_ptr<lcio::LCRelationNavigator> GetRelations(lcio::LCEvent* evt, std::string RelName);

  /* helper function to look up subdetector IDs by name from the DD4hep geometry.
   * Each entry is 0 (i.e. "not used") with a warning printed if no subdetector with that name exists.
   */
  std::vector<int> GetSubDetIDs(dd4hep::Detector& detector, const std::vector<std::string>& detNames);

  // helper function to identify ILD@FCC-ee detector models
  bool IsFCCeeModel(dd4hep::Detector& detector) const;

  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name{};

  /** Input track relations name for refitting.
   */
  std::string _input_track_rel_name{};

  /** refitted track collection name.
   */
  std::string _output_track_col_name{};

  /** Output track relations name for refitting.
   */
  std::string _output_track_rel_name{};

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem{nullptr};

  bool _MSOn{};
  bool _ElossOn{};
  bool _SmoothOn{};

  float _initialTrackError_d0{};
  float _initialTrackError_phi0{};
  float _initialTrackError_omega{};
  float _initialTrackError_z0{};
  float _initialTrackError_tanL{};
  float _maxChi2PerHit{};
  double _mass{};

  int _n_run{};
  int _n_evt{};

  int _initialTrackState{};
  int _fitDirection{};

  std::string _trkSystemName{};

  float _bField{};

  /** subdetector IDs contributing to each hit-count category, looked up by name
   *  from the DD4hep geometry
   */
  std::vector<int> _vxdIDs{};
  std::vector<int> _ftdIDs{};
  std::vector<int> _sitIDs{};
  std::vector<int> _tpcIDs{};
  std::vector<int> _setIDs{};
};

#endif
