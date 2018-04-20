#ifndef ClonesAndSplitTracksFinder_h
#define ClonesAndSplitTracksFinder_h 1

#include <marlin/Processor.h>

#include <UTIL/BitField64.h>

#include <EVENT/Track.h>
#include <IMPL/TrackImpl.h>

#include <cfloat>

namespace MarlinTrk {
class IMarlinTrkSystem;
class IMarlinTrack;
}

class ClonesAndSplitTracksFinder : public marlin::Processor {

public:
  virtual marlin::Processor *newProcessor() { return new ClonesAndSplitTracksFinder; }

  ClonesAndSplitTracksFinder();
  ClonesAndSplitTracksFinder(const ClonesAndSplitTracksFinder &) = delete;
  ClonesAndSplitTracksFinder &operator=(const ClonesAndSplitTracksFinder &) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init();

  // Called at the beginning of every run
  virtual void processRunHeader(lcio::LCRunHeader *run);

  // Run over each event - the main algorithm
  virtual void processEvent(lcio::LCEvent *evt);

  // Run at the end of each event
  virtual void check(lcio::LCEvent *evt);

  // Called at the very end for cleanup, histogram saving, etc.
  virtual void end();

protected:

  // Checks for overlapping hits
  int overlappingHits(const Track*, const Track*);
  
  // Picks up the best track between two clones (based on chi2 and length requirements)
  void bestInClones(Track*, Track*, int, Track*&);

  // Service function to set the information from a Track* object to a TrackImpl* object
  void fromTrackToTrackImpl(const Track*, TrackImpl*&);

  // Merges hits from two tracks in one and fits it
  void mergeAndFit(Track*, Track*, Track*&);

  // Removes doubles (from clone treatments and track merging) and filters multiple connections (clones and mergeable tracks treated differently)
  void filterClonesAndMergedTracks(std::multimap<int,std::pair<int,Track*>>&, LCCollection*&, TrackVec&, bool);

  lcio::LCCollection *GetCollection(lcio::LCEvent *evt, std::string colName);

  std::string _input_track_col_name;
  std::string _output_track_col_name;

  MarlinTrk::IMarlinTrkSystem *_trksystem = nullptr;

  int _n_run = -1;
  int _n_evt = -1;

  bool _MSOn = true;
  bool _ElossOn = true;
  bool _SmoothOn = false;
  double _magneticField = 0.0;
  bool _extrapolateForward = true;

  double _maxDeltaTheta = 0.0, _maxDeltaPhi = 0.0, _maxDeltaPt = 0.0;

  bool _clonesOnly = true;

  // Track fit parameters
  double _initialTrackError_d0;
  double _initialTrackError_phi0;
  double _initialTrackError_omega;
  double _initialTrackError_z0;
  double _initialTrackError_tanL;
  double _maxChi2perHit;

  std::shared_ptr<UTIL::BitField64> _encoder{};

};

bool sort_by_r(EVENT::TrackerHit*, EVENT::TrackerHit*);

#endif
