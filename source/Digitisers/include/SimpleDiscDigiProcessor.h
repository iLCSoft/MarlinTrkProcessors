/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef SimpleDiscDigiProcessor_h
#define SimpleDiscDigiProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

using namespace lcio ;
using namespace marlin ;

class FTDLayer ;


/** ======= SimpleDiscDigiProcessor ========== <br>
 * Produces a TrackerHit collection from SimTrackerHit collection. <br> 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in (x,y) plane according to the specified point resolution. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits <br>
 * @param CollectionName The name of input collection of the SimTrackerHits <br>
 * (default name FTDCollection) <br>
 * @param OutputCollectionName The name of output collection of TrackerHits <br>
 * (default name FTDTrackerHits) <br> 
 * @param PointResolution Point resolution in (x,y) for the planar detectors (in mm) <br>
 * (default value 0.01) <br>
 * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio <br>
 * (default value ILDDetID::FTD) <br>
 * @param keepHitsFromDeltas Whether to include hits from secondary particles. Mokka also stores hits from secondary particles. 
 * Due to the way they are created by SensitiveDetector class in Mokka they are
 * not on the exact detector measurement  surface as the surface of a hit is calculated as (entrypoint+exitpoint)/2.
 * Secondary particles however come to existence within the silicon and will die there too. So their entry point
 * is not really the beginning of the silicon but already within. The same is true for the exit point.
 * The middle between them is therefore usually not at the exact z-postion of the detector measurement plane. 
 * So only if keepHitsFromDeltas is set to true, will hits with a z-position not corresponding to the measurement plane produce TrackerHits.
 * <br>
 *
 */


class SimpleDiscDigiProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleDiscDigiProcessor ; }
  
  
  SimpleDiscDigiProcessor() ;
  
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
  
  
 protected:

  /** Input collection name.
   */
  std::string _inColName ;
  std::string _outColName ;

  int _sub_det_id ;

  int _nRun ;
  int _nEvt ;
 
  float _pointReso;
  bool _keepHitsFromDeltas;
  std::vector< double > _FTDZCoordinate;
  
  bool hasCorrectZPos ( double z );
  
  // gsl random number generator
  gsl_rng * r ;


} ;

#endif



