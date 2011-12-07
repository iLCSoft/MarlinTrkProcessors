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

namespace EVENT{
  class SimTrackerHit;
}

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
 * (default value lcio::ILDDetID::FTD) <br>
 * @param keepHitsFromDeltas Whether to include hits from secondary particles. Mokka also stores hits from secondary particles. 
 * Due to the way they are created by SensitiveDetector class in Mokka they are
 * not on the exact detector measurement  surface as the surface of a hit is calculated as (entrypoint+exitpoint)/2.
 * Secondary particles however come to existence within the silicon and will die there too. So their entry point
 * is not really the beginning of the silicon but already within. The same is true for the exit point.
 * The middle between them is therefore usually not at the exact z-postion of the detector measurement plane. 
 * So only if keepHitsFromDeltas is set to true, will hits with a z-position not corresponding to the measurement plane produce TrackerHits.
 * <br>
 * (default value false) <br>
 *
 * @param PetalsPerDisk The number of petals per disk. (Although there are no real petals implemented, hits get a CellID0 corresponding
 * to the petal they would be on.) A value of 1 means there is only one petal spanning the whole disk.
 * The petals start from the x-axis with the beginning of the first petal aligned with the axis.  <br>
 * (default value 1) <br>
 *
 * @param SensorsPerPetal The number of sensors per Petal. A petal will be split up into n Sensors, equally distributed. So if the values is 1,
 * there will be one sensor filling the whole petal. If the value is 2, there will be two sensors, splitting the petal at (Rmin-Rmax)/2.
 * As this digitiser simulates a PixelDisk, there there are no sensors at the back. <br>
 * (default value 1) <br> 
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
  
  void process_hits_loi( LCEvent* evt, LCCollection* STHcol );

  void process_hits_new( LCEvent* evt, LCCollection* STHcol );

  /** Input collection name.
   */
  std::string _inColName ;
  std::string _outColName ;
  std::string _outRelColName ;
  
  int _sub_det_id ;
  
  int _nRun ;
  int _nEvt ;
  
  bool _SimHits_encoded_with_cellID;

  bool _use_FTDLayerLayout_from_GEAR;

  float _pointReso;
  bool _keepHitsFromDeltas;

  std::vector< double > _FTDZCoordinate;
  std::vector< double > _diskInnerRadius;
  std::vector< double > _diskOuterRadius;
  
  bool hasCorrectZPos ( SimTrackerHit* hit );
  
  int getPetalNumber ( int layer , double x , double y );
  int getSensorNumber ( int layer , double x , double y );
  
  int _petalsPerDisk;
  int _sensorsPerPetal;
  
  // gsl random number generator
  gsl_rng * r ;
  
  
} ;

#endif



