#ifndef SimplePlanarDigiProcessor_h
#define SimplePlanarDigiProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>

#include <gsl/gsl_rng.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>

//class TGeoVolume;

using namespace lcio ;
using namespace marlin ;


/** ======= SimplePlanarDigiProcessor ========== <br>
 * Produces VTX TrackerHit collection from SimTrackerHit collections. <br> 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions. 
 * There are two different levels of resolution applied to the layers corresponding to Inner and Outer 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collections of SimTrackerHits in vertex detector <br>
 * <h4>Output</h4>
 * Processor produces collection of digitized TrackerHits in the vertex detector<br>
 * @param VTXCollectionName The name of input collection of VTX SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param VTXHitCollection The name of output collection of digitized VTX TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param PointResolutionRPhi_Inner Point resolution perpendicular to the ladder for the Inner layers (in mm) <br>
 * (default value 0.004) <br>
 * @param PointResolutionZ_Inner Point resolution along the ladder for the Inner layers (in mm) <br>
 * (default value 0.004) <br>
 * @param PointResolutionRPhi_Outer Point resolution perpendicular to the ladder for the Outer layers (in mm) <br>
 * (default value 0.01) <br>
 * @param PointResolutionZ_Outer Point resolution along the ladder for the Outer layers (in mm) <br>
 * (default value 0.01) <br>
 * @param Last_Inner_Layer Layer Number counting from 0 of the last Inner Layer <br>
 * (default value 6) <br>
 * @param Ladder_Number_encoded_in_cellID ladder number has been encoded in the cellID <br>
 * (default value false) <br>
 * <br>
 * 
 */
class SimplePlanarDigiProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimplePlanarDigiProcessor ; }
  
  
  SimplePlanarDigiProcessor() ;
  
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

  // find phi in correct range, taken from gear::VXDParameters
  double correctPhiRange( double Phi ) const ;  
  
  //  void fillMapTVolumeChildren(TGeoVolume* volume) ;
  
 protected:

  std::string _colNameVTX ;

  std::string _outColNameVTX ;

  int _nRun ;
  int _nEvt ;

  float _pointResoRPhi, _pointResoRPhi_Inner, _pointResoRPhi_Outer ;
  float _pointResoZ, _pointResoZ_Inner, _pointResoZ_Outer ;

  int _last_inner_layer;

  bool _ladder_Number_encoded_in_cellID;

  gsl_rng * _rng ;


} ;

#endif



