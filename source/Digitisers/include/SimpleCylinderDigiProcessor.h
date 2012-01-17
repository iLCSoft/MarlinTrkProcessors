#ifndef SimpleCylinderDigiProcessor_h
#define SimpleCylinderDigiProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>

#include <gsl/gsl_rng.h>


using namespace lcio ;
using namespace marlin ;


/** ======= SimpleCylinderDigiProcessor ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The geometry should be supplied via parameters in a Gear File
 * Containing the following named parameters:
 *
 *     SensitiveLayerInnerRadius
 *     SensitiveLayerThickness
 *     LayerHalfLength
 *
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in R-Phi and Z according to the specified point resolutions. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name SITCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name SITTrackerHits) <br>
 * @param PointResolutionRPhi_Inner Point resolution perpendicular to the ladder (in mm) <br>
 * (default value 0.004) <br>
 * @param PointResolutionZ_Inner Point resolution along the ladder (in mm) <br>
 * (default value 0.004) <br>
 * @param HitsEncodedWithCellID Mokka has encoded the hits in the cellID0 according to UTIL/ILDConf.h <br>
 * (default value true) <br>
 * @param Sub_Detector_ID ID of Sub-Detector using UTIL/ILDConf.h from lcio only used if HitsEncodedWithCellID false <br>
 * (default value lcio::ILDDetID::SIT) <br>
 * @param GearParametersName Name of the Gear Parameters to be used <br>
 * (default value SIT_Simple) <br>
 * <br>
 * 
 */
class SimpleCylinderDigiProcessor : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new SimpleCylinderDigiProcessor ; }
  
  
  SimpleCylinderDigiProcessor() ;
  
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
  
  std::string _inColName ;
  
  std::string _outColName ;
  std::string _outRelColName ;
 
  std::string _gearParametersName;
  
  int _sub_det_id ;
  
  int _nRun ;
  int _nEvt ;
  
  float _pointResoRPhi ;
  float _pointResoZ ;
  
  bool _hits_encoded_with_cellID;
  
  gsl_rng* _rng ;
  
  
} ;

#endif



