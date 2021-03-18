#ifndef SplitCollectionByPolarAngle_h
#define SplitCollectionByPolarAngle_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <map>
#include <vector>

#include <EVENT/LCCollection.h>

using namespace lcio;
using namespace marlin;

/**  Hit filtering processor for marlin.
 * 
 * @param TrackerHitCollectionName Name of the input hit collection
 * @param SplitHitCollection Base name of the output hit collections
 * 
 * @author F. Meloni, DESY; R. Simoniello, CERN
 * @version $Id: SplitCollectionByPolarAngle.h,v 0.1 2020-09-27 11:24:21 fmeloni Exp $ 
 */

class SplitCollectionByPolarAngle : public Processor
{

public:
  virtual Processor *newProcessor() { return new SplitCollectionByPolarAngle; }

  SplitCollectionByPolarAngle();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  // Call to get collections
  void getCollection(LCCollection *&, std::string, LCEvent *);

protected:
  // Collection names for (in/out)put
  std::string m_inputHitCollection = "";
  std::string m_outputHitCollection = "";

  int _nRun{};
  int _nEvt{};
};

#endif
