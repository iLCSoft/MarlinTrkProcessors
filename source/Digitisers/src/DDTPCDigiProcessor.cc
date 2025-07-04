/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "DDTPCDigiProcessor.h"

#include "FixedPadSizeDiskLayout.h"
#include "TPCModularEndplate.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include <gsl/gsl_randist.h>

#include "Circle.h"
#include "LCCylinder.h"
#include "SimpleHelix.h"
#include "constants.h"
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>

// stl exception handler
#include "constants.h"
#include "voxel.h"
#include <stdexcept>

//
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>

// --- DD4hep ---
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

using namespace lcio;
using namespace marlin;
using namespace constants;

#ifdef MARLIN_USE_AIDA
using namespace AIDA;
#endif

DDTPCDigiProcessor aDDTPCDigiProcessor;

bool compare_phi(Voxel_tpc* a, Voxel_tpc* b) { return (a->getPhiIndex() < b->getPhiIndex()); }

bool compare_z(Voxel_tpc* a, Voxel_tpc* b) { return (a->getZIndex() < b->getZIndex()); }

DDTPCDigiProcessor::~DDTPCDigiProcessor() { delete _tpcEP; }

DDTPCDigiProcessor::DDTPCDigiProcessor() : Processor("DDTPCDigiProcessor") {
  // modify processor description
  _description =
      "Produces TPC TrackerHit collection from SimTrackerHit collection, smeared in RPhi and Z. A search is made for "
      "adjacent hits on a pad row, if they are closer in z and r-phi than the steering parameters _doubleHitResRPhi "
      "(default value 2.0 mm) and _doubleHitResZ (default value 5.0 mm) they are considered to overlap. Clusters of "
      "hits smaller than _maxMerge (default value 3) are merged into a single tracker hit, with the position given as "
      "the average poision of the hits in phi and in z. Clusters which have _maxMerge hits or more are determined to "
      "be identifiable as multiple hits, and are not added to the tracker hit collection. This of course means that "
      "good hits caught up in a cluster of background hits will be lossed.";

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::SIMTRACKERHIT, "TPCPadRowHitCollectionName",
                          "Name of the default pad-row based SimTrackerHit collection", _padRowHitColName,
                          std::string("TPCCollection"));

  registerInputCollection(
      LCIO::SIMTRACKERHIT, "TPCSpacePointCollectionName",
      "Name of the additional space point collection which provides additional guide hits between pad row centers.",
      _spacePointColName, std::string("TPCSpacePointCollection"));

  registerInputCollection(LCIO::SIMTRACKERHIT, "TPCLowPtCollectionName",
                          "Name of the LowPt SimTrackerHit collection Produced by Mokka TPC Driver TPC0X",
                          _lowPtHitscolName, std::string("TPCLowPtCollection"));

  registerOutputCollection(LCIO::TRACKERHIT, "TPCTrackerHitsCol", "Name of the Output TrackerHit collection",
                           _TPCTrackerHitsCol, std::string("TPCTrackerHits"));

  registerOutputCollection(LCIO::LCRELATION, "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection", _outRelColName,
                           std::string("TPCTrackerHitRelations"));

  registerProcessorParameter("UseRawHitsToStoreSimhitPointer",
                             "Store the pointer to the SimTrackerHits in RawHits (deprecated) ",
                             _use_raw_hits_to_store_simhit_pointer, bool(false));

  registerProcessorParameter("PointResolutionPadPhi", "Pad Phi Resolution constant in TPC", _pointResoPadPhi,
                             (float)0.900);

  registerProcessorParameter("RejectCellID0", "whether or not to use hits without proper cell ID (pad row)",
                             _rejectCellID0, (int)1);

  registerProcessorParameter("PointResolutionRPhi", "R-Phi Resolution constant in TPC", _pointResoRPhi0, (float)0.050);

  registerProcessorParameter("DiffusionCoeffRPhi", "R-Phi Diffusion Coefficent in TPC", _diffRPhi, (float)0.025);

  registerProcessorParameter("N_eff", "Number of Effective electrons per pad in TPC", _nEff, (int)22);

  registerProcessorParameter("PointResolutionZ", "TPC Z Resolution Coefficent independent of diffusion", _pointResoZ0,
                             (float)0.4);

  registerProcessorParameter("DiffusionCoeffZ", "Z Diffusion Coefficent in TPC", _diffZ, (float)0.08);

  registerProcessorParameter("HitSortingBinningZ", "Defines spatial slice in Z", _binningZ, (float)5.0);

  registerProcessorParameter("HitSortingBinningRPhi", "Defines spatial slice in RP", _binningRPhi, (float)2.0);

  registerProcessorParameter("DoubleHitResolutionZ", "Defines the minimum distance for two seperable hits in Z",
                             _doubleHitResZ, (float)5.0);

  registerProcessorParameter("DoubleHitResolutionRPhi", "Defines the minimum distance for two seperable hits in RPhi",
                             _doubleHitResRPhi, (float)2.0);

  registerProcessorParameter("MaxClusterSizeForMerge",
                             "Defines the maximum number of adjacent hits which can be merged", _maxMerge, (int)3);

  // fg: these are the numbers for the large ILD TPC
  IntVec tpcEPModNumExample = {14, 18, 23, 28, 32, 37, 42, 46};

  registerProcessorParameter("TPCEndPlateModuleNumbers", "Number of modules in the rings of the TPC endplate",
                             _tpcEndPlateModuleNumbers, tpcEPModNumExample);

  // fg: these are just guestimates - need to get correct phi0 values for the TPC modules
  FloatVec tpcEPModPhi0Example = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07};

  registerProcessorParameter("TPCEndPlateModulePhi0s", "Phi0s of modules in the rings of the TPC endplate",
                             _tpcEndPlateModulePhi0s, tpcEPModPhi0Example);

  registerProcessorParameter("TPCEndPlateModuleGapPhi",
                             "Gap size in mm of the gaps between the endplace modules in Phi", _tpcEndPlateModuleGapPhi,
                             (float)1.);

  registerProcessorParameter("TPCEndPlateModuleGapR", "Gap size in mm of the gaps between the endplace modules in R",
                             _tpcEndPlateModuleGapR, (float)1.);
}

void DDTPCDigiProcessor::init() {
  // From GNU documentation:
  // A replacement for the standard terminate_handler which prints
  // more information about the terminating exception (if any) on stderr. Call ...
  // std::set_terminate (__gnu_cxx::__verbose_terminate_handler);

#ifdef DIGIPLOTS
  /// Hook an AIDA implementation -----------------------------------------------

  // First create a pointer to the "IAnalysisFactory" of a specific AIDA
  // implementation. This factory can then be used to produce all other
  // factories.
  _AF = AIDA_createAnalysisFactory();

  // Create a ITreeFactory. -----------------------------------------------------
  // A ITree can be used to store AIDA objects in memory or on disk.

  _TRF = _AF->createTreeFactory();

  /// Create a ITree object which is bound to a file. ---------------------------
  // You must always create a "ITree" object to create any other factory.
  /*
   * Creates a new tree and associates it with a store.
   * The store is assumed to be read/write.
   * The store will be created if it does not exist.
   * @param storeName The name of the store, if empty (""), the tree is
   *                  created in memory and therefore will not be associated
   *                  with a file.
   * @param storeType Implementation specific string, may control store type
   * @param readOnly If true the store is opened readonly, an exception if it
   *                 does not exist
   * @param createNew If false the file must exist, if true the file will be
   *                  created
   * @param options Other options, currently are not specified
   */
  // ITree * ITreeFactory::create(const std::string & storeName,
  //                              const std::string & storeType = "",
  //                              bool readOnly = false,
  //                              bool createNew = false,
  //                              const std::string & options = "") ;

  _TREE = _TRF->create("DDTPCDigi.root", "root", false, true);

  /// Create an IHistogramFactory which is bound to the tree "*_TREE". -----------

  /*
   * Create an IHistogramFactory.
   * @param tree The ITree which created histograms will be associated to.
   * @return     The IHistogramFactory.
   */
  // IHistogramFactory * IAnalysisFactory::createHistogramFactory(ITree & tree);

  _HF = _AF->createHistogramFactory(*_TREE);

  _TREE->mkdir("Histograms");

  /*
   * Create a IHistogram1D.
   * @param path      The path of the created IHistogram. The path can either
   *                  be a relative or full path.
   *                  ("/folder1/folder2/dataName" and
   *                  "../folder/dataName" are valid paths).
   *                  All the directories in the path must exist. The
   *                  characther `/` cannot be used in names; it is only
   *                  used to delimit directories within paths.
   * @param title     The title of the IHistogram1D.
   * @param nBins     The number of bins of the x axis.
   * @param lowerEdge The lower edge of the x axis.
   * @param upperEdge The upper edge of the x axis.
   * @param options   The options for the IHistogram1D. The default is "".
   *                  "type=efficiency" for an efficiency IHistogram1D.
   * @return          The newly created IHistogram1D.
   */

  _phiDiffHisto = _HF->createHistogram1D("Histograms/phi_diff", "Calculated Phi - Track Phi", 201, -0.05, 0.05);

  _thetaDiffHisto = _HF->createHistogram1D("Histograms/theta_diff", "Calculated Theta - Track Theta", 201, -0.05, 0.05);

  _phiRelHisto = _HF->createHistogram1D("Histograms/padPhi", "Phi Relative to the Pad", 201, 0.0, 6.3);

  _thetaRelHisto = _HF->createHistogram1D("Histograms/padtheta", "Theta Relative to the pad", 201, 0.0, 6.3);

  _rPhiDiffHisto = _HF->createHistogram1D("Histograms/rPhiDiff", "rPhi_rec - rPhi_sim", 201, -1.0, 1.0);

  _zDiffHisto = _HF->createHistogram1D("Histograms/zDiff", "Z_rec - Z_sim", 201, -1.0, 1.0);

  _zPullHisto = _HF->createHistogram1D("Histograms/zPull", "(z_rec - z_sim) / Sigma_z", 201, -10.0, 10.0);

  _phiDistHisto = _HF->createHistogram1D("Histograms/phiDist", "phi_rec - Phi_sim", 201, -1.0, 1.0);

  _rPhiPullHisto =
      _HF->createHistogram1D("Histograms/rPhiPull", "(rPhi_rec - rPhi_sim) / Sigma_rPhi", 201, -10.0, 10.0);

  _zSigmaVsZHisto =
      _HF->createHistogram2D("Histograms/zSigmaVsZ", "z Sigma vs Z ", 3000, 0.0, 3000.0, 201, -0.20, 5.20);

  _zSigmaHisto = _HF->createHistogram1D("Histograms/zSigma", "z Sigma ", 201, -0.20, 5.20);

  _rPhiSigmaHisto = _HF->createHistogram1D("Histograms/rPhiSigma", "rPhi Sigma", 201, -0.20, 0.20);

  _radiusCheckHisto = _HF->createHistogram1D("Histograms/radiusCheck",
                                             "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight)", 201, -0.20, 0.20);

  _ResidualsRPhiHisto = _HF->createHistogram1D("Histograms/ResidualsRPhi", "MC Track Phi - Hit Phi", 50, -0.001, 0.001);

  _NSimTPCHitsHisto = _HF->createHistogram1D("Histograms/SimTPCHits", "Number of SimTPC Hits", 100, 0.0, 1000000.0);

  _NBackgroundSimTPCHitsHisto = _HF->createHistogram1D("Histograms/NBackgroundSimTPCHits",
                                                       "Number of Background SimTPC Hits", 100, 0.0, 1000000.0);

  _NPhysicsSimTPCHitsHisto =
      _HF->createHistogram1D("Histograms/NPhysicsSimTPCHits", "Number of PhysicsSimTPC Hits", 100, 0.0, 100000.0);

  _NPhysicsAbove02GeVSimTPCHitsHisto = _HF->createHistogram1D(
      "Histograms/NPhysicsAbove02GeVTPCHits", "Number of PhysicsSimTPC Hits above 0.2GeV pt", 100, 0.0, 100000.0);

  _NPhysicsAbove1GeVSimTPCHitsHisto = _HF->createHistogram1D(
      "Histograms/NPhysicsAbove1GeVPtTPCHits", "Number of PhysicsSimTPC Hits above 1.0 GeV pt", 100, 0.0, 100000.0);

  _NRecTPCHitsHisto = _HF->createHistogram1D("Histograms/NRecTPCHits", "Number of Rec TPC Hits", 50, 0.0, 100000.0);

  _NLostPhysicsTPCHitsHisto =
      _HF->createHistogram1D("Histograms/NLostPhysicsTPCHits", "Number of PhysicsSimTPC Hits Lost", 100, 0.0, 5000.0);

  _NLostPhysicsAbove02GeVPtTPCHitsHisto =
      _HF->createHistogram1D("Histograms/NLostPhysicsAbove02GeVPtTPCHits",
                             "Number of PhysicsSimTPC Hits Lost above 0.2 GeV pt", 100, 0.0, 5000.0);

  _NLostPhysicsAbove1GeVPtTPCHitsHisto =
      _HF->createHistogram1D("Histograms/NLostPhysicsAbove1GeVPtTPCHits",
                             "Number of PhysicsSimTPC Hits Lost above 1.0 GeV pt", 100, 0.0, 1000.0);

  _NRevomedHitsHisto =
      _HF->createHistogram1D("Histograms/NRevomedHits", "Number of Removed TPC hits", 100, 0.0, 1000000.0);

  _NKeptPhysicsTPCHitsHistoPercent = _HF->createHistogram1D("Histograms/NKeptPhysicsTPCHitsPercent",
                                                            "Number of PhysicsSimTPC Hits Kept", 303, 0.0, 1.01);

  _NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent =
      _HF->createHistogram1D("Histograms/NKeptPhysicsAbove02GeVPtTPCHitsPercent",
                             "Number of PhysicsSimTPC Hits Kept above 0.2 GeV pt", 303, 0.0, 1.01);

  _NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent =
      _HF->createHistogram1D("Histograms/NKeptPhysicsAbove1GeVPtTPCHitsPercent",
                             "Number of PhysicsSimTPC Hits Kept above 1.0 GeV pt", 303, 0.0, 1.01);

#endif

  //--- get the geometry data from dd4hep
  dd4hep::Detector& theDet = dd4hep::Detector::getInstance();
  dd4hep::DetElement tpcDE = theDet.detector("TPC");
  _tpc = tpcDE.extension<dd4hep::rec::FixedPadSizeTPCData>();

  // fill the data for the TPC endplate
  _tpcEP = new TPCModularEndplate(_tpc);

  if (_tpcEndPlateModuleNumbers.size() != _tpcEndPlateModulePhi0s.size()) {
    throw Exception(" DDTPCDigiProcessor: parameters tpcEndPlateModuleNumbers and tpcEndPlateModulePhi0s dont have the "
                    "same number of elements ( module rings ) !! ");
  }

  for (unsigned i = 0, N = _tpcEndPlateModuleNumbers.size(); i < N; ++i) {
    _tpcEP->addModuleRing(_tpcEndPlateModuleNumbers[i], _tpcEndPlateModulePhi0s[i]);
  }
  _tpcEP->initialize();

  // -----

  streamlog_out(DEBUG6) << " initialized TPC geometry from TPCData: " << *_tpc << std::endl;

  double bfieldV[3];
  theDet.field().magneticField({0., 0., 0.}, bfieldV);
  _bField = bfieldV[2] / dd4hep::tesla;
  //----

  printParameters();

  // intialise random number generator
  _random = gsl_rng_alloc(gsl_rng_ranlxs2);
  marlin::Global::EVENTSEEDER->registerProcessor(this);

  _cellid_encoder = 0;
  _nRun = 0;
  _nEvt = 0;
}

void DDTPCDigiProcessor::processRunHeader(LCRunHeader*) { _nRun++; }

void DDTPCDigiProcessor::processEvent(LCEvent* evt) {
  gsl_rng_set(_random, marlin::Global::EVENTSEEDER->getSeed(this));
  streamlog_out(DEBUG) << "seed set to " << marlin::Global::EVENTSEEDER->getSeed(this) << " for event number "
                       << evt->getEventNumber() << std::endl;

  int numberOfVoxelsCreated(0);

  _NSimTPCHits = 0;
  _NBackgroundSimTPCHits = 0;
  _NPhysicsSimTPCHits = 0;
  _NPhysicsAbove02GeVSimTPCHits = 0;
  _NPhysicsAbove1GeVSimTPCHits = 0;
  _NRecTPCHits = 0;

  _NLostPhysicsTPCHits = 0;
  _NLostPhysicsAbove02GeVPtTPCHits = 0;
  _NLostPhysicsAbove1GeVPtTPCHits = 0;
  _NRevomedHits = 0;

  static bool firstEvent = true;
  _tpcHitMap.clear();
  _tpcRowHits.clear();

  streamlog_out(DEBUG8) << "  =========  processing event " << std::setw(9) << evt->getEventNumber() << " run "
                        << std::setw(9) << evt->getRunNumber() << "  ========= " << endl;

  if (firstEvent == true) {
    if (!_use_raw_hits_to_store_simhit_pointer) {
      streamlog_out(DEBUG4) << "The relations to SimTrackerHits are now stored in relation collection "
                            << _outRelColName
                            << "\n SimTrackerHits are no longer stored in RawTrackerHits. Enable this deprecated "
                               "feature by setting UseRawHitsToStoreSimhitPointer to true in steering file."
                            << std::endl;

    } else {
      streamlog_out(DEBUG4)
          << "SimTrackerHits will be stored in RawTrackerHits. This is a deprecated please use the relations stored in "
          << _outRelColName << std::endl;
    }
  }

  firstEvent = false;

  _padWidth = _tpc->padWidth / dd4hep::mm;
  // set size of row_hits to hold (n_rows) vectors
  _tpcRowHits.resize(_tpc->maxRow);

  // created the collection which will be written out
  _trkhitVec = new LCCollectionVec(LCIO::TRACKERHIT);
  // relations from created trackerhits to the SimTrackerHits that caused them
  auto hitSimHitNav = UTIL::LCRelationNavigator(LCIO::TRACKERHIT, LCIO::SIMTRACKERHIT);

  _cellid_encoder = new CellIDEncoder<TrackerHitImpl>(lcio::LCTrackerCellID::encoding_string(), _trkhitVec);

  // first deal with the pad-row based hits from Mokka
  LCCollection* STHcol = 0;
  try {
    STHcol = evt->getCollection(_padRowHitColName);
  } catch (DataNotAvailableException& e) {
  }

  const FixedPadSizeDiskLayout padLayout(_tpc);
  const double TPCPadPlaneRMin = _tpc->rMinReadout / dd4hep::mm;
  const double TPCPadPlaneRMax = _tpc->rMaxReadout / dd4hep::mm;

  float edep0 = 0.0;
  if (STHcol != 0) {
    int n_sim_hits = STHcol->getNumberOfElements();

    LCFlagImpl colFlag(STHcol->getFlag());

    _NSimTPCHits = n_sim_hits;

    streamlog_out(DEBUG4) << "number of Pad-Row based SimHits = " << n_sim_hits << std::endl;

    // make sure that all the pointers are initialise to NULL
    _mcp = NULL;
    _previousMCP = NULL;
    _nextMCP = NULL;
    _nMinus2MCP = NULL;
    _nPlus2MCP = NULL;

    _SimTHit = NULL;
    _previousSimTHit = NULL;
    _nextSimTHit = NULL;
    _nMinus2SimHit = NULL;
    _nPlus2SimHit = NULL;

    // loop over all the pad row based sim hits
    for (int i = 0; i < n_sim_hits; i++) {
      // this will used for nominaml smearing for very low pt rubish, so set it to zero initially
      double ptSqrdMC = 0;

      _SimTHit = dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));

      float edep;
      double padPhi(0.0);
      double padTheta(0.0);

      streamlog_out(DEBUG3) << "processing hit " << i << std::endl;
      streamlog_out(DEBUG3) << " address = " << _SimTHit << " x = " << _SimTHit->getPosition()[0]
                            << " y = " << _SimTHit->getPosition()[1] << " z = " << _SimTHit->getPosition()[2]
                            << std::endl;

      CLHEP::Hep3Vector thisPoint(_SimTHit->getPosition()[0], _SimTHit->getPosition()[1], _SimTHit->getPosition()[2]);
      double padheight = _tpc->padHeight / dd4hep::mm;

      // conversion constant. r = pt / (FCT*_bField)
      const double FCT = 2.99792458E-4;

      _mcp = _SimTHit->getMCParticle();

      // increase the counters for the different classification of simhits
      if (_mcp) {
        // get the pt of the MCParticle, this will used later to uses nominal smearing for low momentum rubish
        const double* momentumMC = _mcp->getMomentum();
        ptSqrdMC = momentumMC[0] * momentumMC[0] + momentumMC[1] * momentumMC[1];

        streamlog_out(DEBUG3) << " mcp address = " << _mcp << " px = " << momentumMC[0] << " py = " << momentumMC[1]
                              << " pz = " << momentumMC[2] << std::endl;

        // SJA:FIXME: the fact that it is a physics hit relies on the fact that for overlay
        // the pointer to the mcp is set to NULL. This distinction may not always be true ...
        ++_NPhysicsSimTPCHits;
        if (ptSqrdMC > (0.2 * 0.2))
          ++_NPhysicsAbove02GeVSimTPCHits;
        if (ptSqrdMC > 1.0)
          ++_NPhysicsAbove1GeVSimTPCHits;

#ifdef DIGIPLOTS
        if (_mcp)
          plotHelixHitResidual(_mcp, &thisPoint);
#endif
      } else {
        ++_NBackgroundSimTPCHits;
      }

      // if the hits contain the momentum of the particle use this to calculate the angles relative to the pad
      if (colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
        const float* mcpMomentum = _SimTHit->getMomentum();

        CLHEP::Hep3Vector mom(mcpMomentum[0], mcpMomentum[1], mcpMomentum[2]);

        // const double pt = mom.perp();
        // const double radius = pt / (FCT*_bField);

        // const double tanLambda = mom.z()/pt;

        padPhi = fabs(thisPoint.deltaPhi(mom));
        padTheta = mom.theta();

      }

      else { // LCIO::THBIT_MOMENTUM not set

        // as the momentum vector is not available from the hits use triplets of
        // hits to fit a circle and calculate theta and phi relative to the pad

        if (!_mcp || (sqrt(ptSqrdMC) / (FCT * _bField)) < (padheight / (0.1 * twopi))) {
          // if the hit has no record of it MCParticle then there is no way to know if this hit has consecutive hits
          // from the same MCParticle so just set nominal values theta=phi=90 here make a cut for particles which will
          // suffer more than a 10 percent change in phi over the distance of the pad R > padheight/(0.1*2PI) in both
          // cases set the angles to 90 degrees
          padTheta = twopi / 4.0;
          padPhi = twopi / 4.0;
        } else {
          // if there is at least one more hit after this one, set the pointer to the MCParticle for the next hit
          if (i < (n_sim_hits - 1)) {
            _nextSimTHit = dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i + 1));
            _nextMCP = _nextSimTHit->getMCParticle();
          } else { // set make sure that the pointers are set back to NULL so that the comparisons later hold
            _nextSimTHit = NULL;
            _nextMCP = NULL;
          }
          // if there is at least two more hits after this one, set the pointer to the MCParticle for the next but one
          // hit
          if (i < (n_sim_hits - 2)) {
            _nPlus2SimHit = dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i + 2));
            _nPlus2MCP = _nPlus2SimHit->getMCParticle();
          } else { // set make sure that the pointers are set back to NULL so that the comparisons later hold
            _nPlus2SimHit = NULL;
            _nPlus2MCP = NULL;
          }

          if (_mcp == _previousMCP && _mcp == _nextMCP) { // middle hit of 3 from the same MCParticle

            CLHEP::Hep3Vector precedingPoint(_previousSimTHit->getPosition()[0], _previousSimTHit->getPosition()[1],
                                             _previousSimTHit->getPosition()[2]);
            CLHEP::Hep3Vector followingPoint(_nextSimTHit->getPosition()[0], _nextSimTHit->getPosition()[1],
                                             _nextSimTHit->getPosition()[2]);

            streamlog_out(DEBUG3) << "address of _previousSimTHit = " << _previousSimTHit
                                  << " x = " << _previousSimTHit->getPosition()[0]
                                  << " y = " << _previousSimTHit->getPosition()[1]
                                  << " z = " << _previousSimTHit->getPosition()[2] << std::endl;

            streamlog_out(DEBUG4) << "address of _nextSimTHit = " << _nextSimTHit
                                  << " x = " << _nextSimTHit->getPosition()[0]
                                  << " y = " << _nextSimTHit->getPosition()[1]
                                  << " z = " << _nextSimTHit->getPosition()[2] << std::endl;

            // get phi and theta using functions defined below
            padPhi = getPadPhi(&thisPoint, &precedingPoint, &thisPoint, &followingPoint);
            padTheta = getPadTheta(&precedingPoint, &thisPoint, &followingPoint);

          } else if (_mcp == _nextMCP && _mcp == _nPlus2MCP) { // first  hit of 3 from the same MCParticle

            CLHEP::Hep3Vector followingPoint(_nextSimTHit->getPosition()[0], _nextSimTHit->getPosition()[1],
                                             _nextSimTHit->getPosition()[2]);
            CLHEP::Hep3Vector nPlus2Point(_nPlus2SimHit->getPosition()[0], _nPlus2SimHit->getPosition()[1],
                                          _nPlus2SimHit->getPosition()[2]);

            // get phi and theta using functions defined below
            padPhi = getPadPhi(&thisPoint, &thisPoint, &followingPoint, &nPlus2Point);
            padTheta = getPadTheta(&thisPoint, &followingPoint, &nPlus2Point);

          } else if (_mcp == _previousMCP && _mcp == _nMinus2MCP) { // last   hit of 3 from the same MCParticle

            CLHEP::Hep3Vector nMinus2Point(_nMinus2SimHit->getPosition()[0], _nMinus2SimHit->getPosition()[1],
                                           _nMinus2SimHit->getPosition()[2]);
            CLHEP::Hep3Vector precedingPoint(_previousSimTHit->getPosition()[0], _previousSimTHit->getPosition()[1],
                                             _previousSimTHit->getPosition()[2]);

            // get phi and theta using functions defined below
            padPhi = getPadPhi(&thisPoint, &nMinus2Point, &precedingPoint, &thisPoint);
            padTheta = getPadTheta(&nMinus2Point, &precedingPoint, &thisPoint);

          } else { // the hit is isolated as either a single hit, or a pair of hits, from a single MCParticle
            padTheta = twopi / 4.0;
            padPhi = twopi / 4.0;
          }
        }

#ifdef DIGIPLOTS
        if (colFlag.bitSet(LCIO::THBIT_MOMENTUM)) {
          const float* mcpMomentum = _SimTHit->getMomentum();

          CLHEP::Hep3Vector mom(mcpMomentum[0], mcpMomentum[1], mcpMomentum[2]);

          double trackPhi = mom.phi();

          if (trackPhi < 0.0)
            trackPhi = trackPhi + twopi;
          if (trackPhi > twopi)
            trackPhi = trackPhi - twopi;
          if (trackPhi > twopi / 2.0)
            trackPhi = trackPhi - twopi / 2.0;

          double localPhi = thisPoint.phi() - padPhi;

          _phiRelHisto->fill(padPhi);
          _phiDiffHisto->fill((fabs(localPhi - trackPhi)) / trackPhi);
          _thetaRelHisto->fill(padTheta);
          _thetaDiffHisto->fill((sin(padTheta) - sin(mom.theta())) / sin(mom.theta()));

          streamlog_out(DEBUG3) << "track Phi = " << trackPhi * (360.0 / twopi) << endl;
          streamlog_out(DEBUG3) << "localPhi = " << localPhi * (360.0 / twopi) << endl;
          streamlog_out(DEBUG3) << "pad Phi = " << padPhi * (360.0 / twopi) << endl;
          streamlog_out(DEBUG3) << "pad Phi from track mom = " << (thisPoint.phi() - trackPhi) * (360.0 / twopi)
                                << endl;
          streamlog_out(DEBUG3) << "padTheta = " << padTheta * (360.0 / twopi) << endl;
          streamlog_out(DEBUG3) << "padTheta from track mom = " << mom.theta() * (360.0 / twopi) << endl;
        }
#endif
      }

      //      int pad = padLayout.getNearestPad(thisPoint.perp(),thisPoint.phi());
      int layerNumber = _SimTHit->getCellID0();

      if (_rejectCellID0 && (layerNumber < 1)) {
        continue;
      }

      edep = _SimTHit->getEDep();

      // Calculate Point Resolutions according to Ron's Formula

      // sigma_{RPhi}^2 = sigma_0^2 + Cd^2/N_{eff} * L_{drift}

      // sigma_0^2 = (50micron)^2 + (900micron*sin(phi))^2
      // Cd^2/N_{eff}} = 25^2/(22/sin(theta)*h/6mm)
      // Cd = 25 ( microns / cm^(1/2) )
      // (this is for B=4T, h is the pad height = pad-row pitch in mm,
      // theta is the polar angle)

      // sigma_{z}^2 = (400microns)^2 + L_{drift}cm * (80micron/sqrt(cm))^2

      double aReso =
          _pointResoRPhi0 * _pointResoRPhi0 + (_pointResoPadPhi * _pointResoPadPhi * sin(padPhi) * sin(padPhi));
      double driftLength = _tpc->driftLength / dd4hep::mm - (fabs(thisPoint.z()));

      if (driftLength < 0) {
        streamlog_out(DEBUG3) << " DDTPCDigiProcessor : Warning! driftLength < 0 " << driftLength
                              << " --> wrong data in dd4hep::rec::FixedPadSizeTPCData ? " << std::endl;
        streamlog_out(DEBUG3) << "Setting driftLength to 0.1" << std::endl;
        streamlog_out(DEBUG3) << "_tpc->driftLength/dd4hep::mm = " << _tpc->driftLength / dd4hep::mm << std::endl;
        driftLength = 0.10;
      }

      padheight = _tpc->padHeight / dd4hep::mm;

      // double bReso = ( (_diffRPhi * _diffRPhi) / _nEff ) * sin(padTheta) * ( 6.0 / (padheight) )  * ( 4.0 / _bField
      // ) ;
      //  formula with new quadratic B-field correction term
      double bReso =
          ((_diffRPhi * _diffRPhi) / _nEff) * sin(padTheta) * (6.0 / (padheight)) * ((4.0 * 4.0) / (_bField * _bField));

      double tpcRPhiRes = sqrt(aReso + bReso * (driftLength / 10.0)); // driftLength in cm

      double tpcZRes =
          sqrt((_pointResoZ0 * _pointResoZ0) + (_diffZ * _diffZ) * (driftLength / 10.0)); // driftLength in cm

      int padIndex = padLayout.getNearestPad(thisPoint.perp(), thisPoint.phi());

      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);

      int NBinsZ = (int)((2.0 * _tpc->driftLength / dd4hep::mm) / _binningZ);
      int iZHit = (int)((float)NBinsZ * (_tpc->driftLength / dd4hep::mm + thisPoint.z()) /
                        (2.0 * _tpc->driftLength / dd4hep::mm));

      if (iZHit < 0)
        iZHit = 0;
      if (iZHit > NBinsZ)
        iZHit = NBinsZ;

      // make sure that the hit lies at the middle of the pad ring
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);

      if ((thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax)) {
        streamlog_out(DEBUG3) << "Hit R not in TPC " << endl;
        streamlog_out(DEBUG3) << "R = " << thisPoint.perp() << endl;
        streamlog_out(DEBUG3) << "the tpc InnerRadius = " << TPCPadPlaneRMin << endl;
        streamlog_out(DEBUG3) << "the tpc OuterRadius = " << TPCPadPlaneRMax << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue;
      }

      if ((fabs(thisPoint.z()) > _tpc->driftLength / dd4hep::mm)) {
        streamlog_out(DEBUG3) << "Hit Z not in TPC " << endl;
        streamlog_out(DEBUG3) << "Z = " << thisPoint.z() << endl;
        streamlog_out(DEBUG3) << "the tpc Max Z = " << _tpc->driftLength / dd4hep::mm << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue;
      }

      // get energy deposit of this row
      edep = _SimTHit->getEDep();

      // create a tpc voxel hit and store it for this row
      Voxel_tpc* atpcVoxel = new Voxel_tpc(iRowHit, iPhiHit, iZHit, thisPoint, edep, tpcRPhiRes, tpcZRes);

      _tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;

      // store the simhit pointer for this tpcvoxel hit in the hit map
      _tpcHitMap[atpcVoxel] = _SimTHit;

      // move the pointers on
      _nMinus2MCP = _previousMCP;
      _previousMCP = _mcp;
      _nMinus2SimHit = _previousSimTHit;
      _previousSimTHit = _SimTHit;
    }
  }

  // now process the LowPt collection
  LCCollection* STHcolLowPt = 0;
  try {
    STHcolLowPt = evt->getCollection(_lowPtHitscolName);
  } catch (DataNotAvailableException& e) {
  }

  if (STHcolLowPt != NULL) {
    int n_sim_hitsLowPt = STHcolLowPt->getNumberOfElements();

    _NBackgroundSimTPCHits += n_sim_hitsLowPt;
    _NSimTPCHits += n_sim_hitsLowPt;

    streamlog_out(DEBUG4) << "number of LowPt hits:" << n_sim_hitsLowPt << std::endl;

    // loop over the LowPt hit collection
    for (int i = 0; i < n_sim_hitsLowPt; i++) {
      _SimTHit = dynamic_cast<SimTrackerHit*>(STHcolLowPt->getElementAt(i));

      CLHEP::Hep3Vector thisPoint(_SimTHit->getPosition()[0], _SimTHit->getPosition()[1], _SimTHit->getPosition()[2]);

      int NBinsZ = (int)((2.0 * _tpc->driftLength / dd4hep::mm) / _binningZ);

      if ((thisPoint.perp() < TPCPadPlaneRMin) || (thisPoint.perp() > TPCPadPlaneRMax)) {
        streamlog_out(DEBUG3) << "Hit R not in TPC " << endl;
        streamlog_out(DEBUG3) << "R = " << thisPoint.perp() << endl;
        streamlog_out(DEBUG3) << "the tpc InnerRadius = " << TPCPadPlaneRMin << endl;
        streamlog_out(DEBUG3) << "the tpc OuterRadius = " << TPCPadPlaneRMax << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue;
      }

      if ((fabs(thisPoint.z()) > _tpc->driftLength / dd4hep::mm)) {
        streamlog_out(DEBUG3) << "Hit Z not in TPC " << endl;
        streamlog_out(DEBUG3) << "Z = " << thisPoint.z() << endl;
        streamlog_out(DEBUG3) << "the tpc Max Z = " << _tpc->driftLength / dd4hep::mm << endl;
        streamlog_out(DEBUG3) << "Hit Dropped " << endl;
        continue;
      }

      int padIndex = padLayout.getNearestPad(thisPoint.perp(), thisPoint.phi());

      int iRowHit = padLayout.getRowNumber(padIndex);
      int iPhiHit = padLayout.getPadNumber(padIndex);
      int iZHit = (int)((float)NBinsZ * (_tpc->driftLength / dd4hep::mm + thisPoint.z()) /
                        (2.0 * _tpc->driftLength / dd4hep::mm));

      // shift the hit in r-phi to the nearest pad-row centre
      thisPoint.setPerp(padLayout.getPadCenter(padIndex)[0]);

      // set the resolutions to the pads to digital like values
      double tpcRPhiRes = _padWidth;
      double tpcZRes = _binningZ;

      // get energy deposit of this hit
      edep0 = _SimTHit->getEDep();

      // create a tpc voxel hit for this simhit and store it for this tpc pad row
      Voxel_tpc* atpcVoxel = new Voxel_tpc(iRowHit, iPhiHit, iZHit, thisPoint, edep0, tpcRPhiRes, tpcZRes);

      _tpcRowHits.at(iRowHit).push_back(atpcVoxel);
      ++numberOfVoxelsCreated;

      // store the simhit pointer for this voxel hit in a map
      _tpcHitMap[atpcVoxel] = _SimTHit;
    }
  }

  int number_of_adjacent_hits(0);

  streamlog_out(DEBUG4) << "finished looping over simhits, number of voxels = " << numberOfVoxelsCreated << endl;

  int numberOfhitsTreated(0);

  vector<Voxel_tpc*> row_hits;

  // loop over the tpc rows containing hits and check for merged hits
  for (unsigned int i = 0; i < _tpcRowHits.size(); ++i) {
    row_hits = _tpcRowHits.at(i);
    std::sort(row_hits.begin(), row_hits.end(), compare_phi);

    // double loop over the hits in this row
    for (unsigned int j = 0; j < row_hits.size(); ++j) {
      ++numberOfhitsTreated;

      for (unsigned int k = j + 1; k < row_hits.size(); ++k) {
        if (row_hits[k]->getPhiIndex() >
            (row_hits[j]->getPhiIndex()) + 2) { // SJA:FIXME: here we need an OR to catch the wrap around
          break;                                // only compare hits in adjacent phi bins
        }

        // look to see if the two hit occupy the same pad in phi or if not whether they are within the r-phi double hit
        // resolution
        else if (row_hits[k]->getPhiIndex() == row_hits[j]->getPhiIndex() ||
                 ((fabs(row_hits[k]->getHep3Vector().deltaPhi(row_hits[j]->getHep3Vector()))) * row_hits[j]->getR()) <
                     _doubleHitResRPhi) {
          // if neighboring in phi then compare z
          map<Voxel_tpc*, SimTrackerHit*>::iterator it;

          SimTrackerHit* Hit1 = NULL;
          SimTrackerHit* Hit2 = NULL;

          // search of the simhit pointers in the tpchit map
          it = _tpcHitMap.find(row_hits[j]);
          if (it != _tpcHitMap.end()) {
            Hit1 = it->second; // hit found
          }

          it = _tpcHitMap.find(row_hits[k]);
          if (it != _tpcHitMap.end()) {
            Hit2 = it->second; // hit found
          }

          double pathlengthZ1(0.0);
          double pathlengthZ2(0.0);

          if (Hit1 && Hit2) { // if both sim hits were found

            // check if the track momentum has been stored for the hits
            bool momentum_set = true;

            if (STHcol != NULL) {
              LCFlagImpl colFlag(STHcol->getFlag());
              momentum_set = momentum_set && colFlag.bitSet(LCIO::THBIT_MOMENTUM);
            }

            if (STHcolLowPt != NULL) {
              LCFlagImpl colFlag(STHcolLowPt->getFlag());
              momentum_set = momentum_set && colFlag.bitSet(LCIO::THBIT_MOMENTUM);
            }

            if (momentum_set) {
              const float* Momentum1 = Hit1->getMomentum();
              const float* Momentum2 = Hit2->getMomentum();

              CLHEP::Hep3Vector mom1(Momentum1[0], Momentum1[1], Momentum1[2]);
              CLHEP::Hep3Vector mom2(Momentum2[0], Momentum2[1], Momentum2[2]);

              pathlengthZ1 = fabs(Hit1->getPathLength() * mom1.cosTheta());
              pathlengthZ2 = fabs(Hit2->getPathLength() * mom2.cosTheta());
            } else {
              pathlengthZ1 = _doubleHitResZ; // assume the worst i.e. that the track is moving in z
              pathlengthZ2 = _doubleHitResZ; // assume the worst i.e. that the track is moving in z
            }

            double dZ = fabs(row_hits[j]->getZ() - row_hits[k]->getZ());

            double spacial_coverage = 0.5 * (pathlengthZ1 + pathlengthZ2) + _binningZ;

            if ((dZ - spacial_coverage) < _doubleHitResZ) {
              row_hits[j]->setAdjacent(row_hits[k]);
              row_hits[k]->setAdjacent(row_hits[j]);
              ++number_of_adjacent_hits;
            }
          } else {
            streamlog_out(DEBUG3) << "Hit1=" << Hit1 << "Hit2=" << Hit2 << endl;
          }
        }
      }
    }

    // now all hits have been checked for adjacent hits, go throught and write out the hits or merge

    for (unsigned int j = 0; j < row_hits.size(); ++j) {
      Voxel_tpc* seed_hit = row_hits[j];

      if (seed_hit->IsMerged() || seed_hit->IsClusterHit()) {
        continue;
      }

      if (seed_hit->getNumberOfAdjacent() == 0) { // no adjacent hits so smear and write to hit collection
        writeVoxelToHit(seed_hit, hitSimHitNav);
      }

      else if (seed_hit->getNumberOfAdjacent() <
               (_maxMerge)) { // potential 3-hit cluster, can use simple average merge.

        vector<Voxel_tpc*>* hitsToMerge = new vector<Voxel_tpc*>;

        int clusterSize = seed_hit->clusterFind(hitsToMerge);

        if (clusterSize <= _maxMerge) { // merge cluster
          seed_hit->setIsMerged();
          writeMergedVoxelsToHit(hitsToMerge, hitSimHitNav);
        }
        delete hitsToMerge;
      }
    }
  }

  int numberOfHits(0);
  // count up the number of hits merged or lost
  for (unsigned int i = 0; i < _tpcRowHits.size(); ++i) {
    row_hits = _tpcRowHits.at(i);
    for (unsigned int j = 0; j < row_hits.size(); ++j) {
      numberOfHits++;
      Voxel_tpc* seed_hit = row_hits[j];
      if (seed_hit->IsMerged() || seed_hit->IsClusterHit() || seed_hit->getNumberOfAdjacent() > _maxMerge) {
        ++_NRevomedHits;
        _mcp = (_tpcHitMap[seed_hit])->getMCParticle();
        if (_mcp != NULL) {
          ++_NLostPhysicsTPCHits;
          const double* mom = _mcp->getMomentum();
          double ptSQRD = mom[0] * mom[0] + mom[1] * mom[1];
          if (ptSQRD > (0.2 * 0.2))
            ++_NLostPhysicsAbove02GeVPtTPCHits;
          if (ptSQRD > 1.0)
            ++_NLostPhysicsAbove1GeVPtTPCHits;
        }
      }
    }
  }

  streamlog_out(DEBUG4) << "the number of adjacent hits is " << number_of_adjacent_hits << "  _doubleHitResZ "
                        << _doubleHitResZ << endl;
  streamlog_out(DEBUG4) << "number of rec_hits = " << _NRecTPCHits << endl;
  streamlog_out(DEBUG4) << "finished row hits " << numberOfHits << " " << numberOfhitsTreated << endl;

  // set the parameters to decode the type information in the collection
  // for the time being this has to be done manually
  // in the future we should provide a more convenient mechanism to
  // decode this sort of meta information

  StringVec typeNames;
  IntVec typeValues;
  typeNames.push_back(LCIO::TPCHIT);
  typeValues.push_back(1);
  _trkhitVec->parameters().setValues("TrackerHitTypeNames", typeNames);
  _trkhitVec->parameters().setValues("TrackerHitTypeValues", typeValues);

  // add the collection to the event
  evt->addCollection(_trkhitVec, _TPCTrackerHitsCol);
  auto relCol = hitSimHitNav.createLCCollection();
  evt->addCollection(relCol, _outRelColName);

  // delete voxels
  for (unsigned int i = 0; i < _tpcRowHits.size(); ++i) {
    vector<Voxel_tpc*>* current_row = &_tpcRowHits.at(i);
    for (unsigned int j = 0; j < current_row->size(); ++j) {
      delete current_row->at(j);
    }
  }

#ifdef DIGIPLOTS
  _NSimTPCHitsHisto->fill(_NSimTPCHits);
  _NBackgroundSimTPCHitsHisto->fill(_NBackgroundSimTPCHits);
  _NPhysicsSimTPCHitsHisto->fill(_NPhysicsSimTPCHits);
  _NPhysicsAbove02GeVSimTPCHitsHisto->fill(_NPhysicsAbove02GeVSimTPCHits);
  _NPhysicsAbove1GeVSimTPCHitsHisto->fill(_NPhysicsAbove1GeVSimTPCHits);
  _NRecTPCHitsHisto->fill(_NRecTPCHits);

  _NLostPhysicsTPCHitsHisto->fill(_NLostPhysicsTPCHits);
  _NLostPhysicsAbove02GeVPtTPCHitsHisto->fill(_NLostPhysicsAbove02GeVPtTPCHits);
  _NLostPhysicsAbove1GeVPtTPCHitsHisto->fill(_NLostPhysicsAbove1GeVPtTPCHits);
  _NRevomedHitsHisto->fill(_NRevomedHits);

  _NKeptPhysicsTPCHitsHistoPercent->fill((float)(_NPhysicsSimTPCHits - _NLostPhysicsTPCHits) /
                                         (float)_NPhysicsSimTPCHits);
  _NKeptPhysicsAbove02GeVPtTPCHitsHistoPercent->fill(
      (float)(_NPhysicsAbove02GeVSimTPCHits - _NLostPhysicsAbove02GeVPtTPCHits) / (float)_NPhysicsAbove02GeVSimTPCHits);
  _NKeptPhysicsAbove1GeVPtTPCHitsHistoPercent->fill(
      (float)(_NPhysicsAbove1GeVSimTPCHits - _NLostPhysicsAbove1GeVPtTPCHits) / (float)_NPhysicsAbove1GeVSimTPCHits);
#endif

  streamlog_out(DEBUG4) << "_NSimTPCHits = " << _NSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NBackgroundSimTPCHits = " << _NBackgroundSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NPhysicsSimTPCHits = " << _NPhysicsSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NPhysicsAbove02GeVSimTPCHits = " << _NPhysicsAbove02GeVSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NPhysicsAbove1GeVSimTPCHits = " << _NPhysicsAbove1GeVSimTPCHits << endl;
  streamlog_out(DEBUG4) << "_NRecTPCHits = " << _NRecTPCHits << endl;
  streamlog_out(DEBUG4) << "_NLostPhysicsTPCHits = " << _NLostPhysicsTPCHits << endl;
  streamlog_out(DEBUG4) << "_NLostPhysicsAbove02GeVPtTPCHits = " << _NLostPhysicsAbove02GeVPtTPCHits << endl;
  streamlog_out(DEBUG4) << "_NLostPhysicsAbove1GeVPtTPCHits = " << _NLostPhysicsAbove1GeVPtTPCHits << endl;
  streamlog_out(DEBUG4) << "_NRevomedHits = " << _NRevomedHits << endl;

  _nEvt++;
  // Clear the maps and the end of the event.
  _tpcHitMap.clear();
  _tpcRowHits.clear();

  delete _cellid_encoder;
}

void DDTPCDigiProcessor::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void DDTPCDigiProcessor::end() {
#ifdef DIGIPLOTS
  _TREE->commit();
  _TREE->cd("/Histograms");
  _TREE->ls("..");

  _TREE->close();
  streamlog_out(MESSAGE) << "DIGICHECKPLOTS Finished" << endl;
#endif

  gsl_rng_free(_random);
  streamlog_out(MESSAGE) << "DDTPCDigiProcessor::end()  " << name() << " processed " << _nEvt << " events in " << _nRun
                         << " runs " << endl;
  //
}

void DDTPCDigiProcessor::writeVoxelToHit(Voxel_tpc* aVoxel, UTIL::LCRelationNavigator& hitSimHitNav) {
  Voxel_tpc* seed_hit = aVoxel;

  //  if( seed_hit->getRowIndex() > 5 ) return ;

  // store hit variables
  TrackerHitImpl* trkHit = new TrackerHitImpl;
  // now the hit pos has to be smeared

  double tpcRPhiRes = seed_hit->getRPhiRes();
  double tpcZRes = seed_hit->getZRes();

  CLHEP::Hep3Vector point(seed_hit->getX(), seed_hit->getY(), seed_hit->getZ());

  //-------- fg: remove the hit if it lies within a module boundary in phi --------

  dd4hep::rec::Vector3D hitPos(seed_hit->getX(), seed_hit->getY(), seed_hit->getZ());

  double hitGapDist = _tpcEP->computeDistanceRPhi(hitPos);

  if (hitGapDist < _tpcEndPlateModuleGapPhi / 2.) {
    streamlog_out(DEBUG2) << " removing hit in endplate module gap : " << hitPos << " - distance : " << hitGapDist
                          << std::endl;

    return;
  }

  // fg: here we could add additional effects, e.g larger smearing of hits in the
  //     vicinity of the gaps ( need extra parameters )
  //     or larger smearing in the pad row next to a module ring boundary in r
  //     ...
  //------------------------------------------------------------------------------

  double unsmearedPhi = point.phi();

  double randrp = gsl_ran_gaussian(_random, tpcRPhiRes);
  double randz = gsl_ran_gaussian(_random, tpcZRes);

  point.setPhi(point.phi() + randrp / point.perp());
  point.setZ(point.z() + randz);

  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if (fabs(point.z()) > _tpc->driftLength / dd4hep::mm)
    point.setZ((fabs(point.z()) / point.z()) * _tpc->driftLength / dd4hep::mm);

#if 0 // fg: it turns out that is more correct to allow reconstructed hits to be outside of the drift volume due to the
      // large uncertainties
  // make sure the hit is not smeared beyond the cathode:
  double dzCathode = _tpc->zMinReadout/dd4hep::mm ;
  if( point.z() > 0. && hitPos.z() < 0. ){
    point.setZ( -dzCathode ) ;
  }
  if( point.z() < 0. && hitPos.z() > 0. ) {
    point.setZ( dzCathode ) ;
  }
#endif

  double pos[3] = {point.x(), point.y(), point.z()};
  trkHit->setPosition(pos);
  trkHit->setEDep(seed_hit->getEDep());
  //  trkHit->setType( 500 );

  //  int side = lcio::ILDDetID::barrel ;
  //
  //  if( pos[2] < 0.0 ) side = 1 ;

  (*_cellid_encoder)[lcio::LCTrackerCellID::subdet()] = lcio::ILDDetID::TPC;
  (*_cellid_encoder)[lcio::LCTrackerCellID::layer()] = seed_hit->getRowIndex();
  (*_cellid_encoder)[lcio::LCTrackerCellID::module()] = 0;
  (*_cellid_encoder)[lcio::LCTrackerCellID::side()] = (pos[2] < 0 ? -1 : 1);

  _cellid_encoder->setCellID(trkHit);

  // check values for inf and nan
  if (std::isnan(unsmearedPhi) || std::isinf(unsmearedPhi) || std::isnan(tpcRPhiRes) || std::isinf(tpcRPhiRes)) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n"
             << "unsmearedPhi = " << unsmearedPhi << " tpcRPhiRes = " << tpcRPhiRes << "\n";
    throw Exception(errorMsg.str());
  }

  // For no error in R
  float covMat[TRKHITNCOVMATRIX] = {float(sin(unsmearedPhi) * sin(unsmearedPhi) * tpcRPhiRes * tpcRPhiRes),
                                    float(-cos(unsmearedPhi) * sin(unsmearedPhi) * tpcRPhiRes * tpcRPhiRes),
                                    float(cos(unsmearedPhi) * cos(unsmearedPhi) * tpcRPhiRes * tpcRPhiRes),
                                    float(0.),
                                    float(0.),
                                    float(tpcZRes * tpcZRes)};

  trkHit->setCovMatrix(covMat);

  if (_tpcHitMap[seed_hit] == NULL) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n"
             << "SimTracker Pointer is NULL throwing exception\n"
             << "\n";
    throw Exception(errorMsg.str());
  }

  if (pos[0] * pos[0] + pos[1] * pos[1] > 0.0) {
    //    push back the SimTHit for this TrackerHit

    if (_use_raw_hits_to_store_simhit_pointer) {
      trkHit->rawHits().push_back(_tpcHitMap[seed_hit]);
    }

    hitSimHitNav.addRelation(trkHit, _tpcHitMap[seed_hit], 1.0);

    _trkhitVec->addElement(trkHit);
    _NRecTPCHits++;
  }

#ifdef DIGIPLOTS
  SimTrackerHit* theSimHit = _tpcHitMap[seed_hit];
  double rSimSqrd = theSimHit->getPosition()[0] * theSimHit->getPosition()[0] +
                    theSimHit->getPosition()[1] * theSimHit->getPosition()[1];

  double phiSim = atan2(theSimHit->getPosition()[1], theSimHit->getPosition()[0]);

  double rPhiDiff = (point.phi() - phiSim) * sqrt(rSimSqrd);
  double rPhiPull =
      ((point.phi() - phiSim) * sqrt(rSimSqrd)) / (sqrt((covMat[2]) / (cos(point.phi()) * cos(point.phi()))));

  double zDiff = point.getZ() - theSimHit->getPosition()[2];
  double zPull = zDiff / sqrt(covMat[5]);

  _rPhiDiffHisto->fill(rPhiDiff);
  _rPhiPullHisto->fill(rPhiPull);
  _phiDistHisto->fill(point.phi() - phiSim);
  _zDiffHisto->fill(zDiff);
  _zPullHisto->fill(zPull);

  _zSigmaVsZHisto->fill(seed_hit->getZ(), sqrt(covMat[5]));
  _rPhiSigmaHisto->fill(sqrt((covMat[2]) / (cos(point.phi()) * cos(point.phi()))));
  _zSigmaHisto->fill(sqrt(covMat[5]));
#endif
}

void DDTPCDigiProcessor::writeMergedVoxelsToHit(vector<Voxel_tpc*>* hitsToMerge,
                                                UTIL::LCRelationNavigator& hitSimHitNav) {
  TrackerHitImpl* trkHit = new TrackerHitImpl;

  double sumZ = 0;
  double sumPhi = 0;
  double sumEDep = 0;
  //  double R = 0;
  double lastR = 0;

  unsigned number_of_hits_to_merge = hitsToMerge->size();

  for (unsigned int ihitCluster = 0; ihitCluster < number_of_hits_to_merge; ++ihitCluster) {
    sumZ += hitsToMerge->at(ihitCluster)->getZ();
    sumPhi += hitsToMerge->at(ihitCluster)->getPhi();
    sumEDep += hitsToMerge->at(ihitCluster)->getEDep();
    hitsToMerge->at(ihitCluster)->setIsMerged();
    lastR = hitsToMerge->at(ihitCluster)->getR();

    if (_use_raw_hits_to_store_simhit_pointer) {
      trkHit->rawHits().push_back(_tpcHitMap[hitsToMerge->at(ihitCluster)]);
    }

    hitSimHitNav.addRelation(trkHit, _tpcHitMap[hitsToMerge->at(ihitCluster)], float(1.0 / number_of_hits_to_merge));
  }

  double avgZ = sumZ / (hitsToMerge->size());
  double avgPhi = sumPhi / (hitsToMerge->size());

  // set deposit energy as mean of merged hits
  sumEDep = sumEDep / (double)number_of_hits_to_merge;

  CLHEP::Hep3Vector* mergedPoint = new CLHEP::Hep3Vector(1.0, 1.0, 1.0);
  mergedPoint->setPerp(lastR);
  mergedPoint->setPhi(avgPhi);
  mergedPoint->setZ(avgZ);

  // store hit variables

  // first the hit pos has to be smeared------------------------------------------------

  // FIXME: which errors should we use for smearing the merged hits ?
  //        this might be a bit large ....
  double tpcRPhiRes = _padWidth;
  double tpcZRes = _binningZ;

  CLHEP::Hep3Vector point(mergedPoint->getX(), mergedPoint->getY(), mergedPoint->getZ());

  //  double unsmearedPhi = point.phi();

  double randrp = gsl_ran_gaussian(_random, tpcRPhiRes);
  double randz = gsl_ran_gaussian(_random, tpcZRes);

  point.setPhi(point.phi() + randrp / point.perp());
  point.setZ(point.z() + randz);

  // make sure the hit is not smeared beyond the TPC Max DriftLength
  if (fabs(point.z()) > _tpc->driftLength / dd4hep::mm)
    point.setZ((fabs(point.z()) / point.z()) * _tpc->driftLength / dd4hep::mm);

#if 0 // fg: it turns out that is more correct to allow reconstructed hits to be outside of the drift volume due to the
      // large uncertainties
  // make sure the hit is not smeared onto the other side of the cathode:
  double dzCathode = _tpc->zMinReadout/dd4hep::mm ;
  if( point.z() > 0. && mergedPoint->getZ() < 0. ){
    point.setZ( -dzCathode ) ;
  }
  if( point.z() < 0. && mergedPoint->getZ() > 0. ){
    point.setZ( dzCathode ) ;
  }
#endif

  double pos[3] = {point.x(), point.y(), point.z()};

  //---------------------------------------------------------------------------------
  trkHit->setPosition(pos);
  trkHit->setEDep(sumEDep);
  //  trkHit->setType( 500 );

  FixedPadSizeDiskLayout padLayout(_tpc);
  int padIndex = padLayout.getNearestPad(mergedPoint->perp(), mergedPoint->phi());
  int row = padLayout.getRowNumber(padIndex);

  (*_cellid_encoder)[lcio::LCTrackerCellID::subdet()] = lcio::ILDDetID::TPC;
  (*_cellid_encoder)[lcio::LCTrackerCellID::layer()] = row;
  (*_cellid_encoder)[lcio::LCTrackerCellID::module()] = 0;
  (*_cellid_encoder)[lcio::LCTrackerCellID::side()] = (pos[2] < 0 ? -1 : 1);

  _cellid_encoder->setCellID(trkHit);

  double phi = mergedPoint->getPhi();

  // check values for inf and nan
  if (std::isnan(phi) || std::isinf(phi) || std::isnan(tpcRPhiRes) || std::isinf(tpcRPhiRes)) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n"
             << "phi = " << phi << " tpcRPhiRes = " << tpcRPhiRes << "\n";
    throw Exception(errorMsg.str());
  }

  // For no error in R
  float covMat[TRKHITNCOVMATRIX] = {float(sin(phi) * sin(phi) * tpcRPhiRes * tpcRPhiRes),
                                    float(-cos(phi) * sin(phi) * tpcRPhiRes * tpcRPhiRes),
                                    float(cos(phi) * cos(phi) * tpcRPhiRes * tpcRPhiRes),
                                    float(0.),
                                    float(0.),
                                    float(tpcZRes * tpcZRes)};

  trkHit->setCovMatrix(covMat);

  //  if(pos[0]*pos[0]+pos[1]*pos[1]>0.0){
  _trkhitVec->addElement(trkHit);
  ++_nRechits;
  //  } else {
  //    delete trkHit;
  //  }

  delete mergedPoint;
}

#ifdef DIGIPLOTS
void DDTPCDigiProcessor::plotHelixHitResidual(MCParticle* mcp, CLHEP::Hep3Vector* thisPoint) {
  const double FCT = 2.99792458E-4;
  double charge = mcp->getCharge();
  const double* mom = mcp->getMomentum();
  double pt = sqrt(mom[0] * mom[0] + mom[1] * mom[1]);
  double radius = pt / (FCT * _bField);
  double tanLambda = mom[2] / pt;
  double omega = charge / radius;

  if (pt > 1.0) {
    // FIXME SJA: this is only valid for tracks from the IP and should be done correctly for non prompt tracks
    double Z0 = 0.;
    double D0 = 0.;

    LCVector3D refPoint(0., 0., 0);

    SimpleHelix* helix = new SimpleHelix(D0, atan2(mom[1], mom[0]), omega, Z0, tanLambda, refPoint);

    // an almost "infinite" cylinder in z
    LCVector3D startCylinder(0., 0., -1000000.0);
    LCVector3D endCylinder(0., 0., 1000000.0);
    bool endplane = true;

    LCCylinder cylinder(startCylinder, endCylinder, thisPoint->perp(), endplane);

    bool pointExists = false;

    double pathlength = helix->getIntersectionWithCylinder(cylinder, pointExists);

    LCErrorMatrix* errors = new LCErrorMatrix();

    if (pointExists) {
      LCVector3D intersection = helix->getPosition(pathlength, errors);

      double intersectionPhi = atan2(intersection[1], intersection[0]);
      double residualRPhi = ((intersectionPhi - thisPoint->phi()));
      _ResidualsRPhiHisto->fill(residualRPhi);
    }

    delete errors;
    delete helix;

    int row = padLayout.getRowNumber(padLayout.getNearestPad(thisPoint->perp(), thisPoint->phi()));
    int pad = padLayout.getNearestPad(thisPoint->perp(), thisPoint->phi());

    double rHit_diff = thisPoint->perp() - TPCPadPlaneRMin - ((row + 0.5) * _tpc->padHeight / dd4hep::mm);

    _radiusCheckHisto->fill(rHit_diff);

    //      streamlog_out(MESSAGE) << "$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$" << endl;
    //      streamlog_out(MESSAGE) << "thisPoint->perp() = " << thisPoint->perp() << endl;
    //      streamlog_out(MESSAGE) << "TPC Sensitive rMin = " << padLayout.getPlaneExtent()[0] << endl;
    //      streamlog_out(MESSAGE) << "Row number + 0.5 = " <<  row + 0.5 << endl;
    //      streamlog_out(MESSAGE) << "Pad Height = " <<  padLayout.getPadHeight(pad) << endl;
    //      streamlog_out(MESSAGE) << "Row Height = " <<   padLayout.getRowHeight(row) << endl;
    //      streamlog_out(MESSAGE) << "R_hit - TPC Rmin - ((RowIndex + 0.5 )* padheight) = " << rHit_diff << endl;
  }
  return;
}
#endif

double DDTPCDigiProcessor::getPadPhi(CLHEP::Hep3Vector* thisPoint, CLHEP::Hep3Vector* firstPoint,
                                     CLHEP::Hep3Vector* middlePoint, CLHEP::Hep3Vector* lastPoint) {
  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(), firstPoint->y());
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(), middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(), lastPoint->y());

  // check that the points are not the same, at least at the level of a tenth of a micron
  if ((fabs(firstPointRPhi.x() - middlePointRPhi.x()) < 1.e-05 &&
       fabs(firstPointRPhi.y() - middlePointRPhi.y()) < 1.e-05) ||
      (fabs(middlePointRPhi.x() - lastPointRPhi.x()) < 1.e-05 &&
       fabs(middlePointRPhi.y() - lastPointRPhi.y()) < 1.e-05) ||
      (fabs(firstPointRPhi.x() - lastPointRPhi.x()) < 1.e-05 &&
       fabs(firstPointRPhi.y() - lastPointRPhi.y()) < 1.e-05)) {
    streamlog_out(WARNING)
        << " DDTPCDigiProcessor::getPadPhi "
        << "2 of the 3 SimTracker hits passed to Circle Fit are the same hit taking pad phi as PI/2\n"
        << " firstPoint->x() " << firstPoint->x() << " firstPoint->y() " << firstPoint->y() << " firstPoint->z() "
        << firstPoint->z() << " middlePoint->x() " << middlePoint->x() << " middlePoint->y() " << middlePoint->y()
        << " middlePoint->z() " << middlePoint->z() << " lastPoint->x() " << lastPoint->x() << " lastPoint->y() "
        << lastPoint->y() << " lastPoint.z() " << lastPoint->z() << std::endl;

    return twopi / 4.0;
  }

  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);

  double localPhi =
      atan2((thisPoint->y() - theCircle.GetCenter()->y()), (thisPoint->x() - theCircle.GetCenter()->x())) +
      (twopi / 4.0);

  if (localPhi > twopi)
    localPhi = localPhi - twopi;
  if (localPhi < 0.0)
    localPhi = localPhi + twopi;
  if (localPhi > twopi / 2.0)
    localPhi = localPhi - twopi / 2.0;

  double pointPhi = thisPoint->phi();

  if (pointPhi > twopi)
    pointPhi = pointPhi - twopi;
  if (pointPhi < 0.0)
    pointPhi = pointPhi + twopi;
  if (pointPhi > twopi / 2.0)
    pointPhi = pointPhi - twopi / 2.0;

  double padPhi = fabs(pointPhi - localPhi);

  // check that the value returned is reasonable
  if (std::isnan(padPhi) || std::isinf(padPhi)) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n"
             << "padPhi = " << padPhi << "\n";
    throw Exception(errorMsg.str());
  }

  return padPhi;
}

double DDTPCDigiProcessor::getPadTheta(CLHEP::Hep3Vector* firstPoint, CLHEP::Hep3Vector* middlePoint,
                                       CLHEP::Hep3Vector* lastPoint) {
  // Calculate thetaPad for current hit
  CLHEP::Hep2Vector firstPointRPhi(firstPoint->x(), firstPoint->y());
  CLHEP::Hep2Vector middlePointRPhi(middlePoint->x(), middlePoint->y());
  CLHEP::Hep2Vector lastPointRPhi(lastPoint->x(), lastPoint->y());

  // check that the points are not the same, at least at the level of a tenth of a micron
  if ((fabs(firstPointRPhi.x() - middlePointRPhi.x()) < 1.e-05 &&
       fabs(firstPointRPhi.y() - middlePointRPhi.y()) < 1.e-05) ||
      (fabs(middlePointRPhi.x() - lastPointRPhi.x()) < 1.e-05 &&
       fabs(middlePointRPhi.y() - lastPointRPhi.y()) < 1.e-05) ||
      (fabs(firstPointRPhi.x() - lastPointRPhi.x()) < 1.e-05 &&
       fabs(firstPointRPhi.y() - lastPointRPhi.y()) < 1.e-05)) {
    streamlog_out(WARNING)
        << " DDTPCDigiProcessor::getPadTheta "
        << "2 of the 3 SimTracker hits passed to Circle Fit are the same hit taking pad phi as PI/2\n"
        << " firstPoint->x() " << firstPoint->x() << " firstPoint->y() " << firstPoint->y() << " firstPoint->z() "
        << firstPoint->z() << " middlePoint->x() " << middlePoint->x() << " middlePoint->y() " << middlePoint->y()
        << " middlePoint->z() " << middlePoint->z() << " lastPoint->x() " << lastPoint->x() << " lastPoint->y() "
        << lastPoint->y() << " lastPoint.z() " << lastPoint->z() << std::endl;

    return twopi / 4.0;
  }

  Circle theCircle(&firstPointRPhi, &middlePointRPhi, &lastPointRPhi);

  double deltaPhi = firstPoint->deltaPhi(*lastPoint);

  double pathlength = fabs(deltaPhi) * theCircle.GetRadius();

  double padTheta = atan(pathlength / fabs(lastPoint->z() - firstPoint->z()));

  double pathlength1 =
      2.0 *
      asin((sqrt((middlePointRPhi.x() - firstPointRPhi.x()) * (middlePointRPhi.x() - firstPointRPhi.x()) +
                 (middlePointRPhi.y() - firstPointRPhi.y()) * (middlePointRPhi.y() - firstPointRPhi.y())) /
            2.0) /
           theCircle.GetRadius()) *
      theCircle.GetRadius();

  double pathlength2 =
      ((sqrt((lastPointRPhi.x() - middlePointRPhi.x()) * (lastPointRPhi.x() - middlePointRPhi.x()) +
             (lastPointRPhi.y() - middlePointRPhi.y()) * (lastPointRPhi.y() - middlePointRPhi.y())) /
        2.0) /
           theCircle.GetRadius() >=
       1.0)
          ? (2.0 * asin(1.0) * theCircle.GetRadius())
          : (2.0 *
             asin((sqrt((lastPointRPhi.x() - middlePointRPhi.x()) * (lastPointRPhi.x() - middlePointRPhi.x()) +
                        (lastPointRPhi.y() - middlePointRPhi.y()) * (lastPointRPhi.y() - middlePointRPhi.y())) /
                   2.0) /
                  theCircle.GetRadius()) *
             theCircle.GetRadius());

  padTheta = atan((fabs(pathlength1 + pathlength2)) / (fabs(lastPoint->z() - firstPoint->z())));

  // check that the value returned is reasonable
  if (std::isnan(padTheta) || std::isinf(padTheta)) {
    std::stringstream errorMsg;
    errorMsg << "\nProcessor: DDTPCDigiProcessor \n"
             << "padTheta = " << padTheta << "\n";
    throw Exception(errorMsg.str());
  }

  return padTheta;
}
