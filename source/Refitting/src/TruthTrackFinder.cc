/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "TruthTrackFinder.h"
#include "DDRec/API/IDDecoder.h"

//#include "DD4hep/LCDD.h"
//#include "DD4hep/VolumeManager.h"

#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/Global.h"

#include "CLHEP/Vector/TwoVector.h"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace DD4hep ;
using namespace AIDA ;

TruthTrackFinder aTruthTrackFinder;

/*

 This code performs a true pattern recognition by looping over all MC particles and adding all hits
 associated to them onto a track. This is then fitted and output.
 
*/

TruthTrackFinder::TruthTrackFinder() : Processor("TruthTrackFinder") {
	
	// modify processor description
	_description = "TruthTrackFinder builds tracks out of all hits associated to an MC particle" ;

	
  // Input collections - mc particles, tracker hits and the relationships between them
	
	std::vector<std::string> inputTrackerHitCollections;
	inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));

  registerInputCollections( LCIO::TRACKERHITPLANE,
                          "TrackerHitCollectionNames" ,
                          "Name of the TrackerHit input collections"  ,
                          m_inputTrackerHitCollections ,
                          inputTrackerHitCollections ) ;

//	std::vector<std::string> inputSimTrackerHitCollections;
//	inputSimTrackerHitCollections.push_back(std::string("VXDCollection"));
//
//  registerInputCollections( LCIO::SIMTRACKERHIT,
//                          "SimTrackerHitCollectionNames" ,
//                          "Name of the SimTrackerHit input collections"  ,
//                          m_inputSimTrackerHitCollections ,
//                          inputSimTrackerHitCollections ) ;
	
	std::vector<std::string> inputTrackerHitRelationCollections;
	inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));

  registerInputCollections( LCIO::LCRELATION,
                          "SimTrackerHitRelCollectionNames",
                          "Name of TrackerHit SimTrackHit relation collections",
                          m_inputTrackerHitRelationCollections,
                          inputTrackerHitRelationCollections );

  registerInputCollection( LCIO::MCPARTICLE,
                          "MCParticleCollectionName",
                          "Name of the MCParticle input collection",
                          m_inputParticleCollection,
                          std::string("MCParticle"));

  // Output collections - tracks and relations
  registerOutputCollection( LCIO::TRACK,
                           "SiTrackCollectionName",
                           "Silicon track Collection Name",
                           m_outputTrackCollection,
                           std::string("SiTracks"));

  registerOutputCollection( LCIO::LCRELATION,
                           "SiTrackRelationCollectionName",
                           "Silicon track particle relation Collection Name",
                           m_outputTrackRelationCollection,
                           std::string("SiTrackRelations"));
	
}

bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2){
	double radius1 = sqrt((hit1->getPosition()[0])*(hit1->getPosition()[0]) + (hit1->getPosition()[1])*(hit1->getPosition()[1]));
	double radius2 = sqrt((hit2->getPosition()[0])*(hit2->getPosition()[0]) + (hit2->getPosition()[1])*(hit2->getPosition()[1]));
 return (radius1 < radius2);
}

void TruthTrackFinder::init() {
	
	// Print the initial parameters
	printParameters() ;
	
	// Reset counters
	m_runNumber = 0 ;
	m_eventNumber = 0 ;
	m_fitFails = 0;
  
  // Set up the track fit factory
  const gear::GearMgr* fakeGear = 0;
  trackFactory =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , fakeGear , "" ) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useQMS,        true) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::usedEdx,       true) ;
  trackFactory->setOption( MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing,  false) ;
  trackFactory->init() ;

  // Put default values for track fitting
  m_initialTrackError_d0 = 1.e6;
  m_initialTrackError_phi0 = 1.e2;
  m_initialTrackError_omega = 1.e-4;
  m_initialTrackError_z0 = 1.e6;
  m_initialTrackError_tanL = 1.e2;
  m_maxChi2perHit = 1.e2;
  
  // Get the magnetic field
  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
  m_magneticField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)

	// Register this process
	Global::EVENTSEEDER->registerProcessor(this);
		
}


void TruthTrackFinder::processRunHeader( LCRunHeader* run) {
	m_runNumber++ ;
}

void TruthTrackFinder::processEvent( LCEvent* evt ) {
	
  // Get the collection of MC particles
  LCCollection* particleCollection = 0 ;
  getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;

	// Make objects to hold all of the tracker hit, simulated hit and relation collections
	std::vector<LCCollection*> trackerHitCollections;
	std::vector<LCCollection*> trackerHitRelationCollections;
	std::vector<LCRelationNavigator*> relations;
	
	// Loop over each input collection and get the data
	for(unsigned int collection=0; collection<m_inputTrackerHitCollections.size();collection++){
		
		// Get the collection of tracker hits
		LCCollection* trackerHitCollection = 0 ;
		getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt); if(trackerHitCollection == 0) return;
		trackerHitCollections.push_back(trackerHitCollection);

	  // Get the collection of tracker hit relations
	  LCCollection* trackerHitRelationCollection = 0 ;
	  getCollection(trackerHitRelationCollection, m_inputTrackerHitRelationCollections[collection], evt); if(trackerHitRelationCollection == 0) return;
		trackerHitRelationCollections.push_back(trackerHitRelationCollection);
  
	  // Create the relations navigator
	  LCRelationNavigator* relation = new LCRelationNavigator( trackerHitRelationCollection );
		relations.push_back(relation);
  
	}
	
  // Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec( LCIO::TRACK )  ;
	
  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0) ;
  trkFlag.setBit( LCIO::TRBIT_HITS ) ;
  trackCollection->setFlag( trkFlag.getFlag()  ) ;
  
  // Make the output particle to track relation collection
  LCCollectionVec* trackRelationCollection = new LCCollectionVec( LCIO::LCRELATION )  ;

 	/* 
   Now for each MC particle we want the list of hits belonging to it. The most 
   efficient way is to loop over all hits once, and store the pointers in a 
   map, with the key a pointer to the MC particle. We can then loop over each
   MC particle at the end and get all of the hits, before making a track.
  */
  
  // Make the container
  std::map<MCParticle*, std::vector<TrackerHit*> > particleHits;

	// Loop over all input collections
	for(unsigned int collection=0; collection<trackerHitCollections.size();collection++){
  	// Loop over tracker hits
	  int nHits = trackerHitCollections[collection]->getNumberOfElements();
	  for(int itHit=0;itHit<nHits;itHit++){
        
	    // Get the hit
	    TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
    
	    // Get the related simulated hit(s)
	    const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects( hit );

	    // Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
	    SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));

	    // Get the particle belonging to that hit
	    MCParticle* particle = simHit->getMCParticle();

	    // Push back the element into the container
	    particleHits[particle].push_back(hit);

	  }
	}
	
  // Now loop over all particles and get the list of hits
  int nParticles = particleCollection->getNumberOfElements();
  for(int itP=0;itP<nParticles;itP++){
		
    // Get the particle
    MCParticle* mcParticle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itP) ) ;
    
    // Get the vector of hits from the container
    if(particleHits.count(mcParticle) == 0) continue;
    std::vector<TrackerHit*> trackHits = particleHits[mcParticle];

    // Only make tracks with 3 or more hits
    if(trackHits.size() < 3) continue;
		
		// Sort the hits from smaller to larger radius
		std::sort(trackHits.begin(),trackHits.end(),sort_by_radius);

    /*
     Fit - this gets complicated. 
     Need to pass a series of objects, including some initial states,
     covariance matrix etc. Set these up, then call the fit.
    */
		
		// Make the track object and relations object
		LCRelationImpl* relationTrack = new LCRelationImpl;
		TrackImpl* track = new TrackImpl ;

    // First, for some reason there are 2 track objects, one which gets saved and one which is used for fitting. Don't ask...
    MarlinTrk::IMarlinTrack* marlinTrack = trackFactory->createTrack();

    // Save a vector of the hits to be used (why is this not attached to the track directly?? MarlinTrkUtils to be updated?)
    EVENT::TrackerHitVec trackfitHits;
    for(unsigned int itTrackHit=0;itTrackHit<trackHits.size();itTrackHit++){
      trackfitHits.push_back(trackHits[itTrackHit]);
    }
		
    // Make an initial covariance matrix with very broad default values
    EVENT::FloatVec covMatrix (15,0); // Size 15, filled with 0s
    covMatrix[0]  = ( m_initialTrackError_d0    ); //sigma_d0^2
    covMatrix[2]  = ( m_initialTrackError_phi0  ); //sigma_phi0^2
    covMatrix[5]  = ( m_initialTrackError_omega ); //sigma_omega^2
    covMatrix[9]  = ( m_initialTrackError_z0    ); //sigma_z0^2
    covMatrix[14] = ( m_initialTrackError_tanL  ); //sigma_tanl^2

    // Try to fit
    int fitError = 0;
//    try {
      fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackfitHits, track, MarlinTrk::IMarlinTrack::forward, covMatrix, m_magneticField, m_maxChi2perHit);
//    } catch (...) {
//      throw;
//    }
		
		// Check track quality - if fit fails chi2 will be 0
		if(track->getChi2() == 0.){delete track; delete relationTrack; delete marlinTrack; m_fitFails++; continue;}
		
    // Add hit information - this is just a fudge for the moment, since we only use vertex hits. Should do for each subdetector once enabled
    track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
    track->subdetectorHitNumbers()[ 2 * lcio::ILDDetID::VXD - 2 ] = trackHits.size();

    // Push back to the output container
    trackCollection->addElement(track);
		
    // Make the particle to track link
    relationTrack->setFrom(track);
    relationTrack->setTo(mcParticle);
    relationTrack->setWeight(1.0);
    trackRelationCollection->addElement(relationTrack);

		delete marlinTrack;
  }
  
  // Save the output track collection
  evt->addCollection( trackCollection , m_outputTrackCollection ) ;
  // Save the output particle to track relation collection
  evt->addCollection( trackRelationCollection , m_outputTrackRelationCollection ) ;

	// Increment the event number
	m_eventNumber++ ;
	for(unsigned int collection=0; collection<trackerHitCollections.size();collection++) delete relations[collection];
	
}

void TruthTrackFinder::check( LCEvent * evt ) {
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TruthTrackFinder::end(){

	streamlog_out(MESSAGE) << " end()  " << name()
	<< " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
	<< std::endl ;
	
//	std::cout<<"Fit fails: "<<m_fitFails<<std::endl;
	
}

void TruthTrackFinder::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    std::cout<<"- cannot get collections !!"<<std::endl;
		std::cout << "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}


  
  
  
  
  
  
  
