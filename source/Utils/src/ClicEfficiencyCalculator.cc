/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ClicEfficiencyCalculator.h"
#include "DDRec/API/IDDecoder.h"

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

ClicEfficiencyCalculator aClicEfficiencyCalculator;

/*

 Efficiency calculator for tracking. WARNING: Only works with 1 input collection at present
 
*/

ClicEfficiencyCalculator::ClicEfficiencyCalculator() : Processor("ClicEfficiencyCalculator") {
	
	// modify processor description
	_description = "ClicEfficiencyCalculator calculates the tracking efficiency and makes some performance plots" ;

	// Input collections
	registerInputCollection( LCIO::TRACK,
													"TrackCollectionName",
													"Track collection name",
													m_inputTrackCollection,
													std::string("SiTracks"));
	
	registerInputCollection( LCIO::LCRELATION,
													"TrackRelationCollectionName",
													"Track relation collection name",
													m_inputTrackRelationCollection,
													std::string("SiTrackRelations"));
	
	registerInputCollection( LCIO::MCPARTICLE,
													"MCParticleCollectionName",
													"Name of the MCParticle input collection",
													m_inputParticleCollection,
													std::string("MCParticle"));
	
	// All tracker hit collections
	std::vector<std::string> inputTrackerHitCollections;
	inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
	
	registerInputCollections( LCIO::TRACKERHITPLANE,
													 "TrackerHitCollectionNames" ,
													 "Name of the TrackerHit input collections"  ,
													 m_inputTrackerHitCollections ,
													 inputTrackerHitCollections ) ;
	
	// All tracker hit relation collections
	std::vector<std::string> inputTrackerHitRelationCollections;
	inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
	
	registerInputCollections( LCIO::LCRELATION,
													 "SimTrackerHitRelCollectionNames",
													 "Name of TrackerHit relation collections",
													 m_inputTrackerHitRelationCollections,
													 inputTrackerHitRelationCollections );
	
}


void ClicEfficiencyCalculator::init() {
	
	// Print the initial parameters
	printParameters() ;
	
	// Reset counters
	m_runNumber = 0 ;
	m_eventNumber = 0 ;
	
	// Get the magnetic field
	DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
	const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
	double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
	lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
	m_magneticField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
	
	// Define the purity cut - should be made configurable
	m_purity = 0.5;
	
	// Register this process
	Global::EVENTSEEDER->registerProcessor(this);
	
}


void ClicEfficiencyCalculator::processRunHeader( LCRunHeader* run) {
	++m_runNumber ;
}

void ClicEfficiencyCalculator::processEvent( LCEvent* evt ) {
	
	// First pick up all of the collections that will be used - tracks, MCparticles, hits from relevent subdetectors - and their relations
	
	// Get the collection of tracks
	LCCollection* trackCollection = 0 ;
	getCollection(trackCollection, m_inputTrackCollection, evt); if(trackCollection == 0) return;
	
	// Get the collection of MC particles
	LCCollection* particleCollection = 0 ;
	getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;
	
	// Get the collection of track relations to MC particles
	LCCollection* trackRelationCollection = 0 ;
	getCollection(trackRelationCollection, m_inputTrackRelationCollection, evt); if(trackRelationCollection == 0) return;
	
	// Create the relations navigator
	LCRelationNavigator* trackRelation = new LCRelationNavigator( trackRelationCollection );
	
	// Make objects to hold all of the tracker hit and relation collections
	std::vector<LCCollection*> trackerHitCollections;
	std::vector<LCCollection*> trackerHitRelationCollections;
	std::vector<LCRelationNavigator*> relations;
	
	// Loop over each input collection and get the data
	for(unsigned int collection=0; collection<m_inputTrackerHitCollections.size();collection++){
		
		// Get the collection of tracker hits
		LCCollection* trackerHitCollection = 0 ;
		getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt); if(trackerHitCollection == 0) return;
		trackerHitCollections.push_back(trackerHitCollection);
		
//		CellIDDecoder<TrackerHitPlane> cellid_decoder( STHcol) ;
//		int layer  = cellid_decoder( simTHit )["layer"];

		// Get the collection of tracker hit relations
		LCCollection* trackerHitRelationCollection = 0 ;
		getCollection(trackerHitRelationCollection, m_inputTrackerHitRelationCollections[collection], evt); if(trackerHitRelationCollection == 0) return;
		trackerHitRelationCollections.push_back(trackerHitRelationCollection);
  
		// Create the relations navigator
		LCRelationNavigator* relation = new LCRelationNavigator( trackerHitRelationCollection );
		relations.push_back(relation);
  
	}

	std::map<MCParticle*,int> particleTracks;
	
	/*
		Look at all tracks that were reconstructed and calculate their purity. This can be used to point to all
	 	reconstructed MC particles. Then loop over MC particles and for each that was not reconstructed check the
	  particle properties. Use these to define efficiency etc. For this we also need the list of hits belonging 
	  to each particle.
  */
	
	// Loop over all tracks
	int nTracks = trackCollection->getNumberOfElements();
	for(int itTrack=0;itTrack<nTracks;itTrack++){
		
		// Get the track
		Track* track = dynamic_cast<Track*>( trackCollection->getElementAt(itTrack) ) ;
		
		// Get the hits
		const TrackerHitVec& hitVector = track->getTrackerHits();
		
		// Some storage to keep track of which particles have hits on this track
		std::vector<MCParticle*> trackParticles;
		std::map<MCParticle*,int> trackParticleHits;
		
		// Loop over hits and check the MC particle that they correspond to, and which subdetector they are on
		int nHits = hitVector.size();
		for(int itHit=0;itHit<nHits;itHit++){
			// Get the tracker hit
			TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitVector.at(itHit));
			// Get the simulated hit
			const LCObjectVec& simHitVector = relations[0]->getRelatedToObjects( hit ); // FIX! this only lets you use 1 input collection
			// Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
			SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
			// Get the particle belonging to that hit
			MCParticle* particle = simHit->getMCParticle();
			// Check if this particle already has hits on this track
			if( std::find(trackParticles.begin(), trackParticles.end(), particle) == trackParticles.end() ) trackParticles.push_back(particle);
			// Increment the number of counts
			trackParticleHits[particle]++;
		}
		
		// Check the purity of the track, treat ghosts
		int nTrackParticles = trackParticles.size(); int maxHits=0;
		MCParticle* associatedParticle = 0;
		for(int itParticle=0; itParticle<nTrackParticles; itParticle++){
			int nParticleHits = trackParticleHits[trackParticles[itParticle]];
			if(nParticleHits>maxHits){ maxHits=nParticleHits; associatedParticle = trackParticles[itParticle]; }
		}
		double purity = (double)maxHits/(double)nHits;
		if(purity < m_purity) continue; // do something with ghosts here

		// Now have a track which is associated to a particle
		particleTracks[associatedParticle]++;

	}
	
	
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
	
	/*
	   Now there is a map containing the number of tracks associated to each MC particle, and a 
	   map of all hits on each particle. Loop over all MC particles, decide if they are 
	   reconstructable or not, and check if they were reconstructed more than once (clones).
	 */
	
	// Loop over particles
	int nParticles = particleCollection->getNumberOfElements();
	for(int itParticle=0;itParticle<nParticles;itParticle++){
		
		// Get the particle
		MCParticle* particle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itParticle) ) ;
		
		// Check if it was reconstructed
		if(particleTracks.count(particle)){
			// Assumption: if particle was reconstructed then it is reconstructable!
			m_particles["all"]++; m_particles["reconstructable"]++;
			m_reconstructedParticles["all"]++;
			// Check if clones were produced (1 particle, more than 1 track)
			if(particleTracks[particle] > 1) m_reconstructedParticles["clones"]+=(particleTracks[particle]-1);
			continue;
		}
		
		// Now decide on criteria for different particle types/classifications
		m_particles["all"]++; // all particles
		
		// No hits in the input collections
		if(particleHits.count(particle) == 0) continue; // No hits in the input collections

		// Exclude particles which are not "real"
		if(particle->vertexIsNotEndpointOfParent()) continue;
		
		// Only make tracks with 3 or more hits
		std::vector<TrackerHit*> trackHits = particleHits[particle];
		if(trackHits.size() < 3) continue;

		m_particles["reconstructable"]++; // reconstructable particles
	
		
	}
	
	// Increment the event number
	m_eventNumber++ ;
	
}

void ClicEfficiencyCalculator::check( LCEvent * evt ) {
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ClicEfficiencyCalculator::end(){

	streamlog_out(MESSAGE) << " end()  " << name()
	<< " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
	<< std::endl ;
	
	// Calculate efficiency results
	std::cout<<"Reconstructable particle efficiency: "<<100.*m_reconstructedParticles["all"]/m_particles["reconstructable"]<<" % ("<<m_reconstructedParticles["all"]<<"/"<<m_particles["reconstructable"]<<")"<<std::endl;
	
}

void ClicEfficiencyCalculator::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
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


  
  
  
  
  
  
  
