#include "FilterConeHits.h"
#include <iostream>
#include <cmath>
#include <set>

#include <EVENT/MCParticle.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

#include "marlin/AIDAProcessor.h"


using namespace lcio ;
using namespace marlin ;


FilterConeHits aFilterConeHits ;


FilterConeHits::FilterConeHits() : Processor("FilterConeHits") {

  // --- Processor description:

  _description = "FilterConeHits selects tracker hits in a cone opened around a MC particle direction";

    
  // --- Processor parameters:
  
  registerProcessorParameter("MCParticleCollection",
			     "Name of the MCParticle collection",
			     m_inputMCParticlesCollName,
			     std::string("MCParticle") );

  registerProcessorParameter("TrackerHitInputCollections",
			     "Name of the tracker hit input collections",
			     m_inputTrackerHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerSimHitInputCollections",
			     "Name of the tracker simhit input collections",
			     m_inputTrackerSimHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerHitInputRelations",
			     "Name of the tracker hit relation collections",
			     m_inputTrackerHitRelNames,
			     {} );

  registerProcessorParameter("TrackerHitOutputCollections",
			     "Name of the tracker hit output collections",
			     m_outputTrackerHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerSimHitOutputCollections",
			     "Name of the tracker simhit output collections",
			     m_outputTrackerSimHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerHitOutputRelations",
			     "Name of the tracker hit relation collections",
			     m_outputTrackerHitRelNames,
			     {} );

  registerProcessorParameter( "DeltaRCut" ,
			      "Maximum angular distance between the hits and the particle direction" ,
			      m_deltaRCut,
			      double(1.) );

  registerProcessorParameter( "FillHistograms",
			      "Flag to fill the diagnostic histograms",
			      m_fillHistos,
			      false );

    
}



void FilterConeHits::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  // --- Print the processor parameters:
  printParameters() ;
  
  // --- Initialize the run and event counters:
  _nRun = 0 ;
  _nEvt = 0 ;

  // --- Initialize the AIDAProcessor and book the diagnostic histograms: 

  AIDAProcessor::histogramFactory(this);

  m_deltaR  = new TH1F("m_deltaR", "#DeltaR(hit,MCparticle);#Delta R",
		       1000, 0., 1.);

}


void FilterConeHits::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;

} 



void FilterConeHits::processEvent( LCEvent * evt ) { 

  // --- Check whether the number of input and output collections match

  if ( m_inputTrackerHitsCollNames.size() != m_inputTrackerSimHitsCollNames.size() ||
       m_inputTrackerHitsCollNames.size() != m_inputTrackerHitRelNames.size()  ){

    std::stringstream err_msg;
    err_msg << "Mismatch between the reco and sim hits input collections"
	    << std::endl ;

    throw EVENT::Exception( err_msg.str() ) ;

  }

  if ( m_outputTrackerHitsCollNames.size() != m_outputTrackerSimHitsCollNames.size() ||
       m_outputTrackerHitsCollNames.size() != m_outputTrackerHitRelNames.size()  ){

    std::stringstream err_msg;
    err_msg << "Mismatch between the reco and sim hits output collections"
	    << std::endl ;

    throw EVENT::Exception( err_msg.str() ) ;

  }
  

  // --- Get the MC particles collection:

  LCCollection* m_inputMCParticles = NULL;
  try {
    m_inputMCParticles = evt->getCollection( m_inputMCParticlesCollName );
  }
  catch( lcio::DataNotAvailableException& e ) {
    streamlog_out(WARNING) << m_inputMCParticlesCollName << " collection not available" << std::endl;
    return;
  }


  // --- Get the input hit collections and create the output collections:

  const unsigned int nTrackerHitCol = m_inputTrackerHitsCollNames.size();
  std::vector<LCCollection*> inputHitColls(nTrackerHitCol);
  std::vector<LCCollection*> inputSimHitColls(nTrackerHitCol);
  std::vector<LCCollection*> inputHitRels(nTrackerHitCol);

  std::vector<LCCollectionVec*> outputTrackerHitColls(nTrackerHitCol);
  std::vector<LCCollectionVec*> outputTrackerSimHitColls(nTrackerHitCol);
  std::vector<LCCollectionVec*> outputTrackerHitRels(nTrackerHitCol);

  for (unsigned int icol=0; icol<nTrackerHitCol ; ++icol) {

    // get the reco hits
    try {
      inputHitColls[icol] = evt->getCollection(m_inputTrackerHitsCollNames[icol]);
    }
    catch( lcio::DataNotAvailableException& e ) {
      streamlog_out(WARNING) << m_inputTrackerHitsCollNames[icol]
			     << " collection not available" << std::endl;
      continue;
    }

    // get the sim hits
    try {
      inputSimHitColls[icol] = evt->getCollection(m_inputTrackerSimHitsCollNames[icol]);
    }
    catch( lcio::DataNotAvailableException& e ) {
      streamlog_out(WARNING) << m_inputTrackerSimHitsCollNames[icol]
			     << " collection not available" << std::endl;
      continue;
    }

    // get the reco-sim relations
    try {
      inputHitRels[icol] = evt->getCollection(m_inputTrackerHitRelNames[icol]);
    }
    catch( lcio::DataNotAvailableException& e ) {
      streamlog_out(WARNING) << m_inputTrackerHitRelNames[icol]
			     << " collection not available" << std::endl;
      continue;
    }

    // reco hit output collection
    std::string encoderString = inputHitColls[icol]->getParameters().getStringVal( "CellIDEncoding" );
    outputTrackerHitColls[icol] = new LCCollectionVec( inputHitColls[icol]->getTypeName() );
    outputTrackerHitColls[icol]->parameters().setValue( "CellIDEncoding", encoderString );
    LCFlagImpl lcFlag(inputHitColls[icol]->getFlag());
    outputTrackerHitColls[icol]->setFlag(lcFlag.getFlag());
    
    // sim hit output collections
    outputTrackerSimHitColls[icol] = new LCCollectionVec( inputSimHitColls[icol]->getTypeName() );
    outputTrackerSimHitColls[icol]->parameters().setValue( "CellIDEncoding", encoderString );
    LCFlagImpl lcFlag_sim(inputSimHitColls[icol]->getFlag());
    outputTrackerSimHitColls[icol]->setFlag(lcFlag_sim.getFlag());

    // reco-sim relation output collections
    outputTrackerHitRels[icol] = new LCCollectionVec( inputHitRels[icol]->getTypeName() );
    LCFlagImpl lcFlag_rel(inputHitRels[icol]->getFlag()) ;
    outputTrackerHitRels[icol]->setFlag( lcFlag_rel.getFlag() ) ;
    
  }


  // --- Loop over the MC particles:
  
  std::vector<std::set<int> > hits_to_save(nTrackerHitCol);
  
  for (int ipart=0; ipart<m_inputMCParticles->getNumberOfElements(); ++ipart){

    MCParticle* part = dynamic_cast<MCParticle*>( m_inputMCParticles->getElementAt(ipart) );
    
    if ( part->getGeneratorStatus() != 1 ) continue;

    double part_pt  = sqrt(part->getMomentum()[0]*part->getMomentum()[0]+
			   part->getMomentum()[1]*part->getMomentum()[1]);
    double part_phi = atan2(part->getMomentum()[1],part->getMomentum()[0]);
    double part_theta = atan2(part_pt,part->getMomentum()[2]);


    // --- Loop over the tracker hits and selects hits inside the cone:

    for (unsigned int icol=0; icol<inputHitColls.size(); ++icol){

      LCCollection* hit_col  =  inputHitColls[icol];
      if( !hit_col ) continue ;
      
      for (int ihit=0; ihit<hit_col->getNumberOfElements(); ++ihit){

	TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hit_col->getElementAt(ihit));

    	double hit_rho = sqrt( (hit->getPosition()[0]-part->getVertex()[0])*
			       (hit->getPosition()[0]-part->getVertex()[0])+
			       (hit->getPosition()[1]-part->getVertex()[1])*
			       (hit->getPosition()[1]-part->getVertex()[1]) );
    	double hit_phi = atan2( hit->getPosition()[1]-part->getVertex()[1],
				hit->getPosition()[0]-part->getVertex()[0] );
    	double hit_theta = atan2( hit_rho,hit->getPosition()[2]-part->getVertex()[2]);
    	
	double deltaR = sqrt( (hit_phi-part_phi)*(hit_phi-part_phi) +
			      (hit_theta-part_theta)*(hit_theta-part_theta) );

	if ( m_fillHistos )
	  m_deltaR->Fill(deltaR);

	
	if ( deltaR < m_deltaRCut )
	  hits_to_save[icol].insert(ihit);
	
      } // ihit loop
       
    } // icol loop

  } // ipart loop


  // --- Add the filtered hits to the output collections:
  
  for (unsigned int icol=0; icol<inputHitColls.size(); ++icol){

    for ( auto& ihit: hits_to_save[icol] ){

      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(inputHitColls[icol]->getElementAt(ihit));
      TrackerHitPlaneImpl* hit_new = new TrackerHitPlaneImpl();

      hit_new->setCellID0(hit->getCellID0());
      hit_new->setCellID1(hit->getCellID1());
      hit_new->setType(hit->getType());
      hit_new->setPosition(hit->getPosition());
      hit_new->setU(hit->getU());         
      hit_new->setV(hit->getV());     
      hit_new->setdU(hit->getdU());
      hit_new->setdV(hit->getdV());
      hit_new->setEDep(hit->getEDep());      
      hit_new->setEDepError(hit->getEDepError()); 
      hit_new->setTime(hit->getTime());
      hit_new->setQuality(hit->getQuality());   

      outputTrackerHitColls[icol]->addElement( hit_new );


      LCRelation* rel = dynamic_cast<LCRelation*>(inputHitRels[icol]->getElementAt(ihit));

      
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(rel->getTo());
      SimTrackerHitImpl* simhit_new = new SimTrackerHitImpl();

      simhit_new->setCellID0(simhit->getCellID0());
      simhit_new->setCellID1(simhit->getCellID1());
      simhit_new->setPosition(simhit->getPosition());
      simhit_new->setEDep(simhit->getEDep());
      simhit_new->setTime(simhit->getTime());
      simhit_new->setMCParticle(simhit->getMCParticle());
      simhit_new->setMomentum(simhit->getMomentum());
      simhit_new->setPathLength(simhit->getPathLength());
      simhit_new->setQuality(simhit->getQuality());
      simhit_new->setOverlay(simhit->isOverlay());
      simhit_new->setProducedBySecondary(simhit->isProducedBySecondary());

      outputTrackerSimHitColls[icol]->addElement( simhit_new );

      
      LCRelationImpl* rel_new = new LCRelationImpl();
      
      rel_new->setFrom(hit_new);
      rel_new->setTo(simhit_new);
      rel_new->setWeight(rel->getWeight());

      outputTrackerHitRels[icol]->addElement( rel_new );

    } // ihit loop

    streamlog_out( MESSAGE ) << " " << hits_to_save[icol].size() << " hits added to the collections: "
			     << m_outputTrackerHitsCollNames[icol] << ", "
			     << m_outputTrackerSimHitsCollNames[icol] << ", "
			     << m_outputTrackerHitRelNames[icol] << std::endl;

    evt->addCollection( outputTrackerHitColls[icol], m_outputTrackerHitsCollNames[icol] ) ;
    evt->addCollection( outputTrackerSimHitColls[icol], m_outputTrackerSimHitsCollNames[icol] ) ;
    evt->addCollection( outputTrackerHitRels[icol], m_outputTrackerHitRelNames[icol] ) ;

    streamlog_out( DEBUG5 ) << " output collection " << m_outputTrackerHitsCollNames[icol] << " of type "
			    << outputTrackerHitColls[icol]->getTypeName() << " added to the event \n"
			    << " output collection " << m_outputTrackerSimHitsCollNames[icol] << " of type "
			    << outputTrackerSimHitColls[icol]->getTypeName() << " added to the event \n"
			    << " output collection " << m_outputTrackerHitRelNames[icol] << " of type "
			    << outputTrackerHitRels[icol]->getTypeName() << " added to the event  "
			    << std::endl ;

  } // icol loop

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;
  
  _nEvt ++ ;

}



void FilterConeHits::check( LCEvent * evt ) { 
}



void FilterConeHits::end(){ 

  std::cout << "FilterConeHits::end()  " << name() 
   	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;

}
