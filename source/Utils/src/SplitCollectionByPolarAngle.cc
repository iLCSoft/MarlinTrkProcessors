#include "SplitCollectionByPolarAngle.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>

#include "TMath.h"
#include "TVector3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

SplitCollectionByPolarAngle aSplitCollectionByPolarAngle;

SplitCollectionByPolarAngle::SplitCollectionByPolarAngle() : Processor("SplitCollectionByPolarAngle")
{

    // Modify processor description
    _description = "SplitCollectionByPolarAngle splits hits in theta bins";

    // Input collection
    registerProcessorParameter("TrackerHitCollectionName",
                               "Name of the TrackerHit input collection",
                               m_inputHitCollection,
                               std::string("InputCollection"));

    // Output collection
    registerProcessorParameter("SplitHitCollection",
                               "Split hits from tracker",
                               m_outputHitCollection,
                               std::string("SplitCollection"));
}

void SplitCollectionByPolarAngle::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void SplitCollectionByPolarAngle::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void SplitCollectionByPolarAngle::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of tracker hits
    LCCollection *trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputHitCollection, evt);

    std::string encoderString = trackerHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<TrackerHitPlane> myCellIDEncoding(encoderString);

    // Make the output collections
    LCCollectionVec *SplitHitsCollection030 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollection030->setSubset(true);
    SplitHitsCollection030->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollection3050 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollection3050->setSubset(true);
    SplitHitsCollection3050->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollection5070 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollection5070->setSubset(true);
    SplitHitsCollection5070->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollection7090 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollection7090->setSubset(true);
    SplitHitsCollection7090->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollectionN030 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollectionN030->setSubset(true);
    SplitHitsCollectionN030->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollectionN3050 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollectionN3050->setSubset(true);
    SplitHitsCollectionN3050->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollectionN5070 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollectionN5070->setSubset(true);
    SplitHitsCollectionN5070->parameters().setValue("CellIDEncoding", encoderString);
    LCCollectionVec *SplitHitsCollectionN7090 = new LCCollectionVec(trackerHitCollection->getTypeName());
    SplitHitsCollectionN7090->setSubset(true);
    SplitHitsCollectionN7090->parameters().setValue("CellIDEncoding", encoderString);

    int nHits = trackerHitCollection->getNumberOfElements();

    // Loop over tracker hits
    for (size_t itHit = 0; itHit < nHits; itHit++)
    {
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));
        TVector3 pos(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
        double deg_theta = pos.Theta() * 180. / TMath::Pi();

        if (deg_theta >= 150.)
        {
            SplitHitsCollectionN030->addElement(hit);
        }
        else if (deg_theta < 150. && deg_theta >= 130.)
        {
            SplitHitsCollectionN3050->addElement(hit);
        }
        else if (deg_theta < 130. && deg_theta >= 110.)
        {
            SplitHitsCollectionN5070->addElement(hit);
        }
        else if (deg_theta < 110. && deg_theta >= 90.)
        {
            SplitHitsCollectionN7090->addElement(hit);
        }
        else if (deg_theta < 90. && deg_theta >= 70.)
        {
            SplitHitsCollection7090->addElement(hit);
        }
        else if (deg_theta < 70. && deg_theta >= 50.)
        {
               SplitHitsCollection5070->addElement(hit);
        }
        else if (deg_theta < 50. && deg_theta >= 30.)
        {
               SplitHitsCollection3050->addElement(hit);
        }
        else
        {
            SplitHitsCollection030->addElement(hit);
        }
    }

    // Store the filtered hit collections
    evt->addCollection(SplitHitsCollection030, m_outputHitCollection + "030");
    evt->addCollection(SplitHitsCollection3050, m_outputHitCollection + "3050");
    evt->addCollection(SplitHitsCollection5070, m_outputHitCollection + "5070");
    evt->addCollection(SplitHitsCollection7090, m_outputHitCollection + "7090");
    evt->addCollection(SplitHitsCollectionN030, m_outputHitCollection + "N030");
    evt->addCollection(SplitHitsCollectionN3050, m_outputHitCollection + "N3050");
    evt->addCollection(SplitHitsCollectionN5070, m_outputHitCollection + "N5070");
    evt->addCollection(SplitHitsCollectionN7090, m_outputHitCollection + "N7090");

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;
    
    _nEvt++;
}

void SplitCollectionByPolarAngle::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SplitCollectionByPolarAngle::end()
{

    //   std::cout << "SplitCollectionByPolarAngle::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void SplitCollectionByPolarAngle::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
{
    try
    {
        collection = evt->getCollection(collectionName);
    }
    catch (DataNotAvailableException &e)
    {
        streamlog_out(DEBUG5) << "- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
        return;
    }
    return;
}
