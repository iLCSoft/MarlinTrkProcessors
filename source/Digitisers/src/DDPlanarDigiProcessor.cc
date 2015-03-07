/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DDPlanarDigiProcessor.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>


#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "AIDA/AIDA.h"

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

//using namespace DD4hep ;

DDPlanarDigiProcessor aDDPlanarDigiProcessor ;

DDPlanarDigiProcessor::DDPlanarDigiProcessor() : Processor("DDPlanarDigiProcessor") {
  
  // modify processor description
  _description = "DDPlanarDigiProcessor creates TrackerHits from SimTrackerHits, smearing them according to the input parameters."
    "The geoemtry of the surface is taken from the DDRec::Surface asscociated to the hit via the cellID" ;
  
  
  // register steering parameters: name, description, class-variable, default value
  
  FloatVec resUEx ;
  resUEx.push_back( 0.0040 ) ;
  
  registerProcessorParameter( "ResolutionU" ,
                              "resolution in direction of u - either one per layer or one for all layers "  ,
                              _resU ,
                              resUEx) ;
  
  FloatVec resVEx ;
  resVEx.push_back( 0.0040 ) ;

  registerProcessorParameter( "ResolutionV" , 
                              "resolution in direction of v - either one per layer or one for all layers " ,
                             _resV ,
                              resVEx );

  registerProcessorParameter( "IsStrip",
                              "whether hits are 1D strip hits",
                              _isStrip,
                              bool(false) );
  
  
  registerProcessorParameter( "Sub_Detector_ID" , 
                             "ID of Sub-Detector from UTIL/ILDConf.h from lcio. Either for VXD, SIT or SET" ,
                             _sub_det_id ,
                             int(lcio::ILDDetID::VXD));
  
  
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                          "SimTrackHitCollectionName" , 
                          "Name of the Input SimTrackerHit collection"  ,
                          _inColName ,
                          std::string("VXDCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHITPLANE,
                           "TrackerHitCollectionName" , 
                           "Name of the TrackerHit output collection"  ,
                           _outColName ,
                           std::string("VTXTrackerHits") ) ;
  
  registerOutputCollection(LCIO::LCRELATION,
                           "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection",
                           _outRelColName,
                           std::string("VTXTrackerHitRelations"));
  
  
  // setup the list of supported detectors
  
  
}

enum {
  hu = 0,
  hv,
  hSize 
} ;

void DDPlanarDigiProcessor::init() { 
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);


  _h.resize( hSize ) ;

  Global::EVENTSEEDER->registerProcessor(this);

  
  if( _resU.size() !=  _resV.size() ) {
    
    std::stringstream ss ;
    ss << name() << "::init() - Inconsistent number of resolutions given for U and V coordinate: " 
       << "ResolutionU  :" <<   _resU.size() << " != ResolutionV : " <<  _resV.size() ;

    throw EVENT::Exception( ss.str() ) ;
  }

  //============ populate surface map =============================

  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  DD4hep::Geometry::DetElement world = lcdd.world() ;
  
  DD4hep::Geometry::DetElement::Children detectors = world.children() ;

  lcio::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;


  streamlog_out( DEBUG3 ) << " DDPlanarDigiProcessor::init(): searching for detector with id = " << _sub_det_id << std::endl ;

  for ( DD4hep::Geometry::DetElement::Children::const_iterator it=detectors.begin() ; it != detectors.end() ; ++it ){
    
    DD4hep::Geometry::DetElement det = (*it).second ;
    
    if( det.id() == _sub_det_id ){
        
      DD4hep::DDRec::SurfaceManager ds( det ) ;
    
      const DD4hep::DDRec::SurfaceList& detSL = ds.surfaceList() ;

      for( DD4hep::DDRec::SurfaceList::const_iterator it = detSL.begin() ; it != detSL.end() ; ++it ){

        DD4hep::DDRec::Surface* surf =  *it ;
        
        encoder.setValue( surf->id() ) ;
        streamlog_out( DEBUG5 ) << " ++++  found surface : " << encoder.valueString() << std::endl ;

        _map[ surf->id() ] = surf ;
      }

      break ;
      //-----
    }
  }

  streamlog_out( DEBUG3 ) << " DDPlanarDigiProcessor::init(): found " << _map.size() << " surfaces for detector with id:" <<  _sub_det_id << std::endl ;

  if( _map.empty() ) {   std::stringstream err  ; err << " Could not find any surfaces for detector with id: " << _sub_det_id << " in processor " 
                                                      <<   this->name() ;
    throw Exception( err.str() ) ;
  }
  
}


void DDPlanarDigiProcessor::processRunHeader( LCRunHeader* run) { 
  ++_nRun ;
} 

void DDPlanarDigiProcessor::processEvent( LCEvent * evt ) { 


  if( isFirstEvent() ) {

    streamlog_out( MESSAGE ) << " *** DDPlanarDigiProcessor::processEvent():   creating histograms for first event ! " << std::endl ;


    //AIDA::ICloud1D* hMCPEnergy  = 
    AIDAProcessor::histogramFactory(this) ; //->createHistogram1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ;

    _h[ hu ] = new TH1F( "hu" , "smearing u" , 50, -5. , +5. );
    _h[ hv ] = new TH1F( "hv" , "smearing v" , 50, -5. , +5. );
  }

  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG4 ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;
  
  LCCollection* STHcol = 0 ;
  try{
    STHcol = evt->getCollection( _inColName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG4) << "Collection " << _inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }
  
  if( STHcol != 0 ){    
    
    unsigned nCreatedHits=0;
    unsigned nDismissedHits=0;
    
    LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHITPLANE )  ;
    
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( lcio::ILDCellID0::encoder_string , trkhitVec ) ;

    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;
    
    CellIDDecoder<SimTrackerHit> cellid_decoder( STHcol) ;
    
    
    int nSimHits = STHcol->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << " processing collection " << _inColName  << " with " <<  nSimHits  << " hits ... " << std::endl ;
    
    for(int i=0; i< nSimHits; ++i){
      
      SimTrackerHit* simTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
      
      const int cellID0 = simTHit->getCellID0() ;
  
      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************
      
      const DD4hep::DDRec::Surface* surf = _map[ cellID0 ] ;
      
      if( ! surf ){    std::stringstream err ; err << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " 
                                                   <<   cellid_decoder( simTHit ).valueString()  ;
        throw Exception ( err.str() ) ;
      }

      int layer  = cellid_decoder( simTHit )["layer"];

      DDSurfaces::Vector3D oldPos( simTHit->getPosition()[0], simTHit->getPosition()[1], simTHit->getPosition()[2] );
      
      DDSurfaces::Vector3D newPos ;

     //************************************************************
      // Check if Hit is inside senstive 
      //************************************************************
      
      if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {
        
        streamlog_out( DEBUG3 ) << "  hit at " << oldPos 
                                << " " << cellid_decoder( simTHit).valueString() 
                                << " is not on surface " 
                                << *surf  
                                << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                                << std::endl;        
        
        ++nDismissedHits;
        
        continue; 
      }
      
      //**************************************************************************
      // Try to smear the hit but ensure the hit is inside the sensitive region
      //**************************************************************************
      
      DDSurfaces::Vector3D u = surf->u() ;
      DDSurfaces::Vector3D v = surf->v() ;
      
      // get local coordinates on surface
      DDSurfaces::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
      double uL = lv[0] / dd4hep::mm ;
      double vL = lv[1] / dd4hep::mm ;

      bool accept_hit = false ;
      unsigned  tries   =  0 ;              
      static const unsigned MaxTries = 10 ; 
      
      float resU = ( _resU.size() > 1 ?   _resU.at(  layer )     : _resU.at(0)   )  ;
      float resV = ( _resV.size() > 1 ?   _resV.at(  layer )     : _resV.at(0)   )  ; 

      while( tries < MaxTries ) {
        
        if( tries > 0 ) streamlog_out(DEBUG0) << "retry smearing for " <<  cellid_decoder( simTHit ).valueString() << " : retries " << tries << std::endl;
        
        double uSmear  = gsl_ran_gaussian( _rng, resU ) ;
        double vSmear  = gsl_ran_gaussian( _rng, resV ) ;
        
        // DDSurfaces::Vector3D newPosTmp = oldPos +  uSmear * u ;  
        // if( ! _isStrip )  newPosTmp = newPosTmp +  vSmear * v ;  
        
        
        
        DDSurfaces::Vector3D newPosTmp = 1./dd4hep::mm  * ( ! _isStrip  ? 
                                                            surf->localToGlobal( DDSurfaces::Vector2D (  ( uL + uSmear ) * dd4hep::mm, ( vL + vSmear )  *dd4hep::mm ) )  :
                                                            surf->localToGlobal( DDSurfaces::Vector2D (  ( uL + uSmear ) * dd4hep::mm,          0.                  ) ) 
                                                            ) ;

        streamlog_out( DEBUG1 ) << " hit at    : " << oldPos 
                                << " smeared to: " << newPosTmp
                                << " uL: " << uL 
                                << " vL: " << vL 
                                << " uSmear: " << uSmear
                                << " vSmear: " << vSmear
                                << std::endl ;


        if ( surf->insideBounds( dd4hep::mm * newPosTmp ) ) {    
          
          accept_hit = true ;
          newPos     = newPosTmp ;

          _h[hu]->Fill(  uSmear / resU ) ; 
          _h[hv]->Fill(  vSmear / resV ) ; 


          break;  

        } else { 
          
          streamlog_out( DEBUG1 ) << "  hit at " << newPosTmp 
                                  << " " << cellid_decoder( simTHit).valueString() 
                                  << " is not on surface " 
                                  << " distance: " << surf->distance( dd4hep::mm * newPosTmp ) 
                                  << std::endl;        
        }
        
        ++tries;
       }
      
      
      if( accept_hit == false ) {
        streamlog_out(DEBUG4) << "hit could not be smeared within ladder after " << MaxTries << "  tries: hit dropped"  << std::endl;
        ++nDismissedHits;
        continue; 
      } 
      
      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************
      

      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;
                  
      const int cellID1 = simTHit->getCellID1() ;
      trkHit->setCellID0( cellID0 ) ;
      trkHit->setCellID1( cellID1 ) ;
      
      trkHit->setPosition( newPos.const_array()  ) ;

      trkHit->setEDep( simTHit->getEDep() ) ;

      float u_direction[2] ;
      u_direction[0] = u.theta();
      u_direction[1] = u.phi();
      
      float v_direction[2] ;
      v_direction[0] = v.theta();
      v_direction[1] = v.phi();
      
      streamlog_out(DEBUG0)  << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1] 
                             << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
                             << std::endl ;

      trkHit->setU( u_direction ) ;
      trkHit->setV( v_direction ) ;
      
      trkHit->setdU( resU ) ;

      if( _isStrip ) trkHit->setdV( 0 ); // no error in v direction for strip hits as there is no meesurement information in v direction
      else trkHit->setdV( resV ) ;
      

      if( _isStrip ){
        trkHit->setType( UTIL::set_bit( trkHit->getType() ,  UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) ) ;
      }

      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************    
         
      LCRelationImpl* rel = new LCRelationImpl;

      rel->setFrom (trkHit);
      rel->setTo (simTHit);
      rel->setWeight( 1.0 );
      relCol->addElement(rel);

      
      //**************************************************************************
      // Add hit to collection
      //**************************************************************************    
      
      trkhitVec->addElement( trkHit ) ; 
      
      ++nCreatedHits;
      
      streamlog_out(DEBUG3) << "-------------------------------------------------------" << std::endl;
      
    }      
    
    
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
    
    evt->addCollection( trkhitVec , _outColName ) ;
    evt->addCollection( relCol , _outRelColName ) ;
    
    streamlog_out(DEBUG4) << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed as not on sensitive element\n";
    
  }
  _nEvt ++ ;
}



void DDPlanarDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DDPlanarDigiProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
  << " processed " << _nEvt << " events in " << _nRun << " runs "
  << std::endl ;
  
}
