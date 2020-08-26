#include "FilterDoubleLayerHits.h"
#include <iostream>
#include <cstdlib>
#include <string.h>

#include <IMPL/LCCollectionVec.h>
// #include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
// #include <EVENT/SimCalorimeterHit.h>
// #include <EVENT/CalorimeterHit.h>

#include "marlin/AIDAProcessor.h"

#include <UTIL/CellIDDecoder.h>
// #include <UTIL/BitField64.h>
// #include "UTIL/LCTrackerConf.h"
// #include <UTIL/LCTOOLS.h>

#include "DD4hep/DD4hepUnits.h"

#include <TH1F.h>
#include <TH2F.h>


using namespace lcio ;
using namespace marlin ;


FilterDoubleLayerHits aFilterDoubleLayerHits ;

FilterDoubleLayerHits::FilterDoubleLayerHits() : Processor("FilterDoubleLayerHits") ,
			   _nRun(0), _nEvt(0) {

  // modify processor description
  _description = "remove hits in double layers if they don't have a corresponding close-by hit in the other sublayer" ;


  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter( "InputCollection" ,
                             "Name of the input collection with hits" ,
                             _inColName ,
                             std::string("VXDTrackerHitPlanes")
                             );

  registerProcessorParameter( "OutputCollection" ,
                             "Name of the output collection with filtered hits" ,
                             _outColName ,
                             std::string("VXDTrackerHitPlanes_DLFiltered")
                             );

  registerProcessorParameter( "SubDetectorName" ,
                             "Name of dub detector" ,
                             _subDetName ,
                             std::string("Vertex") );

  registerProcessorParameter( "FillHistograms" ,
                             "Whether to fill diagnostic histograms" ,
                             _fillHistos ,
                             false );

  StringVec dlCutConfigsEx ;
  dlCutConfigsEx.push_back("0") ;
  dlCutConfigsEx.push_back("1") ;
  dlCutConfigsEx.push_back("0.5") ;
  dlCutConfigsEx.push_back("0.05") ;

  registerProcessorParameter("DoubleLayerCuts" ,
			     "Layer IDs and cuts to be applied: <layer 0> <layer 1> <dX> <dTheta>" ,
			     _dlCutConfigs ,
			     dlCutConfigsEx
			     );


}


dd4hep::rec::Vector2D FilterDoubleLayerHits::globalToLocal(long int cellID, const dd4hep::rec::Vector3D& posGlobal) {
  // Finding the surface corresponding to  the cellID
  dd4hep::rec::SurfaceMap::const_iterator surfIt = _map->find( cellID );
  if( surfIt == _map->end() ){
    std::stringstream err;
    err << " FilterDoubleLayerHits::processEvent(): no surface found for cellID: " << cellID;
    throw Exception ( err.str() );
  }
  const dd4hep::rec::ISurface* surf = surfIt->second;
  // Converting global position to local
  return surf->globalToLocal(  posGlobal ) ;
}


void FilterDoubleLayerHits::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  printParameters() ;

  // Extracting double-layer cut configurations
  _dlCuts.resize( _dlCutConfigs.size() / 4 ) ;

  unsigned i=0,index=0 ;
  while( i < _dlCutConfigs.size() ){
    _dlCuts[index].layer0      = std::atoi( _dlCutConfigs[ i++ ].c_str() ) ;
    _dlCuts[index].layer1      = std::atoi( _dlCutConfigs[ i++ ].c_str() ) ;
    _dlCuts[index].dX_max      = std::atof( _dlCutConfigs[ i++ ].c_str() ) ;
    _dlCuts[index].dTheta_max  = std::atof( _dlCutConfigs[ i++ ].c_str() ) / 1e3 ;  // converting mrad -> rad
    ++index ;
  }

  // Get the surface map from the SurfaceManager
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>();
  dd4hep::DetElement det = theDetector.detector( _subDetName );

  _map = surfMan.map( det.name() );

  if( !_map ) {
    std::stringstream err;
    err << " Could not find surface map for detector: " << _subDetName << " in SurfaceManager";
    throw Exception( err.str() );
  }

  // Booking diagnostic histograms for each configured cut
  AIDAProcessor::histogramFactory(this);
  char hname[100];
  for(size_t iCut=0, nCuts=_dlCuts.size() ; iCut<nCuts ; ++iCut){
    const DoubleLayerCut& cut = _dlCuts.at(iCut);
    if (_fillHistos) {
      sprintf(hname, "h_dX_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[ std::string(hname) ] = new TH1F( hname , ";dX [mm]; Hit pairs", 5000, -25, 25 );
      sprintf(hname, "h2_dX_dPhi_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[ std::string(hname) ] = new TH2F( hname , ";dX [mm]; d#phi [mrad]", 500, -5, 5, 200, -20, 20);
      sprintf(hname, "h_dTheta_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[ std::string(hname) ] = new TH1F( hname , ";d#theta [mrad]; Hit pairs", 5000, -10, 10 );
      sprintf(hname, "h_dPhi_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[ std::string(hname) ] = new TH1F( hname , ";d#phi [mrad]; Hit pairs", 5000, -100, 100 );
      sprintf(hname, "h_posU_layer%d", cut.layer0);
      _histos[ std::string(hname) ] = new TH1F( hname , ";U [mm]; Hits", 1000, -50, 50 );
      sprintf(hname, "h_posV_layer%d", cut.layer0);
      _histos[ std::string(hname) ] = new TH1F( hname , ";V [mm]; Hits", 1000, -200, 200 );
    }
    // Printing the configured cut
    streamlog_out( DEBUG5 ) <<  iCut << ". layers: " << cut.layer0 << " >> " << cut.layer1 << ";  dX: "
    << cut.dX_max << " mm;  dTheta: " << cut.dTheta_max << " rad" << std::endl;
  }

  _nRun = 0 ;
  _nEvt = 0 ;

}


void FilterDoubleLayerHits::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
}




void FilterDoubleLayerHits::processEvent( LCEvent * evt ) {


  LCCollection* col(0);

  try {
    col = evt->getCollection( _inColName );
  }
  catch( lcio::DataNotAvailableException& e) {
    streamlog_out( WARNING ) <<  " input collection not in event : " << _inColName << "   - nothing to do    !!! " << std::endl;
    return;
  }


  std::string encoderString = col->getParameters().getStringVal( "CellIDEncoding" );
  UTIL::CellIDDecoder<TrackerHitPlane> decoder( encoderString ) ;

  //---- create the output collection
  LCCollectionVec* colOut = new LCCollectionVec( col->getTypeName() );
  colOut->setSubset( true ) ;
  colOut->parameters().setValue( "CellIDEncoding", encoderString );

  // Set acceptance flags for all hits to FALSE
  const size_t nHit = col->getNumberOfElements();
  memset(&_hitAccepted, false, nHit);

  //---- loop over hits
  for(size_t iHit = 0; iHit < nHit ; iHit++){

    // Skipping hits that are already accepted
    if (_hitAccepted[iHit]) continue;

    TrackerHitPlane* h = (TrackerHitPlane*)col->getElementAt( iHit );

    unsigned int layerID = decoder(h)["layer"];
    unsigned int sideID = decoder(h)["side"];
    unsigned int ladderID = decoder(h)["module"];
    unsigned int moduleID = decoder(h)["sensor"];
    streamlog_out( DEBUG0 ) << " Checking 1st hit at layer: " << layerID << "  ladder: " << ladderID << "  module: " << moduleID <<  std::endl ;


    // Checking if the hit is at the inner double layer to be filtered
    const DoubleLayerCut* dlCut(0);
    for(int iCut=0, nCuts=_dlCuts.size(); iCut<nCuts; ++iCut){
      const DoubleLayerCut& cut = _dlCuts.at(iCut);
      if( ( layerID != cut.layer0 ) && ( layerID != cut.layer1 ) ) {
        _hitAccepted[iHit] = true;
        continue;
      }
      dlCut = &cut;
      break;
    }

    // Accepting hit immediately if it's not affected by the double-layer cuts
    if (!dlCut) {
      _hitAccepted[iHit] = true;
      continue;
    }

    // Skipping the hit from the first pass if it belongs to the outer sublayer of the cut
    if (layerID == dlCut->layer1) continue;

    // Getting local and global hit positions
    dd4hep::rec::Vector3D posGlobal( h->getPosition()[0], h->getPosition()[1], h->getPosition()[2] );
    dd4hep::rec::Vector2D posLocal = globalToLocal( decoder( h ).getValue(), posGlobal );

    // Filling diagnostic histograms
    if (_fillHistos) {
      char hname[100];
      sprintf(hname, "h_posU_layer%d", layerID);
      _histos[hname]->Fill(posLocal.u());
      sprintf(hname, "h_posV_layer%d", layerID);
      _histos[hname]->Fill(posLocal.v());
    }

    // Looking for the compliment hits in the 2nd sublayer
    size_t nCompatibleHits(0);
    unsigned int cutLayerID = layerID == dlCut->layer0 ? dlCut->layer1 : dlCut->layer0;
    for(size_t iHit2 = 0; iHit2 < nHit; iHit2++){
      if (iHit2 == iHit) continue;
      TrackerHitPlane* h2 = (TrackerHitPlane*)col->getElementAt( iHit2 );

      // Ensuring the hit is on the proper layer
      unsigned int layerID2 = decoder(h2)["layer"];
      if (layerID2 != cutLayerID) continue;

      // Ensuring the hit is on the same side, ladder and module as the first hit
      if (decoder(h2)["side"] != sideID) continue;
      if (decoder(h2)["module"] != ladderID) continue;
      if (decoder(h2)["sensor"] != moduleID) continue;

      // Getting the local and global hit positions
      dd4hep::rec::Vector3D posGlobal2( h2->getPosition()[0], h2->getPosition()[1], h2->getPosition()[2] );
      dd4hep::rec::Vector2D posLocal2 = globalToLocal( decoder( h2 ).getValue(), posGlobal2 );

      // Checking whether hit is close enough to the 1st one
      double dX = posLocal2.u() - posLocal.u();
      double dTheta = posGlobal2.theta() - posGlobal.theta();
      double dPhi = posGlobal2.phi() - posGlobal.phi();
      streamlog_out( DEBUG0 ) << " Checking 2nd hit at layer: " << layerID2 << ";  dX: " << dX << ";  dTheta: " <<  dTheta << std::endl;

      // Filling diagnostic histograms
      if (_fillHistos) {
        char hname[100];
        sprintf(hname, "h_dX_layers_%d_%d", layerID, layerID2);
        _histos[hname]->Fill(dX);
        sprintf(hname, "h2_dX_dPhi_layers_%d_%d", layerID, layerID2);
        _histos[hname]->Fill(dX, dPhi*1e3);
        sprintf(hname, "h_dTheta_layers_%d_%d", layerID, layerID2);
        _histos[hname]->Fill(dTheta*1e3);
        sprintf(hname, "h_dPhi_layers_%d_%d", layerID, layerID2);
        _histos[hname]->Fill(dPhi*1e3);
      }

      // Skipping the hit if it's outside the cut window
      if (fabs(dX) > dlCut->dX_max) continue;
      if (fabs(dTheta) > dlCut->dTheta_max) continue;

      nCompatibleHits++;
      _hitAccepted[iHit2] = true;
      streamlog_out( DEBUG5 ) << " Accepted 2nd hit at layer: " << layerID2 << ";  dX: " << dX << ";  dTheta: " <<  dTheta << std::endl;
    }

    // Accepting the first hit if it has at least one compatible pair
    if (nCompatibleHits > 0) {
      _hitAccepted[iHit] = true;
      streamlog_out( DEBUG5 ) << " Accepted 1st hit at layer: " << layerID << std::endl;
    }
  }

  // Adding accepted hits to the output collection
  size_t nHitsAccepted(0);
  for(size_t iHit = 0; iHit < nHit; iHit++){
    if (!_hitAccepted[iHit]) continue;
  	colOut->addElement( col->getElementAt( iHit ) );
    nHitsAccepted++;
  }
	streamlog_out( MESSAGE ) << " " << nHitsAccepted << " hits added to collection: " << _outColName << std::endl;


  evt->addCollection(  colOut , _outColName ) ;
  streamlog_out( DEBUG5 ) << " output collection " << _outColName << " of type " <<  col->getTypeName() << " added to the event  " << std::endl ;

  _nEvt ++ ;
}



void FilterDoubleLayerHits::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void FilterDoubleLayerHits::end(){

  streamlog_out( MESSAGE ) << "FilterDoubleLayerHits::end()  " << name()
			   << " processed " << _nEvt << " events in " << _nRun << " runs "
			   << std::endl ;

}


