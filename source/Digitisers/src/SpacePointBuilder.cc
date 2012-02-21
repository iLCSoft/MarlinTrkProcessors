#include "SpacePointBuilder.h"

#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/TrackerHitPlaneImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCRelationImpl.h"
#include "UTIL/ILDConf.h"
#include "UTIL/LCRelationNavigator.h"


#include "marlin/VerbosityLevels.h"

#include "gear/GEAR.h"
#include "gear/GearParameters.h"
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"
#include "gear/ZPlanarParameters.h"
#include "gear/ZPlanarLayerLayout.h"
#include "marlin/Global.h"
#include "UTIL/ILDConf.h"

//FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
#include "MarlinTrk/Factory.h"

#include "gear/gearsurf/MeasurementSurfaceStore.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/ICoordinateSystem.h"
#include "gear/gearsurf/CartesianCoordinateSystem.h"

#include <cmath>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace CLHEP;

SpacePointBuilder aSpacePointBuilder ;


SpacePointBuilder::SpacePointBuilder() : Processor("SpacePointBuilder") {

   // modify processor description
   _description = "SpacePointBuilder combine si-strip measurements into 3D spacepoints (1TrackerHitPlanar+1TrackHitPlanar = 1 TrackerHit), that can be used by reconstruction" ;


   // register steering parameters: name, description, class-variable, default value
   registerInputCollection(LCIO::TRACKERHIT,
                           "TrackerHitCollection",
                           "TrackerHitCollection",
                           _TrackerHitCollection,
                           std::string("FTDTrackerHits")); 

   registerInputCollection(LCIO::LCRELATION,
                           "TrackerHitSimHitRelCollection",
                           "The name of the input collection of the relations of the TrackerHits to SimHits",
                           _TrackerHitSimHitRelCollection,
                           std::string("FTDTrackerHitRelations")); 
   
   registerOutputCollection(LCIO::TRACKERHIT,
                            "SpacePointsCollection",
                            "SpacePointsCollection",
                            _SpacePointsCollection,
                            std::string("FTDSpacePoints"));

   registerOutputCollection(LCIO::LCRELATION,
                            "SimHitSpacePointRelCollection",
                            "Name of the SpacePoint SimTrackerHit relation collection",
                            _relColName,
                            std::string("FTDSimHitSpacepointRelations"));
    
   
}




void SpacePointBuilder::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;


  //FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with KalTest/KalDet
  MarlinTrk::IMarlinTrkSystem* trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "KalTest" , marlin::Global::GEAR , "" ) ;
  
  
  if( trksystem == 0 ) {
    
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest" )  ) ;
    
  }
  
  trksystem->init() ;  
  
  //FIXME:SJA gear surface store has now been filled so we can dispose of the MarlinTrkSystem
  delete trksystem;

  
}


void SpacePointBuilder::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 



void SpacePointBuilder::processEvent( LCEvent * evt ) { 


  LCCollection* col = evt->getCollection( _TrackerHitCollection ) ;
  LCRelationNavigator* nav = new LCRelationNavigator(evt->getCollection( _TrackerHitSimHitRelCollection ));

  if( col != NULL ){
    
    
    unsigned createdSpacePoints = 0;
    unsigned rawStripHits = 0;
    unsigned possibleSpacePoints = 0;
    _nOutOfBoundary = 0;
    _nStripsTooParallel = 0;
    _nPlanesNotParallel = 0;
    
    
    LCCollectionVec * spCol = new LCCollectionVec(LCIO::TRACKERHIT);    // output spacepoint collection
    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);    // outpur relation collection
    
    // to store the weights
    LCFlagImpl lcFlag(0) ;
    lcFlag.setBit( LCIO::LCREL_WEIGHTED ) ;
    relCol->setFlag( lcFlag.getFlag()  ) ;
    
    
    unsigned nHits = col->getNumberOfElements()  ;
    
    streamlog_out( DEBUG4 ) << "Number of hits: " << nHits <<"\n";
    
    //store hits in map according to their CellID0
    std::map< int , std::vector< TrackerHitPlane* > > map_cellID0_hits;
    std::map< int , std::vector< TrackerHitPlane* > >::iterator it;
    
    for( unsigned i=0; i<nHits; i++){
      
      TrackerHitPlane* trkHit = dynamic_cast<TrackerHitPlane*>( col->getElementAt( i ) );
      if( trkHit != NULL) map_cellID0_hits[ trkHit->getCellID0() ].push_back( trkHit ); 
      
    }
    
    
    
    // now loop over all CellID0s
    for( it= map_cellID0_hits.begin(); it!= map_cellID0_hits.end(); it++ ){
     
      
      rawStripHits += it->second.size();
      
      std::vector< TrackerHitPlane* > hitsFront = it->second;
      
      int cellID0 = it->first;
     
      //get the CellID0s at the back of this sensor
      std::vector< int > cellID0sBack = getCellID0sAtBack( cellID0 );
     
      
      for( unsigned i=0; i< cellID0sBack.size(); i++ ){ 
        
        
        int cellID0Back = cellID0sBack[i];
        std::vector< TrackerHitPlane* > hitsBack = map_cellID0_hits[ cellID0Back ];
        
        streamlog_out( DEBUG3 ) << "strips: CellID0 " << cellID0  << " " << getCellID0Info( cellID0 )  << "(" << hitsFront.size()
          << " hits) <---> CellID0 " << cellID0Back << getCellID0Info( cellID0Back )
          << "(" << hitsBack.size() << " hits)\n"
          << "--> " << hitsFront.size() * hitsBack.size() << " possible combinations\n";
        
        possibleSpacePoints += hitsFront.size() * hitsBack.size();
        
        
        // Now iterate over all combinations and store those that make sense
        for( unsigned i=0; i<hitsFront.size(); i++ ){
          
          TrackerHitPlane* hitFront = hitsFront[i];
          
          for( unsigned j=0; j<hitsBack.size(); j++ ){
            
            
            TrackerHitPlane* hitBack = hitsBack[j];
            
            
            TrackerHitImpl* spacePoint = createSpacePoint( hitFront, hitBack );
            if ( spacePoint == NULL ) continue;
            
            
            CellIDEncoder<TrackerHitImpl> cellid_encoder( ILDCellID0::encoder_string , spCol );
            cellid_encoder.setValue( cellID0 ); //give the new hit, the CellID0 of the front hit
            cellid_encoder.setCellID( spacePoint ) ;
            
            // store the hits it's composed of:
            spacePoint->rawHits().push_back( hitFront );
            spacePoint->rawHits().push_back( hitBack );
            
            spacePoint->setType( UTIL::set_bit( spacePoint->getType() ,  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT ) ) ;
           
            spCol->addElement( spacePoint ) ; 
            
            createdSpacePoints++;
            
            
            ///////////////////////////////
            // make the relations
            const LCObjectVec& simHitsFront = nav->getRelatedToObjects( hitFront );
            
            if( simHitsFront.size() == 1 ){
              
              SimTrackerHit* simHit = dynamic_cast< SimTrackerHit* >( simHitsFront[0] );
              
              if( simHit != NULL ){
                LCRelationImpl* rel = new LCRelationImpl;
                rel->setFrom (spacePoint);
                rel->setTo  (simHit);
                rel->setWeight( 0.5 );
                relCol->addElement(rel);
              }
            }
            
            const LCObjectVec& rawObjectsBack = nav->getRelatedToObjects( hitBack );
            
            if( rawObjectsBack.size() == 1 ){
              
              SimTrackerHit* simHit = dynamic_cast< SimTrackerHit* >( rawObjectsBack[0] );
              
              if( simHit != NULL ){
                LCRelationImpl* rel = new LCRelationImpl;
                rel->setFrom (spacePoint);
                rel->setTo  (simHit);
                rel->setWeight( 0.5 );
                relCol->addElement(rel);
              }
            }
            //////////////////////////////////
            
            
          }
          
        }
        
      }
      
    }
    
    evt->addCollection( spCol, _SpacePointsCollection);
    evt->addCollection( relCol , _relColName ) ;
    
    streamlog_out( DEBUG4 )<< "\nCreated " << createdSpacePoints
      << " space points ( raw strip hits: " << rawStripHits << ")\n";
    
    streamlog_out( DEBUG3 ) << "  There were " << rawStripHits << " strip hits available, giving " 
      << possibleSpacePoints << " possible space points\n";
    
    streamlog_out( DEBUG3 ) << "  " << _nStripsTooParallel << " space points couldn't be created, because the strips were too parallel\n";
    streamlog_out( DEBUG3 ) << "  " << _nPlanesNotParallel << " space points couldn't be created, because the planes of the measurement surfaces where not parallel enough\n";
    streamlog_out( DEBUG3 ) << "  " << _nOutOfBoundary     << " space points couldn't be created, because the result was outside the sensor boundary\n"; 
    
    
    streamlog_out( DEBUG4 ) << "\n";
    
  }


  _nEvt ++ ;
  
  delete nav;
  
}





void SpacePointBuilder::check( LCEvent * evt ) {}


void SpacePointBuilder::end(){
   
   
}

TrackerHitImpl* SpacePointBuilder::createSpacePoint( TrackerHitPlane* a , TrackerHitPlane* b ){
  
  const double* p1 = a->getPosition();
  double x1 = p1[0];
  double y1 = p1[1];
  double z1 = p1[2];
  Hep3Vector P1( x1,y1,z1 );
  
  
  
  gear::MeasurementSurface const* ms1 = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( a->getCellID0() );
  gear::CartesianCoordinateSystem* ccs1 = dynamic_cast< gear::CartesianCoordinateSystem* >( ms1->getCoordinateSystem() );
  CLHEP::Hep3Vector W1 = ccs1->getLocalZAxis(); // the vector W of the local coordinate system the measurement surface has
  CLHEP::Hep3Vector V1 = ccs1->getLocalYAxis(); // the vector W of the local coordinate system the measurement surface has
  
  const double* p2 = b->getPosition();
  double x2 = p2[0];
  double y2 = p2[1];
  double z2 = p2[2];
  Hep3Vector P2( x2,y2,z2 );
  
  gear::MeasurementSurface const* ms2 = Global::GEAR->getMeasurementSurfaceStore().GetMeasurementSurface( b->getCellID0() );
  gear::CartesianCoordinateSystem* ccs2 = dynamic_cast< gear::CartesianCoordinateSystem* >( ms2->getCoordinateSystem() );
  CLHEP::Hep3Vector W2 = ccs2->getLocalZAxis(); // the vector W of the local coordinate system the measurement surface has
  CLHEP::Hep3Vector V2 = ccs2->getLocalYAxis(); // the vector W of the local coordinate system the measurement surface has
  
  
  streamlog_out( DEBUG2 ) << "\t ( " << x1 << " " << y1 << " " << z1 << " ) <--> ( " << x2 << " " << y2 << " " << z2 << " )\n";
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // First: check if the two measurement surfaces are parallel (i.e. the w are parallel or antiparallel)
  double angle = fabs(W2.angle(W1));
  double angleMax = 1.*M_PI/180.;
  if(( angle > angleMax )&&( angle < M_PI-angleMax )){
    
    _nPlanesNotParallel++;
    streamlog_out( DEBUG2 ) << "\tThe planes of the measurement surfaces are not parallel enough, the angle between the W vectors is " << angle
    << " where the angle has to be smaller than " << angleMax << " or bigger than " << M_PI-angleMax << "\n\n";
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Next: check if the angle between the strips is not 0
  angle = fabs(V2.angle(V1));
  double angleMin= 1.*M_PI/180.;
  if(( angle < angleMin )||( angle > M_PI-angleMin )){
    
    _nStripsTooParallel++;
    streamlog_out( DEBUG2 ) << "\tThe strips (V vectors) of the measurement surfaces are too parallel, the angle between the V vectors is " << angle
    << " where the angle has to be between " << angleMax << " or bigger than " << M_PI-angleMin << "\n\n";
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
    
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Next we want to calculate the crossing point.
  
  Hep3Vector point;
  calculatePointBetweenTwoLines( P1, V1, P2, V2, point );
  
  streamlog_out( DEBUG2 ) << "\tPosition of space point (global) : ( " << point.x() << " " << point.y() << " " << point.z() << " )\n";
 
  
  
  // Check if the new hit is within the boundaries
  CLHEP::Hep3Vector localPoint1 = ccs1->getLocalPoint(point);
  localPoint1.setZ( 0. ); // we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  
  CLHEP::Hep3Vector localPoint2 = ccs2->getLocalPoint(point);
  localPoint2.setZ( 0. ); // we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  
  
  if( !ms1->isLocalInBoundary( localPoint1 ) ){
    
    _nOutOfBoundary++;
    streamlog_out( DEBUG2 ) << "\tFirst hit is out of boundary: local coordinates are ( " 
    << localPoint1.x() << " " << localPoint1.y() << " " << localPoint1.z() << " )\n\n";
    
    return NULL;
    
  }
  if( !ms2->isLocalInBoundary( localPoint2 ) ){
    
    _nOutOfBoundary++;
    streamlog_out( DEBUG2 ) << "\tSecond hit is out of boundary: local coordinates are ( " 
    << localPoint2.x() << " " << localPoint2.y() << " " << localPoint2.z() << " )\n\n";
    
    return NULL;
    
  }
  
  
  //Create the new TrackerHit
  TrackerHitImpl* spacePoint = new TrackerHitImpl();
  
  double pos[3] = {point.x(), point.y(), point.z() };
  spacePoint->setPosition(  pos  ) ;
  
  streamlog_out( DEBUG2 ) << "\tHit accepted\n\n";
  
  return spacePoint;
  
}
/*
TrackerHitImpl* SpacePointBuilder::createSpacePointOld( TrackerHitPlane* a , TrackerHitPlane* b ){
  
  streamlog_out( DEBUG2 ) << "\t OLD OLD OLD OLD\n";
  
  
  const double* p1 = a->getPosition();
  double x1 = p1[0];
  double y1 = p1[1];
  double z1 = p1[2];
  const float* v1 = a->getV();
  float ex1 = cos( v1[1] ) * sin( v1[0] ); 
  float ey1 = sin( v1[1] ) * sin( v1[0] );
  
  const double* p2 = b->getPosition();
  double x2 = p2[0];
  double y2 = p2[1];
  double z2 = p2[2];
  const float* v2 = b->getV();
  float ex2 = cos( v2[1] ) * sin( v2[0] ); 
  float ey2 = sin( v2[1] ) * sin( v2[0] );
  
  streamlog_out( DEBUG2 ) << "\t ( " << x1 << " " << y1 << " " << z1 << " ) <--> ( " << x2 << " " << y2 << " " << z2 << " )\n";
  
  double x=0.;
  double y=0.;
  
  if ( calculateXingPoint( x1, y1, ex1, ey1, x2, y2, ex2, ey2, x, y ) != 0 ){
    
    _nStripsTooParallel++;
    streamlog_out( DEBUG2 ) << "\tStrips too parallel\n\n";
    return NULL; //calculate the xing point and if that fails don't create a spacepoint
  
  }
  
  double z= (z1 + z2)/2.;
  
  streamlog_out( DEBUG2 ) << "\tPosition of space point (global) : ( " << x << " " << y << " " << z << " )\n";
  
  // Check if the new hit is within the boundary
  CLHEP::Hep3Vector globalPoint(x,y,z);
  //  gear::MeasurementSurface* ms = gear::MeasurementSurfaceStore::Instance().GetMeasurementSurface( a->getCellID0() );
  CLHEP::Hep3Vector localPoint = ms->getCoordinateSystem()->getLocalPoint(globalPoint);
  localPoint.setZ( 0. ); // we set w to 0 so it is in the plane ( we are only interested if u and v are in or out of range, to exclude w from the check it is set to 0)
  if( !ms->isLocalInBoundary( localPoint ) ){
    
    _nOutOfBoundary++;
    streamlog_out( DEBUG2 ) << "\tHit is out of boundary: local coordinates are ( " 
      << localPoint.x() << " " << localPoint.y() << " " << localPoint.z() << " )\n\n";
    
    return NULL;
    
  }
  
  
  //Create the new TrackerHit
  TrackerHitImpl* spacePoint = new TrackerHitImpl();
  
  double pos[3] = {x,y,z};
  spacePoint->setPosition(  pos  ) ;
  
  streamlog_out( DEBUG2 ) << "\tHit accepted\n\n";

  return spacePoint;
  
}
*/

int SpacePointBuilder::calculatePointBetweenTwoLines( const Hep3Vector& P1, const Hep3Vector& V1, const Hep3Vector& P2, const Hep3Vector& V2, Hep3Vector& point ){
  
  // Richgungsvektor normal auf die anderen beiden:
  Hep3Vector n = V1.cross( V2 );
  
  // Now we want to rotate into a coordinate system, where n is parallel to the z axis
  // For this: first set phi to 0
  // then: set theta to 0 (we set phi to 0 first, so we can then rotate arount the y axis)
  HepRotation rot;
  rot.rotateZ( -n.phi() );
  Hep3Vector nPrime = rot * n; //now the phi of nPrime should be 0
  streamlog_out( DEBUG0 ) << "phi of n' = " << nPrime.phi() << " (it should be 0!!!)\n";
  rot.rotateY( -n.theta() );
  nPrime = rot * n;
  streamlog_out( DEBUG0 ) << "phi of n'' = " << nPrime.phi() << " (it should be 0!!!)\n";
  streamlog_out( DEBUG0 ) << "theta of n'' = " << nPrime.theta() <<  " (it should be 0!!!)\n";
  
  // Now rotate all the vectors and points into this coordinatesystem.
  Hep3Vector P1prime = rot * P1;
  Hep3Vector V1prime = rot * V1;
  Hep3Vector P2prime = rot * P2;
  Hep3Vector V2prime = rot * V2;
  
  // What is the gain of rotating into this system?
  // A: 
  double x;
  double y;
  int res = calculateXingPoint( P1prime.x(), P1prime.y(), V1prime.x(), V1prime.y(), P2prime.x(), P2prime.y(), V2prime.x(), V2prime.y(), x, y );
  
  if ( res != 0 ) return 1;
  
  point.setX( x );
  point.setY( y );
  point.setZ( (P1prime.z() + P2prime.z())/2. );
  
  // Now transform back to the global coordinates
  point = rot.inverse() * point;
  
  
  return 0;
  
}


int SpacePointBuilder::calculateXingPoint( double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2, double& x, double& y ){


  float a = (x1*ey1 - y1*ex1) - (x2*ey1 - y2*ex1);
  float b = ex2*ey1 - ex1*ey2;

  const float epsilon = 0.00001;

  if( fabs(b) < epsilon ) return 1; // if b==0 the two directions e1 and e2 are parallel and there is no crossing!

  float t = a/b;

  x = x2 + t*ex2;
  y = y2 + t*ey2;

  return 0;

  

}
 
 
std::vector< int > SpacePointBuilder::getCellID0sAtBack( int cellID0 ){
  
  std::vector< int > back;
  
  UTIL::BitField64  cellID( ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
  int subdet = cellID[ ILDCellID0::subdet ] ;
  
  if( subdet == ILDDetID::FTD ) return getCellID0sAtBackOfFTD( cellID0 );
  if( subdet == ILDDetID::SIT ) return getCellID0sAtBackOfSIT( cellID0 );
  if( subdet == ILDDetID::SET ) return getCellID0sAtBackOfSET( cellID0 );
  
  return back;
  
}


std::vector< int > SpacePointBuilder::getCellID0sAtBackOfFTD( int cellID0 ){
  
  std::vector< int > back;
  
  const gear::FTDLayerLayout& ftdLayers = Global::GEAR->getFTDParameters().getFTDLayerLayout();

  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
//   int side   = cellID[ ILDCellID0::side ];
//   int module = cellID[ ILDCellID0::module ];
  int sensor = cellID[ ILDCellID0::sensor ];
  int layer  = cellID[ ILDCellID0::layer ];
  
  
  //check if sensor is in front
  if(( ftdLayers.isDoubleSided( layer ) ) && ( sensor <= ftdLayers.getNSensors( layer ) / 2 ) ){
   
    cellID[ ILDCellID0::sensor ] = sensor + ftdLayers.getNSensors( layer ) / 2; 
    // it is assumed (according to current gear and mokka), that sensors 1 until n/2 will be on front
    // and sensor n/2 + 1 until n are at the back
    // so the sensor x, will have sensor x+n/2 at the back
    
    back.push_back( cellID.lowWord() );
    
  }
  
  return back;
  
  
}

std::vector< int > SpacePointBuilder::getCellID0sAtBackOfSIT( int cellID0 ){
  
  std::vector< int > back;
  
//   const gear::ZPlanarLayerLayout& sitLayout = Global::GEAR->getSITParameters().getZPlanarLayerLayout();
  
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
//   int side   = cellID[ ILDCellID0::side ];
//   int module = cellID[ ILDCellID0::module ];
//   int sensor = cellID[ ILDCellID0::sensor ];
  int layer  = cellID[ ILDCellID0::layer ];
  
  
  //check if sensor is in front
  if( layer%2 == 0 ){ // even layers are front sensors
    
    cellID[ ILDCellID0::layer ] = layer + 1; 
    // it is assumed that the even layers are the front layers
    // and the following odd ones the back layers
    
    back.push_back( cellID.lowWord() );
    
  }
  
  return back;

  
}

std::vector< int > SpacePointBuilder::getCellID0sAtBackOfSET( int cellID0 ){
  
  std::vector< int > back;
  
//   const gear::ZPlanarLayerLayout& setLayout = Global::GEAR->getSETParameters().getZPlanarLayerLayout();
  
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );
  
  //   int side   = cellID[ ILDCellID0::side ];
  //   int module = cellID[ ILDCellID0::module ];
  //   int sensor = cellID[ ILDCellID0::sensor ];
  int layer  = cellID[ ILDCellID0::layer ];
  
  
  //check if sensor is in front
  if( layer%2 == 0 ){ // even layers are front sensors
    
    cellID[ ILDCellID0::layer ] = layer + 1; 
    // it is assumed that the even layers are the front layers
    // and the following odd ones the back layers
    
    back.push_back( cellID.lowWord() );
    
  }
  
  return back;
  
  
}






std::string SpacePointBuilder::getCellID0Info( int cellID0 ){

  std::stringstream s;
  
  //find out layer, module, sensor
  UTIL::BitField64  cellID( ILDCellID0::encoder_string );
  cellID.setValue( cellID0 );

  int subdet = cellID[ ILDCellID0::subdet ] ;
  int side   = cellID[ ILDCellID0::side ];
  int module = cellID[ ILDCellID0::module ];
  int sensor = cellID[ ILDCellID0::sensor ];
  int layer  = cellID[ ILDCellID0::layer ];
  
  s << "(su" << subdet << ",si" << side << ",la" << layer << ",mo" << module << ",se" << sensor << ")";
  
  return s.str();
  
}



