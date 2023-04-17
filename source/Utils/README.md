### FilterDoubleLayerHits
A filter processor that uses double-layer sensor arrangement for selecting hits that are close-by on the two sublayers.
It takes a collection of type `TrackerHitPlane` as input and returns a copy of it without the hits that have a pair satisfying any of the provided cuts.

A single cut defines an inner and outer `layer ID` for the pair of hits, maximum allowed distance in the local `U` coordinate and maximum allowed distance in the global `Theta` coordinate.
Hits from layers not included in the cuts are kept in the output collection.

#### Example configuration:

```xml
<processor name="FilterDL_VXDB" type="FilterDoubleLayerHits">
  <!-- Name of the corresponding detector in the geometry -->
  <parameter name="SubDetectorName" type="string" value="Vertex" />
  <!-- Name of the input hit collection -->
  <parameter name="InputCollection" type="string" value="VXDBTrackerHits" />
  <!-- Name of the output filtered hit collection -->
  <parameter name="OutputCollection" type="string" value="VXDBTrackerHits_DL" />
  <!-- Configuration of the maximum angular distance between a pair of hits in a double layer -->
  <!-- 4 numbers per double-layer: <inner layer ID>  <outer layer ID>  <dPhi max [mrad]>  <dTheta max [mrad]> -->
  <parameter name="DoubleLayerCuts" type="StringVec">
      0 1 0.6 0.35
      2 3 0.6 0.33
      4 5 0.5 0.27
      6 7 0.4 0.21
  </parameter>
  <!-- Maximum time difference between 2 hits to form a stub candidate -->
  <parameter name="DeltaTimeMax" type="float"> 0.18 </parameter>
  <!-- Whether to fill diagnostic histograms about affected hits -->
  <parameter name="FillHistograms" type="bool" value="true" />
  <!-- Verbosity level ("DEBUG0-9,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT") -->
  <parameter name="Verbosity" type="string"> MESSAGE </parameter>
</processor>
```

### FilterConeHits
A filter processor for particle gun samples overlaid to the beam-induced background. The processor selects and saves the tracker hits that are included in a cone around the MC particle trajectory along with the corresponding sim hits and the reco-sim relations.

#### Example configuration:

```xml
<processor name="FilterCone" type="FilterConeHits">
  <!-- Name of the input MC particle collection -->
  <parameter name="MCParticleCollection" type="string" value="MCParticle" />
  <!-- Name of the input hit collections -->
  <parameter name="TrackerHitInputCollections" type="StringVec">
    VBTrackerHits
    VETrackerHits
    IBTrackerHits
    IETrackerHits
    OBTrackerHits
    OETrackerHits
  </parameter>
  <!-- Name of the input simhit collections -->
  <parameter name="TrackerSimHitInputCollections" type="StringVec">
    VertexBarrelCollection
    VertexEndcapCollection
    InnerTrackerBarrelCollection
    InnerTrackerEndcapCollection
    OuterTrackerBarrelCollection
    OuterTrackerEndcapCollection
  </parameter>
  <!-- Name of the input hit relation collections -->
  <parameter name="TrackerHitInputRelations" type="StringVec">
    VBTrackerHitsRelations
    VETrackerHitsRelations
    IBTrackerHitsRelations
    IETrackerHitsRelations
    OBTrackerHitsRelations
    OETrackerHitsRelations
  </parameter>
  <!-- Name of the output filtered hit collection -->
  <parameter name="TrackerHitOutputCollections" type="StringVec">
    VBTrackerHits_Cone
    VETrackerHits_Cone
    IBTrackerHits_Cone
    IETrackerHits_Cone
    OBTrackerHits_Cone
    OETrackerHits_Cone
  </parameter>
  <!-- Name of the output filtered simhit collection -->
  <parameter name="TrackerSimHitOutputCollections" type="StringVec">
    VBTrackerSimHits_Cone
    VETrackerSimHits_Cone
    IBTrackerSimHits_Cone
    IETrackerSimHits_Cone
    OBTrackerSimHits_Cone
    OETrackerSimHits_Cone
  </parameter>
  <!-- Name of the output filtered hit relation collection -->
  <parameter name="TrackerHitOutputRelations" type="StringVec">
    VBTrackerRel_Cone
    VETrackerRel_Cone
    IBTrackerRel_Cone
    IETrackerRel_Cone
    OBTrackerRel_Cone
    OETrackerRel_Cone
  </parameter>
  <!-- Angular distance between the hits and the particle direction -->
  <parameter name="DeltaRCut" type="float" value="0.05" />
  <parameter name="FillHistograms" type="bool" value="true" />
  <parameter name="Verbosity" type="string"> MESSAGE0 </parameter>
</processor>
```


### SplitCollectionByPolarAngle
A filter processor for particle gun samples overlaid to the beam-induced background. The processor selects and saves the tracker hits with a polar angle within a certain range, defined by the lower and upper limits, along with the corresponding sim hits and the reco-sim relations.

#### Example configuration:

```xml
<processor name="SplitByPolarAngle" type="SplitCollectionByPolarAngle">
  <!-- Name of the input hit collections -->
  <parameter name="TrackerHitInputCollections" type="StringVec">
    VBTrackerHits
    VETrackerHits
    IBTrackerHits
    IETrackerHits
    OBTrackerHits
    OETrackerHits
    </parameter>
  <!-- Name of the input simhit collections -->
  <parameter name="TrackerSimHitInputCollections" type="StringVec">
    VertexBarrelCollection
    VertexEndcapCollection
    InnerTrackerBarrelCollection
    InnerTrackerEndcapCollection
    OuterTrackerBarrelCollection
    OuterTrackerEndcapCollection
    </parameter>
  <!-- Name of the input hit relation collections -->
  <parameter name="TrackerHitInputRelations" type="StringVec">
    VBTrackerHitsRelations
    VETrackerHitsRelations
    IBTrackerHitsRelations
    IETrackerHitsRelations
    OBTrackerHitsRelations
    OETrackerHitsRelations
  </parameter>
  <!-- Name of the output filtered hit collection -->
  <parameter name="TrackerHitOutputCollections" type="StringVec">
    VBTrackerHits_SplittedByPolarAngle
    VETrackerHits_SplittedByPolarAngle
    IBTrackerHits_SplittedByPolarAngle
    IETrackerHits_SplittedByPolarAngle
    OBTrackerHits_SplittedByPolarAngle
    OETrackerHits_SplittedByPolarAngle
  </parameter>
  <!-- Name of the output filtered simhit collection -->
  <parameter name="TrackerSimHitOutputCollections" type="StringVec">
    VBTrackerSimHits_SplittedByPolarAngle
    VETrackerSimHits_SplittedByPolarAngle
    IBTrackerSimHits_SplittedByPolarAngle
    IETrackerSimHits_SplittedByPolarAngle
    OBTrackerSimHits_SplittedByPolarAngle
    OETrackerSimHits_SplittedByPolarAngle
  </parameter>
  <!-- Name of the output filtered hit relation collection -->
  <parameter name="TrackerHitOutputRelations" type="StringVec">
    VBTrackerRel_SplittedByPolarAngle
    VETrackerRel_SplittedByPolarAngle
    IBTrackerRel_SplittedByPolarAngle
    IETrackerRel_SplittedByPolarAngle
    OBTrackerRel_SplittedByPolarAngle
    OETrackerRel_SplittedByPolarAngle
  </parameter>
  <!-- Lower limit and upper limit on the hit polar angle [deg] -->
    <parameter name="PolarAngleLowerLimit" type="double" value="50.0" />
    <parameter name="PolarAngleUpperLimit" type="double" value="130.0" />
    <parameter name="FillHistograms" type="bool" value="true" />
    <parameter name="Verbosity" type="string"> MESSAGE0 </parameter>
 </processor>
```
