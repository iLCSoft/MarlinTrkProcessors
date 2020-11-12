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
  <!-- Configuration of the maximum allowed dU and dTheta between a pair of hits at the inner and outer layer -->
  <!-- 4 numbers per double-layer: <inner layer ID>  <outer layer ID>  <dU max [mm]>  <dTheta max [mrad]> -->
  <parameter name="DoubleLayerCuts" type="StringVec">
      0 1 1.0 0.6
      2 3 1.0 0.33
      4 5 1.0 0.27
      6 7 1.0 0.21
  </parameter>
  <!-- Whether to fill diagnostic histograms about affected hits -->
  <parameter name="FillHistograms" type="bool" value="true" />
  <!-- Verbosity level ("DEBUG0-9,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT") -->
  <parameter name="Verbosity" type="string"> MESSAGE </parameter>
</processor>
```
