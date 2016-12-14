# MarlinTrkProcessors

## A collection of Tracking Relelated Processors Based on MarlinTrk

Authors:

- S. Aplin, DESY
- F. Gaede, DESY
- Y.Voutsinas, DESY

## License and Copyright

Copyright (C) 2011-2016 DESY, Deutsches Elektronen Synchrotron

MarlinTrkProcessors is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see [http://www.gnu.org/licenses](http://www.gnu.org/licenses).



## Changelog - Release Notes 

###  v02-04-01 

- FPCCDSiliconTracking_MarlinTrk.cc:
	- patch from K.Fujii for macos (use std::isnan)

  
### v02-04
  

  F.Gaede
   - add parameter ParticleMass to RefitProcessor
      - can fit with different mass hypotheses
   

  Y. Voutsinas 
   - added DDSpacePointBuilder :  using DDRec rather than GEAR 
      -> to become the new default

   - improved CellsAutomatonMV


  R. Simoniello
   - cleaned up ExtrToTracker and DDCellsAutomatonMV ( uss DDRec rather than gear, ...)



  
### v02-03
  

  R. Simoniello
 
 - Add procedure to handle with spiral tracks. Fit direction and fit initialisation configurable by steering macro (default fit direction: backward)
  - improved speed, map of hits on each element introduced, find neighbours of element ID and extrapolate to them (to be improved to some problem in the geo datat structure)	

  F. Gaede
  - made compatible with c++11
	- removed -ansi -pedantic -Wno-long-long
	- fixed narrowing in initializer lists

  Y. Voutsinas
  - Minivector creation from 1D SIT hits added 


### v02-02
  
  R.Simoniello
   - RefitProcessor
     - for each subdetctor X, store nhits in fit in lcio::ILDDetID::X -1 and nhits 
       in track in lcio::ILDDetIDX -2. Structure needed for Pandora.
   - file added for correct comparison for float and double: fpcompare.h 
   -DDCellsAutomatonMV/ ExtrToTracker
     - criteria based on nhits ans QI for CA vertex tracks, remove the possibility of geo from gear file (option GEO1 in compilation) 
       for extrapolated tracks

  F.Gaede
  - DDPlanarDigiProcessor:
   - introduced parameter ForceHitsOntoSurface:
     - if true hits are projected onto the planar surface
       if they are putside the bounds
     - no check on distance is applied !
    - store the resolution from the length of the wafer for 1D hits
      ( could be treated as 2d measurements )

  - RefitProcessor:
     - no longer write incomplete tracks (w/o) track state at calo
     - rely on MarlinTrkTrack's implementation of initialise() if
       no already fitted TrackState is used



  
### v02-01
  

   R.Simoniello
   - updated Refitting/src/ExtrToTracker.cc
     - extrapolation also in endcap subdetectors + CA selecting longest (nhits) tracks

   F.Gaede
   - handle empty collections gracefully ( for lcgeo/ddsim ) 
     in Digitisers/src/PlanarDigiProcessor.cc

   - use abstract ISurface


  
### v02-00  
  

   - adapted to  DD4hep::Surfaces and DDKalTest


   - DDCellsAutomatonMV: adapted CA to DD4hep (R.Simoniello)

   - modified RefitProcessor (F.Gaede)
    - InitialTrackState: allows to select one of the track states (or use hits) 
    - FitDirection: allow to fit forwards or backwards (default)
    -> can use for CLIC SiBarrel tracker to refit TrackState AtIP 
       in forward direction

   - adapted to new MarlinTrk::Factory: no longer delete the IMarlinTrkSystem pointer
     at the end of the processing for all fitters

   - ExtrToTracker: add new processor that extracts an existing track to another 
     tracking subdetector (R.Simoniello)

   - added DDPlanarDigiProcessor which uses the DDRec::Surface
     for smearing the hits (no gear needed) (F.Gaede)

   - add option to use DDKalTest to all fitters (TrackSubsetProcessor, FullLDCTracking_MarlinTrk,
      TrackSubsetProcessor, RefittProcessor)    (F.Gaede)
      - set parameter "TrackSystemName" to "DDKalTest"

   - PlanarDigiProcessor: add new utility processor that allows to split a collection 
    of hits into several collections, based on the layerID
    of the hits (works for all lcio hit types) (F.Gaede)

   - ExtrToSIT: added processor for refitting and propagating track from TPC or VXD to SIT (Y.Voutsinas)

   - reverted changes for dEdx calculationin FullLDCMarlinTrk: moved to new processor (M.Kurata)


   - new processors: 
      Refitting/src/TruthTrackFinder.cc (D.Hynds)
         alternative track cheating
      Utils/src/ClicEfficiencyCalculator.cc (D.Hynds)
         tracking efficiency for CLIC (all silicon)

  
### v01-11  
  
 -  added processor for VXD tracking using a cellular automaton algorithm based on mini - vectors


  
### v01-10  
  
 -  added FPCCDTracking code by Tatsuya Mori

### v01-09-01  

  - updated calling attrinbutes of ced_hit_ID to new CED ### v01-09

  
### v01-09
  
  - FullLDCTracking_MarlinTrk
    - Use the Chi2 probability to remove badly fitted silicon tracks so that the don't pollute the TPC Si track merging.
    - Reject combinations of TPC and Silicon Tracks if more that a certain number of Silicon Hits, default 2, get rejected by the fit of the combined track.

  
### v01-08 
  

 - General 
   - SimplePlanarTestDigiProcessor renamed as PlanarDigiProcessor
   - Fixed signed unsigned warning.
   - Fix memory leaks.

 - PlanarDigiProcessor
   - Made processor parameters ResolutionU(V) vectors
   	 in order to allow for different point resolution values per layer
   	 ( needed for the VXD in the DBD)

 - FullLDCTracking_MarlinTrk
   - Added SET hits. 
   - Make sure hits flagged as not to be used in the fit are not included.
   - ForceTPCSegmentsMerging off by default
   - Ensure that trying to add a false hits, cannot derail the whole track creation. 
     Test number of Outliers using MaxAllowedPercentageOfOutliersForTrackCombination. 
   - When not using ForceTPCSegmentsMerging, make sure that TPC tracks which do not get 
   	 merged with Silicon Tracks have their hits set to setUsedInFit.


  
### v01-07 
  

  General:

   - Debug output has been made more consistent throughout.	

  FullLDCTracking_MarlinTrk
  
   - Use existing tracks parameters when refitting. 
   - Updated for new lcio TrackState copy constructor taking const reference. 
     
  TruthTracker	

   - Updated for new lcio TrackState copy constructor taking const reference. 
   - Corrected Helix orientation for pre-fit.


  RefitProcessor

   - Corrected Helix orientation for pre-fit.			
   - Protect against missing truth relations.

  
### v01-06-01 
  

  General:			

   - RIADA dependency added

  FullLDCTracking_MarlinTrk:

 - Increased debug output
 - Use fit from TPC to initialise fit when doing SiTrack TPCTrack comparisons.
 	Added debug output which can we removed at a later stage
 - Use Tracks from the TPC which are composed only of the innermost segment, which contains the fit needed 
 	for matching to the Silicon Tracks. The Silicon Track, TPC Track, and the segments are added to the final track, 
	which is then itself added to the Final Track collection. 
 - Added type bits for the final Track collection. 

 SiliconTracking_MarlinTrk:
    
- Make sure that float[6] is used for getPointOnCircle. 
- Added diagnostic histograms. 
- Check error return of fast fitter. 
- Use iopt 2 instead of 3, as newtonian part of fast fitter is not performing well.     
- Added type bits for the final Track collection. Increased debug for FTD tracking.
 
  SpacePointBuilder:
   
    - Added extra check for valid line intersection when creating space points. This fixes problems from very low pt tracks.


  
### v01-06 
   

  General:			 
  - Removed use of sort predicate. Use sorted list of std pairs instead.
  - Use _maxChi2PerHit in fits. 

  FullLDCTracking_MarlinTrk:	 
  - Corrected uninitialised covariance matrix in combine tracks method.
  - Fixed sorting hits in r2 during TPCTrack SiTracks combination.
  - Improved hit and geometry handling. 
  - Removed unused parameters. 
  - Fixed problem in allocating left over silicon hits, track fit was incorrect after adding a hit. 
  - Corrected point_res_rphi for FTD SpacePoints. 
  - Removed unused code.

  TruthTracker:	 
 - Added new mode to fit tracks iteratively, off by default. Allow the CED event display to be used during track creation.	 
 - For the iterative mode split Tracks now have the individual segments added to a separate collection. 
 - These segments are used to produce a combined track which uses the innermost and outermost segments to 
 - provided the track states at the IP, FirstHit, LastHit and Calo respectively. 
 - The combined track is added to the main track collection, the NDF and chi2 are set with the values from the innermost segment.
 - The NDF and chi2 values for the other segments can be retrieved from the tracks added to the composite track.  
	 
 -	Hits with no MCParticle are discarded.

  SiliconTracking_MarlinTrk:	 
  - Added protection against being overwhelmed by a very large number of hits within one sector. 
  - The default is set to 100. For now an ERROR message is printed and the sector is dropped.
  - Quality parameter added to output collection set to "Poor"
  - This will be removed when the problem is more properly addressed.
  - Make sure sort predicate for chi2 fulfils strict weak ordering.
  - Fixed problem where using return instead of continue caused loss of all fitted tracks when single track fit failed.
  - Improved hit and geometry handling. 
  - Instrumented with CEDEventDisplay to improve debugging. 
  - Removed unused parameters.

  SpacePointBuilder:		 
  - Corrected SpacePoint creation bug. Endpoint of rear strip swapped with endpoint of front strip. Added debug output for checking whether a strip is created from a real combination or ghost combination.

  TrackSubsetProcessor:		 
  - Add hit numbers and additional trackstates. Use MarlinTrk utilities. 

  
### v01-05 
   

 General:			 Added dependencies on KiTrack and KiTrackMarlin to MarlinTrkProcessors.
 
 Refitting:			 Added the TrackSubsetProcessor, to get the best subset from a redundant set of tracks. 
 				 TruthTracker restructured to use new convenience methods from MarlinTrk to produce 
				 cannonical trackstates. The diagonal elements of the Initial error matrix can now be set individually.
				 Fixed r_comparison sort problem.

 SiliconTracking_MarlinTrk	 Use convenience functions for track fitting from MarlinTrk. Provide steering of initial covariance matrix,  				 
& FullLDCTracking_MarlinTrk	 and max delta chi2 for adding hits. 
  				 Removed the creation of relations for MC <-> Tracks.

 SpacePointBuilder:		 Corrected bias in position by using vertex constraint	

  
### v01-04 
   

  SimplePlanarDigiProcessor:	 When calling getLadderNumber the "ladderNumber" was given as input instead of the "layerNumber" (typo).
  				 In getLadderNumber some hits were not assigned to the right layer due to rounding errors 
				 (they got INT_MAX which caused an exception when writing it to the bit field). Introduced a (hard coded) 1 nm rounding 
				 acceptance and all hits were found to be in their layer. In the case this should fail additionally put in that 
				 the layer with the closest distance to the layer centre is returned.
  				 Fixes from M. Killenberg 
   				
  SimplePlanarTestDigiProcessor: Updated to be able to work for all planar sub detectors, VXD, SIT, SET and FTD.
  				 Will very likely supperseed SimplePlanarDigiProcessor in the next release. 

  SiliconTracking_MarlinTrk &	 Added Strip Hits for FTD. 
  FullLDCTracking_MarlinTrk:     FTD restored to track finding and fitting.
 				 

  RefitProcessor:		 Updated to deal with recent addition of Strip Hits in FTD.

  TruthTracker:			 Allow Helix parameters for fit to be take from MCParticle using UseMCParticleParametersFotInitOfFit == true. 
  				 Allow steering of values used for the initial diagonal elements of the trackfit via InitialTrackErrors.

  
### v01-03 
   

  General:			Make use of MeasurementSurface classes from GEAR, as well as new UTIL::BitSet32 and 
   				UTIL::ILDTrkHitTypeBit from ILDConf (LCIO). 			

  SiliconTracking_MarlinTrk: 	Adapted to use composite Space-points composed from strip hits. 
     			    	Added new parameter NHitsChi2 which controls the maximal number of hits for which a track 
				with n hits is aways better than one with n-1 hits.
				General Cleanup and merger from changes from CLIC CDR. 				

  FullLDCTracking_MarlinTrk:	Adapted to use composite Space-points composed from strip hits.
	
  SimplePlanarDigiProcessor:	Adapted to have option to produce strip hits.

  TruthTracker:			Flexible collection input and able to construct track including composite Space-points.

  SpacePointBuilder:		New - Builds 3D composite Space-points from strip hits.



  
### v01-02 
   

SiliconTracking_MarlinTrk &	First version which is able to reconstruct LOI data.
 FullLDCTracking_MarlinTrk:     Fixed situation where, although more that 3 hits are excepted for fitting, 
 				less than 3 actually are included in the fit, by the kalman filter.

 SimpleCylinderDigiProcessor:	New - Uses simple gaussian smearing for cylinders.

 TruthTracker:			Fixed bug in SimTrackerHitSortPredicate.

 SimpleDiscDigiProcessor:	For LOI data the z coordinate is now set to that of the disk and the tolerance for dz has been increased to 10 microns.


  
### v01-01 
   

 SiliconTracking_MarlinTrk &   Modified to use getHitsInFit from MarlinTrk to get the correct TrackState for the first and last hits. 
 FullLDCTracking_MarlinTrk:    SiliconTracking has had the work around applied for the optimization problem of compare_r as applied 
                               to FullLDCTracking. 

  
### v01-00 
   

  SimplePlanarDigiProcessor:  
  			      Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 			      The plannar geometry should be either VXD, SIT or SET described using ZPlannarLayout
 			      The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 			      of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions.

  SimpleDiscDigiProcessor: 
  			      Produces a TrackerHit collection from SimTrackerHit collection. Intended for use with the FTD.
			      The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
			      of SimTrackerHits in (x,y) plane according to the specified point resolution. 
	
  SiliconTracking_MarlinTrk:  Processor performing stand-alone pattern recognition in the VTX, FTD and SIT.
  			      Adaption of the of the original SiliconTracking from A. Raspereza, to use MarlinTrk.

  FullLDCTracking_MarlinTrk:  Processor performing track finding procedure in the entire ILD detector by linking track segments found 
                              by the SiliconTracking module in the silicon detectors and by Clupatra in the TPC. 
			      Adaption of the of the original SiliconTracking from A. Raspereza, to use MarlinTrk.

  RefitProcessor:             Track Refitter processor for marlin. Refits an input track collection using MarlinTrk, 
                              producing a new collection of tracks. Primarily to demonstrate the use of MarlinTrk.

  TruthTracker:     	      Track creation based on MC truth. Uses Relations between TrackerHits and SimTrackerHits to produce a track 
  			      collection. It can also fit the tracks using MarlinTrk or fill the track collection with the helix parameters
                              of the MCParticle.			   
