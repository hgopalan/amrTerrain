# amrTerrain
This python utility provides a workflow to create input files for immersed forcing terrain simulation in AMR-Wind. It can also be adopted for ERF in future. 

## DEM Files 

There are two methods available to input the DEM files within the python scripts:

### Elevation Package 

THis is done automatically within the code based on the user input. The size of the domain is restricted to 100 kms. 

It is possible to use the he script stitchtif.py in src/backend to download multiple individual tiles. Once the tiles are downloaded, the tiles can be stitched using following commands: 

gdalbuildvrt terrain.vrt *.tif 

gdalwarp -of GTiff  terrain.tif


### User-provided Tiff File 

This option is useful for generating tiff files when the domain size is larger than 100kms. The tiff file can be generated from: https://portal.opentopography.org/raster?opentopoID=OTSRTM.082016.4326.1



### Sample Workflow 

The input files are created by running python src/backend/backendinterface.py pathtoyamlfile 

A sample yaml file is discussed below below. 

#### Solver type: Currently only amr-wind. May include ERF in future. 
solver: "amrWind"

Currently only amrWind is supported. ERF support will be added in the future. 

#### Directories
caseParent: "/Users/hgopalan/Documents/P101_AMR-Wind/Data/tempGUI"

caseFolder: "wfip2_25km_rans_domain_noprecursor_flowturning"

These variables provide the parent directory and name of case folder for writing the files. 

#### 4 Case types: precursor, terrain, terrain_noprecursor and turbine 
caseType: "terrain_noprecursor"

The recommended mode of execution is terrain_noprecursor. The other case types are added for backward compatibility with the initial immersed forcing developed and may be depcrecated in the future. 

#### Placeholders 
caseInitial: "amr"

domainType: "center"

These are currently placeholders and will be updated in the future. 

#### Center point around which to generate terrain 
centerLat: 45.63374

centerLon: -120.66047

This is the center point of the domain around which the tiff files is extracted and the domain is created. 

#### Distances measured in meters from center lat to chop the tiff file 
west: 12500

east: 12500

south: 12500

north: 12500

This is the span of the domain from the center latitutde and longitude. 

#### Smoothing from terrain to flat surface in meters 
westSlope: 2500

eastSlope: 2500

northSlope: 2500 

southSlope: 2500 

This is the distance within which the terrain is smooted to a flat surface. 

#### Span for the flat surface in meters 
westFlat: 2500

eastFlat: 2500

northFlat: 2500 

southFlat: 2500

This is the span of the flat surface which is used to create a fringe zone near inflow-outflow boundaries. 


#### Meshing 
#### Horizontal cell size in meters 
cellSize: 128

The cell size in the horizontal direction. 

#### Vertical cell size is cellSize/verticalAR 
verticalAR: 4

Aspect ratio to compute dz. AR of 1 increases computational cost and it is recommended to use 4. 

#### Size to write the terrain file in m 
terrainSize: 32

The size of the terrain when writing the terrain file to be read in AMR wind. It is recommended to use a terrainSize to be smaller than the finest horizontal grid size. 

#### Turbulence Model: RANS/LES 
turbulenceModel: "RANS"

AMR-Wind has support for several LES and RANS models. However, only the non-linear LES model and the one-equation K model for RANS has terrain support. 

#### Placeholder for initial wind if no met mast is specified 
windX: 10.0

windY: 0.0 

windZ: 0.0 

This value will be overwritten with met-mast values. The input is required in the AMR-Wind file and it will be overwritten. 

#### Physics values used for calculation 
refTemperature: 300.0

refRoughness: 0.1

refHeatflux: 0.0 

RefTemperature is required for the Bouissnesq terms. Currently, only a single roughness can be specified for terrain. Will be replaced with non-uniform roughness in future. Heat flux is specified at the surface for accounting for the stratification effects. This will be modified by non uniform heat flux in future. 

#### Monin-Obukhov length required for RANS/LES simulations 
molLength: 1e10

As an alterntive to specifying heat-flux, Monin Obukhov length can be specified. This is the commonly used method in commercial codes. When mollength is specified, refHeatflux is ignored. This method is still under development. 

#### If case type is turbine all turbines within US included 
turbineMarkType: "database"

The turbineMarkType can be used with the caseType:"turbine" to automatically include all the turbines in the domain within continental united states. This feature may be depcrecated in future and replaced with a new option. 

#### Writes a STL file of the terrain 
writeTerrain: true

Writes the terrain tiff file as a STL. 

#### Runs the 1-D rans solver to create initialconditions 
#### Mandatory for terrain_noprecursor to avoid solver divergence 
rans1D: true

It is recommended to run the 1-D solver for all cases. The solver creates a vertical profile generated from the single column RANS model. The output overwrites the uniform initial conditions in AMR-Wind. 


#### Errors between met-mast wind and the 1-D solver computed 
#### met-mast height winds 
allowedError: 0.05

The 1-D solver is run to convergence iteratively until the allowed wind speed error between the specified met-mast wind and the 1-D solver is less than the allowedError. 

#### Write terrain aligned sampling to compute speed-up maps 
writeTerrainSampling: true 

verticalLevels: [10,80,100,200]

Writes terrain-aligned speed-up maps at the specified vertical levels. The output needs to be postprocessed manually. A script may be added in future to do it. 

#### Refinment Regions 
refinementRegions: ["roi"]

refinementMinX: [-3000]

refinementMaxX: [3000]

refinementMinY: [-3000]

refinementMaxY: [3000]

heightAboveTerrain: [200]

refinementLevels: [1]

Specifies multiple refinment regions based on the user requirements. Each refinement region automatically creates a terrain adaptive mesh refinement within the bounds.  

#### Location of met-masts 
metMastNames: ["mast1","mast2","mast3"]

metMastLatLon: false 

metMastX: [-500,0,500]

metMastY: [0,0,0]

metMastRadius: [500,500,500]

metMastHeight: [100,100,100]

metMastRefinementLevel: [3,3,3]

metMastLineSampling: true 

A cylindrical refinement region is created around each mast. The metMastHeight is automatically computed above the terrain level. The metMastLineSampling writes out postprocessing results for each met-mast. Each met-mast refinement automatically creates a terrain adaptive mesh refinement within the bounds. 

#### Use one of the met-masts as driving wind 

metMastWind: [12,2]

metMastHeight: 100

For a specified wind direction, use the metMast closest to the domain boundary to provide the reference wind. 


#### AMR Refinement 

refineTerrain: false

refineLow: [-3000,-3000,300]

refineHigh: [3000,3000,600]

refineTerrainMaxLevel: 4

This set of variables provides additional AMR refinement to terrain in addition to refinementRegions and MetMastRegions. This can be useful when the upstream region of certain masts have high curvature and we want to capture it. 

#### Guess Geostrophic Wind 
initialUG: 17.1748 

initialVG: -4.10299

The 1-D solver generates the wind speed profile by computing the geostrophic wind required to generate the wind speed value at the met mast height. The process runs multiple iterations of the 1-D solver. A closer guess speeds up the computation. This value is optional and does not have to be included. 

#### Comment out the line below and change path to use user-specified tiff file 
useTiff: "/Users/hgopalan/Documents/P101_AMR-Wind/Data/tempGUI/output_SRTMGL1.tif"

The user-specified tif file should be used when the domain size is greater than 100 kms. 
