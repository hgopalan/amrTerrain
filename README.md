# amrGUI
A python backend and front-end (hopefully) to setup cases in AMR-Wind with Terrain and generate postprocessing to obtain acceleration maps, etc. 

The code is executed by running python src/backend/backendinterface.py pathtoyamlfile 

Check the tests/yaml_files/wfip2_25km_rans_domain_noprecursor.yaml for detailed information of the various inputs. 

The default mode of operation for the code is to download the SRTM 1 degree data and construct the terrain. However, this can be avoided by specifying a user-defined tiff file. It is also possible to use the script stitchtif.py in src/backend. This will download the individual tiles. Once the tiles are saved run the following commands 

gdalbuildvrt terrain.vrt *.tif 

gdalwarp -of GTiff  terrain.tif

This terrain.tif file can be used inside the yaml file to avoid the download step for areas > 100 km as the elevation package will throw an error message. 
