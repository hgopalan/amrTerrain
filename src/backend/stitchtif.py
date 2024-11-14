# This is experimental code. Need some cleaning in future. 
from terrain import SRTM


centerLat= 45.63374
centerLon=-120.66047
dx=0.5
latmin=centerLat-4.0
latmax=centerLat+4.0
longmin=centerLon-4
longmax=centerLon+4
filename=1
product = 'SRTM1'
while (latmin<=latmax):
    while (longmin<=longmax):
        srtm_bounds = west, south, east, north = (longmin,latmin,longmin+0.5,latmin+0.5)
        srtm_output="/Users/hgopalan/Documents/P101_AMR-Wind/Codes/amrGUI/"+str(filename)+".tif"
        srtm = SRTM(srtm_bounds, fpath=srtm_output, product=product)
        srtm.download()
        filename=filename+1
        longmin=longmin+0.5
        print(latmin,longmin)
    latmin=latmin+0.5
    longmin=centerLon-4
