'''
!----------------------------------------------------------------!
AMR - Wind back-end for cases with terrain                       !
!----------------------------------------------------------------!
'''

import yaml
from pathlib import Path
from terrain import SRTM

class amrBackend():
    def __init__(self,yamlFile):
        self.yamlFilePath = Path(yamlFile)
        self.yamlFile=yaml.safe_load(self.yamlFilePath.open())
        self.createCase()

    def createCase(self):
        self.caseParent=self.yamlFile['caseParent']
        self.caseName=self.yamlFile['caseFolder']
        self.caseType=self.yamlFile['caseType']
        self.caseInitial=self.yamlFile['caseInitial']
        caseDir=Path(self.caseParent,self.caseName)
        self.caseDir=caseDir.as_posix()
        caseDir.mkdir(parents=True,exist_ok=True)
        if(self.caseType=="terrain"):
            casePrecursor=Path(self.caseParent,self.caseName,"precursor")
            self.casePrecursor=casePrecursor.as_posix()
            casePrecursor.mkdir(parents=True,exist_ok=True)
            caseTerrain=Path(self.caseParent,self.caseName,"terrain")
            self.caseTerrain=caseTerrain.as_posix()
            caseTerrain.mkdir(parents=True,exist_ok=True)
        self.caseDomainType=self.yamlFile['domainType']
        if(self.caseDomainType=="corners"):
            pass
            self.caseNorth=self.yamlFile['north']
            self.caseSouth=self.yamlFile['south']
            self.caseEast=self.yamlFile['east']
            self.caseWest=self.yamlFile['west']
        else:
            self.caseCenterLat=self.yamlFile["centerLat"]
            self.caseCenterLon=self.yamlFile["centerLon"]
            self.caseNorth=self.caseCenterLat+self.yamlFile['north']
            self.caseSouth=self.caseCenterLat-self.yamlFile['south']
            self.caseEast=self.caseCenterLon+self.yamlFile['east']
            self.caseWest=self.caseCenterLon-self.yamlFile['west']

    def createDomain(self):
        bounds = self.caseWest, self.caseSouth, self.caseEast, self.caseNorth 
        dx = dy = 90.
        product = 'SRTM1' 
        srtm_output = Path(self.caseParent,self.caseName,"terrain.tif")
        srtm = SRTM(bounds,fpath=srtm_output.as_posix(),product=product)
        srtm.download()
        x1,x2,x3 = srtm.to_terrain(dx,dy)
        cornerCut=12
        self.caseX1=x1[cornerCut:x1.shape[0]-cornerCut,cornerCut:x1.shape[1]-cornerCut]
        self.caseX2=x2[cornerCut:x2.shape[0]-cornerCut,cornerCut:x2.shape[1]-cornerCut]
        self.caseX3=x3[cornerCut:x3.shape[0]-cornerCut,cornerCut:x3.shape[1]-cornerCut]
        self.caseLatList=0*self.caseX1
        self.caseLonList=0*self.caseX2
        for i in range(0,self.caseX1.shape[0]):
            for j in range(0,self.caseX1.shape[1]):
                self.caseLatList[i,j],self.caseLonList[i,j]=srtm.to_latlon(self.caseX1[i,j],self.caseX2[i,j])
        # import pyvista as pv
        # import numpy as np 
        # x1=self.caseX1.flatten(order='F')
        # x2=self.caseX2.flatten(order='F')
        # x3=self.caseX3.flatten(order='F')
        # data=np.column_stack([x1,x2,x3])
        # print(data.shape)
        # mesh=pv.PolyData(data)
        # mesh['elevation']=data[:,2]
        # surf = mesh.delaunay_2d()
        # pl = pv.Plotter()
        # pl.add_mesh(surf)
        # pl.set_scale(zscale=5)
        # pl.view_xy()
        # pl.show_axes()
        # pl.show()
    
    def createAMRFiles(self):
        self.amrPrecursorFile=Path(self.caseParent,self.caseName,"precursor","precursor.inp").open("w")
        self.amrPrecursorFile.write("# Generating the precursor file\n")
        self.amrTerrainFile=Path(self.caseParent,self.caseName,"terrain","terrain.inp").open("w")
        self.amrTerrainFile.write("# Generating the terrain file\n")

    def closeAMRFiles(self):
        self.amrPrecursorFile.close()
        self.amrTerrainFile.close()





from sys import argv 
amrRef=amrBackend(argv[1])
amrRef.createDomain()
amrRef.createAMRFiles()