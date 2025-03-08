'''
!----------------------------------------------------------------!
AMR - Wind back-end for cases with terrain                       !
!----------------------------------------------------------------!
'''

import yaml
from pathlib import Path
from terrain import SRTM
import numpy as np 
import warnings 
import utm
from scipy.interpolate import NearestNDInterpolator


class amrBackend():
    def __init__(self,yamlFile):
        self.yamlFilePath = Path(yamlFile)
        self.yamlFile=yaml.safe_load(self.yamlFilePath.open())
        self.caseCellSize=96  
        self.caseverticalAR=4
        self.turbulence_model='RANS'
        self.case_end_time=7200
        self.plotOutput=1800
        self.restartOutput=1800
        # Rearranging YAML file reading 
        self.readVariables()
        self.createCase()
        self.setCaseType()
    
    def readVariables(self):
        # 1 - Folders 
        self.caseParent=self.yamlFile['caseParent']
        self.caseName=self.yamlFile['caseFolder']
        # 2 - Case Type 
        try:
            self.caseType=self.yamlFile['caseType']
        except:
            self.caseType="terrain_noprecursor"
        # 3 - Farm (even for precursor)
        self.caseCenterLat=self.yamlFile["centerLat"]
        self.caseCenterLon=self.yamlFile["centerLon"]
        # Method 1 - Assuming a circular farm but amr-wind uses rectangular grid so using inscribed circle 
        try:
            farmRadius=self.yamlFile['farmRadius']
            self.caseNorth=farmRadius
            self.caseSouth=farmRadius
            self.caseEast=farmRadius
            self.caseWest=farmRadius
        except:
            # Explicitly setting size of farms. Useful when the farm is longer in one direction 
            self.caseNorth=self.yamlFile['north']
            self.caseSouth=self.yamlFile['south']
            self.caseEast=self.yamlFile['east']
            self.caseWest=self.yamlFile['west']  
        # Currently a dummy variable. Will be used with WRF data 
        try:
            self.refHeight=self.yamlFile["refHeight"]
        except:
            self.refHeight=2184
        # 4 - Define Fringe Regions 
        # Slopes from the terrain to flat surface at boundaries. 
        # Default value of 5% is used if not specified 
        try:
            self.caseNorthSlope=self.yamlFile['northSlope']
        except:
            self.caseNorthSlope=0.05*(self.caseNorth + self.caseSouth)
        try:
            self.caseSouthSlope=self.yamlFile['southSlope']
        except:
            self.caseSouthSlope=0.05*(self.caseNorth + self.caseSouth)
        try:
            self.caseEastSlope=self.yamlFile['eastSlope']
        except:
            self.caseEastSlope=0.05*(self.caseEast+self.caseWest)
        try:
            self.caseWestSlope=self.yamlFile['westSlope']
        except:
            self.caseWestSlope=0.05*(self.caseEast+self.caseWest)
        # Length of the flat region to apply the forcing and fringe boundary conditions 
        # Default value of 5% is used if not specified 
        try:
            self.caseNorthFlat=self.yamlFile['northFlat']
        except:
            self.caseNorthFlat=0.05*(self.caseNorth + self.caseSouth)
        try:
            self.caseSouthFlat=self.yamlFile['southFlat']
        except:
            self.caseSouthFlat=0.05*(self.caseNorth + self.caseSouth)
        try:
            self.caseEastFlat=self.yamlFile['eastFlat']
        except:
            self.caseEastFlat=0.05*(self.caseEast+self.caseWest)
        try:
            self.caseWestFlat=self.yamlFile['westFlat']
        except:
            self.caseWestFlat=0.05*(self.caseEast+self.caseWest)
        # 5 Define Mesh Sizing 
        try:
            self.caseCellSize=self.yamlFile['cellSize']
        except:
            self.caseCellSize=96.0
        try:
            self.caseverticalAR=self.yamlFile['verticalAR']
        except:
            self.caseverticalAR=4
        # 6 Define Turbulence Model 
        try:
            self.turbulence_model=self.yamlFile['turbulenceModel']
        except:
            self.turbulence_model='RANS'
        try:
            self.rans_1d=self.yamlFile["rans1D"]
            # 1-D solver is mandatory for the RANS model while optional for LES 
            if(self.turbulence_model=='RANS'):
                self.rans_1d=True 
        except:
            self.rans_1d=True
        # 7 Define Physics 
        try:
            self.refTemperature=self.yamlFile["refTemperature"]
        except:
            self.refTemperature=300.0
        try:
            self.refRoughness=float(self.yamlFile["refRoughness"])
        except:
            self.refRoughness=0.1
        # When we are sweeping roughness changes 
        try:
            wind=self.yamlFile["metMastWind"]
            self.metMastHeight=self.yamlFile["metMastHeight"]
            self.metMastWind=[wind[0],wind[1]]
            #print(self.metMastHeight)
            #print(self.metMastWind)
        except:
            # Set defaults if not specified 
            self.metMastWind=[10,0]
            self.metMastHeight=[100] 
        # 8 Read if we are doing AEP 
        try:
            self.sweep_angle_increment=self.yamlFile["sweepAngle"]
        except:
            self.setSweep=False 
        else:
            self.setSweep=True 
        # 9 Time Stepping 
        try:
            self.timeMethod=self.yamlFile['timeMethod']
        except:
            self.timeMethod="end_time"
        if(self.timeMethod=="step"):
            self.timeSteps=self.yamlFile["numOfSteps"]
        else:
            try:
                self.case_end_time=self.yamlFile["endTime"]
            except:
                pass
        try:
            self.plotOutput=self.yamlFile['plotOutput']
        except:
            pass
        # 10 - Special FF handling 
        try:
            self.fastBoxes=self.yamlFile['fastBoxes']      
        except:
            self.fastBoxes=False
        # Met Mast Driving 
        try:
            self.metmast_horizontal_radius=self.yamlFile["metmast_horizontal_radius"]
        except:
            self.metmast_horizontal_radius=500.0
        try:
            self.metmast_vertical_radius=self.yamlFile["metmast_vertical_radius"]
        except:
            self.metmast_vertical_radius=5.0
        try:
            self.metmast_damping_radius=self.yamlFile["metmast_damping_radius"]
        except:
            self.metmast_damping_radius=100.0

    def createCase(self):
        caseDir=Path(self.caseParent,self.caseName)
        self.caseDir=caseDir.as_posix()
        caseDir.mkdir(parents=True,exist_ok=True)

    def setCaseType(self):
        if(self.caseType=="terrain_noprecursor"):
            pass
        else:
            casePrecursor=Path(self.caseParent,self.caseName,"precursor")
            self.casePrecursor=casePrecursor.as_posix()
            casePrecursor.mkdir(parents=True,exist_ok=True)
        if(self.caseType=="terrain" or self.caseType=="terrain_noprecursor"):
            caseTerrain=Path(self.caseParent,self.caseName,"terrain")
            self.caseTerrain=caseTerrain.as_posix()
            caseTerrain.mkdir(parents=True,exist_ok=True)
        elif(self.caseType=="terrainTurbine"):
            caseTerrain=Path(self.caseParent,self.caseName,"terrainTurbine")
            self.caseTerrain=caseTerrain.as_posix()
            caseTerrain.mkdir(parents=True,exist_ok=True)  
    

    def createDomain(self):
        try:
            self.write_stl=self.yamlFile["writeTerrain"]
        except:
            self.write_stl=False
        try:
            self.terrainSTL=self.yamlFile["terrainSTL"]
        except:
            import SRTM_to_STL_example as converter
            try:
                self.usetiff=self.yamlFile['useTiff']
            except:
                self.usetiff=' '
                westlon=-0.5
                eastlon=0.5
                southlat=-0.5
                northlat=0.5
            else:
                import rasterio
                dataset = rasterio.open(self.usetiff)
                westlon=dataset.bounds[0]-self.caseCenterLon
                eastlon=dataset.bounds[2]-self.caseCenterLon
                southlat=dataset.bounds[1]-self.caseCenterLat
                northlat=dataset.bounds[3]-self.caseCenterLat
            try:
                self.xref,self.yref,self.zRef,self.srtm,self.zone_number=converter.SRTM_Converter(Path(self.caseParent,self.caseName).as_posix(),self.caseCenterLat,self.caseCenterLon,self.refHeight, \
                                                            self.caseWest,self.caseEast,self.caseSouth,self.caseNorth, \
                                                            self.caseWestSlope,self.caseEastSlope,self.caseSouthSlope,self.caseNorthSlope, \
                                                            self.caseWestFlat,self.caseEastFlat,self.caseSouthFlat,self.caseNorthFlat, self.usetiff, \
                                                                self.write_stl,westlon,eastlon,southlat,northlat)
                                                                #-3,3,-1.5,1.5)
            except:
                print("Cannot connect to internet to download file")
                exit(-1)
        else:
            warnings.warn("Turbine or Met Mast Locations should be in STL coordinates. Not Lat/Lon")
            self.xref=0
            self.yref=0
            self.zRef=0
        stlFile=Path(self.caseParent,self.caseName,"terrain.vtk").as_posix()
        import pyvista as pv 
        mesh=pv.read(stlFile)
        x1=mesh.points[:,0]
        x2=mesh.points[:,1]
        x3=mesh.points[:,2]
        for i in range(0,len(x3)):
            # Water measurements can go negative 
            x3[i]=max(x3[i],0.0)
        self.terrainX1=x1[:]
        self.terrainX2=x2[:]
        self.terrainX3=x3[:]
        if(not self.write_stl):
            Path(self.caseParent,self.caseName,"terrain.vtk").unlink()
        try:
            roughnessFile=self.yamlFile["roughnessFile"]
        except:
            pass
        else:
            self.makeRoughness()
    
    def createAMRFiles(self):
        if(self.caseType=="terrain_noprecursor"):
            pass
        else:
            self.amrPrecursorFile=Path(self.caseParent,self.caseName,"precursor","precursor.inp").open("w")
            self.amrPrecursorFile.write("# Generating the precursor file\n")
            self.createPrecursorFiles()
        # Sweep angles for multiple case creation 
        if(self.setSweep):
            import shutil
            reference_angle=[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340]
            true_angle=[270,250,230,210,190,170,150,130,110,90,70,50,30,10,350,330,310,290]
            angle=0
            while(angle<360):
                print("Sweeping:",angle)
                self.sweep_angle=angle 
                self.amrTerrainFile=Path(self.caseParent,self.caseName,"terrain","terrain.inp").open("w")
                self.amrTerrainFile.write("# Generating the terrain file\n")
                self.createTerrainFiles("terrain")
                self.closeAMRFiles() 
                index=np.where(np.array(reference_angle)==angle)
                #print(true_angle[index[0][0]])
                destination=Path(self.caseParent,self.caseName,"terrain_"+str(int(true_angle[index[0][0]])))
                try:
                    shutil.move(Path(self.caseParent,self.caseName,"terrain").as_posix(),destination)
                except:
                    shutil.rmtree(destination)
                    shutil.move(Path(self.caseParent,self.caseName,"terrain").as_posix(),destination)
                caseTerrain=Path(self.caseParent,self.caseName,"terrain")
                angle=angle+self.sweep_angle_increment
                # 2025-03-06 Need to update roughness for angle 

                if(angle<360):
                    self.caseTerrain=caseTerrain.as_posix()
                    caseTerrain.mkdir(parents=True,exist_ok=True)
                    self.yamlFile=yaml.safe_load(self.yamlFilePath.open())
        else:
            if(self.caseType=='terrain' or self.caseType=="terrain_noprecursor"):
                self.amrTerrainFile=Path(self.caseParent,self.caseName,"terrain","terrain.inp").open("w")
                self.amrTerrainFile.write("# Generating the terrain file\n")
                self.createTerrainFiles("terrain")
                self.closeAMRFiles()    
            elif(self.caseType=='terrainTurbine'):
                self.amrTerrainFile=Path(self.caseParent,self.caseName,"terrainTurbine","terrainTurbine.inp").open("w")
                self.amrTerrainFile.write("# Generating the terrain file\n")
                self.createTerrainFiles("terrainTurbine")
                self.closeAMRFiles()
                # Turbine files are not properly written
                self.createTurbineFiles()
                self.plotTurbines()

    def createPrecursorFiles(self):
        #self.refHeatFlux=0.0
        if(self.caseType=="terrain_noprecursor"):
            return 
        #print("Creating precursor")
        self.createAMRGeometry(self.amrPrecursorFile,1)
        # Variables missing
        if(self.rans_1d):
            self.createAMR1dSolver()
        self.createAMRGrid(self.amrPrecursorFile)
        self.createAMRTime(self.amrPrecursorFile)
        self.createSolverInfo(self.amrPrecursorFile)
        self.createAMRTransport(self.amrPrecursorFile)
        self.createAMRTurbulence(self.amrPrecursorFile)
        self.createAMRABLData(self.amrPrecursorFile,0,1) 
        if(self.turbulence_model=="RANS"):
            self.createAMRSourceTerm(self.amrPrecursorFile,1)    
        else:
            self.createAMRSourceTerm(self.amrPrecursorFile) 
        self.createAMRBC(self.amrPrecursorFile)
        self.createAMRTolerance(self.amrPrecursorFile)
        #print(" Done creating precursor")

    def createTerrainFiles(self,folder):
        #print("Creating Terrain Files")
        self.createAMRGeometry(self.amrTerrainFile,-1)
        if(self.rans_1d and self.caseType=="terrain_noprecursor"):
            self.createAMR1dSolver()
        self.createAMRGrid(self.amrTerrainFile)
        self.createAMRTime(self.amrTerrainFile,1)
        if(self.caseType=='terrainTurbine'):
            self.createSolverInfo(self.amrTerrainFile,1,1)
        else:
            self.createSolverInfo(self.amrTerrainFile,1)
        self.createAMRTransport(self.amrTerrainFile)
        self.createAMRTurbulence(self.amrTerrainFile)
        self.createAMRABLData(self.amrTerrainFile,1,0)  
        if(self.caseType=='terrainTurbine'):
            self.createAMRSourceTerm(self.amrTerrainFile,-1,1,1) 
        else:
            self.createAMRSourceTerm(self.amrTerrainFile,-1,1)    
        self.createAMRBC(self.amrTerrainFile,1)
        self.createAMRTolerance(self.amrTerrainFile,1)
        self.writeTerrainData(folder)
        self.writeRestart(self.amrTerrainFile)
        if(self.rans_1d==True):
            self.fillrans1dinfo(self.amrTerrainFile,1)
        if(self.caseType=="terrain_noprecursor"):
            pass
        else:
            self.createAMRPrecursorSampling(self.amrPrecursorFile)
        # Creating refinement regions 
        self.writeRefinementRegions(self.amrTerrainFile)
        self.writeAccelerationMaps(self.amrTerrainFile)
        # Terrain Monitoring and Refining 
        #self.metMastRefinement(self.amrTerrainFile)
        #self.metMastMonitoring(self.amrTerrainFile)

    def createAMRGeometry(self,target,periodic=-1):
        target.write("# Geometry\n")
        minX=np.amin(self.terrainX1)
        minY=np.amin(self.terrainX2)
        minZ=np.amin(self.terrainX3)
        maxX=np.amax(self.terrainX1)
        maxY=np.amax(self.terrainX2)
        self.terrainZMax=np.amax(self.terrainX3)
        ranges=np.arange(512,8192,64)
        idx = np.argmin(np.abs(ranges - self.terrainZMax))
        if(ranges[idx]>self.terrainZMax):
            self.terrainZMax=ranges[idx]
        else:
            self.terrainZMax=ranges[idx+1]
        self.ABLHeight=max(self.terrainZMax,2048)
        self.RDLHeight=max(self.terrainZMax,2048)
        self.maxZ=self.terrainZMax+self.ABLHeight+self.RDLHeight
        #print(self.terrainZMax,self.ABLHeight,self.RDLHeight,self.maxZ)
        target.write("%-50s = %g %g %g \n"%("geometry.prob_lo",minX,minY,minZ))
        target.write("%-50s = %g %g %g \n"%("geometry.prob_hi",maxX,maxY,self.maxZ))
        if(periodic==1):
            target.write("%-50s = 1 1 0\n"%("geometry.is_periodic"))
        else:
            target.write("%-50s = 0 0 0\n"%("geometry.is_periodic"))

    def createAMRGrid(self,target):
        nx=int((np.amax(self.terrainX1)-np.amin(self.terrainX1))/self.caseCellSize)
        while (nx%8 !=0):
            nx=nx+1
        ny=int((np.amax(self.terrainX2)-np.amin(self.terrainX2))/self.caseCellSize)
        while (ny%8 !=0):
            ny=ny+1
        nz=self.caseverticalAR*int(self.maxZ/self.caseCellSize)
        while (nz%8 !=0):
            nz=nz+1
        #print("dx,dz",self.caseCellSize,self.maxZ/nz)
        target.write("# Grid \n")
        target.write("%-50s = %g %g %g\n"%("amr.n_cell",nx,ny,nz))
        #target.write("%-50s = 0\n"%("amr.max_level"))

    def createAMRTime(self,target,blanking=-1):
        if(self.caseType=="terrain_noprecursor"):
            if(self.turbulence_model=="RANS"):
                self.case_end_time=600
                self.plotOutput=600
                self.restartOutput=-1
            else:
                self.case_end_time=7200
                self.plotOutput=1800
                self.restartOutput=-1
        else:
            if(self.turbulence_model=="RANS"):
                if(blanking==1):
                    self.case_end_time=2400
                else:
                    self.case_end_time=2405
                self.plotOutput=600
                self.restartOutput=600
            else:
                if(blanking==1):
                    self.case_end_time=7200
                else:
                    self.case_end_time=7205
                self.plotOutput=1800
                self.restartOutput=1800
        if(self.timeMethod=="step"):
            target.write("%-50s = -1\n"%("time.stop_time"))
            target.write("%-50s = %g\n"%("time.max_step",self.timeSteps))
        else:
            target.write("%-50s = %g\n"%("time.stop_time",self.case_end_time))
        target.write("%-50s = 1.0\n"%("time.initial_dt"))
        try: 
            if(target==self.amrPrecursorFile):
                target.write("%-50s = 1\n"%("time.fixed_dt"))
            else:
                target.write("%-50s = -1\n"%("time.fixed_dt"))
        except:
            target.write("%-50s = -1\n"%("time.fixed_dt"))
        if(self.turbulence_model=="RANS"):
            target.write("%-50s = 0.5\n"%("time.cfl"))  
        else:
            target.write("%-50s = 0.9\n"%("time.cfl"))
        if(self.timeMethod=="step"):
            target.write('%-50s = %g\n'%("time.plot_interval",self.plotOutput))
            target.write("%-50s = %g\n"%("time.checkpoint_interval",self.restartOutput))
        else:
            target.write('%-50s = %g\n'%("time.plot_time_interval",self.plotOutput))
            target.write("%-50s = %g\n"%("time.checkpoint_time_interval",self.restartOutput))
        # Writing io 
        target.write("%-50s = false \n"%("io.output_default_variables"))
        if(self.turbulence_model=="RANS"):
            target.write("%-50s = velocity temperature mu_turb tke terrainz0 \n"%("io.outputs"))
        else:
            target.write("%-50s = velocity temperature mu_turb terrainz0 \n"%("io.outputs"))
        if(blanking==1):
            target.write("%-50s = terrain_blank terrain_drag \n"%("io.int_outputs"))

    def createSolverInfo(self,target,terrain=-1,turbine=-1):
        self.caseWindspeedX=self.metMastWind[0]
        self.caseWindspeedY=self.metMastWind[1]
        self.caseWindspeedZ=0.0
        try:
            self.caseWindspeedX=self.geostropicX
        except:
            pass 
        try:
            self.caseWindspeedY=self.geostropicY
        except:
            pass
        target.write("# incflo \n")
        if(terrain==1 and turbine==1):
            if(self.fastBoxes):
                target.write("%-50s = ABL TerrainDrag \n"%("incflo.physics"))
            else:
                target.write("%-50s = ABL TerrainDrag Actuator\n"%("incflo.physics"))
        elif(terrain==1):
            target.write("%-50s = ABL TerrainDrag\n"%("incflo.physics"))            
        elif(turbine==1):
            target.write("%-50s = ABL Actuator\n"%("incflo.physics"))
        else:
            target.write("%-50s = ABL\n"%("incflo.physics"))
        target.write("%-50s = 1.225\n"%("incflo.density "))
        target.write("%-50s = 0.  0. -9.81 \n"%("incflo.gravity "))
        target.write("%-50s = %g %g %g \n"%("incflo.velocity",self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        target.write("%-50s = 0\n"%("incflo.verbose"))
        target.write("%-50s = 8\n"%("incflo.initial_iterations"))
        target.write("%-50s = true\n"%("incflo.do_initial_proj"))
        target.write("%-50s = true\n"%("incflo.constant_density"))
        target.write("%-50s = true\n"%("incflo.use_godunov "))
        if(self.turbulence_model=="RANS"):
            target.write('%-50s = "ppm"\n'%("incflo.godunov_type"))
        else:
            target.write('%-50s = "weno_z"\n'%("incflo.godunov_type "))
        target.write("%-50s = 2\n"%('incflo.diffusion_type'))

    def createAMRTransport(self,target):
        target.write("# transport equation parameters \n")
        target.write("%-50s = ConstTransport\n"%("transport.model"))
        target.write("%-50s = 1e-5\n"%("transport.viscosity"))
        target.write("%-50s = 0.7\n"%("transport.laminar_prandtl"))
        target.write("%-50s = 0.333\n"%("transport.turbulent_prandtl "))     

    def createAMRTurbulence(self,target):
        # Default option is to do LES 
        # RANS capability added if requested 
        if(self.turbulence_model=="RANS"):
            target.write("# turbulence equation parameters\n")
            target.write("%-50s = KLAxell\n"%("turbulence.model"))
            target.write("%-50s = KransAxell\n"%("TKE.source_terms"))
        else:
            target.write("# turbulence equation parameters \n")
            target.write("%-50s = Kosovic\n"%("turbulence.model"))
            target.write("%-50s = -1e30\n"%("Kosovic.refMOL"))

    def createAMRABLData(self,target,iomode=-1,fluctuations=1):
        #self.refHeatFlux=0.0
        target.write("# Atmospheric boundary layer\n")
        if(fluctuations==1 and (not self.turbulence_model=='RANS')):
            UPeriod=int((np.amax(self.terrainX1)-np.amin(self.terrainX1))/200)
            VPeriod=int((np.amax(self.terrainX1)-np.amin(self.terrainX1))/200)
            target.write("%-50s = %g\n"%("ABL.Uperiods",UPeriod))
            target.write("%-50s = %g\n"%("ABL.Vperiods",VPeriod))
            target.write("%-50s = 50.0\n"%("ABL.cutoff_height"))
            target.write("%-50s = 1.0\n"%("ABL.deltaU"))
            target.write("%-50s = 1.0\n"%("ABL.deltaV")) 
            target.write("%-50s = 50.0\n"%("ABL.perturb_ref_height"))
            target.write("%-50s = true\n"%("ABL.perturb_velocity "))
            target.write("%-50s = false\n"%("ABL.perturb_temperature "))
        else:
            target.write("%-50s = false\n"%("ABL.perturb_velocity"))
            target.write("%-50s = false\n"%("ABL.perturb_temperature"))
        target.write("%-50s = .41\n"%("ABL.kappa"))
        target.write("%-50s = 2\n"%("ABL.normal_direction"))
        target.write("%-50s = %g\n"%("ABL.reference_temperature",self.refTemperature))
        target.write("%-50s = netcdf\n"%("ABL.stats_output_format "))
        # Update roughness for sweep angles 
        try:
            int(self.sweep_angle)
        except:
            xdir=self.xref
            ydir=self.yref
        else:
            farmRadius=self.yamlFile['farmRadius']
            if(self.sweep_angle<=90):
                xdir= self.xref-farmRadius* np.cos(self.sweep_angle*np.pi/180)
                ydir= self.yref-farmRadius* np.sin(self.sweep_angle*np.pi/180)
            elif(self.sweep_angle>90 and self.sweep_angle<=180):
                xdir= self.xref+farmRadius* np.cos((180-self.sweep_angle)*np.pi/180)
                ydir= self.yref- farmRadius* np.sin((180-self.sweep_angle)*np.pi/180)
            elif(self.sweep_angle>180 and self.sweep_angle<=270):
                xdir= self.xref+ farmRadius* np.cos((self.sweep_angle-180)*np.pi/180)
                ydir= self.yref+ farmRadius* np.sin((self.sweep_angle-180)*np.pi/180)
            else:
                xdir= self.xref- farmRadius* np.cos((360-self.sweep_angle)*np.pi/180)
                ydir= self.yref+ farmRadius* np.sin((360-self.sweep_angle)*np.pi/180)   
        self.refRoughness=self.roughness_interp(xdir-self.xref,ydir-self.yref)
        target.write("%-50s = %g\n"%("ABL.surface_roughness_z0 ",self.refRoughness))
        # Write Heights 
        if(not self.rans_1d):
            inversionHeight=round(self.terrainZMax,-3)+1000
            inversionLayerThickness=round(self.terrainZMax,-3)+1000+100
            lapseRate=0.003 
            target.write("%-50s = %g %g %g %g %g \n"%("ABL.temperature_heights",0.0,round(self.terrainZMax,-3),inversionHeight,inversionLayerThickness,self.maxZ)) 
            TRef=300
            target.write("%-50s = %g %g %g %g %g "%("ABL.temperature_values ",TRef,TRef,TRef,TRef+5,TRef+5+lapseRate*(self.maxZ-inversionLayerThickness)))    
            target.write("\n")
        target.write("%-50s = local\n"%("ABL.wall_shear_stress_type"))
        try:
            mol_length=self.yamlFile["molLength"]
        except:
            mol_length=-1e30
        try:
            target.write("%-50s = %g\n"%("ABL.surface_temp_flux",self.Qh))
        except:
            target.write("%-50s = %g\n"%("ABL.surface_temp_flux",0.0))
        target.write("%-50s = %g\n"%("ABL.mol_length",float(mol_length)))
        target.write("%-50s = %s\n"%("ABL.wall_het_model",'"mol"'))
        if(self.caseType=="terrain_noprecursor"):
            target.write("%-50s = true\n"%("ABL.horizontal_sponge_temp"))
            target.write("%-50s = true\n"%("ABL.horizontal_sponge_tke"))
        if(iomode==0):
            target.write('%-50s = "bndry_files"\n'%("ABL.bndry_file "))
            target.write("%-50s = 1\n"%("ABL.bndry_write_frequency"))
            target.write("%-50s = 0\n"%("ABL.bndry_io_mode"))
            target.write("%-50s = xlo xhi ylo yhi \n"%("ABL.bndry_planes")) 
            # if(self.caseWindspeedX>=0 and self.caseWindspeedY>=0):
            #     target.write("%-50s = xlo ylo \n"%("ABL.bndry_planes"))     
            # elif(self.caseWindspeedX>=0 and self.caseWindspeedY<0):
            #     target.write("%-50s = xlo yhi \n"%("ABL.bndry_planes"))      
            # if(self.caseWindspeedX<0 and self.caseWindspeedY>=0):
            #     target.write("%-50s = xhi ylo \n"%("ABL.bndry_planes"))     
            # elif(self.caseWindspeedX<0 and self.caseWindspeedY<0):
            #     target.write("%-50s = xhi yhi \n"%("ABL.bndry_planes "))    
            startTime=600
            target.write("%-50s = %g\n"%("ABL.bndry_output_start_time",startTime))
            if(self.turbulence_model=="RANS"):
                target.write("%-50s = velocity temperature tke\n"%("ABL.bndry_var_names"))
            else:
                target.write("%-50s = velocity temperature\n"%("ABL.bndry_var_names"))
            target.write("%-50s = native\n"%("ABL.bndry_output_format"))
        elif(iomode==1 and (not self.caseType=="terrain_noprecursor")):
            target.write('%-50s = "../precursor/bndry_files"\n'%("ABL.bndry_file"))
            target.write("%-50s = 1\n"%("ABL.bndry_io_mode"))
            if(self.turbulence_model=="RANS"):
                target.write("%-50s = velocity temperature tke\n"%("ABL.bndry_var_names"))
            else:
                target.write("%-50s = velocity temperature\n"%("ABL.bndry_var_names"))
            target.write("%-50s = native\n"%("ABL.bndry_output_format"))


    def createAMRSourceTerm(self,target,sponge=-1,terrain=-1,turbine=-1):
        target.write("# Source\n")
        if((sponge==1 and self.turbulence_model=="RANS") or (self.caseType=="terrain_noprecursor")):
            #forcingterms="WindSpongeForcing ABLMeanBoussinesq BoussinesqBuoyancy  "
            try:
                molLength=self.yamlFile["molLength"]
            except:
                forcingterms="VelocityFreeAtmosphereForcing ABLMeanBoussinesq BoussinesqBuoyancy "            
            else:
                if(abs(float(molLength))>5000):
                    forcingterms="VelocityFreeAtmosphereForcing ABLMeanBoussinesq BoussinesqBuoyancy "  
                else:
                    forcingterms="ABLMeanBoussinesq BoussinesqBuoyancy " 
        elif(self.caseType=="precursor"):
            forcingterms="VelocityFreeAtmosphereForcing BoussinesqBuoyancy"
        else:
            forcingterms="ABLMeanBoussinesq BoussinesqBuoyancy RayleighDamping "
        try: 
            self.includeCoriolis=self.yamlFile["includeCoriolis"]
        except:
            forcingterms=forcingterms+" CoriolisForcing "
        else:
            if(self.includeCoriolis):
                forcingterms=forcingterms+" CoriolisForcing "
        try:
            self.forcingHeight=self.yamlFile["forcingHeight"]
        except:
            forcingterms=forcingterms+" GeostrophicForcing "
        # try:
        #     metMastRegions=self.yamlFile["metMastNames"]
        # except:
        #     pass
        # else:
        #     forcingterms=forcingterms+" MetMastForcing "
        #     target.write('%-50s = "metmast.info"\n'%("ABL.metmast_1dprofile_file"))
        if(not self.turbulence_model=="RANS"):
            forcingterms=forcingterms+" NonLinearSGSTerm "
        if(terrain==1 or turbine==1):
            forcingterms=forcingterms+" DragForcing "
        if(turbine==1):
            forcingterms=forcingterms+" ActuatorForcing "
        target.write("%-50s = %s\n"%("ICNS.source_terms",forcingterms))
        try:
            self.forcingHeight=self.yamlFile["forcingHeight"]
        except:
            target.write("%-50s = %g %g %g\n"%("GeostrophicForcing.geostrophic_wind",self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        else:
            target.write("%-50s = %g \n"%("ABLForcing.abl_forcing_height",self.forcingHeight))
        if(terrain==1 or turbine==1):
            target.write('%-50s = "terrain.amrwind.new" \n'%("TerrainDrag.terrain_file"))
            if((self.turbulence_model=="RANS" and sponge==1) or (self.caseType=="terrain_noprecursor")):
                target.write("%-50s = TemperatureFreeAtmosphereForcing  DragTempForcing\n"%("Temperature.source_terms"))
            else:
                target.write("%-50s = DragTempForcing\n"%("Temperature.source_terms"))
        try:
            roughnessFile=self.yamlFile["roughnessFile"]
        except:
            pass
        else:
            target.write('%-50s = "roughness.amrwind.new" \n'%("TerrainDrag.roughness_file"))
        if(self.caseType=="terrain_noprecursor"):
            target.write("%-50s = 1\n"%("DragForcing.sponge_west"))
            target.write("%-50s = 1\n"%("DragForcing.sponge_east"))
            target.write("%-50s = 1\n"%("DragForcing.sponge_north"))
            target.write("%-50s = 1\n"%("DragForcing.sponge_south"))
            target.write("%-50s = -%g\n"%("DragForcing.sponge_distance_west",self.caseWestFlat))
            target.write("%-50s = %g\n"%("DragForcing.sponge_distance_east",self.caseEastFlat))
            target.write("%-50s = -%g\n"%("DragForcing.sponge_distance_south",self.caseSouthFlat))
            target.write("%-50s = %g\n"%("DragForcing.sponge_distance_north",self.caseNorthFlat))

        else:
            target.write("%-50s = 0\n"%("DragForcing.sponge_west"))
            target.write("%-50s = 0\n"%("DragForcing.sponge_east"))
            target.write("%-50s = 0\n"%("DragForcing.sponge_north"))
            target.write("%-50s = 0\n"%("DragForcing.sponge_south"))
        try:
            target.write("%-50s = %g\n"%("DragTempForcing.soil_temperature",self.soil_temperature))
        except:
            target.write("%-50s = %g\n"%("DragTempForcing.soil_temperature",300.0))
        target.write("%-50s = 0 0 1\n"%("RayleighDamping.force_coord_directions"))
        target.write("%-50s = %g\n"%("BoussinesqBuoyancy.reference_temperature",self.refTemperature))
        target.write("%-50s = %g\n"%("BoussinesqBuoyancy.thermal_expansion_coeff",1.0/self.refTemperature))
        #if(self.includeCoriolis):
        # Write the coriolis term for geostrophic forcing term 
        target.write("%-50s = 1.0 0.0 0.0 \n"%("CoriolisForcing.east_vector"))
        target.write("%-50s = 0.0 1.0 0.0 \n"%("CoriolisForcing.north_vector"))
        target.write("%-50s = %g \n"%("CoriolisForcing.latitude",self.caseCenterLat))
        target.write("%-50s = %g %g %g\n"%("RayleighDamping.reference_velocity",self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        startRayleigh=self.maxZ-self.RDLHeight
        target.write("%-50s = %g\n"%("RayleighDamping.length_sloped_damping",512))
        target.write("%-50s = %g\n"%("RayleighDamping.length_complete_damping",self.maxZ-startRayleigh-512))
        target.write("%-50s = 20.0\n"%("RayleighDamping.time_scale"))     
        target.write("%-50s = %g\n"%("ABL.metmast_horizontal_radius",self.metmast_horizontal_radius))
        target.write("%-50s = %g\n"%("ABL.metmast_vertical_radius",self.metmast_vertical_radius))
        target.write("%-50s = %g\n"%("ABL.metmast_damping_radius",self.metmast_damping_radius))

    def createAMRBC(self,target,inflowOutflow=-1):
        target.write("# BC \n")
        boundaries=["xlo","xhi","ylo","yhi"]
        if(inflowOutflow==1 and (not self.caseType=="terrain_noprecursor")):
            for boundary in boundaries:
                target.write('%-50s = "mass_inflow_outflow"\n'%(boundary+".type "))
                target.write("%-50s = 1.225\n"%(boundary+".density"))
                target.write("%-50s = 300\n"%(boundary+".temperature"))
                if(self.turbulence_model=="RANS"):
                    target.write("%-50s = 0.1 \n"%(boundary+".tke"))
        if(self.caseType=="terrain_noprecursor"):
            for boundary in boundaries:
                target.write('%-50s = "pressure_outflow"\n'%(boundary+".type "))
        #if(inflowOutflow==1):
            # if(self.caseWindspeedX>=0):
            #     target.write('%-50s = "mass_inflow"\n'%("xlo.type "))
            #     target.write("%-50s = 1.225\n"%("xlo.density"))
            #     target.write("%-50s = 300\n"%("xlo.temperature"))
            #     if(self.turbulence_model=="RANS"):
            #         target.write("%-50s = 0.1 \n"%("xlo.tke"))
            #     target.write('%-50s = "pressure_outflow"\n'%("xhi.type"))
            # else:
            #     target.write('%-50s = "mass_inflow"\n'%("xhi.type"))
            #     target.write("%-50s = 1.225\n"%("xhi.density"))
            #     target.write("%-50s = 300\n"%("xhi.temperature "))
            #     if(self.turbulence_model=="RANS"):
            #         target.write("%-50s = 0.1 \n"%("xhi.tke"))
            #     target.write('%-50s = "pressure_outflow"\n'%("xlo.type"))   
            # if(self.caseWindspeedY>=0):
            #     target.write('%-50s = "mass_inflow"\n'%("ylo.type"))
            #     target.write("%-50s = 1.225\n"%("ylo.density"))
            #     target.write("%-50s = 300\n"%("ylo.temperature "))
            #     if(self.turbulence_model=="RANS"):
            #         target.write("%-50s = 0.1 \n"%("ylo.tke"))
            #     target.write('%-50s = "pressure_outflow"\n'%("yhi.type"))
            # else:
            #     target.write('%-50s = "mass_inflow"\n'%("yhi.type"))
            #     target.write("%-50s = 1.225\n"%("yhi.density"))
            #     target.write("%-50s = 300\n"%("yhi.temperature "))
            #     if(self.turbulence_model=="RANS"):
            #         target.write("%-50s = 0.1 \n"%("yhi.tke"))
            #     target.write('%-50s = "pressure_outflow"\n'%("ylo.type"))  
        target.write('%-50s = "slip_wall"\n'%("zhi.type"))
        target.write('%-50s = "fixed_gradient"\n'%("zhi.temperature_type"))
        target.write("%-50s =  0.003\n"%("zhi.temperature"))
        target.write('%-50s = "wall_model"\n'%("zlo.type"))

    def createAMRTolerance(self,target,modify=-1):
        #if(modify==1):
        if(self.caseverticalAR==3 or self.caseverticalAR==4):
            self.smoothing=8
        elif(self.caseverticalAR>4 and self.caseverticalAR<=8):
            self.smoothing=32
        elif(self.caseverticalAR>8 and self.caseverticalAR<=16):
            self.smoothing=64
        if(self.caseverticalAR>=3):
            target.write("%-50s = %g \n"%("mac_proj.num_pre_smooth",self.smoothing))
            target.write("%-50s = %g \n"%("mac_proj.num_post_smooth",self.smoothing))
        if(modify==1):
            target.write("%-50s = -1 \n"%("mac_proj.mg_rtol"))
            target.write("%-50s = 1e-4 \n"%("mac_proj.mg_atol"))
        else:
            target.write("%-50s = -1 \n"%("mac_proj.mg_rtol"))
            target.write("%-50s = 1e-4 \n"%("mac_proj.mg_atol"))
        target.write("%-50s = 25 \n"%("mac_proj.maxiter "))
        target.write("%-50s = 4\n"%("mac_proj.fmg_maxiter"))
        if(self.caseverticalAR>=3):
            target.write("%-50s = %g \n"%("nodal_proj.num_pre_smooth",self.smoothing))
            target.write("%-50s = %g \n"%("nodal_proj.num_post_smooth",self.smoothing))
        if(modify==1):
            target.write("%-50s = -1 \n"%("nodal_proj.mg_rtol"))
            target.write("%-50s = 1e-4 \n"%("nodal_proj.mg_atol "))                
        else:
            target.write("%-50s = -1 \n"%("nodal_proj.mg_rtol"))
            target.write("%-50s = 1e-4 \n"%("nodal_proj.mg_atol"))
        target.write("%-50s = 25 \n"%("nodal_proj.maxiter"))  
        target.write("%-50s = 4\n"%("nodal_proj.fmg_maxiter")) 
        target.write("%-50s = -1 \n"%("diffusion.mg_rtol"))
        target.write("%-50s = 1e-4 \n"%("diffusion.mg_atol "))
        target.write("%-50s = -1 \n"%("temperature_diffusion.mg_rtol"))
        target.write("%-50s = 1e-4 \n"%("temperature_diffusion.mg_atol"))
        target.write("%-50s = -1 \n"%("tke_diffusion.mg_rtol"))
        target.write("%-50s = 1e-4 \n"%("tke_diffusion.mg_atol"))

    
    def createAMR1dSolver(self):
        try:
            mol_length=float(self.yamlFile["molLength"])
        except:
            mol_length=-1e30
        try:
            allowed_error=self.yamlFile["allowedError"]
        except:
            allowed_error=0.25
        try:
            zheight=self.yamlFile["ransDomainTop"]
        except:
            zheight=self.terrainZMax+self.ABLHeight
        wind=self.metMastWind
        dz=8.0
        npts=int(zheight/dz)
        num_of_steps=30000
        tolerance=1e-3
        if(mol_length>0 and mol_length<500):
            num_of_steps=50000
            tolerance=1e-4
        # Update roughness for sweep angles 
        try:
            int(self.sweep_angle)
        except:
            xdir=self.xref
            ydir=self.yref
        else:
            farmRadius=self.yamlFile['farmRadius']
            if(self.sweep_angle<=90):
                xdir= self.xref-farmRadius* np.cos(self.sweep_angle*np.pi/180)
                ydir= self.yref-farmRadius* np.sin(self.sweep_angle*np.pi/180)
            elif(self.sweep_angle>90 and self.sweep_angle<=180):
                xdir= self.xref+farmRadius* np.cos((180-self.sweep_angle)*np.pi/180)
                ydir= self.yref- farmRadius* np.sin((180-self.sweep_angle)*np.pi/180)
            elif(self.sweep_angle>180 and self.sweep_angle<=270):
                xdir= self.xref+ farmRadius* np.cos((self.sweep_angle-180)*np.pi/180)
                ydir= self.yref+ farmRadius* np.sin((self.sweep_angle-180)*np.pi/180)
            else:
                xdir= self.xref- farmRadius* np.cos((360-self.sweep_angle)*np.pi/180)
                ydir= self.yref+ farmRadius* np.sin((360-self.sweep_angle)*np.pi/180)   
        self.refRoughness=self.roughness_interp(xdir-self.xref,ydir-self.yref)
        #print(xdir,ydir,self.refRoughness)
        roughness_length=self.refRoughness
        terrain_ht=0
        coriolis=self.caseCenterLat
        if(abs(float(mol_length))>500):
            inv_height=np.amax(self.terrainX3)+1500
        elif(mol_length<0):
            inv_height=843
        inv_width=0
        inv_strength=0
        lapse_rate=0.003
        heat_flux_mode=4
        # Change the wind if we are sweeping 
        try:
            int(self.sweep_angle)
        except:
            pass
        else:
            windspeed=10.0
            wind[0]=windspeed*np.cos(self.sweep_angle*np.pi/180)
            wind[1]=windspeed*np.sin(self.sweep_angle*np.pi/180)
            #print(wind)
        #initial_ug=wind[0]
        #if(self.caseCenterLat>0):
        #    initial_vg=max(wind[0],wind[1])*np.sin(self.caseCenterLat)
        #else:
        #    initial_vg=-max(wind[0],wind[1])*np.sin(self.caseCenterLat)
        # Setting approximate Geostrophic Wind 
        # https://www.researchgate.net/publication/230284417_Optimal_turbine_spacing_in_fully_developed_wind_farm_boundary_layers
        # M=np.sqrt(wind[0]**2+wind[1]**2)
        # ustar=M*0.41/np.log(self.metMastHeight/self.refRoughness)
        # f=abs(2*7.27e-5*np.sin(self.caseCenterLat*np.pi/180))
        # if(self.caseCenterLat>0):
        #     initial_ug=ustar/0.41*np.log(ustar/(f*self.refRoughness))-4.5
        #     initial_vg= -11.25*ustar
        # else:
        #     initial_ug=-ustar/0.41*np.log(ustar/(f*self.refRoughness))-4.5
        #     initial_vg= 11.25*ustar  
        # print(ustar,initial_ug,initial_vg)
        # 
        try:
            initial_ug=self.yamlFile["initialUG"]
            initial_vg=self.yamlFile["initialVG"]
            #print("User-Specified:",initial_ug,initial_vg)
        except:
            M=np.sqrt(wind[0]**2+wind[1]**2)
            wt=4e-3*M
            zi=1000.0
            f=2*7.27e-5*np.sin(self.caseCenterLat*np.pi/180)
            #print(wt/(f*zi))
            if(self.caseCenterLat<10):
                initial_ug=wind[0]
                initial_vg=wind[1]
            else:
                initial_ug=wind[0]+wt/(f*zi)*wind[1]
                initial_vg=wind[1]-wt/(f*zi)*wind[0]
        include_ti=False
        initial_ug,initial_vg,z0=self.generate_profile(allowed_error,self.metMastHeight,self.metMastWind,npts,zheight,roughness_length,terrain_ht, \
                        coriolis,inv_height,inv_width,inv_strength,lapse_rate,heat_flux_mode,mol_length,num_of_steps,tolerance, \
                            initial_ug,initial_vg,include_ti)
        try:
            self.fillrans1dinfo(self.amrPrecursorFile)
        except:
            self.fillrans1dinfo(' ')
        self.geostropicX=initial_ug
        self.geostropicY=initial_vg

    def fillrans1dinfo(self,target,sponge=-1):
        if(target==' '):
            pass
        else:
            stringtowrite="ABL.initial_wind_profile"
            target.write("%-50s = true\n"%(stringtowrite))
            stringtowrite="ABL.rans_1dprofile_file"
            target.write('%-50s = "rans_1d.info" \n'%(stringtowrite))
            zstart=self.terrainZMax+self.ABLHeight
            if(self.turbulence_model=="RANS"):
                stringtowrite="ABL.meso_sponge_start "
                target.write('%-50s = %g \n'%(stringtowrite,zstart))
            else:
                stringtowrite="ABL.meso_sponge_start "
                target.write('%-50s = %g \n'%(stringtowrite,zstart))

        # Write for AMR-Wind 
        data=np.genfromtxt(Path(self.caseParent,self.caseName,"1dSolverOutput.info").as_posix())
        zvals=data[:,0]
        uvals=data[:,1]
        vvals=data[:,2]
        wvals=data[:,3]
        tempvals=data[:,4]
        tkevals=data[:,5]
        if(target==' '):
            pass
        else:
            for i in range(0,len(zvals)):
                if(i==0):
                    stringtowrite="ABL.temperature_heights"
                    target.write("%-50s = %g"%(stringtowrite,zvals[i]))
                else:
                    target.write("  %g "%(zvals[i]))
            target.write(" %g \n"%(self.maxZ))
            for i in range(0,len(zvals)):
                if(i==0):
                    stringtowrite="ABL.temperature_values"
                    target.write("%-50s = %g"%(stringtowrite,tempvals[i]))
                else:
                    target.write("  %g "%(tempvals[i]))
            target.write(" %g \n"%(tempvals[i]))
        try:
            newtarget=open(Path(self.caseParent,self.caseName,"precursor","rans_1d.info").as_posix(),"w")
        except:
            pass
        else:
            for i in range(0,len(zvals)):
                if(zvals[i]>2048):
                    tkevalue=tkevals[i]
                    fixValue=True
                    break
            for i in range(0,len(zvals)):
                if(zvals[i]>2048 and fixValue):
                    newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevalue))
                else:
                    newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevals[i]))
            newtarget.write("%g %g %g 0 %g\n"%(self.maxZ,uvals[i],vvals[i],tkevalue))
            newtarget.close()  
        try:
            newtarget=open(Path(self.caseParent,self.caseName,"terrain","rans_1d.info").as_posix(),"w")
        except:
            pass 
        else:
            for i in range(0,len(zvals)):
                if(zvals[i]>2048):
                    tkevalue=tkevals[i]
                    fixValue=True
                    break
            for i in range(0,len(zvals)):
                if(zvals[i]>2048 and fixValue):
                    newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevalue))
                else:
                    newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevals[i]))
            newtarget.write("%g %g %g 0 %g\n"%(self.maxZ,uvals[i],vvals[i],tkevalue))
            newtarget.close()  
        try:
            newtarget=open(Path(self.caseParent,self.caseName,"terrainTurbine","rans_1d.info").as_posix(),"w")
        except:
            pass 
        else:
            tkevalue=1e-10
            for i in range(0,len(zvals)):
                if(zvals[i]>2048):
                    tkevalue=tkevals[i]
                    fixValue=True
                    break
            for i in range(0,len(zvals)):
                if(zvals[i]>2048 and fixValue):
                    newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevalue))
                else:
                    newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevals[i]))
            newtarget.write("%g %g %g 0 %g\n"%(self.maxZ,uvals[i],vvals[i],tkevalue))
            newtarget.close()    

    # Modify this function 
    def generate_profile(self,allowed_error,metMastHeight,metMastWind,npts,zheight,roughness_length,terrain_ht, \
                    coriolis,inv_height,inv_width,inv_strength,lapse_rate,heat_flux_mode,mol_length,num_of_steps,tolerance, \
                        initial_ug,initial_vg,include_ti=False):
        # Generate Geostrphic Wind for Terrain Height  Consistent Wind Speed 
        # Initial Guess of Geostropic Wind 
        # Coarse grid run to identify close wind speed 
        # A good initial guess reduces number of while loop iterations 
        print("Running 1-D Solver - Coarse Grid ")
        print("%-15s|%-15s|%-15s|%-15s"%("Geostrophic Wind ","Reference Wind "," CFD Wind"," Error"))
        residualx=100
        residualy=100 
        # Initialize Grid Npts, Height, Roughness Length and Terrain Height for IB 
        from amr1DSolver import amr1dSolver
        pathToWrite=Path(self.caseParent,self.caseName,"1dSolverOutput.info").as_posix()
        # Coarse Run 
        try:
            zheight=self.yamlFile["ransDomainTop"]
        except:
            zheight=self.terrainZMax+self.ABLHeight
        dz=32.0
        npts=int(zheight/dz)
        amr1D=amr1dSolver(npts,zheight,roughness_length,terrain_ht,pathToWrite)
        ug=[initial_ug,initial_vg]
        import time 
        start=time.time()
        while (residualx>allowed_error or residualy>allowed_error):
            # Initialize Phyiscs: ux,uy,T,tke, ustar, pblh (can leave tke ustar and pblh to default value)
            amr1D.initialize_physics(ug[0],ug[1],300,0.4,0.4,1000)
            # Coriolis 
            amr1D.initialize_coriolis(coriolis)
            # Temperature profile inversion height, width of strong inversion rate, strong inversion strength 
            # and lapse rate (can leave at default)
            amr1D.temperature_inversion(inv_height,inv_width,inv_strength,lapse_rate)
            # Heat Flux model: 1 - Heat flux specified ; 2 - Surface temperature specified 
            # 3 - Heating or cooling rate specified ; 4 - Monin Obukhov length specified 
            amr1D.heat_flux_model(heat_flux_mode,mol_length)
            # Simulation Iteration, Convergence Tolerance 
            try:
                forceMetMast=self.yamlFile["forceMetMast"]
            except:
                forceMetMast=False 
            if(forceMetMast):
                amr1D.run_simulation(num_of_steps,tolerance,True,self.metMastHeight,self.metMastWind[0],self.metMastWind[1])
            else:
                amr1D.run_simulation(num_of_steps,tolerance)
            # Error calculation 
            # Wind speed at the metMast 
            z=amr1D.z
            ux=amr1D.ux
            uy=amr1D.uy
            # Correcting met-mast height 
            # Find the z terrain location for data
            error=10000 
            for j in range(0,len(self.terrainX1)):
                residual=np.sqrt((self.terrainX1[j])**2+(self.terrainX2[j])**2)
                if(residual<error):
                    error=residual
                    xterrain=self.terrainX1[j]
                    yterrain=self.terrainX2[j]
                    zterrainmin=self.terrainX3[j]
            #print(xterrain,yterrain,zterrainmin)
            try:
                int(self.sweep_angle)
                if(len(metMastHeight)==1):
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight+500,z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight+500,z,uy)
                else:
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight[0]+500,z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight[0]+500,z,uy)
            except:
                if(len(metMastHeight)==1):
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight,z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight,z,uy)
                else:
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight[0],z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight[0],z,uy)
            #print(met_mast_cfd_ux,met_mast_cfd_uy)
            errx=met_mast_cfd_ux-metMastWind[0]
            erry=met_mast_cfd_uy-metMastWind[1]
            print("%-8s %-7s|%-8s %-7s|%-8s %-7s|%-8s %-7s"%(round(ug[0],2),round(ug[1],2), \
                                                            round(metMastWind[0],2),round(metMastWind[1],2), \
                                                            round(met_mast_cfd_ux,2),round(met_mast_cfd_uy,2),\
                                                            round(errx,2),round(erry,2)))
                                                            
            #print("Met Mast Wind:[%g %g]"%(metMastWind[0],metMastWind[1]))
            #print("Specified Geostrophic Wind: [%g %g]"%(ug[0],ug[1]))
            #print("CFD Met Mast Wind and Error:[%g %g]  [%g %g]"%(met_mast_cfd_ux,met_mast_cfd_uy,met_mast_cfd_ux-metMastWind[0],met_mast_cfd_uy-metMastWind[1]))
            tke=np.interp(metMastHeight,z,amr1D.tke)
            M=np.sqrt(met_mast_cfd_ux**2+met_mast_cfd_uy**2)
            #print("TI:",np.sqrt(2.0/3.0*tke)/M*100)
            residualx=abs(met_mast_cfd_ux-metMastWind[0])
            residualy=abs(met_mast_cfd_uy-metMastWind[1])
            # Reduce only the higher error to  speed-up 
            if(residualx<allowed_error and residualy<allowed_error):
                print("Coarse grid converged")
                pass
            #elif(residualx>residualy):
            else:
                if(metMastWind[0]>0):
                    if(met_mast_cfd_ux>metMastWind[0]):
                        ug[0]=ug[0]-max(0.5*residualx,allowed_error)
                    else:
                        ug[0]=ug[0]+max(0.5*residualx,allowed_error)
                else:
                    if(met_mast_cfd_ux<metMastWind[0]):
                        ug[0]=ug[0]+max(0.5*residualx,allowed_error)
                    else:
                        ug[0]=ug[0]-max(0.5*residualx,allowed_error)
            #else:
                if(metMastWind[1]>0):
                    if(met_mast_cfd_uy>metMastWind[1]):
                        ug[1]=ug[1]-max(0.5*residualy,allowed_error)
                    else:
                        ug[1]=ug[1]+max(0.5*residualy,allowed_error)
                else:
                    if(met_mast_cfd_uy<metMastWind[1]):
                        ug[1]=ug[1]+max(0.5*residualy,allowed_error)
                    else:
                        ug[1]=ug[1]-max(0.5*residualy,allowed_error)
            ux,uy=amr1D.return_windspeed()
            amr1D.reinitialize_windspeed(ux,uy)
            end = time.time()
        print("Coarse grid time:",end-start)
        z0=roughness_length
        print("Running 1-D Solver - Fine Grid ")
        # Fine Run 
        # A good initial guess reduces number of while loop iterations 
        # print("Fine Run")
        residualx=100
        residualy=100 
        # Initialize Grid Npts, Height, Roughness Length and Terrain Height for IB 
        pathToWrite=Path(self.caseParent,self.caseName,"1dSolverOutput.info").as_posix()
        ux,uy=amr1D.return_windspeed()
        # Need to interpolate 
        # Coarse Run 
        start=time.time()
        try:
            zheight=self.yamlFile["ransDomainTop"]
        except:
            zheight=self.terrainZMax+self.ABLHeight
        dz=16.0
        npts=int(zheight/dz)
        znew=np.linspace(0,zheight,npts)
        amr1D=amr1dSolver(npts,zheight,roughness_length,terrain_ht,pathToWrite)
        uxnew=0*znew
        uynew=0*znew
        for i in range(0,len(znew)):
            uxnew[i]=np.interp(znew[i],z,ux)
            uynew[i]=np.interp(znew[i],z,uy)
        amr1D.reinitialize_windspeed(uxnew,uynew)
        start=time.time()
        while (residualx>allowed_error or residualy>allowed_error):
            # Initialize Phyiscs: ux,uy,T,tke, ustar, pblh (can leave tke ustar and pblh to default value)
            amr1D.initialize_physics(ug[0],ug[1],300,0.4,0.4,100)
            # Coriolis 
            amr1D.initialize_coriolis(coriolis)
            # Temperature profile inversion height, width of strong inversion rate, strong inversion strength 
            # and lapse rate (can leave at default)
            amr1D.temperature_inversion(inv_height,inv_width,inv_strength,lapse_rate)
            # Heat Flux model: 1 - Heat flux specified ; 2 - Surface temperature specified 
            # 3 - Heating or cooling rate specified ; 4 - Monin Obukhov length specified 
            amr1D.heat_flux_model(heat_flux_mode,mol_length)
            # Simulation Iteration, Convergence Tolerance 
            try:
                forceMetMast=self.yamlFile["forceMetMast"]
            except:
                forceMetMast=False 
            if(forceMetMast):
                amr1D.run_simulation(num_of_steps,tolerance,True,self.metMastHeight,self.metMastWind[0],self.metMastWind[1])
            else:
                amr1D.run_simulation(num_of_steps,tolerance)
            # Error calculation 
            # Wind speed at the metMast 
            z=amr1D.z
            ux=amr1D.ux
            uy=amr1D.uy
            try:
                int(self.sweep_angle)
                if(len(metMastHeight)==1):
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight+500,z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight+500,z,uy)
                else:
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight[0]+500,z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight[0]+500,z,uy)
            except:
                if(len(metMastHeight)==1):
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight,z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight,z,uy)
                else:
                    met_mast_cfd_ux=np.interp(zterrainmin+metMastHeight[0],z,ux)
                    met_mast_cfd_uy=np.interp(zterrainmin+metMastHeight[0],z,uy)
            #met_mast_cfd_ux=np.interp(metMastHeight,z,ux)
            #met_mast_cfd_uy=np.interp(metMastHeight,z,uy)
            errx=met_mast_cfd_ux-metMastWind[0]
            erry=met_mast_cfd_uy-metMastWind[1]
            print("%-8s %-7s|%-8s %-7s|%-8s %-7s|%-8s %-7s"%(round(ug[0],2),round(ug[1],2), \
                                                            round(metMastWind[0],2),round(metMastWind[1],2), \
                                                            round(met_mast_cfd_ux,2),round(met_mast_cfd_uy,2),\
                                                            round(errx,2),round(erry,2)))
            # print("Met Mast Wind:[%g %g]"%(metMastWind[0],metMastWind[1]))
            # print("Specified Geostrophic Wind: [%g %g]"%(ug[0],ug[1]))
            # print("CFD Met Mast Wind and Error:[%g %g]  [%g %g]"%(met_mast_cfd_ux,met_mast_cfd_uy,met_mast_cfd_ux-metMastWind[0],met_mast_cfd_uy-metMastWind[1]))
            tke=np.interp(metMastHeight,z,amr1D.tke)
            M=np.sqrt(met_mast_cfd_ux**2+met_mast_cfd_uy**2)
            #print("TI:",np.sqrt(2.0/3.0*tke)/M*100)
            residualx=abs(met_mast_cfd_ux-metMastWind[0])
            residualy=abs(met_mast_cfd_uy-metMastWind[1])
            # Reduce only the higher error to  speed-up 
            if(residualx<allowed_error and residualy<allowed_error):
                print("Fine grid converged")
                pass
            elif(residualx>residualy):
                if(metMastWind[0]>0):
                    if(met_mast_cfd_ux>metMastWind[0]):
                        ug[0]=ug[0]-max(0.5*residualx,allowed_error)
                    else:
                        ug[0]=ug[0]+max(0.5*residualx,allowed_error)
                else:
                    if(met_mast_cfd_ux<metMastWind[0]):
                        ug[0]=ug[0]+max(0.5*residualx,allowed_error)
                    else:
                        ug[0]=ug[0]-max(0.5*residualx,allowed_error)
            else:
                if(metMastWind[1]>0):
                    if(met_mast_cfd_uy>metMastWind[1]):
                        ug[1]=ug[1]-max(0.5*residualy,allowed_error)
                    else:
                        ug[1]=ug[1]+max(0.5*residualy,allowed_error)
                else:
                    if(met_mast_cfd_uy<metMastWind[1]):
                        ug[1]=ug[1]+max(0.5*residualy,allowed_error)
                    else:
                        ug[1]=ug[1]-max(0.5*residualy,allowed_error)
            ux,uy=amr1D.return_windspeed()
            amr1D.reinitialize_windspeed(ux,uy)
            end = time.time()
        z0=roughness_length
        print("Fine grid time:",end-start)
        self.Qh=amr1D.qh()
        self.soil_temperature=amr1D.temperature[0]
        return ug[0],ug[1],z0


    def createAMRPrecursorSampling(self,target):
        pass
        # x1=self.terrainX1.flatten(order='F')
        # x2=self.terrainX2.flatten(order='F')
        # x3=self.terrainX3.flatten(order='F')
        # target.write("# Cloud \n")
        # target.write("sampling.labels \t\t\t = height80m \n")
        # target.write("sampling.height80m.type \t\t\t = ProbeSampler \n")
        # target.write('sampling.height80m.probe_location_file \t\t\t = "height80m.txt" \n')
        # height80File=Path(self.caseParent,self.caseName,"precursor","height80m.txt").open("w")
        # height80File.write("%g\n"%(len(self.smoothTerrainX1)))
        # for i in range(0,len(self.terrainX1)):
        #     height80File.write("%g %g %g \n"%(self.smoothTerrainX1[i],self.smoothTerrainX2[i],self.smoothTerrainX3[i]+80))
        # height80File.close()

    def writeRefinementRegions(self,target):
        try: 
            refinementRegions=self.yamlFile["refinementRegions"]
        except:
            refinementRegions=[] 
        try:
            metMastRegions=self.yamlFile["metMastNames"]
        except:
            metMastRegions=[]
        if(len(refinementRegions)>0):
            try:
                refinementMinx=self.yamlFile["refinementMinX"]
            except:
                warnings.warn("Missing minimum X values. No refinements written")
                return 
            try:
                refinementMaxx=self.yamlFile["refinementMaxX"]
            except:
                warnings.warn("Missing maximum X values. No refinements written")
                return 
            try:
                refinementMiny=self.yamlFile["refinementMinY"]
            except:
                warnings.warn("Missing minimum Y values. No refinements written")
                return 
            try:
                refinementMaxy=self.yamlFile["refinementMaxY"]
            except:
                warnings.warn("Missing maximum Y values. No refinements written")
                return 
            try:
                refinementMaxz=self.yamlFile["heightAboveTerrain"]
            except:
                warnings.warn("Missing maximum Z values. No refinements written")
                return 
            try:
                refinementLevels=self.yamlFile["refinementLevels"]
            except:
                warnings.warn("No refinement levels specified")
                return 
        # Field refinement 
        try:
            fieldRefinement=self.yamlFile['refineTerrain']
        except:
            fieldRefinement=False
        target.write("# tagging\n")
        for i in range(0,len(refinementRegions)):
            if(fieldRefinement and i==0):
                target.write("%-50s = f1 %s %s "%("tagging.labels",refinementRegions[i],refinementRegions[i]+"terrain"))
            elif(i==0):
                target.write("%-50s = %s %s "%("tagging.labels",refinementRegions[i],refinementRegions[i]+"terrain"))
            else:
                target.write(" %s %s "%(refinementRegions[i],refinementRegions[i]+"terrain"))
        for i in range(0,len(metMastRegions)):
            if(i==0 and len(refinementRegions)==0 and fieldRefinement==True):
                target.write("%-50s = f1 %s %s "%("tagging.labels",metMastRegions[i],metMastRegions[i]+"terrain"))
            elif(i==0):
                target.write("%-50s = %s %s "%("tagging.labels",metMastRegions[i],metMastRegions[i]+"terrain"))
            else:
                target.write(" %s %s "%(metMastRegions[i],metMastRegions[i]+"terrain"))           
        target.write("\n")
        if(fieldRefinement):
            target.write("%-50s = FieldRefinement\n"%("tagging.f1.type"))
            target.write("%-50s = terrain_blank\n"%("tagging.f1.field_name"))
            target.write("%-50s = 0.1 0.1 0.1 0.1 0.1 0.1\n"%("tagging.f1.grad_error"))
            try:
                box_lo=self.yamlFile["refineLow"]
                box_hi=self.yamlFile["refineHigh"]
            except: 
                pass
            else:
                minTerrain=1e30
                maxTerrain=-1e30
                for j in range(0,len(self.terrainX1)):
                    if(self.terrainX1[j]>box_lo[0] and self.terrainX1[j]<box_hi[0] \
                       and self.terrainX2[j]>box_lo[1] and self.terrainX2[j]<box_hi[1]):
                        if(self.terrainX3[j]<minTerrain):
                            minTerrain=self.terrainX3[j]
                        if(self.terrainX3[j]>maxTerrain):
                            maxTerrain=self.terrainX3[j]
                #print("Terrain Hee Haw",minTerrain,maxTerrain)
                #print("MinMax",np.amin(self.terrainX3),np.amax(self.terrainX3))
                # target.write("%-50s = %g %g %g \n"%("tagging.f1.box_lo",box_lo[0],box_lo[1],box_lo[2]))
                # target.write("%-50s = %g %g %g \n"%("tagging.f1.box_hi",box_hi[0],box_hi[1],box_hi[2]))
                target.write("%-50s = %g %g %g \n"%("tagging.f1.box_lo",box_lo[0],box_lo[1],minTerrain+self.caseCellSize/self.caseverticalAR))
                target.write("%-50s = %g %g %g \n"%("tagging.f1.box_hi",box_hi[0],box_hi[1],maxTerrain+self.caseCellSize))
        for i in range(0,len(refinementRegions)):
            # Find the z terrain location for data
            error=10000 
            for j in range(0,len(self.terrainX1)):
                residual=np.sqrt((self.terrainX1[j]-refinementMinx[i])**2+(self.terrainX2[j]-refinementMiny[i])**2)
                if(residual<error):
                    error=residual
                    zterrainmin=self.terrainX3[j]
            error=10000 
            for j in range(0,len(self.terrainX1)):
                residual=np.sqrt((self.terrainX1[j]-refinementMaxx[i])**2+(self.terrainX2[j]-refinementMaxy[i])**2)
                if(residual<error):
                    error=residual
                    zterrainmax=self.terrainX3[j]
            zstart=min(zterrainmin,zterrainmax)-100
            zlength=abs(zterrainmin-zterrainmax)+refinementMaxz[i]
            #print(zstart,zlength)
            taggingstring="tagging."+refinementRegions[i]+".type"
            target.write("%-50s = GeometryRefinement\n"%(taggingstring))
            taggingstring="tagging."+refinementRegions[i]+".shapes"
            target.write("%-50s = object%g\n"%(taggingstring,i))
            taggingstring="tagging."+refinementRegions[i]+".min_level"          
            target.write("%-50s = %g\n"%(taggingstring,0))
            taggingstring="tagging."+refinementRegions[i]+".max_level"          
            target.write("%-50s = %g\n"%(taggingstring,max(refinementLevels[i]-1,0)))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".type"
            target.write("%-50s = box\n"%(taggingstring))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".origin"
            target.write("%-50s = %g %g %g\n"%(taggingstring,refinementMinx[i],refinementMiny[i],zstart))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".xaxis"
            target.write("%-50s = %g 0 0\n"%(taggingstring,
                                             refinementMaxx[i]-refinementMinx[i]))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".yaxis"
            target.write("%-50s = 0 %g 0\n"%(taggingstring,
                                             refinementMaxy[i]-refinementMiny[i]))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".zaxis"
            target.write("%-50s = 0 0 %g\n"%(taggingstring,zlength))
            taggingstring="tagging."+refinementRegions[i]+"terrain"+".type"
            target.write("%-50s = FieldRefinement\n"%(taggingstring))
            taggingstring="tagging."+refinementRegions[i]+"terrain"+".field_name"
            target.write("%-50s = terrain_blank\n"%(taggingstring))
            taggingstring="tagging."+refinementRegions[i]+"terrain"+".grad_error"
            target.write("%-50s = 0.1 0.1 0.1 0.1 0.1 0.1\n"%(taggingstring))
            taggingstring="tagging."+refinementRegions[i]+"terrain"+".box_lo"
            target.write("%-50s = %g %g %g \n"%(taggingstring,refinementMinx[i],refinementMiny[i],zstart))
            taggingstring="tagging."+refinementRegions[i]+"terrain"+".box_hi"
            target.write("%-50s = %g %g %g \n"%(taggingstring,refinementMaxx[i],refinementMaxy[i],zstart+abs(zterrainmin-zterrainmax)+200))
        
        if(len(metMastRegions)==0):
            pass
        else:
            try:
                metMastLatLon=self.yamlFile["metMastLatLon"]
            except:
                metMastLatLon=False 
            if(not metMastLatLon):
                try:
                    self.metMastX=self.yamlFile["self.metMastX"]
                except:
                    warnings.warn("Missing X values. No refinements written")
                    return    
                try:
                    self.metMastY=self.yamlFile["self.metMastY"]
                except:
                    warnings.warn("Missing Y values. No refinements written")
                    return        
            else: 
                try:
                    metMastLon=self.yamlFile["metMastLon"]
                except:
                    warnings.warn("Missing Longitude values. No refinements written")
                    return    
                try:
                    metMastLat=self.yamlFile["metMastLat"]
                    #print(metMastLat)
                except:
                    warnings.warn("Missing Latitude values. No refinements written")
                    return        
                self.metMastX=metMastLat
                self.metMastY=metMastLon
                for i in range(0,len(metMastRegions)):
                    self.metMastX[i],self.metMastY[i]=self.srtm.to_xy(metMastLat[i],metMastLon[i])
                    self.metMastX[i]=self.metMastX[i]-self.xref
                    self.metMastY[i]=self.metMastY[i]-self.yref
                #print(self.metMastX,self.metMastY)
            try:
                metMastRadius=self.yamlFile["metMastRadius"]
            except:
                metMastRadius=[]
                metMastRadius.append(500.0)
                for i in range(1,len(metMastRegions)):
                    metMastRadius.append(500.0)
            try:
                metMastRefinementLevel=self.yamlFile["metMastRefinementLevel"]
            except: 
                metMastRefinementLevel=[]
                metMastRefinementLevel.append(3)
                for i in range(1,len(metMastRegions)):
                    metMastRefinementLevel.append(3)
            try:
                metMastHeight=self.yamlFile["metMastHeight"]
            except:
                metMastHeight=200
            for i in range(0,len(metMastRegions)):
                error=10000 
                for j in range(0,len(self.terrainX1)):
                    residual=np.sqrt((self.terrainX1[j]-self.metMastX[i])**2+(self.terrainX2[j]-self.metMastY[i])**2)
                    if(residual<error):
                        error=residual
                        zterrainmin=self.terrainX3[j]
                taggingstring="tagging."+metMastRegions[i]+".type"
                target.write("%-50s = GeometryRefinement\n"%(taggingstring))
                taggingstring="tagging."+metMastRegions[i]+".shapes"
                target.write("%-50s = metmastobject%g\n"%(taggingstring,i))
                taggingstring="tagging."+metMastRegions[i]+".min_level"          
                target.write("%-50s = %g\n"%(taggingstring,0))
                taggingstring="tagging."+metMastRegions[i]+".max_level"          
                target.write("%-50s = %g\n"%(taggingstring,max(metMastRefinementLevel[i]-1,0)))
                taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".type"
                target.write("%-50s = cylinder \n"%(taggingstring))
                taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".start"
                target.write("%-50s = %g %g %g\n"%(taggingstring, \
                            self.metMastX[i],self.metMastY[i],zterrainmin))
                taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".end"
                target.write("%-50s = %g %g %g\n"%(taggingstring, \
                            self.metMastX[i],self.metMastY[i],zterrainmin+metMastHeight[i]+100))
                taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".outer_radius"
                target.write("%-50s = %g\n"%(taggingstring,metMastRadius[i]))
                taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".inner_radius"
                target.write("%-50s = %g\n"%(taggingstring,0.0))
                taggingstring="tagging."+metMastRegions[i]+"terrain"+".type"
                target.write("%-50s = FieldRefinement\n"%(taggingstring))
                taggingstring="tagging."+metMastRegions[i]+"terrain"+".field_name"
                target.write("%-50s = terrain_blank\n"%(taggingstring))
                taggingstring="tagging."+metMastRegions[i]+"terrain"+".grad_error"
                target.write("%-50s = 0.1 0.1 0.1 0.1 0.1 0.1\n"%(taggingstring))
                xmin=self.metMastX[i]-metMastRadius[i]
                xmax=self.metMastX[i]+metMastRadius[i]
                ymin=self.metMastY[i]-metMastRadius[i]
                ymax=self.metMastY[i]+metMastRadius[i]
                taggingstring="tagging."+metMastRegions[i]+"terrain"+".box_lo"
                target.write("%-50s = %g %g %g \n"%(taggingstring,xmin,ymin,zterrainmin-50))
                taggingstring="tagging."+metMastRegions[i]+"terrain"+".box_hi"
                target.write("%-50s = %g %g %g \n"%(taggingstring,xmax,ymax,zterrainmin+50))
        level=0
        try:
            level=max(level,max(refinementLevels))
        except:
            level=0
        try:
            level=max(level,max(metMastRefinementLevel))
        except:
            pass
        try:
            fieldRefinement=self.yamlFile['refineTerrain']
        except:
            fieldRefinement=False
        if(fieldRefinement):
            level=max(level,1)
        try:
            refineTerrainMaxLevel=self.yamlFile["refineTerrainMaxLevel"]
            level=max(level,refineTerrainMaxLevel)
        except:
            pass
        stringtowrite="amr.max_level "
        target.write("%-50s = %d\n"%(stringtowrite,level))
        if(self.caseType=="terrain_noprecursor"):
            pass
        else:
            self.amrPrecursorFile.write("%-50s = %d\n"%(stringtowrite,0))
        # Write for AEP
        newtarget=Path(self.caseParent,self.caseName,"utm.info").open("w")
        i=0
        for j in range(0,len(self.terrainX1)):
            residual=np.sqrt((self.terrainX1[j]-self.metMastX[i])**2+(self.terrainX2[j]-self.metMastY[i])**2)
            if(residual<error):
                error=residual
                zterrainmin=self.terrainX3[j]
        newtarget.write("MetMast%g,%g,%g,%g\n"%(i,self.metMastX[i],self.metMastY[i],zterrainmin))
        newtarget.write("utmreference,%g,%g,0\n"%(self.xref,self.yref))
        newtarget.close()
        # Write the metmast file 
        newtarget=Path(self.caseParent,self.caseName,"terrain","metmast.info").open("w")
        newtarget.write("%g %g %g %g %g 0.0 %g\n"%(self.metMastX[0],self.metMastY[0],self.metMastHeight[0], \
                                                self.caseWindspeedX,self.caseWindspeedY,self.refTemperature))
        #print("$$",zterrainmin)
        newtarget.close()



    def writeAccelerationMaps(self,target):
        try:
            writeTerrainSampling=self.yamlFile["writeTerrainSampling"]
        except: 
            writeTerrainSampling=True 
        if(writeTerrainSampling):
            # Write a different sampling file for each level since native processing is messy 
            try:
                offsets=self.yamlFile["verticalLevels"]
            except:
                offsets=[10,50,100,150,200,250]
            target.write("# postprocessing\n")
            for i in range(0,len(offsets)):
                if(i==0):
                    target.write("%-50s = %s "%("incflo.post_processing","terrain"+str(offsets[i])))
                else:
                    target.write(" %s "%("terrain"+str(offsets[i])))
                try:
                    roiLat=self.yamlFile["roiLat"]
                except:
                    pass
                else:
                    target.write(" %s "%("roi"+str(offsets[i])))

        metMastRegions=[]
        try:
            metMastLineSampling=self.yamlFile["metMastLineSampling"]
        except:
            metMastLineSampling=False
        else:
            if(metMastLineSampling):
                try:
                    metMastRegions=self.yamlFile["metMastNames"]
                except:
                    pass
                else:
                    for i in range(0,len(metMastRegions)):
                        target.write(" %s "%(metMastRegions[i]))
        target.write("\n")
        if(writeTerrainSampling):
            for levels in offsets:
                target.write("%-50s = velocity temperature \n"%("terrain"+str(levels)+".fields"))
                target.write('%-50s = "native"\n'%("terrain"+str(levels)+".output_format"))
                target.write("%-50s = 100\n"%("terrain"+str(levels)+".output_frequency"))
                target.write("%-50s = %s \n"%("terrain"+str(levels)+".labels","terrain"+str(levels)))
                samplingentity="terrain"+str(levels)+".terrain"+str(levels)+".type"
                target.write("%-50s = ProbeSampler\n"%(samplingentity))
                samplingentity="terrain"+str(levels)+".terrain"+str(levels)+".probe_location_file"
                target.write("%-50s = %s\n"%(samplingentity,'"terrain.csv"'))
                samplingentity="terrain"+str(levels)+".terrain"+str(levels)+".offset_vector"
                target.write("%-50s = 0 0 1\n"%(samplingentity))
                samplingentity="terrain"+str(levels)+".terrain"+str(levels)+".offsets"
                target.write("%-50s = %g \n"%(samplingentity,levels))
                try:
                    roiLat=self.yamlFile["roiLat"]
                except:
                    pass
                else:
                    target.write("%-50s = velocity temperature \n"%("roi"+str(levels)+".fields"))
                    target.write('%-50s = "native"\n'%("roi"+str(levels)+".output_format"))
                    target.write("%-50s = 100\n"%("roi"+str(levels)+".output_frequency"))
                    target.write("%-50s = %s \n"%("roi"+str(levels)+".labels","roi"+str(levels)))
                    samplingentity="roi"+str(levels)+".roi"+str(levels)+".type"
                    target.write("%-50s = ProbeSampler\n"%(samplingentity))
                    samplingentity="roi"+str(levels)+".roi"+str(levels)+".probe_location_file"
                    target.write("%-50s = %s\n"%(samplingentity,'"region_of_interest.csv"'))
                    samplingentity="roi"+str(levels)+".roi"+str(levels)+".offset_vector"
                    target.write("%-50s = 0 0 1\n"%(samplingentity))
                    samplingentity="roi"+str(levels)+".roi"+str(levels)+".offsets"
                    target.write("%-50s = %g \n"%(samplingentity,levels))
        if(metMastLineSampling and len(metMastRegions)>0):
            for i in range(0,len(metMastRegions)):
                target.write("%-50s = velocity temperature \n"%(str(metMastRegions[i])+".fields"))
                target.write('%-50s = "native"\n'%(str(metMastRegions[i])+".output_format"))
                target.write("%-50s = 60\n"%(str(metMastRegions[i])+".output_frequency"))
                target.write("%-50s = %s \n"%(str(metMastRegions[i])+".labels",str(metMastRegions[i])))
                samplingentity=str(metMastRegions[i])+"."+str(metMastRegions[i])+".type"
                target.write("%-50s = LineSampler\n"%(samplingentity))
                samplingentity=str(metMastRegions[i])+"."+str(metMastRegions[i])+".num_points"
                target.write("%-50s = 50\n"%(samplingentity))
                samplingentity=str(metMastRegions[i])+"."+str(metMastRegions[i])+".start"
                zstart=0
                # Find z from terrain 
                error=10000
                for ii in range(0,len(self.terrainX1)):
                    residual=np.sqrt((self.metMastX[i]-self.terrainX1[ii])**2+(self.metMastY[i]-self.terrainX2[ii])**2)
                    if(residual<error):
                        error=residual 
                        zstart=self.terrainX3[ii]
                zstart=zstart-20
                target.write("%-50s = %g %g %g\n"%(samplingentity,self.metMastX[i],self.metMastY[i],zstart))
                samplingentity=str(metMastRegions[i])+"."+str(metMastRegions[i])+".end"
                target.write("%-50s = %g %g %g\n"%(samplingentity,self.metMastX[i],self.metMastY[i],zstart+200))



    def closeAMRFiles(self):
        try:
            self.amrPrecursorFile.close()
        except: 
            pass 
        try:
            self.amrTerrainFile.close()
        except:
            pass

    def writeTerrainData(self,folder):
        #("Writing Terrain Data")
        x1=self.terrainX1.flatten(order='F')
        x2=self.terrainX2.flatten(order='F')
        x3=self.terrainX3.flatten(order='F')  
        #x=np.arange(np.amin(x1),np.amax(x1),self.caseCellSize)
        #y=np.arange(np.amin(x2),np.amax(x2),self.caseCellSize)
        #x=np.arange(np.amin(x1),np.amax(x1),self.caseCellSize)
        #y=np.arange(np.amin(x2),np.amax(x2),self.caseCellSize)
        x=np.arange(np.amin(x1),np.amax(x1),30)
        y=np.arange(np.amin(x2),np.amax(x2),30)
        from scipy.interpolate import NearestNDInterpolator
        self.interp = NearestNDInterpolator(list(zip(x1,x2)),x3) 
        xterrain,yterrain=np.meshgrid(x,y)
        zterrain = self.interp(xterrain,yterrain)    
        import matplotlib.pylab as plt 
        plt.contourf(xterrain,yterrain,zterrain)
        x1=xterrain.flatten(order='F')
        x2=yterrain.flatten(order='F')
        x3=zterrain.flatten(order='F')
        #print("Shape:",xterrain.shape)
        # target=Path(self.caseParent,self.caseName,folder,"terrain.amrwind.new").open("w")
        # target.write("%d\n"%(xterrain.shape[0]))
        # target.write("%d\n"%(yterrain.shape[0]))
        # xvals=np.linspace(min(x1),max(x1),xterrain.shape[0])
        # for i in range(0,len(xvals)):
        #     target.write("%g\n"%(xvals[i]))
        # yvals=np.linspace(min(x2),max(x2),xterrain.shape[1]) 
        # for i in range(0,len(yvals)):
        #     target.write("%g\n"%(yvals[i]))
        # for i in range(0,xterrain.shape[0]):
        #     for j in range(0,xterrain.shape[1]):
        #         target.write("%g\n"%(zterrain[i,j]))
        # target.close()
        target=Path(self.caseParent,self.caseName,folder,"terrain.amrwind").open("w")
        for i in range(0,len(x1)):
             target.write("%g %g %g\n"%(x1[i],x2[i],x3[i]))
        target.close()

        data = np.loadtxt(Path(self.caseParent,self.caseName,folder,"terrain.amrwind").as_posix())
        x = np.unique(data[:, 0])
        y = np.unique(data[:, 1])
        z = data[:, 2]
        assert len(z) == (len(x) * len(y))
        with open(Path(self.caseParent,self.caseName,folder,"terrain.amrwind.new").as_posix(), "w") as f:
            f.write(f"{len(x)}\n")
            f.write(f"{len(y)}\n")
            x.tofile(f, sep="\n")
            f.write(f"\n")
            y.tofile(f, sep="\n")
            f.write(f"\n")
            z.tofile(f, sep="\n")
        if(self.write_stl):
            data=np.column_stack([x1,x2,x3])
            import pyvista as pv
            mesh=pv.PolyData(data)
            mesh['elevation']=data[:,2]
            mesh.save(Path(self.caseParent,self.caseName,folder,"terrainPoints.vtk").as_posix())
        # Coarser Terrain file 
        x=np.arange(np.amin(x1),np.amax(x1),self.caseCellSize)
        y=np.arange(np.amin(x2),np.amax(x2),self.caseCellSize)
        xterrain,yterrain=np.meshgrid(x,y)
        zterrain = self.interp(xterrain,yterrain)  
        x1=xterrain.flatten(order='F')
        x2=yterrain.flatten(order='F')
        x3=zterrain.flatten(order='F')
        target=Path(self.caseParent,self.caseName,folder,"terrain.csv").open("w")
        metMastLat=self.yamlFile["metMastLat"]
        metMastLon=self.yamlFile["metMastLon"]
        target.write("%d \n"%(len(x1)+len(metMastLat)))
        for i in range(0,len(metMastLat)):
            tempx,tempy=self.srtm.to_xy(metMastLat[i],metMastLon[i])
            tempx-=self.xref
            tempy-=self.yref
            error=10000 
            for j in range(0,len(self.terrainX1)):
                residual=np.sqrt((self.terrainX1[j]-tempx)**2+(self.terrainX2[j]-tempy)**2)
                if(residual<error):
                    error=residual
                    zterrainmin=self.terrainX3[j]
            target.write("%g %g %g\n"%(tempx,tempy,zterrainmin))
        for i in range(0,len(x1)):
             target.write("%g %g %g\n"%(x1[i],x2[i],x3[i]))
        target.close()
        # Finer Wind Farm 
        try:
            roi_lat=self.yamlFile["roiLat"]
            roi_lon=self.yamlFile["roiLon"]
        except:
            pass
        else:
            target=Path(self.caseParent,self.caseName,folder,"region_of_interest.csv").open("w")
            lowerx,lowery= self.srtm.to_xy(roi_lat[0],roi_lon[0])
            lowerx-=self.xref
            lowery-=self.yref
            upperx,uppery= self.srtm.to_xy(roi_lat[1],roi_lon[1])
            upperx-=self.xref
            uppery-=self.yref
            roi_range_x=np.arange(lowerx,upperx,0.25*self.caseCellSize)
            roi_range_y=np.arange(lowery,uppery,0.25*self.caseCellSize)
            xroi,yroi=np.meshgrid(roi_range_x,roi_range_y)
            zroi = self.interp(xroi,yroi)    
            x1=xroi.flatten(order='F')
            x2=yroi.flatten(order='F')
            x3=zroi.flatten(order='F')
            target.write("%d \n"%(len(x1)))
            for i in range(0,len(x1)):
                target.write("%g %g %g\n"%(x1[i],x2[i],x3[i]))
            target.close()
        if(self.write_stl):
            try:
                roi_lat=self.yamlFile["roiLat"]
                roi_lon=self.yamlFile["roiLon"]
            except:
                pass
            else:
                data=np.column_stack([x1,x2,x3])
                import pyvista as pv
                mesh=pv.PolyData(data)
                mesh['elevation']=data[:,2]
                mesh.save(Path(self.caseParent,self.caseName,folder,"roiPoints.vtk").as_posix())


               # Roughness file 
        try:
            roughnessFile=self.yamlFile["roughnessFile"]
        except:
            pass
        else:           
            #print("Opening file",Path(self.caseParent,self.caseName,folder,"roughness.amrwind").as_posix())
            data=np.genfromtxt(Path(self.caseParent,self.caseName,"roughness.amrwind").as_posix())
            dx=np.amax(x)-np.amin(x)
            nx=int(dx/300.0)
            xtemp = np.linspace(np.amin(x),np.amax(x),nx)
            dy=np.amax(y)-np.amin(y)
            ny=int(dx/300.0)
            ytemp = np.linspace(np.amin(y),np.amax(y),ny)
            ztemp=[]
            for i in range (0,len(xtemp)):
                for j in range(0,len(ytemp)):
                    z0=self.roughness_interp(xtemp[i],ytemp[j])
                    ztemp.append(z0)
            #z = data[:, 2]
            #print(latlist[0:10])
            #print(y)
            #print(len(x),len(y),len(z))
            assert len(ztemp) == (len(xtemp) * len(ytemp))
            ztemp=np.array(ztemp)
            # for i in range(0,len(x)):
            #     xtemp,ytemp,_,_=utm.from_latlon(y[i],x[i],force_zone_number=self.zone_number)
            #     x[i]=xtemp-self.xref
            #     if(self.caseCenterLat<0):
            #         ytemp=ytemp-10000000 
            #     y[i]=ytemp-self.yref
            #print("Writing file",Path(self.caseParent,self.caseName,folder,"roughness.amrwind.new").as_posix())
            with open(Path(self.caseParent,self.caseName,folder,"roughness.amrwind.new").as_posix(), "w") as f:
                f.write(f"{len(xtemp)}\n")
                f.write(f"{len(ytemp)}\n")
                xtemp.tofile(f, sep="\n")
                f.write(f"\n")
                ytemp.tofile(f, sep="\n")
                f.write(f"\n")
                ztemp.tofile(f, sep="\n")
        #xterrain,yterrain=np.meshgrid(x,y)
        # for i in range(0,xterrain.shape[0]):
        #     for j in range(0,xterrain.shape[1]):
        #         error=100000000
        #         for k in range(0,len(x1)):
        #             error=np.sqrt((xterrain[i,j]-x1[k])**2+(yterrain[i,j]-x2[k])**2)
        #             if(error<self.caseCellSize):
        #                 print(i,j,xterrain[i,j],yterrain[i,j],x1[k],x2[k],x3[k],error)
        #                 break
        # target=Path(self.caseParent,self.caseName,folder,"terrain.amrwind").open("w")
        # self.smoothTerrainX1=smoothData.points[:,0]
        # self.smoothTerrainX2=smoothData.points[:,1]
        # self.smoothTerrainX3=smoothData.points[:,2]
        # for i in range(0,len(smoothData.points[:,0])):
        #     target.write("%g %g %g\n"%(smoothData.points[i,0],smoothData.points[i,1],smoothData.points[i,2]))
        # target.close()
        # Temporary Roughness Data 
        #print("Writing Roughness Data")
        # self.smoothData=smoothData
        # self.roughnessData=0*smoothData.points[:,0]
        # import geopandas as gpd
        # from shapely import Polygon
        # import pygeohydro as gh
        # offset=0.15
        # GEOM = Polygon(
        #     [
        #         [self.caseWest-offset,self.caseSouth-offset],
        #         [self.caseEast+offset,self.caseSouth-offset],
        #         [self.caseEast+offset,self.caseNorth+offset],
        #         [self.caseWest-offset,self.caseNorth+offset],  
        #         [self.caseWest-offset,self.caseSouth-offset],    
        #     ]
        # )
        # DEF_CRS = 4326
        # ALT_CRS = 3542
        # years = {"cover": [2021]}
        # res = 300
        # geom = gpd.GeoSeries([GEOM], crs=DEF_CRS)
        # lulc = gh.nlcd_bygeom(geom, years=years, resolution=res, crs=ALT_CRS, ssl=False)
        # roughness = gh.overland_roughness(lulc[0].cover_2021)
        # roughNumpy=np.flipud(roughness.to_numpy())
        # lat=np.linspace(self.caseSouth-offset,self.caseNorth+offset,roughNumpy.shape[0])
        # lon=np.linspace(self.caseWest-offset,self.caseEast+offset,roughNumpy.shape[1])
        # self.roughLon,self.roughLat=np.meshgrid(lon,lat)
        # lat=[]
        # lon=[]
        # roughness1D=[]
        # for i in range(0,self.roughLat.shape[0]):
        #     for j in range(0,self.roughLat.shape[1]):
        #         if(not np.isnan(roughNumpy[i,j])):
        #             lat.append(self.roughLat[i,j])
        #             lon.append(self.roughLon[i,j])
        #             roughness1D.append(roughNumpy[i,j])
        # xrougharray=[]
        # yrougharray=[]
        # for i in range(0,len(lat)):
        #     xrough,yrough=self.srtm.to_xy(lat[i],lon[i])
        #     xrougharray.append(xrough)
        #     yrougharray.append(yrough)
        # xrougharray=np.asarray(xrougharray)
        # yrougharray=np.asarray(yrougharray)
        # data=np.column_stack([xrougharray,yrougharray,0*xrougharray])
        # mesh2=pv.PolyData(data)
        # mesh2['roughness']=roughness1D
        # # surf = mesh2.delaunay_2d()
        # # surf.save(Path(self.caseParent,self.caseName,folder,"terrainRoughness.vtk").as_posix())
        # print(self.roughLat.shape,roughNumpy.shape)
        # target=Path(self.caseParent,self.caseName,folder,"terrain.roughness").open("w")
        # for i in range(0,len(roughness1D)):
        #     target.write("%g %g %g\n"%(xrougharray[i],yrougharray[i],roughness1D[i]))
        # target.close()
        #import matplotlib.pylab as plt 
        #plt.contourf(self.roughLon,self.roughLat,roughNumpy)
        #plt.show()
    
    def makeRoughness(self):
        from netCDF4 import Dataset
        roughnessFile=self.yamlFile["roughnessFile"]
        ncfile=Dataset(roughnessFile,'r')
        lat=ncfile.variables['lat'][:]
        lon=ncfile.variables['lon'][:]
        landuseclass=ncfile.variables['lccs_class'][:]
        lowlat=self.caseCenterLat-max(self.caseNorth,self.caseSouth)/110e3-0.2
        highlat=self.caseCenterLat+max(self.caseNorth,self.caseSouth)/110e3+0.2
        lowlon=self.caseCenterLon-max(self.caseEast,self.caseWest)/110e3-0.2
        highlon=self.caseCenterLon+max(self.caseEast,self.caseWest)/110e3+0.2
        lowlatindex=(np.abs(lowlat - lat) <= 0.02).argmax()
        highlatindex=(np.abs(highlat - lat) <= 0.02).argmax()
        lowlonindex=(np.abs(lowlon - lon) <= 0.02).argmax()
        highlonindex=(np.abs(highlon - lon) <= 0.02).argmax()
        latlist=[]
        lonlist=[]
        zval=[]
        roughnesslist=[]
        #print(lat.shape,lon.shape,landuseclass.shape)
        roughness_length=[0.1,0.1,0.1,0.07,0.2,0.3,0.3,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.15,0.3,0.07,0.07,0.07,0.07,0.05,0.05,0.05,0.05,0.05,0.2,0.2,0.2,0.8,0.05,0.05,0.05,1e-4,0.0012]
        lccs_class = [0, 10, 11, 12, 20, 30, 40, 50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 110, 120, 121, 122, 130, 140, 150, 151, 152, 153, 160, 170, 180, 190, 200, 201, 202, 210, 220]
        for j in range(lowlonindex,highlonindex):
            for i in range(lowlatindex,highlatindex,-1):
                latlist.append(lat[i])
                lonlist.append(lon[j])
                zval.append(0)
                findindex=(np.abs(landuseclass[0,i,j] - lccs_class) <= 0.1).argmax()
                z0=roughness_length[findindex]
                roughnesslist.append(z0)
        import pyvista as pv
        data=np.column_stack((lonlist,latlist,zval))
        mesh=pv.PolyData(data)
        mesh['roughness']=roughnesslist
        mesh.save(Path(self.caseParent,self.caseName,'roughness.vtk').as_posix())
        #print("Writing default roughness file")
        target=Path(self.caseParent,self.caseName,"roughness.amrwind").open("w")
        x=[]
        y=[]
        for i in range(0,len(roughnesslist)):
            #target.write("%g %g %g\n"%(lonlist[i],latlist[i],roughnesslist[i])
            xtemp,ytemp,_,_=utm.from_latlon(latlist[i],lonlist[i],force_zone_number=self.zone_number)
            x.append(xtemp-self.xref)
            if(self.caseCenterLat<0):
                ytemp=ytemp-10000000 
            y.append(ytemp-self.yref)
            target.write("%g %g %g\n"%(x[i],y[i],roughnesslist[i]))
        target.close()
        self.roughness_interp = NearestNDInterpolator(list(zip(x, y)), roughnesslist)



    def writeRestart(self,target):
        if(self.caseType=="terrain_noprecursor"):
            pass
        else:
            target.write("#io \n")
            if(not self.turbulence_model=="RANS"):
                target.write('%-50s = "../precursor/chk03600"\n'%("io.restart_file"))
            else:
                target.write('%-50s = "../precursor/chk01200"\n'%("io.restart_file"))

    def createTurbineFiles(self):
        print("Creating Turbines")
        latmin,lonmin=self.srtm.to_latlon(self.xref-abs(np.amin(self.terrainX1)),self.yref-abs(np.amin(self.terrainX2)))
        latmax,lonmax=self.srtm.to_latlon(np.amax(self.terrainX1)+self.xref,np.amax(self.terrainX2)+self.yref)
        print(self.xref,self.yref)
        print(latmin,lonmin)
        print(latmax,lonmax)
        self.turbineMarkType=self.yamlFile['turbineMarkType']
        if(self.turbineMarkType=='database'):
            xmin=latmin+0.01
            ymin=lonmin+0.01
            xmax=latmax-0.01
            ymax=lonmax-0.01
            import pandas as pd
            df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1',usecols=[1,2])   
            data=df.to_numpy()
            lon=data[:,0]
            lat=data[:,1]
            df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1')
            data=df.to_numpy()
            # data=np.genfromtxt("turbine.csv",dtype=str,delimiter=",",invalid_raise = False,skip_header=1)
            location=data[:,0]
            index=0
            self.caseTurbineLat=[]
            self.caseTurbineLon=[]
            self.turbineName=[]
            for i in range(0,len(lat)):
                if(lat[i]>xmin and lat[i]<xmax and lon[i]>ymin and lon[i]<ymax):
                    index=index+1
                    self.caseTurbineLat.append(lat[i])
                    self.caseTurbineLon.append(lon[i]) 
                    self.turbineName.append(location[i])
                    print(location[i],xmin,xmax,lat[i],ymin,ymax,lon[i])
            target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineLabels.info").open("w")
            target.write("Actuator.labels =")
            index=0
            for i in range(0,len(self.caseTurbineLon)):
                    index=index+1
                    target.write(" Turb%g "%(index))
            target.write("\n")
            target.close()
            index=0
            xTerrain=self.terrainX1.flatten(order='F')
            yTerrain=self.terrainX2.flatten(order='F')
            zTerrain=self.terrainX3.flatten(order='F')
            zTurbineLoc=[]
            self.turbineX1=[]
            self.turbineX2=[]
            self.turbineX3=[]
            from scipy.interpolate import NearestNDInterpolator
            self.interp = NearestNDInterpolator(list(zip(xTerrain,yTerrain)),zTerrain)
            self.turbineType=self.yamlFile["turbineType"]
            minx=1000000
            maxx=-1000000
            miny=1000000
            maxy=-1000000
            turbine=1
            for i in range(0,len(self.caseTurbineLon)):
                xturb,yturb=self.srtm.to_xy(self.caseTurbineLat[i],self.caseTurbineLon[i])
                xturb-=self.xref
                yturb-=self.yref
                residual=1000000
                # for k in range(0,len(xTerrain)):
                #     error=np.sqrt((xTerrain[k]-xturb)**2+(yTerrain[k]-yturb)**2)
                #     if(error<residual):
                #         residual=error
                #         zturb=zTerrain[k]
                #         kloc=k
                # zTurbineLoc.append(zturb)
                zturb=self.interp(xturb,yturb)
                self.turbineX1.append(xturb)
                self.turbineX2.append(yturb)
                self.turbineX3.append(zturb)
                index=index+1
                #print(index,xturb,yturb,xTerrain[kloc],yTerrain[kloc],residual)
                if(xturb>-500 and xturb<2000 and yturb>-1000 and yturb<2000):
                    #print(i,len(self.caseTurbineLon),self.turbineName[i],index,xturb,yturb,zturb) 
                    minx=min(minx,xturb)
                    maxx=max(maxx,xturb)
                    miny=min(miny,yturb)
                    maxy=max(maxy,yturb)
                    # box_hr.HighT1_inflow0deg.type         = PlaneSampler
                    # box_hr.HighT1_inflow0deg.num_points   = 46 46
                    # box_hr.HighT1_inflow0deg.origin       = -2586.2500 2436.2500 353.7500
                    # box_hr.HighT1_inflow0deg.axis1        = 112.5000 0.0 0.0
                    # box_hr.HighT1_inflow0deg.axis2        = 0.0 112.5000 0.0
                    # box_hr.HighT1_inflow0deg.normal       = 0.0 0.0 1.0
                    # box_hr.HighT1_inflow0deg.offsets      = 0.0 2.5 5.0 7.5 10.0 12.5 15.0 17.5 20.0 22.5 25.0 27.5 30.0 32.5 35.0 37.5 40.0 \
                    # 42.5 45.0 47.5 50.0 52.5 55.0 57.5 60.0 62.5 65.0 67.5 70.0 72.5 75.0 77.5 80.0 82.5 85.0 87.5 90.0 92.5 95.0 97.5 100.0 \
                    # 102.5 105.0 107.5 110.0 112.5 115.0 117.5 120.0 122.5 125.0 127.5 130.0 132.5 135.0 137.5
                    string_to_use="box_hr.HighT"+str(turbine)
                    print(string_to_use+"_inflow0deg.type         = PlaneSampler")
                    print(string_to_use+"_inflow0deg.num_points   = 20 20")
                    print(string_to_use+"_inflow0deg.origin       = %g %g %g"%(xturb-80,yturb-80,zturb-80))
                    print(string_to_use+"_inflow0deg.axis1        = 160.00 0.0 0.0")
                    print(string_to_use+"_inflow0deg.axis2        = 0.0 160.00 0.0")
                    print(string_to_use+"_inflow0deg.normal       = 0.0 0.0 1.0")
                    print(string_to_use+"_inflow0deg.offsets      = 0.   2.   4.   6.   8.  10.  12.  14.  16.  18.  20. 22.  24.  26.  28.  30.  32.  34.  36.  38.  40.  42. 44.  46.  48.  50.  52.  54.  56.  58.  60.  62.  64. 66.  68.  70.  72.  74.  76.  78.  80.  82.  84.  86. 88.  90.  92.  94.  96.  98. 100. 102. 104. 106. 108. 110. 112. 114. 116. 118. 120. 122. 124. 126. 128. 130. 132. 134. 136. 138. 140. 142. 144. 146. 148. 150. 152. 154. 156. 158. 160.")
                    turbine+=1
                    if(not self.fastBoxes):
                        target=Path(self.caseParent,self.caseName,"terrainTurbine","turbine"+str(i+1)+"Details.info").open("w")
                        if(i==0):
                            self.turbineHeader(self.turbineType,target)
                        self.createDefaultTurbine(self.turbineType,xturb,yturb,zturb,index,target)
                        target.close()
        print("Level 4:",minx,maxx,miny,maxy)
        globalMindistance=100000
        for i in range(0,len(self.turbineX1)):
            distance=100000
            for j in range(0,len(self.turbineX1)):
                tempDistance=np.sqrt((self.turbineX1[i]-self.turbineX1[j])**2+(self.turbineX2[i]-self.turbineX2[j])**2+(self.turbineX3[i]-self.turbineX3[j])**2)
                if(tempDistance<distance and tempDistance>0):
                    distance=tempDistance
            if(distance<globalMindistance):
                globalMindistance=distance
            #print("Turbine %g minimum distance is %g"%(i,distance))
        print("Minimum Distance should be %g"%(globalMindistance))
        if(not self.fastBoxes):
            # Concatenate files 
            tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","tempFile.info").open("w")
            target=Path(self.caseParent,self.caseName,"terrainTurbine","terrainTurbine.inp").open("r")
            tempFile.write(target.read())
            target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineLabels.info").open("r")
            tempFile.write(target.read())
            for i in range(0,len(self.caseTurbineLon)):
                target=Path(self.caseParent,self.caseName,"terrainTurbine","turbine"+str(i+1)+"Details.info").open("r")
                tempFile.write(target.read()) 
            # Does not work below python 3.4
            tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","tempFile.info")
            target=Path(self.caseParent,self.caseName,"terrainTurbine","terrainTurbine.inp")
            tempFile.rename(target.as_posix())
            # Delete temporary files 
            tmpName=Path(self.caseParent,self.caseName,"terrainTurbine","turbineLabels.info")
            tmpName.unlink()
            for i in range(0,len(self.caseTurbineLon)):
                tmpName=Path(self.caseParent,self.caseName,"terrainTurbine","turbine"+str(i+1)+"Details.info")
                tmpName.unlink()
        # Create Refinement Level Labels 
        target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineTaggingLabels.info").open("w")
        target.write("tagging.labels =  g0 g1 g2 \n")
        index=0
        # for i in range(0,len(self.caseTurbineLon)):
        #         target.write(" g%g gg%g ggg%g"%(i,i,i))
        # target.write("\n")
        target.close()
        # Write Refinement Files
        import pyvista as pv
        tempTurbineRadius=max(50,0.4*globalMindistance)
        for i in range(0,len(self.caseTurbineLon)):
            target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineRefinement"+str(i+1)+"Details.info").open("w")
            xstart=self.turbineX1[i]-4*tempTurbineRadius
            ystart=self.turbineX2[i]-4*tempTurbineRadius
            zstart=self.turbineX3[i]-200
            if(i==0):
                target.write("tagging.g0.type \t\t\t = GeometryRefinement\n")   
                target.write("tagging.g1.type \t\t\t = GeometryRefinement\n")  
                target.write("tagging.g2.type \t\t\t = GeometryRefinement\n")  
                for j in range(0,len(self.caseTurbineLon)):
                    if(j==0):
                        target.write("tagging.g0.shapes \t\t\t = b%g "%(j))
                    else:
                        target.write(" b%g  "%(j))
                target.write("\n")
                for j in range(0,len(self.caseTurbineLon)):
                    if(j==0):
                        target.write("tagging.g1.shapes \t\t\t = bb%g "%(j))
                    else:
                        target.write(" bb%g  "%(j))
                target.write("\n")
                for j in range(0,len(self.caseTurbineLon)):
                    if(j==0):
                        target.write("tagging.g2.shapes \t\t\t = bbb%g "%(j))
                    else:
                        target.write(" bbb%g  "%(j))
                target.write("\n")
                target.write("tagging.g0.level \t\t\t = 0\n")   
                target.write("tagging.g1.min_level \t\t\t = 0\n")
                target.write("tagging.g1.max_level \t\t\t = 1\n")
                target.write("tagging.g2.min_level \t\t\t = 0\n")
                target.write("tagging.g2.max_level \t\t\t = 2\n")
            target.write("tagging.g0.b%g.type \t\t\t = box\n"%(i))
            target.write("tagging.g0.b%g.origin = %g %g %g \n"%(i,xstart,ystart,zstart))
            target.write("tagging.g0.b%g.xaxis =  %g %g %g\n"%(i,8*tempTurbineRadius,0,0))
            target.write("tagging.g0.b%g.yaxis =  %g %g %g\n"%(i,0,8*tempTurbineRadius,0))
            target.write("tagging.g0.b%g.zaxis = %g %g %g\n"%(i,0,0,600))
            xstart=self.turbineX1[i]-2*tempTurbineRadius
            ystart=self.turbineX2[i]-2*tempTurbineRadius
            zstart=self.turbineX3[i]-100
            #target.write("tagging.gg%g.type \t\t\t = GeometryRefinement\n"%(i))
            #target.write("tagging.gg%g.shapes \t\t\t = bb%g\n"%(i,i))
            target.write("tagging.g1.bb%g.type \t\t\t = box\n"%(i))
            target.write("tagging.g1.bb%g.origin = %g %g %g \n"%(i,xstart,ystart,zstart))
            target.write("tagging.g1.bb%g.xaxis =  %g %g %g\n"%(i,4*tempTurbineRadius,0,0))
            target.write("tagging.g1.bb%g.yaxis =  %g %g %g\n"%(i,0,4*tempTurbineRadius,0))
            target.write("tagging.g1.bb%g.zaxis = %g %g %g\n"%(i,0,0,400))
            xstart=self.turbineX1[i]-tempTurbineRadius
            ystart=self.turbineX2[i]-tempTurbineRadius
            zstart=self.turbineX3[i]-50
            #target.write("tagging.ggg%g.type \t\t\t = GeometryRefinement\n"%(i))
            #target.write("tagging.ggg%g.shapes \t\t\t = bbb%g\n"%(i,i))
            target.write("tagging.g2.bbb%g.type \t\t\t = box\n"%(i))
            target.write("tagging.g2.bbb%g.origin = %g %g %g \n"%(i,xstart,ystart,zstart))
            target.write("tagging.g2.bbb%g.xaxis =  %g %g %g\n"%(i,2*tempTurbineRadius,0,0))
            target.write("tagging.g2.bbb%g.yaxis =  %g %g %g\n"%(i,0,2*tempTurbineRadius,0))
            target.write("tagging.g2.bbb%g.zaxis = %g %g %g\n"%(i,0,0,300))
        # for i in range(0,len(self.caseTurbineLon)):
        #     target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineRefinement"+str(i+1)+"Details.info").open("w")
        #     xstart=self.turbineX1[i]-4*tempTurbineRadius
        #     ystart=self.turbineX2[i]-4*tempTurbineRadius
        #     zstart=self.turbineX3[i]-200
        #     target.write("tagging.g%g.type \t\t\t = GeometryRefinement\n"%(i))
        #     target.write("tagging.g%g.shapes \t\t\t = b%g bb%g bbb%g\n"%(i,i,i,i))
        #     target.write("tagging.g%g.level \t\t\t = 0\n"%(i))
        #     target.write("tagging.g%g.b%g.type \t\t\t = box\n"%(i,i))
        #     target.write("tagging.g%g.b%g.origin = %g %g %g \n"%(i,i,xstart,ystart,zstart))
        #     target.write("tagging.g%g.b%g.xaxis =  %g %g %g\n"%(i,i,8*tempTurbineRadius,0,0))
        #     target.write("tagging.g%g.b%g.yaxis =  %g %g %g\n"%(i,i,0,8*tempTurbineRadius,0))
        #     target.write("tagging.g%g.b%g.zaxis = %g %g %g\n"%(i,i,0,0,600))
        #     xstart=self.turbineX1[i]-2*tempTurbineRadius
        #     ystart=self.turbineX2[i]-2*tempTurbineRadius
        #     zstart=self.turbineX3[i]-100
        #     #target.write("tagging.gg%g.type \t\t\t = GeometryRefinement\n"%(i))
        #     #target.write("tagging.gg%g.shapes \t\t\t = bb%g\n"%(i,i))
        #     target.write("tagging.g%g.min_level \t\t\t = 0\n"%(i))
        #     target.write("tagging.g%g.max_level \t\t\t = 1\n"%(i))
        #     target.write("tagging.g%g.bb%g.type \t\t\t = box\n"%(i,i))
        #     target.write("tagging.g%g.bb%g.origin = %g %g %g \n"%(i,i,xstart,ystart,zstart))
        #     target.write("tagging.g%g.bb%g.xaxis =  %g %g %g\n"%(i,i,4*tempTurbineRadius,0,0))
        #     target.write("tagging.g%g.bb%g.yaxis =  %g %g %g\n"%(i,i,0,4*tempTurbineRadius,0))
        #     target.write("tagging.g%g.bb%g.zaxis = %g %g %g\n"%(i,i,0,0,400))
        #     xstart=self.turbineX1[i]-tempTurbineRadius
        #     ystart=self.turbineX2[i]-tempTurbineRadius
        #     zstart=self.turbineX3[i]-50
        #     #target.write("tagging.ggg%g.type \t\t\t = GeometryRefinement\n"%(i))
        #     #target.write("tagging.ggg%g.shapes \t\t\t = bbb%g\n"%(i,i))
        #     target.write("tagging.g%g.min_level \t\t\t = 0\n"%(i))
        #     target.write("tagging.g%g.max_level \t\t\t = 2\n"%(i))
        #     target.write("tagging.g%g.bbb%g.type \t\t\t = box\n"%(i,i))
        #     target.write("tagging.g%g.bbb%g.origin = %g %g %g \n"%(i,i,xstart,ystart,zstart))
        #     target.write("tagging.g%g.bbb%g.xaxis =  %g %g %g\n"%(i,i,2*tempTurbineRadius,0,0))
        #     target.write("tagging.g%g.bbb%g.yaxis =  %g %g %g\n"%(i,i,0,2*tempTurbineRadius,0))
        #     target.write("tagging.g%g.bbb%g.zaxis = %g %g %g\n"%(i,i,0,0,300))
            target.close()
            # Writing Bounding Boxes 
            xmin=self.turbineX1[i]-4*tempTurbineRadius
            ymin=self.turbineX2[i]-4*tempTurbineRadius
            xmax=self.turbineX1[i]+4*tempTurbineRadius
            ymax=self.turbineX2[i]+4*tempTurbineRadius
            zmin=self.turbineX3[i]-200
            zmax=self.turbineX3[i]+400
            outerBox=pv.Box(bounds=(xmin,xmax,ymin,ymax,zmin,zmax))
            xmin=self.turbineX1[i]-2*tempTurbineRadius
            ymin=self.turbineX2[i]-2*tempTurbineRadius
            xmax=self.turbineX1[i]+2*tempTurbineRadius
            ymax=self.turbineX2[i]+2*tempTurbineRadius
            zmin=self.turbineX3[i]-100
            zmax=self.turbineX3[i]+300
            middleBox=pv.Box(bounds=(xmin,xmax,ymin,ymax,zmin,zmax))
            xmin=self.turbineX1[i]-tempTurbineRadius
            ymin=self.turbineX2[i]-tempTurbineRadius
            xmax=self.turbineX1[i]+tempTurbineRadius
            ymax=self.turbineX2[i]+tempTurbineRadius
            zmin=self.turbineX3[i]-50
            zmax=self.turbineX3[i]+250
            innerBox=pv.Box(bounds=(xmin,xmax,ymin,ymax,zmin,zmax))
            if(i==0):
                globalOuterBox=outerBox
                globalMiddleBox=middleBox
                globalInnerBox=innerBox
            else:
                localBox=outerBox
                tempbox=globalOuterBox.merge([localBox])
                globalOuterBox=tempbox
                localBox=middleBox
                tempbox=globalMiddleBox.merge([localBox])
                globalMiddleBox=tempbox
                localBox=innerBox
                tempbox=globalInnerBox.merge([localBox])
                globalInnerBox=tempbox
        tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","outerBox.vtk")  
        pv.save_meshio(tempFile.as_posix(),globalOuterBox)    
        tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","middleBox.vtk")  
        pv.save_meshio(tempFile.as_posix(),globalMiddleBox)  
        tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","innerBox.vtk")  
        pv.save_meshio(tempFile.as_posix(),globalInnerBox)  
        # Concatenate files 
        tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","tempFile.info").open("w")
        target=Path(self.caseParent,self.caseName,"terrainTurbine","terrainTurbine.inp").open("r")
        tempFile.write(target.read())
        target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineTaggingLabels.info").open("r")
        tempFile.write(target.read())
        for i in range(0,len(self.caseTurbineLon)):
            target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineRefinement"+str(i+1)+"Details.info").open("r")
            tempFile.write(target.read()) 
        # Does not work below python 3.4
        tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","tempFile.info")
        target=Path(self.caseParent,self.caseName,"terrainTurbine","terrainTurbine.inp")
        tempFile.rename(target.as_posix())
        # Delete temporary files 
        tmpName=Path(self.caseParent,self.caseName,"terrainTurbine","turbineTaggingLabels.info")
        tmpName.unlink()
        for i in range(0,len(self.caseTurbineLon)):
            tmpName=Path(self.caseParent,self.caseName,"terrainTurbine","turbineRefinement"+str(i+1)+"Details.info")
            tmpName.unlink()

    def turbineHeader(self,turbineType,target):
        if(turbineType=="UniformCtDisk"):
            string="Actuator.UniformCtDisk"
        else:
            string="Actuator.JoukowskyDisk"
        target.write(string+".rotor_diameter = 100.0\n")
        target.write(string+".hub_height     = 80.0\n")
        target.write(string+".epsilon        = 5.0 5.0 5.0\n")
        target.write(string+".yaw            = 270.0\n")
        target.write(string+".output_frequency = 100\n")
        target.write(string+".diameters_to_sample = 2.5\n")
        target.write(string+".thrust_coeff  = 1.158093814 0.958191954 0.842281443 0.841410733 0.841410733 0.841410733 0.841410733 0.841410733 0.841410733 0.841410733 ")
        target.write("0.841410733 0.841410733 0.784534886 0.743327664 0.653457742 0.566093507 0.485349168 0.448263929 0.387457153 0.293997037 0.226171155 0.176059266  ") 
        target.write(" 0.13865674 0.110481935 0.089183188 0.072941144 0.060464209 0.05074999 0.043180894 0.037323406\n")
        target.write(string+".wind_speed  = 3.0 3.889649963239854 4.684006996752303 5.377830233987229 5.966542092267928 6.44625847394617 6.8138143922059236 7.066784852 ")
        target.write(" 7.203500851477444 7.22306038896904 7.320786359429763 7.535153078939617 7.864746237154081 8.30739130337076 8.860167873258558 9.519428936578247 10.280 ")
        target.write(" 10.681872976809931 11.13933247768231 12.08928744604103 13.12442240111568 14.237907914913496 15.422397632159566 16.670076738763772 17.972713521 ")
        target.write("19.321713675239476 20.708177009893884 22.122956165519163 23.556716965618207 25.0\n")
        target.write(string+".num_points_r   = 5\n")
        target.write(string+".num_points_t   = 5\n")
        if(not turbineType=="UniformCtDisk"):
            target.write(string+".rpm            = 8 8 8 9.172052607 10.17611854 10.99428938 11.62116715 12.05261594 12.28578923 12.31914861 12.48582322 12.85143216 13.413563 ")
            target.write("14.16850788 15.11128509 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026 ")
            target.write(" 15.16026965 15.16026965 15.16026965 15.16026965\n")
            target.write(string+".num_blades     = 3\n")
            target.write(string+".vortex_core_size = 13.0\n")
            target.write(string+".use_tip_correction = true\n")
            target.write(string+".use_root_correction = true\n")


    def createDefaultTurbine(self,turbineType,xi,yi,zi,index,target):
        string="Actuator.Turb"+str(index)
        target.write("# Turbine %g\n"%(index))
        target.write(string+".type = %s \n"%(turbineType))
        target.write(string+".base_position  = %g %g %g\n"%(xi,yi,zi))


    def plotTurbines(self):
        import pyvista as pv 
        import pandas as pd
        print("Plotting Turbines")
        x1=self.terrainX1.flatten(order='F')
        x2=self.terrainX2.flatten(order='F')
        x3=self.terrainX3.flatten(order='F')
        data=np.column_stack([x1,x2,x3])
        pl = pv.Plotter()
        mesh2=pv.PolyData(data)
        mesh2['elevation']=data[:,2]
        #surf = mesh2.delaunay_2d()
        #pl.add_mesh(surf)
        for i in range(0,len(self.turbineX1)):
                print("Plotting turbine",i+1)
                newDisk=pv.Disc(center=(self.turbineX1[i],self.turbineX2[i],self.turbineX3[i]+80),inner=8,outer=50,normal=(1.0, 0.0,0.0), r_res=1, c_res=24)
                #pl.add_mesh(newDisk)
                if(i==0):
                    globalBox=newDisk
                else:
                    localBox=newDisk
                    tempbox=globalBox.merge([localBox])
                    globalBox=tempbox
        for i in range(0,len(self.turbineX3)):
            self.turbineX3[i]+=80
        data=np.column_stack([self.turbineX1,self.turbineX2,self.turbineX3])
        mesh1=pv.PolyData(data)
        pl.add_mesh(mesh1,render_points_as_spheres=True,point_size=10)      
        tempFile=Path(self.caseParent,self.caseName,"terrainTurbine","turbines.vtk")  
        pv.save_meshio(tempFile.as_posix(),globalBox)
        pl.view_xy()
        pl.show_axes()
        pl.show()


from sys import argv 
amrRef=amrBackend(argv[1])
amrRef.createDomain()
amrRef.createAMRFiles()

