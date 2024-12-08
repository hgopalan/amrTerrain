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


class amrBackend():
    def __init__(self,yamlFile):
        self.yamlFilePath = Path(yamlFile)
        self.yamlFile=yaml.safe_load(self.yamlFilePath.open())
        self.caseCellSize=128 
        self.caseverticalAR=8
        self.turbulence_model='RANS'
        self.case_end_time=7200
        self.plotOutput=1800
        self.restartOutput=1800
        self.createCase()

    def createCase(self):
        self.caseParent=self.yamlFile['caseParent']
        self.caseName=self.yamlFile['caseFolder']
        self.caseType=self.yamlFile['caseType']
        self.caseInitial=self.yamlFile['caseInitial']
        caseDir=Path(self.caseParent,self.caseName)
        self.caseDir=caseDir.as_posix()
        caseDir.mkdir(parents=True,exist_ok=True)
        casePrecursor=Path(self.caseParent,self.caseName,"precursor")
        self.casePrecursor=casePrecursor.as_posix()
        casePrecursor.mkdir(parents=True,exist_ok=True)
        if(self.caseType=="terrain"):
            caseTerrain=Path(self.caseParent,self.caseName,"terrain")
            self.caseTerrain=caseTerrain.as_posix()
            caseTerrain.mkdir(parents=True,exist_ok=True)
        elif(self.caseType=="terrainTurbine"):
            caseTerrain=Path(self.caseParent,self.caseName,"terrainTurbine")
            self.caseTerrain=caseTerrain.as_posix()
            caseTerrain.mkdir(parents=True,exist_ok=True)    
        self.caseNorth=self.yamlFile['north']
        self.caseSouth=self.yamlFile['south']
        self.caseEast=self.yamlFile['east']
        self.caseWest=self.yamlFile['west']
        self.caseCenterLat=self.yamlFile["centerLat"]
        self.caseCenterLon=self.yamlFile["centerLon"]
        self.refHeight=self.yamlFile["refHeight"]
        # Reading optional variables 
        try:
            self.caseCellSize=self.yamlFile['cellSize']
        except:
            self.caseCellSize=128.0
        try:
            self.turbulence_model=self.yamlFile['turbulenceModel']
        except:
            self.turbulence_model='RANS'
        try:
            self.caseverticalAR=self.yamlFile['verticalAR']
        except:
            if(self.turbulence_model=='RANS'):
                self.caseverticalAR=6
            else:
                self.caseverticalAR=4
        try:
            self.caseNorthSlope=self.yamlFile['northSlope']
        except:
            self.caseNorthSlope=3000
        try:
            self.caseSouthSlope=self.yamlFile['southSlope']
        except:
            self.caseSouthSlope=3000
        try:
            self.caseEastSlope=self.yamlFile['eastSlope']
        except:
            self.caseEastSlope=3000
        try:
            self.caseWestSlope=self.yamlFile['westSlope']
        except:
            self.caseWestSlope=3000
        try:
            self.caseNorthFlat=self.yamlFile['northFlat']
        except:
            self.caseNorthFlat=1000
        try:
            self.caseSouthFlat=self.yamlFile['southFlat']
        except:
            self.caseSouthFlat=1000
        try:
            self.caseEastFlat=self.yamlFile['eastFlat']
        except:
            self.caseEastFlat=1000
        try:
            self.caseWestFlat=self.yamlFile['westFlat']
        except:
            self.caseWestFlat=1000
        # Optional
        try:
            self.rans_1d=self.yamlFile["rans1D"]
        except:
            self.rans_1d=False

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
                print(dataset.bounds)
                print(self.caseCenterLon,westlon,eastlon,self.caseCenterLat,southlat,northlat)
                #exit(-1)
            self.xref,self.yref,self.zRef,self.srtm=converter.SRTM_Converter(Path(self.caseParent,self.caseName).as_posix(),self.caseCenterLat,self.caseCenterLon,self.refHeight, \
                                                        self.caseWest,self.caseEast,self.caseSouth,self.caseNorth, \
                                                        self.caseWestSlope,self.caseEastSlope,self.caseSouthSlope,self.caseNorthSlope, \
                                                        self.caseWestFlat,self.caseEastFlat,self.caseSouthFlat,self.caseNorthFlat, self.usetiff, \
                                                            self.write_stl,westlon,eastlon,southlat,northlat)
                                                            #-3,3,-1.5,1.5)
        else:
            warnings.warn("Turbine or Met Mast Locations should be in STL coordinates. Not Lat/Lon")
            self.xref=0
            self.yref=0
            self.zRef=0
        #stlFile=Path(self.caseParent,self.caseName,"terrain.stl").as_posix()
        stlFile=Path(self.caseParent,self.caseName,"terrain.vtk").as_posix()
        import pyvista as pv 
        mesh=pv.read(stlFile)
        x1=mesh.points[:,0]
        x2=mesh.points[:,1]
        x3=mesh.points[:,2]
        for i in range(0,len(x3)):
            x3[i]=max(x3[i],0.0)
        self.terrainX1=x1[:]
        self.terrainX2=x2[:]
        self.terrainX3=x3[:]
        if(not self.write_stl):
            Path(self.caseParent,self.caseName,"terrain.vtk").unlink()
        # else:
        #     print("Writing curvature")
        #     mesh["curvature"]=mesh.curvature()
        #     mesh["terrainHeight"]=x3
        #     print(mesh["terrainHeight"])
        #     Path(self.caseParent,self.caseName,"terrain.vtk").unlink()
        #     mesh.save(stlFile)


    
    def createAMRFiles(self):
        self.amrPrecursorFile=Path(self.caseParent,self.caseName,"precursor","precursor.inp").open("w")
        self.amrPrecursorFile.write("# Generating the precursor file\n")
        self.createPrecursorFiles()
        if(self.caseType=='terrain'):
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
        print("Creating precursor")
        self.createAMRGeometry(self.amrPrecursorFile,1)
        self.createAMRGrid(self.amrPrecursorFile)
        self.createAMRTime(self.amrPrecursorFile)
        # Adding some patching 
        self.refTemperature=self.yamlFile["refTemperature"]
        self.refRoughness=self.yamlFile["refRoughness"]
        self.refHeatFlux=self.yamlFile["refHeatflux"]        
        if(self.rans_1d):
            self.createAMR1dSolver()
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
        print(" Done creating precursor")

    def createTerrainFiles(self,folder):
        print("Creating Terrain Files")
        self.createAMRGeometry(self.amrTerrainFile,-1)
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
        self.fillrans1dinfo(self.amrTerrainFile,1)
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
        # Add 1 km for ABL and 2 km for Rayleigh
        self.ABLHeight=1024
        self.RDLHeight=2048
        if(self.terrainZMax>self.ABLHeight):
            print("Not enough blockage")
            self.ABLHeight=2*self.ABLHeight
        self.maxZ=self.terrainZMax+self.ABLHeight+self.RDLHeight
        print(self.maxZ)
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
        print("dx,dz",self.caseCellSize,self.maxZ/nz)
        target.write("# Grid \n")
        target.write("%-50s = %g %g %g\n"%("amr.n_cell",nx,ny,nz))
        #target.write("%-50s = 0\n"%("amr.max_level"))

    def createAMRTime(self,target,blanking=-1):
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
        if(self.timeMethod=="step"):
            target.write("%-50s = -1\n"%("time.stop_time"))
            target.write("%-50s = %g\n"%("time.max_step",self.timeSteps))
        else:
            target.write("%-50s = %g\n"%("time.stop_time",self.case_end_time))
        target.write("%-50s = 1.0\n"%("time.initial_dt"))
        if(target==self.amrPrecursorFile):
            target.write("%-50s = 1\n"%("time.fixed_dt"))
        else:
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
            target.write("%-50s = velocity temperature mu_turb tke pressure \n"%("io.outputs"))
        else:
            target.write("%-50s = velocity temperature mu_turb pressure \n"%("io.outputs"))
        if(blanking==1):
            target.write("%-50s = terrain_blank terrain_drag \n"%("io.int_outputs"))

    def createSolverInfo(self,target,terrain=-1,turbine=-1):
        self.caseWindspeedX=self.yamlFile['windX']
        self.caseWindspeedY=self.yamlFile['windY']
        self.caseWindspeedZ=self.yamlFile['windZ']   
        try:
            wind=self.yamlFile["metMastWind"]
        except:
            pass
        else:
            self.caseWindspeedX=wind[0]
            self.caseWindspeedY=wind[1]
            self.caseWindspeedZ=0.0
        try:
            self.fastBoxes=self.yamlFile['fastBoxes']      
        except:
            self.fastBoxes=False
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
        self.refTemperature=self.yamlFile["refTemperature"]
        self.refRoughness=self.yamlFile["refRoughness"]
        self.refHeatFlux=self.yamlFile["refHeatflux"]
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
        target.write("%-50s = %g\n"%("ABL.surface_temp_flux",self.refHeatFlux))
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
        elif(iomode==1):
            target.write('%-50s = "../precursor/bndry_files"\n'%("ABL.bndry_file"))
            target.write("%-50s = 1\n"%("ABL.bndry_io_mode"))
            if(self.turbulence_model=="RANS"):
                target.write("%-50s = velocity temperature tke\n"%("ABL.bndry_var_names"))
            else:
                target.write("%-50s = velocity temperature\n"%("ABL.bndry_var_names"))
            target.write("%-50s = native\n"%("ABL.bndry_output_format"))


    def createAMRSourceTerm(self,target,sponge=-1,terrain=-1,turbine=-1):
        target.write("# Source\n")
        if(terrain==1 or turbine==1 or (sponge==1 and self.turbulence_model=="RANS")):
            forcingterms="WindSpongeForcing ABLMeanBoussinesq BoussinesqBuoyancy  "
        else:
            forcingterms="ABLMeanBoussinesq BoussinesqBuoyancy RayleighDamping "
        try: 
            self.includeCoriolis=self.yamlFile["includeCoriolis"]
        except:
            pass
        else:
            forcingterms=forcingterms+" CoriolisForcing "
        try:
            self.forcingHeight=self.yamlFile["forcingHeight"]
        except:
            forcingterms=forcingterms+" GeostrophicForcing "
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
            if(self.turbulence_model=="RANS" and sponge==1):
                target.write("%-50s = TempSpongeForcing  DragTempForcing\n"%("Temperature.source_terms"))
            else:
                target.write("%-50s = DragTempForcing\n"%("Temperature.source_terms"))
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

    def createAMRBC(self,target,inflowOutflow=-1):
        target.write("# BC \n")
        boundaries=["xlo","xhi","ylo","yhi"]
        if(inflowOutflow==1):
            for boundary in boundaries:
                target.write('%-50s = "mass_inflow_outflow"\n'%(boundary+".type "))
                target.write("%-50s = 1.225\n"%(boundary+".density"))
                target.write("%-50s = 300\n"%(boundary+".temperature"))
                if(self.turbulence_model=="RANS"):
                    target.write("%-50s = 0.1 \n"%(boundary+".tke"))
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
            self.smoothing=16
        elif(self.caseverticalAR>8 and self.caseverticalAR<=16):
            self.smoothing=64
        if(self.caseverticalAR>=3):
            target.write("%-50s = %g \n"%("mac_proj.num_pre_smooth",self.smoothing))
            target.write("%-50s = %g \n"%("mac_proj.num_post_smooth",self.smoothing))
        if(modify==1):
            target.write("%-50s = 1.0e-3 \n"%("mac_proj.mg_rtol"))
            target.write("%-50s = 1.0e-4 \n"%("mac_proj.mg_atol"))
        else:
            target.write("%-50s = 1.0e-4 \n"%("mac_proj.mg_rtol"))
            target.write("%-50s = 1.0e-6 \n"%("mac_proj.mg_atol"))
        target.write("%-50s = 25 \n"%("mac_proj.maxiter "))
        target.write("%-50s = 4\n"%("mac_proj.fmg_maxiter"))
        if(self.caseverticalAR>=3):
            target.write("%-50s = %g \n"%("nodal_proj.num_pre_smooth",self.smoothing))
            target.write("%-50s = %g \n"%("nodal_proj.num_post_smooth",self.smoothing))
        if(modify==1):
            target.write("%-50s = 1.0e-3 \n"%("nodal_proj.mg_rtol"))
            target.write("%-50s = 1.0e-4 \n"%("nodal_proj.mg_atol "))                
        else:
            target.write("%-50s = 1.0e-4 \n"%("nodal_proj.mg_rtol"))
            target.write("%-50s = 1.0e-6 \n"%("nodal_proj.mg_atol"))
        target.write("%-50s = 25 \n"%("nodal_proj.maxiter"))  
        target.write("%-50s = 4\n"%("nodal_proj.fmg_maxiter")) 
        target.write("%-50s = 1.0e-4 \n"%("diffusion.mg_rtol"))
        target.write("%-50s = 1.0e-6 \n"%("diffusion.mg_atol "))
        target.write("%-50s = 1.0e-5 \n"%("temperature_diffusion.mg_rtol"))
        target.write("%-50s = 1.0e-6 \n"%("temperature_diffusion.mg_atol"))
        target.write("%-50s = 1.0e-4 \n"%("tke_diffusion.mg_rtol"))
        target.write("%-50s = 1.0e-6 \n"%("tke_diffusion.mg_atol"))

    
    def createAMR1dSolver(self):
        wind=self.yamlFile["metMastWind"]
        mol_length=self.yamlFile["molLength"]
        allowed_error=self.yamlFile["allowedError"]
        self.metMastHeight=self.yamlFile["metMastHeight"]
        self.metMastWind=[wind[0],wind[1]]
        zheight=2048
        dz=8.0
        npts=int(zheight/dz)
        num_of_steps=30000
        tolerance=1e-3
        roughness_length=self.refRoughness
        terrain_ht=0
        coriolis=self.caseCenterLat
        inv_height=1500
        inv_width=0
        inv_strength=0
        lapse_rate=0.003
        heat_flux_mode=4
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
        M=np.sqrt(wind[0]**2+wind[1]**2)
        wt=4e-3*M
        zi=1000.0
        f=2*7.27e-5*np.sin(self.caseCenterLat*np.pi/180)
        print(wt/(f*zi))
        initial_ug=wind[0]+wt/(f*zi)*wind[1]
        initial_vg=wind[1]-wt/(f*zi)*wind[0]
        print(initial_ug,initial_vg)
        include_ti=False
        initial_ug,initial_vg,z0=self.generate_profile(allowed_error,self.metMastHeight,self.metMastWind,npts,zheight,roughness_length,terrain_ht, \
                        coriolis,inv_height,inv_width,inv_strength,lapse_rate,heat_flux_mode,mol_length,num_of_steps,tolerance, \
                            initial_ug,initial_vg,include_ti)
        self.fillrans1dinfo(self.amrPrecursorFile)
        self.geostropicX=initial_ug
        self.geostropicY=initial_vg

    def fillrans1dinfo(self,target,sponge=-1):
        stringtowrite="ABL.initial_wind_profile"
        target.write("%-50s = true\n"%(stringtowrite))
        stringtowrite="ABL.rans_1dprofile_file"
        target.write('%-50s = "rans_1d.info" \n'%(stringtowrite))
        zstart=2000.0
        if(self.turbulence_model=="RANS"):
            stringtowrite="ABL.meso_sponge_start "
            target.write('%-50s = %g \n'%(stringtowrite,zstart))
        else:
            stringtowrite="ABL.meso_sponge_start "
            target.write('%-50s = %g \n'%(stringtowrite,zstart))

        # Write for AMR-Wind 
        data=np.genfromtxt(Path(self.caseParent,self.caseName,"precursor","1dSolverOutput.info").as_posix())
        zvals=data[:,0]
        uvals=data[:,1]
        vvals=data[:,2]
        wvals=data[:,3]
        tempvals=data[:,4]
        tkevals=data[:,5]
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
        newtarget=open(Path(self.caseParent,self.caseName,"precursor","rans_1d.info").as_posix(),"w")
        for i in range(0,len(zvals)):
            newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevals[i]))
        newtarget.write("%g %g %g 0 %g\n"%(self.maxZ,uvals[i],vvals[i],tkevals[i]))
        newtarget.close()
        try:
            newtarget=open(Path(self.caseParent,self.caseName,"terrain","rans_1d.info").as_posix(),"w")
        except:
            pass 
        else:
            for i in range(0,len(zvals)):
                newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevals[i]))
            newtarget.write("%g %g %g 0 %g\n"%(self.maxZ,uvals[i],vvals[i],tkevals[i]))
            newtarget.close()   
        try:
            newtarget=open(Path(self.caseParent,self.caseName,"terrainTurbine","rans_1d.info").as_posix(),"w")
        except:
            pass 
        else:
            for i in range(0,len(zvals)):
                newtarget.write("%g %g %g 0 %g\n"%(zvals[i],uvals[i],vvals[i],tkevals[i]))
            newtarget.write("%g %g %g 0 %g\n"%(self.maxZ,uvals[i],vvals[i],tkevals[i]))
            newtarget.close()    

    # Modify this function 
    def generate_profile(self,allowed_error,metMastHeight,metMastWind,npts,zheight,roughness_length,terrain_ht, \
                    coriolis,inv_height,inv_width,inv_strength,lapse_rate,heat_flux_mode,mol_length,num_of_steps,tolerance, \
                        initial_ug,initial_vg,include_ti=False):
        # Generate Geostrphic Wind for Terrain Height  Consistent Wind Speed 
        # Initial Guess of Geostropic Wind 
        # Coarse grid run to identify close wind speed 
        # A good initial guess reduces number of while loop iterations 
        residualx=100
        residualy=100 
        # Initialize Grid Npts, Height, Roughness Length and Terrain Height for IB 
        from amr1DSolver import amr1dSolver
        pathToWrite=Path(self.caseParent,self.caseName,"precursor","1dSolverOutput.info").as_posix()
        # Coarse Run 
        zheight=2048
        dz=16.0
        npts=int(zheight/dz)
        amr1D=amr1dSolver(npts,zheight,roughness_length,terrain_ht,pathToWrite)
        ug=[initial_ug,initial_vg]
        import time 
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
            met_mast_cfd_ux=np.interp(metMastHeight,z,ux)
            met_mast_cfd_uy=np.interp(metMastHeight,z,uy)
            print("Met Mast Wind:[%g %g]"%(metMastWind[0],metMastWind[1]))
            print("Specified Geostrophic Wind: [%g %g]"%(ug[0],ug[1]))
            print("CFD Met Mast Wind and Error:[%g %g]  [%g %g]"%(met_mast_cfd_ux,met_mast_cfd_uy,met_mast_cfd_ux-metMastWind[0],met_mast_cfd_uy-metMastWind[1]))
            tke=np.interp(metMastHeight,z,amr1D.tke)
            M=np.sqrt(met_mast_cfd_ux**2+met_mast_cfd_uy**2)
            print("TI:",np.sqrt(2.0/3.0*tke)/M*100)
            residualx=abs(met_mast_cfd_ux-metMastWind[0])
            residualy=abs(met_mast_cfd_uy-metMastWind[1])
            # Reduce only the higher error to  speed-up 
            if(residualx<allowed_error and residualy<allowed_error):
                print("Coarse grid converged")
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
        print("Coarse grid time:",end-start)
        z0=roughness_length

        # Fine Run 
        # A good initial guess reduces number of while loop iterations 
        print("Fine Run")
        residualx=100
        residualy=100 
        # Initialize Grid Npts, Height, Roughness Length and Terrain Height for IB 
        pathToWrite=Path(self.caseParent,self.caseName,"precursor","1dSolverOutput.info").as_posix()
        ux,uy=amr1D.return_windspeed()
        # Need to interpolate 
        # Coarse Run 
        start=time.time()
        zheight=2048
        dz=8.0
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
            met_mast_cfd_ux=np.interp(metMastHeight,z,ux)
            met_mast_cfd_uy=np.interp(metMastHeight,z,uy)
            print("Met Mast Wind:[%g %g]"%(metMastWind[0],metMastWind[1]))
            print("Specified Geostrophic Wind: [%g %g]"%(ug[0],ug[1]))
            print("CFD Met Mast Wind and Error:[%g %g]  [%g %g]"%(met_mast_cfd_ux,met_mast_cfd_uy,met_mast_cfd_ux-metMastWind[0],met_mast_cfd_uy-metMastWind[1]))
            tke=np.interp(metMastHeight,z,amr1D.tke)
            M=np.sqrt(met_mast_cfd_ux**2+met_mast_cfd_uy**2)
            print("TI:",np.sqrt(2.0/3.0*tke)/M*100)
            residualx=abs(met_mast_cfd_ux-metMastWind[0])
            residualy=abs(met_mast_cfd_uy-metMastWind[1])
            # Reduce only the higher error to  speed-up 
            if(residualx<allowed_error and residualy<allowed_error):
                print("Coarse grid converged")
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

    def metMastRefinement(self,target):
        # Read Met Mast 
        try:
            latList=self.yamlFile['metMastLat']
        except:
            return 
        lonList=self.yamlFile['metMastLon']
        metMastHeight=self.yamlFile['metMastHeight']
        print(latList,len(latList))
        xref=[]
        yref=[]
        zlower=[]
        zmast=[]
        for i in range(0,len(latList)):
            if(i==0):
                target.write("tagging.labels = metMastGrid1%g metMastGrid2%g metMastGrid3%g "%(i+1,i+1,i+1))
            else:
                target.write(" metMastGrid1%g metMastGrid2%g metMastGrid3%g "%(i+1,i+1,i+1))
            xtemp,ytemp=self.srtm.to_xy(latList[i],lonList[i])
            print(xtemp,ytemp)
            xref.append(xtemp-self.xref)
            yref.append(ytemp-self.yref)
            zmast.append(metMastHeight[i]-self.zRef)
        target.write("\n")
        print(xref,self.xref)
        print(yref,self.yref)
        # Now Write Refinemements 
        import pyvista as pv
        for i in range(0,len(latList)):
            xstart=xref[i]-4000
            ystart=yref[i]-4000
            zstart=self.interp(xref[i],yref[i])
            zdist=zmast[i]-zstart
            print(zmast[i],zstart,zdist,self.zRef)
            target.write("tagging.metMastGrid1%g.type \t\t\t = GeometryRefinement\n"%(i+1))
            target.write("tagging.metMastGrid1%g.shapes \t\t\t = metMast%g\n"%(i+1,i+1))
            target.write("tagging.metMastGrid1%g.level \t\t\t = 0\n"%(i+1))
            target.write("tagging.metMastGrid1%g.metMastGrid1%g.type \t\t\t = box\n"%(i+1,i+1))
            target.write("tagging.metMastGrid1%g.metMastGrid1%g.origin = %g %g %g \n"%(i+1,i+1,xstart,ystart,zstart-64))
            target.write("tagging.metMastGrid1%g.metMastGrid1%g.xaxis =  %g %g %g\n"%(i+1,i+1,8000,0,0))
            target.write("tagging.metMastGrid1%g.metMastGrid1%g.yaxis =  %g %g %g\n"%(i+1,i+1,0,8000,0))
            target.write("tagging.metMastGrid1%g.metMastGrid1%g.zaxis = %g %g %g\n"%(i+1,i+1,0,0,zdist+400))
            xstart=xref[i]-2000
            ystart=yref[i]-2000
            target.write("tagging.metMastGrid2%g.type \t\t\t = GeometryRefinement\n"%(i+1))
            target.write("tagging.metMastGrid2%g.shapes \t\t\t = metMastGrid2%g\n"%(i+1,i+1))
            target.write("tagging.metMastGrid2%g.min_level \t\t\t = 0\n"%(i+1))
            target.write("tagging.metMastGrid2%g.max_level \t\t\t = 1\n"%(i+1))
            target.write("tagging.metMastGrid2%g.metMastGrid2%g.type \t\t\t = box\n"%(i+1,i+1))
            target.write("tagging.metMastGrid2%g.metMastGrid2%g.origin = %g %g %g \n"%(i+1,i+1,xstart,ystart,zstart-32))
            target.write("tagging.metMastGrid2%g.metMastGrid2%g.xaxis =  %g %g %g\n"%(i+1,i+1,4000,0,0))
            target.write("tagging.metMastGrid2%g.metMastGrid2%g.yaxis =  %g %g %g\n"%(i+1,i+1,0,4000,0))
            target.write("tagging.metMastGrid2%g.metMastGrid2%g.zaxis = %g %g %g\n"%(i+1,i+1,0,0,zdist+200))
            xstart=xref[i]-1000
            ystart=yref[i]-1000
            target.write("tagging.metMastGrid3%g.type \t\t\t = GeometryRefinement\n"%(i+1))
            target.write("tagging.metMastGrid3%g.shapes \t\t\t = metMastGrid3%g\n"%(i+1,i+1))
            target.write("tagging.metMastGrid3%g.min_level \t\t\t = 0\n"%(i+1))
            target.write("tagging.metMastGrid3%g.max_level \t\t\t = 2\n"%(i+1))
            target.write("tagging.metMastGrid3%g.metMastGrid3%g.type \t\t\t = box\n"%(i+1,i+1))
            target.write("tagging.metMastGrid3%g.metMastGrid3%g.origin = %g %g %g \n"%(i+1,i+1,xstart,ystart,zstart-16))
            target.write("tagging.metMastGrid3%g.metMastGrid3%g.xaxis =  %g %g %g\n"%(i+1,i+1,2000,0,0))
            target.write("tagging.metMastGrid3%g.metMastGrid3%g.yaxis =  %g %g %g\n"%(i+1,i+1,0,2000,0))
            target.write("tagging.metMastGrid3%g.metMastGrid3%g.zaxis = %g %g %g\n"%(i+1,i+1,0,0,zdist+100))
            # Write Boxes 
            mesh1=pv.Box(bounds=(xref[i]-4000,xref[i]+4000,yref[i]-4000,yref[i]+4000,zstart-64,zstart+zdist+400))
            fileName=Path(self.caseParent,self.caseName,"metMastGrid1"+str(i+1)+".vtk").as_posix()
            mesh1.save(fileName)
            mesh1=pv.Box(bounds=(xref[i]-2000,xref[i]+2000,yref[i]-2000,yref[i]+2000,zstart-32,zstart+zdist+200))
            fileName=Path(self.caseParent,self.caseName,"metMastGrid2"+str(i+1)+".vtk").as_posix()
            mesh1.save(fileName)
            mesh1=pv.Box(bounds=(xref[i]-1000,xref[i]+1000,yref[i]-1000,yref[i]+1000,zstart-16,zstart+zdist+100))
            fileName=Path(self.caseParent,self.caseName,"metMastGrid2"+str(i+1)+".vtk").as_posix()
            mesh1.save(fileName)

    def writeRefinementRegions(self,target):
        try: 
            refinementRegions=self.yamlFile["refinementRegions"]
        except:
            pass 
        try:
            metMastRegions=self.yamlFile["metMastNames"]
        except:
            if(len(metMastRegions)==0):
                return 
        if(len(refinementRegions)==0):
            return 
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
            refinementMinz=self.yamlFile["refinementMinZ"]
        except:
            refinementMinz=0*refinementMiny
        try:
            refinementMaxz=self.yamlFile["refinementMaxZ"]
        except:
            warnings.warn("Missing maximum Z values. No refinements written")
            return 
        try:
            refinementLevels=self.yamlFile["refinementLevels"]
        except:
            warnings.warn("No refinement levels specified")
            return 
        try:
            refinementLatLon=self.yamlFile["refinementLatLon"]
        except:
            refinementLatLon=False 
        if(refinementLatLon):
            warnings.warn("Currently not implemented")
        target.write("# tagging\n")
        for i in range(len(refinementRegions)):
            if(i==0):
                target.write("%-50s = %s "%("tagging.labels",refinementRegions[i]))
            else:
                target.write(" %s "%(refinementRegions[i]))
        for i in range(len(metMastRegions)):
            if(len(refinementRegions)==0):
                target.write("%-50s = %s "%("tagging.labels",metMastRegions[i]))
            else:
                target.write(" %s "%(metMastRegions[i]))           
        target.write("\n")
        for i in range(0,len(refinementRegions)):
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
            target.write("%-50s = %g %g %g\n"%(taggingstring, \
                        refinementMinx[i],refinementMiny[i],refinementMinz[i]))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".xaxis"
            target.write("%-50s = %g 0 0\n"%(taggingstring,
                                             refinementMaxx[i]-refinementMinx[i]))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".yaxis"
            target.write("%-50s = 0 %g 0\n"%(taggingstring,
                                             refinementMaxy[i]-refinementMiny[i]))
            taggingstring="tagging."+refinementRegions[i]+".object"+str(i)+".zaxis"
            target.write("%-50s = 0 0 %g\n"%(taggingstring,
                                             refinementMaxz[i]-refinementMinz[i]))
        
        if(len(metMastRegions)==0):
            return 
        try:
            metMastX=self.yamlFile["metMastX"]
        except:
            warnings.warn("Missing X values. No refinements written")
            return    
        try:
            metMastY=self.yamlFile["metMastY"]
        except:
            warnings.warn("Missing Y values. No refinements written")
            return        
        try:
            metMastRadius=self.yamlFile["metMastRadius"]
        except:
            metMastRadius=500+0*metMastX
        try:
            metMastRefinementLevel=self.yamlFile["metMastRefinementLevel"]
        except:
            metMastRefinementLevel=2+0*metMastX

        for i in range(0,len(metMastRegions)):
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
                        metMastX[i],metMastY[i],0))
            taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".end"
            target.write("%-50s = %g %g %g\n"%(taggingstring, \
                        metMastX[i],metMastY[i],500))
            taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".outer_radius"
            target.write("%-50s = %g\n"%(taggingstring,metMastRadius[i]))
            taggingstring="tagging."+metMastRegions[i]+".metmastobject"+str(i)+".inner_radius"
            target.write("%-50s = %g\n"%(taggingstring,0.0))
        if(len(refinementLevels)>0):
            try:
                refinementLevels=self.yamlFile["refinementLevels"]
            except:
                level=0
            else: 
                level=max(refinementLevels)
        if(len(metMastRegions)>0):
            try:
                metMastRefinementLevel=self.yamlFile["metMastRefinementLevel"]
            except:
                level=max(level,2)
            else:
                level=max(level,max(metMastRefinementLevel))
        stringtowrite="amr.max_level "
        target.write("%-50s = %d\n"%(stringtowrite,level))
        self.amrPrecursorFile.write("%-50s = %d\n"%(stringtowrite,0))



    def writeAccelerationMaps(self,target):
        try:
            writeTerrainSampling=self.yamlFile["writeTerrainSampling"]
        except: 
            return 
        if(writeTerrainSampling):
            target.write("# postprocessing\n")
            target.write("%-50s = %s \n"%("incflo.post_processing","sampling"))
            target.write("%-50s = velocity temperature \n"%("sampling.fields"))
            target.write("%-50s = %s "%("sampling.labels","terrain"))
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
            offsets=self.yamlFile["verticalLevels"]
            samplingentity="sampling.terrain.type"
            target.write("%-50s = ProbeSampler\n"%(samplingentity))
            samplingentity="sampling.terrain.probe_location_file"
            target.write("%-50s = %s\n"%(samplingentity,'"terrain.csv"'))
            samplingentity="sampling.terrain.offset_vector"
            target.write("%-50s = 0 0 1\n"%(samplingentity))
            samplingentity="sampling.terrain.offsets"
            target.write("%-50s ="%(samplingentity))
            for values in offsets:
                target.write(" %g "%(values))
            target.write("\n")
        if(metMastLineSampling and len(metMastRegions)>0):
            metMastX=self.yamlFile["metMastX"]
            metMastY=self.yamlFile["metMastY"]
            for i in range(0,len(metMastRegions)):
                samplingentity="sampling."+str(metMastRegions[i])+".type"
                target.write("%-50s = LineSampler\n"%(samplingentity))
                samplingentity="sampling."+str(metMastRegions[i])+".num_points"
                target.write("%-50s = 50\n"%(samplingentity))
                samplingentity="sampling."+str(metMastRegions[i])+".start"
                zstart=0
                # Find z from terrain 
                error=10000
                for ii in range(0,len(self.terrainX1)):
                    residual=np.sqrt((metMastX[i]-self.terrainX1[ii])**2+(metMastY[i]-self.terrainX2[ii])**2)
                    if(residual<error):
                        zstart=self.terrainX3[ii]
                zstart=zstart-20
                target.write("%-50s = %g %g %g\n"%(samplingentity,metMastX[i],metMastY[i],zstart))
                samplingentity="sampling."+str(metMastRegions[i])+".end"
                target.write("%-50s = %g %g %g\n"%(samplingentity,metMastX[i],metMastY[i],zstart+200))



    def closeAMRFiles(self):
        self.amrPrecursorFile.close()
        try:
            self.amrTerrainFile.close()
        except:
            pass

    def writeTerrainData(self,folder):
        print("Writing Terrain Data")
        x1=self.terrainX1.flatten(order='F')
        x2=self.terrainX2.flatten(order='F')
        x3=self.terrainX3.flatten(order='F')  
        x=np.arange(np.amin(x1),np.amax(x1),self.caseCellSize)
        y=np.arange(np.amin(x2),np.amax(x2),self.caseCellSize)
        #x=np.arange(np.amin(x1),np.amax(x1),self.caseCellSize)
        #y=np.arange(np.amin(x2),np.amax(x2),self.caseCellSize)
        from scipy.interpolate import NearestNDInterpolator
        self.interp = NearestNDInterpolator(list(zip(x1,x2)),x3) 
        xterrain,yterrain=np.meshgrid(x,y)
        zterrain = self.interp(xterrain,yterrain)    
        import matplotlib.pylab as plt 
        plt.contourf(xterrain,yterrain,zterrain)
        x1=xterrain.flatten(order='F')
        x2=yterrain.flatten(order='F')
        x3=zterrain.flatten(order='F')
        print("Shape:",xterrain.shape)
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
        x=np.arange(np.amin(x1),np.amax(x1),2*self.caseCellSize)
        y=np.arange(np.amin(x2),np.amax(x2),2*self.caseCellSize)
        xterrain,yterrain=np.meshgrid(x,y)
        zterrain = self.interp(xterrain,yterrain)  
        x1=xterrain.flatten(order='F')
        x2=yterrain.flatten(order='F')
        x3=zterrain.flatten(order='F')
        target=Path(self.caseParent,self.caseName,folder,"terrain.csv").open("w")
        target.write("%d \n"%(len(x1)))
        for i in range(0,len(x1)):
             target.write("%g %g %g\n"%(x1[i],x2[i],x3[i]))
        target.close()
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
        print("Writing Roughness Data")
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



    def writeRestart(self,target):
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
