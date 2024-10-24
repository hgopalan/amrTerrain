'''
!----------------------------------------------------------------!
AMR - Wind back-end for cases with terrain                       !
!----------------------------------------------------------------!
'''

import yaml
from pathlib import Path
from terrain import SRTM
import numpy as np 


class amrBackend():
    def __init__(self,yamlFile):
        self.yamlFilePath = Path(yamlFile)
        self.yamlFile=yaml.safe_load(self.yamlFilePath.open())
        self.caseCellSize=128 
        self.caseverticalAR=8
        self.turbulence_model='RANS'
        self.case_end_time=14400 
        self.plotOutput=3600
        self.restartOutput=3600
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
                self.caseverticalAR=8
            else:
                self.caseverticalAR=4


    def createDomain(self):
        self.write_stl=self.yamlFile["writeTerrain"]
        import SRTM_to_STL_example as converter
        self.xref,self.yref,self.zRef,self.srtm=converter.SRTM_Converter(Path(self.caseParent,self.caseName).as_posix(),self.caseCenterLat,self.caseCenterLon,self.refHeight, \
                                                    self.caseWest,self.caseEast,self.caseSouth,self.caseNorth,self.write_stl)
        stlFile=Path(self.caseParent,self.caseName,"terrain.stl").as_posix()
        import pyvista as pv 
        mesh=pv.read(stlFile)
        x1=mesh.points[:,0]
        x2=mesh.points[:,1]
        x3=mesh.points[:,2]
        self.terrainX1=x1[:]
        self.terrainX2=x2[:]
        self.terrainX3=x3[:]
        if(not self.write_stl):
            Path(self.caseParent,self.caseName,"terrain.stl").unlink()
    
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
        self.createSolverInfo(self.amrPrecursorFile)
        self.createAMRTransport(self.amrPrecursorFile)
        self.createAMRTurbulence(self.amrPrecursorFile)
        self.createAMRABLData(self.amrPrecursorFile,0,1) 
        self.createAMRSourceTerm(self.amrPrecursorFile)    
        self.createAMRBC(self.amrPrecursorFile)
        self.createAMRTolerance(self.amrPrecursorFile)
        print(" Done creating precursor")

    def createTerrainFiles(self,folder):
        print("Creating Terrain Files")
        self.createAMRGeometry(self.amrTerrainFile,-1)
        self.createAMRGrid(self.amrTerrainFile)
        self.createAMRTime(self.amrTerrainFile)
        if(self.caseType=='terrainTurbine'):
            self.createSolverInfo(self.amrTerrainFile,1,1)
        else:
            self.createSolverInfo(self.amrTerrainFile,1)
        self.createAMRTransport(self.amrTerrainFile)
        self.createAMRTurbulence(self.amrTerrainFile)
        self.createAMRABLData(self.amrTerrainFile,1,0)  
        if(self.caseType=='terrainTurbine'):
            self.createAMRSourceTerm(self.amrTerrainFile,1,1) 
        else:
            self.createAMRSourceTerm(self.amrTerrainFile,1)    
        self.createAMRBC(self.amrTerrainFile,1)
        self.createAMRTolerance(self.amrTerrainFile,1)
        self.writeTerrainData(folder)
        self.writeRestart(self.amrTerrainFile)
        self.createAMRPrecursorSampling(self.amrPrecursorFile)
        # Terrain Monitoring and Refining 
        self.metMastRefinement(self.amrTerrainFile)
        #self.metMastMonitoring(self.amrTerrainFile)

    def createAMRGeometry(self,target,periodic=-1):
        target.write("# Geometry\n")
        minX=np.amin(self.terrainX1)
        minY=np.amin(self.terrainX2)
        minZ=np.amin(self.terrainX3)
        maxX=np.amax(self.terrainX1)
        maxY=np.amax(self.terrainX2)
        self.terrainZMax=np.amax(self.terrainX3)
        # Add 1 km for ABL and 4 km for Rayleigh
        self.ABLHeight=1000
        self.RDLHeight=2000
        self.maxZ=self.terrainZMax+self.ABLHeight+self.RDLHeight
        self.maxZ=round(self.maxZ,-3)
        target.write("geometry.prob_lo \t\t\t = %g %g %g \n"%(minX,minY,minZ))
        target.write("geometry.prob_hi \t\t\t = %g %g %g \n"%(maxX,maxY,self.maxZ))
        if(periodic==1):
            target.write("geometry.is_periodic \t\t\t = 1 1 0\n")
        else:
            target.write("geometry.is_periodic \t\t\t = 0 0 0\n")

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
        target.write("amr.n_cell \t\t\t = %g %g %g\n"%(nx,ny,nz))
        target.write("amr.max_level \t\t\t = 0\n")

    def createAMRTime(self,target):
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
                self.case_end_time=14400
        try:
            self.plotOutput=self.yamlFile['plotOutput']
        except:
            self.plotOutput=-1
        try:
            self.restartOutput=self.yamlFile['restartOutput']
        except:
            self.restartOutput=3600
        if(self.timeMethod=="step"):
            target.write("time.stop_time \t\t\t = -1\n")
            target.write("time.max_step \t\t\t = %g\n"%(self.timeSteps))
        else:
            target.write("time.stop_time \t\t\t = %g\n"%(self.case_end_time))
        target.write("time.initial_dt \t\t\t = 1.0\n")
        target.write("time.fixed_dt \t\t\t = -1\n")
        if(self.turbulence_model=="RANS"):
            target.write("time.cfl \t\t\t = 0.5\n")  
        else:
            target.write("time.cfl \t\t\t = 0.9\n")
        target.write('time.plot_interval \t\t\t = %g\n'%(self.plotOutput))
        target.write("time.checkpoint_interval \t\t\t = %g\n"%(self.restartOutput))

    def createSolverInfo(self,target,terrain=-1,turbine=-1):
        self.caseWindspeedX=self.yamlFile['windX']
        self.caseWindspeedY=self.yamlFile['windY']
        self.caseWindspeedZ=self.yamlFile['windZ']   
        try:
            self.fastBoxes=self.yamlFile['fastBoxes']      
        except:
            self.fastBoxes=False
        target.write("# incflo \n")
        if(terrain==1 and turbine==1):
            if(self.fastBoxes):
                target.write("incflo.physics \t\t\t = ABL TerrainDrag \n")
            else:
                target.write("incflo.physics \t\t\t = ABL TerrainDrag Actuator\n")
        elif(terrain==1):
            target.write("incflo.physics \t\t\t = ABL TerrainDrag\n")            
        elif(turbine==1):
            target.write("incflo.physics \t\t\t = ABL Actuator\n")
        else:
            target.write("incflo.physics \t\t\t = ABL\n")
        target.write("incflo.density \t\t\t = 1.225\n")
        target.write("incflo.gravity \t\t\t = 0.  0. -9.81  # Gravitational force (3D)\n")
        target.write("incflo.velocity \t\t\t = %g %g %g \n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        target.write("incflo.verbose  \t\t\t = 0\n")
        target.write("incflo.initial_iterations \t\t\t = 8\n")
        target.write("incflo.do_initial_proj \t\t\t = true\n")
        target.write("incflo.constant_density \t\t\t = true\n")
        target.write("incflo.use_godunov \t\t\t = true\n")
        print(self.turbulence_model)
        if(self.turbulence_model=="RANS"):
            target.write('incflo.godunov_type \t\t\t = "ppm"\n')
        else:
            target.write('incflo.godunov_type \t\t\t = "weno_z"\n')
        target.write("incflo.diffusion_type \t\t\t = 2\n")

    def createAMRTransport(self,target):
        target.write("# transport equation parameters \n")
        target.write("transport.model \t\t\t = ConstTransport\n")
        target.write("transport.viscosity \t\t\t = 1e-5\n")
        target.write("transport.laminar_prandtl \t\t\t = 0.7\n")
        target.write("transport.turbulent_prandtl \t\t\t = 0.333\n")     

    def createAMRTurbulence(self,target):
        # Default option is to do LES 
        # RANS capability added if requested 
        if(self.turbulence_model=="RANS"):
            target.write("# turbulence equation parameters\n")
            target.write("turbulence.model \t\t\t = OneEqRANS\n")
            target.write("TKE.source_terms \t\t\t = Krans\n")
        else:
            target.write("# turbulence equation parameters \n")
            target.write("turbulence.model \t\t\t = Kosovic\n")
            target.write("Kosovic.refMOL \t\t\t = -1e30\n")

    def createAMRABLData(self,target,iomode=-1,fluctuations=1):
        self.refTemperature=self.yamlFile["refTemperature"]
        self.refRoughness=self.yamlFile["refRoughness"]
        self.refHeatFlux=self.yamlFile["refHeatflux"]
        target.write("# Atmospheric boundary layer\n")
        if(fluctuations==1 and (not self.turbulence_model=='RANS')):
            UPeriod=int((np.amax(self.terrainX1)-np.amin(self.terrainX1))/200)
            VPeriod=int((np.amax(self.terrainX1)-np.amin(self.terrainX1))/200)
            target.write("ABL.Uperiods \t\t\t = %g\n"%(UPeriod))
            target.write("ABL.Vperiods \t\t\t = %g\n"%(VPeriod))
            target.write("ABL.cutoff_height \t\t\t = 50.0\n")
            target.write("ABL.deltaU \t\t\t = 1.0\n")
            target.write("ABL.deltaV \t\t\t = 1.0\n") 
            target.write("ABL.perturb_ref_height \t\t\t = 50.0\n")
            target.write("ABL.perturb_velocity \t\t\t = true\n")
            target.write("ABL.perturb_temperature \t\t\t = false\n")
        else:
            target.write("ABL.perturb_velocity \t\t\t = false\n")
            target.write("ABL.perturb_temperature \t\t\t = false\n")
        target.write("ABL.kappa \t\t\t = .41\n")
        target.write("ABL.normal_direction \t\t\t = 2\n")
        target.write("ABL.reference_temperature \t\t\t = %g\n"%(self.refTemperature))
        target.write("ABL.stats_output_format \t\t\t = netcdf\n")
        target.write("ABL.surface_roughness_z0 \t\t\t = %g\n"%(self.refRoughness))
        # Write Heights 
        inversionHeight=round(self.terrainZMax,-3)+1000
        inversionLayerThickness=round(self.terrainZMax,-3)+1000+100
        lapseRate=0.003 
        target.write("ABL.temperature_heights = %g %g %g %g %g \n"%(0.0,round(self.terrainZMax,-3),inversionHeight,inversionLayerThickness,self.maxZ)) 
        TRef=300
        print(self.maxZ-inversionLayerThickness)
        target.write("ABL.temperature_values  = %g %g %g %g %g "%(TRef,TRef,TRef,TRef+5,TRef+5+lapseRate*(self.maxZ-inversionLayerThickness)))    
        target.write("\n")
        target.write("ABL.wall_shear_stress_type \t\t\t = local\n")
        target.write("ABL.surface_temp_flux \t\t\t = %g\n"%(self.refHeatFlux))
        if(iomode==0):
            target.write('ABL.bndry_file \t\t\t = "bndry_files"\n')
            target.write("ABL.bndry_write_frequency \t\t\t = 1\n")
            target.write("ABL.bndry_io_mode \t\t\t = 0\n")
            if(self.caseWindspeedX>=0 and self.caseWindspeedY>=0):
                target.write("ABL.bndry_planes \t\t\t = xlo ylo \n")     
            elif(self.caseWindspeedX>=0 and self.caseWindspeedY<0):
                target.write("ABL.bndry_planes \t\t\t = xlo yhi \n")      
            if(self.caseWindspeedX<0 and self.caseWindspeedY>=0):
                target.write("ABL.bndry_planes \t\t\t = xhi ylo \n")     
            elif(self.caseWindspeedX<0 and self.caseWindspeedY<0):
                target.write("ABL.bndry_planes \t\t\t = xhi yhi \n")    
            # Guessing a start time far enough 
            M=np.sqrt(self.caseWindspeedX**2+self.caseWindspeedY**2)
            dt=0.9*self.caseCellSize/M
            startTime=max(600,0.25*(self.timeSteps/dt))
            target.write("ABL.bndry_output_start_time \t\t\t = %g\n"%(startTime))
            target.write("ABL.bndry_var_names \t\t\t = velocity temperature\n")
            target.write("ABL.bndry_output_format \t\t\t = native\n")
        elif(iomode==1):
            target.write('ABL.bndry_file \t\t\t = "../precursor/bndry_files"\n')
            target.write("ABL.bndry_io_mode \t\t\t = 1\n")
            target.write("ABL.bndry_var_names \t\t\t = velocity temperature\n")
            target.write("ABL.bndry_output_format \t\t\t = native\n")

    def createAMRSourceTerm(self,target,terrain=-1,turbine=-1):
        target.write("# Source\n")
        #refLat=self.yamlFile["refLat"]
        #refPeriod=self.yamlFile["refPeriod"]
        try: 
            self.includeCoriolis=self.yamlFile["includeCoriolis"]
        except:
            self.includeCoriolis=False 
        if(self.turbulence_model=="RANS"):
            if(self.includeCoriolis):
                try:
                    self.forcingHeight=self.yamlFile["forcingHeight"]
                except: 
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping \n")
                    target.write("GeostrophicForcing.geostrophic_wind \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
                else:
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing ABLForcing RayleighDamping DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing ABLForcing RayleighDamping \n")
                    target.write("ABLForcing.abl_forcing_height \t\t\t = %g \n"%(self.forcingHeight))
            else:
                try:
                    self.forcingHeight=self.yamlFile["forcingHeight"]
                except: 
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping  DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy GeostrophicForcing RayleighDamping DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy GeostrophicForcing RayleighDamping \n")
                    target.write("GeostrophicForcing.geostrophic_wind \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
                else:
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy ABLForcing RayleighDamping DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy ABLForcing RayleighDamping \n")
                    target.write("ABLForcing.abl_forcing_height \t\t\t = %g \n"%(self.forcingHeight))
        else:
            if(self.includeCoriolis):
                try:
                    self.forcingHeight=self.yamlFile["forcingHeight"]
                except: 
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm\n")
                    target.write("GeostrophicForcing.geostrophic_wind \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
                else:
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing ABLForcing RayleighDamping NonLinearSGSTerm DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing ABLForcing RayleighDamping NonLinearSGSTerm\n")
                    target.write("ABLForcing.abl_forcing_height \t\t\t = %g \n"%(self.forcingHeight))
            else:
                try:
                    self.forcingHeight=self.yamlFile["forcingHeight"]
                except: 
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy GeostrophicForcing RayleighDamping NonLinearSGSTerm\n")
                    target.write("GeostrophicForcing.geostrophic_wind \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
                else:
                    if(terrain==1 and turbine==1):
                        if(self.fastBoxes):
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing  \n")
                        else:
                            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing ActuatorForcing \n")
                    elif(terrain==1):    
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy ABLForcing RayleighDamping NonLinearSGSTerm DragForcing\n")
                    else:
                        target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy ABLForcing RayleighDamping NonLinearSGSTerm\n")
                    target.write("ABLForcing.abl_forcing_height \t\t\t = %g \n"%(self.forcingHeight))
        if(terrain==1 or turbine==1):
            target.write("Temperature.source_terms = DragTempForcing\n")
        target.write("RayleighDamping.force_coord_directions= 0 0 1\n")
        target.write("BoussinesqBuoyancy.reference_temperature \t\t\t = %g\n"%(self.refTemperature))
        target.write("BoussinesqBuoyancy.thermal_expansion_coeff \t\t\t = %g\n"%(1.0/self.refTemperature))
        #if(self.includeCoriolis):
        # Write the coriolis term for geostrophic forcing term 
        target.write("CoriolisForcing.east_vector \t\t\t = 1.0 0.0 0.0 \n")
        target.write("CoriolisForcing.north_vector \t\t\t = 0.0 1.0 0.0 \n")
        target.write("CoriolisForcing.latitude \t\t\t = %g \n"%(self.caseCenterLat))
        target.write("RayleighDamping.reference_velocity \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        startRayleigh=self.maxZ-self.RDLHeight
        target.write("RayleighDamping.length_sloped_damping \t\t\t = %g\n"%(500))
        target.write("RayleighDamping.length_complete_damping \t\t\t = %g\n"%(self.maxZ-startRayleigh-500))
        target.write("RayleighDamping.time_scale \t\t\t = 20.0\n")     

    def createAMRBC(self,target,inflowOutflow=-1):
        target.write("# BC \n")
        if(inflowOutflow==1):
            if(self.caseWindspeedX>=0):
                target.write('xlo.type \t\t\t = "mass_inflow"\n')
                target.write("xlo.density \t\t\t = 1.225\n")
                target.write("xlo.temperature \t\t\t = 300\n")
                target.write('xhi.type \t\t\t = "pressure_outflow"\n')
            else:
                target.write('xhi.type \t\t\t = "mass_inflow"\n')
                target.write("xhi.density \t\t\t = 1.225\n")
                target.write("xhi.temperature \t\t\t = 300\n")
                target.write('xlo.type \t\t\t = "pressure_outflow"\n')   
            if(self.caseWindspeedY>=0):
                target.write('ylo.type \t\t\t = "mass_inflow"\n')
                target.write("ylo.density \t\t\t = 1.225\n")
                target.write("ylo.temperature \t\t\t = 300\n")
                target.write('yhi.type \t\t\t = "pressure_outflow"\n')
            else:
                target.write('yhi.type \t\t\t = "mass_inflow"\n')
                target.write("yhi.density \t\t\t = 1.225\n")
                target.write("yhi.temperature \t\t\t = 300\n")
                target.write('ylo.type \t\t\t = "pressure_outflow"\n')  
        target.write('zhi.type \t\t\t = "slip_wall"\n')
        target.write('zhi.temperature_type \t\t\t = "fixed_gradient"\n')
        target.write("zhi.temperature \t\t\t =  0.003\n")
        target.write('zlo.type \t\t\t = "wall_model"\n')

    def createAMRTolerance(self,target,modify=-1):
        #if(modify==1):
        if(self.caseverticalAR==3 or self.caseverticalAR==4):
            self.smoothing=8
        elif(self.caseverticalAR>4 and self.caseverticalAR<=8):
            self.smoothing=32
        elif(self.caseverticalAR>8 and self.caseverticalAR<=16):
            self.smoothing=64
        if(self.caseverticalAR>=3):
            target.write("mac_proj.num_pre_smooth \t\t\t = %g \n"%(self.smoothing))
            target.write("mac_proj.num_post_smooth \t\t\t = %g \n"%(self.smoothing))
        target.write("mac_proj.mg_rtol \t\t\t = 1.0e-4 \n")
        target.write("mac_proj.mg_atol \t\t\t = 1.0e-6 \n")
        target.write("mac_proj.maxiter \t\t\t = 360 \n")
        if(self.caseverticalAR>=3):
            target.write("nodal_proj.num_pre_smooth \t\t\t = %g \n"%(self.smoothing))
            target.write("nodal_proj.num_post_smooth \t\t\t = %g \n"%(self.smoothing))
        target.write("nodal_proj.mg_rtol \t\t\t = 1.0e-4 \n")
        target.write("nodal_proj.mg_atol \t\t\t = 1.0e-6 \n")
        target.write("diffusion.mg_rtol \t\t\t = 1.0e-6 \n")
        target.write("diffusion.mg_atol \t\t\t = 1.0e-8 \n")
        target.write("temperature_diffusion.mg_rtol \t\t\t = 1.0e-6 \n")
        target.write("temperature_diffusion.mg_atol \t\t\t = 1.0e-8 \n")
        target.write("nodal_proj.maxiter \t\t\t = 360 \n")
    
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
        x=np.arange(np.amin(x1),np.amax(x1),0.25*self.caseCellSize)
        y=np.arange(np.amin(x2),np.amax(x2),0.25*self.caseCellSize)
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
        target=Path(self.caseParent,self.caseName,folder,"terrain.amrwind").open("w")
        for i in range(0,len(x1)):
             target.write("%g %g %g\n"%(x1[i],x2[i],x3[i]))
        target.close()
        if(self.write_stl):
            data=np.column_stack([x1,x2,x3])
            import pyvista as pv
            mesh=pv.PolyData(data)
            mesh['elevation']=data[:,2]
            mesh.save(Path(self.caseParent,self.caseName,folder,"terrainPoints.vtk").as_posix())
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
            import math 
            # Guessing a start time far enough 
            M=np.sqrt(self.caseWindspeedX**2+self.caseWindspeedY**2)
            dt=0.9*self.caseCellSize/M
            startTime=0.25*(self.timeSteps/dt)
            startStep=int(startTime/dt)
            print("startStep:",startStep)
            if(startStep<1000):
                startStep=int(math.ceil(startStep/ 100.0)) * 100
            elif(startStep<10000):
                startStep=int(math.ceil(startStep/ 1000.0)) * 1000
            if(startStep<self.restartOutput):
                startStep=self.restartOutput
            else:
                ratio=math.ceil(startStep/self.restartOutput)
                startStep=ratio*self.restartOutput
            if(startStep<1000):
                startStep="00"+str(startStep)
            elif(startStep<10000):
                startStep="0"+str(startStep)
            target.write("#io \n")
            target.write('io.restart_file \t\t\t = "../precursor/chk%s"\n'%(startStep))

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
