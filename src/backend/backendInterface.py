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
        try:
            tmpName=Path(self.caseParent,self.caseName,"terrain.tif")
            tmpName.unlink()
            print("Deleting file")
        except:
            pass 
        srtm_output = Path(self.caseParent,self.caseName,"terrain.tif")
        self.srtm = SRTM(bounds,fpath=srtm_output.as_posix(),product=product)
        self.srtm.download()
        x1,x2,x3 = self.srtm.to_terrain(dx,dy)
        self.purgeCorner=self.yamlFile['purgeCorners']
        cornerCut=self.purgeCorner
        self.terrainX1=x1[cornerCut:x1.shape[0]-cornerCut,cornerCut:x1.shape[1]-cornerCut]
        self.terrainX2=x2[cornerCut:x2.shape[0]-cornerCut,cornerCut:x2.shape[1]-cornerCut]
        self.terrainX3=x3[cornerCut:x3.shape[0]-cornerCut,cornerCut:x3.shape[1]-cornerCut]
        self.caseLatList=0*self.terrainX1
        self.caseLonList=0*self.terrainX2
        print(self.terrainX1.shape)
        meanZ3=np.mean(self.terrainX3.flatten(order='F'))
        print("Mean Height:",meanZ3)
        minX=1.2*np.amin(self.terrainX1)
        minY=1.2*np.amin(self.terrainX2)
        maxX=0.8*np.amax(self.terrainX1)
        maxY=0.8*np.amax(self.terrainX2)
        # Checking for Spurios Zero 
        for i in range(0,self.terrainX1.shape[0]):
            for j in range(0,self.terrainX1.shape[1]):
                if(self.terrainX3[i,j]==0 and  self.terrainX1[i,j]> maxX):
                    self.terrainX3[i,j]=self.terrainX3[i-1,j]
                    print(i,j,self.terrainX3[i,j])
        for i in range(self.terrainX3.shape[0]-1,0,-1):
            for j in range(0,self.terrainX1.shape[1]):
                if(self.terrainX3[i,j]==0 and  self.terrainX1[i,j]< minX):
                    self.terrainX3[i,j]=self.terrainX3[i+1,j]
                    print(i,j,self.terrainX3[i,j])
        print('Done Fixing')
        for i in range(0,self.terrainX1.shape[0]):
            for j in range(0,self.terrainX1.shape[1]):
                self.caseLatList[i,j],self.caseLonList[i,j]=self.srtm.to_latlon(self.terrainX1[i,j],self.terrainX2[i,j])
        # import pyvista as pv
        # import numpy as np 
        # x1=self.terrainX1.flatten(order='F')
        # x2=self.terrainX2.flatten(order='F')
        # x3=self.terrainX3.flatten(order='F')
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
        self.createPrecursorFiles()
        # self.createAMRGeometry(self.amrPrecursorFile,1)
        # self.createAMRGeometry(self.amrTerrainFile,-1)
        # self.createAMRGrid(self.amrPrecursorFile)
        # self.createAMRGrid(self.amrTerrainFile)
        # self.createAMRTime(self.amrPrecursorFile)
        # self.createAMRTime(self.amrTerrainFile)
        # self.createSolverInfo(self.amrPrecursorFile)
        # self.createSolverInfo(self.amrTerrainFile,1)
        # self.createAMRTransport(self.amrPrecursorFile)
        # self.createAMRTransport(self.amrTerrainFile)
        # self.createAMRTurbulence(self.amrPrecursorFile)
        # self.createAMRTurbulence(self.amrTerrainFile)
        # self.createAMRABLData(self.amrPrecursorFile,0,1)
        # self.createAMRABLData(self.amrTerrainFile,1,0)  
        # self.createAMRSourceTerm(self.amrPrecursorFile)
        # self.createAMRSourceTerm(self.amrTerrainFile,1)      
        # self.createAMRBC(self.amrPrecursorFile)
        # self.createAMRBC(self.amrTerrainFile,1)
        # self.createAMRTolerance(self.amrPrecursorFile)
        # self.createAMRTolerance(self.amrTerrainFile)
        # # Write the terrain data 
        # if(self.caseType=='terrain'):
        #     self.writeTerrainData()
        #     self.writeRestart(self.amrTerrainFile)
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


    def createTerrainFiles(self,folder):
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

    def createAMRGeometry(self,target,periodic=-1):
        target.write("# Geometry\n")
        minX=np.amin(self.terrainX1)
        minY=np.amin(self.terrainX2)
        minZ=np.amin(self.terrainX3)
        maxX=np.amax(self.terrainX1)
        maxY=np.amax(self.terrainX2)
        maxZ=max(2*np.amax(self.terrainX3),1600)
        target.write("geometry.prob_lo \t\t\t = %g %g %g \n"%(minX,minY,minZ))
        target.write("geometry.prob_hi \t\t\t = %g %g %g \n"%(maxX,maxY,maxZ))
        if(periodic==1):
            target.write("geometry.is_periodic \t\t\t = 1 1 0\n")
        else:
            target.write("geometry.is_periodic \t\t\t = 0 0 0\n")

    def createAMRGrid(self,target):
        self.caseCellSize=self.yamlFile['cellSize']
        nx=int((np.amax(self.terrainX1)-np.amin(self.terrainX1))/self.caseCellSize)
        while (nx%8 !=0):
            nx=nx+1
        ny=int((np.amax(self.terrainX2)-np.amin(self.terrainX2))/self.caseCellSize)
        while (ny%8 !=0):
            ny=ny+1
        # AMR - Wind cannot handle higher AR 
        nz=2*int((max(2*np.amax(self.terrainX3),1600)-np.amin(self.terrainX3))/self.caseCellSize)
        while (nz%8 !=0):
            nz=nz+1
        target.write("# Grid \n")
        target.write("amr.n_cell \t\t\t = %g %g %g\n"%(nx,ny,nz))
        target.write("amr.max_level \t\t\t = 0\n")

    def createAMRTime(self,target):
        self.timeMethod=self.yamlFile['timeMethod']
        if(self.timeMethod=="step"):
            self.timeSteps=self.yamlFile["numOfSteps"]
        else:
            pass # Need to add parameters for physical time 
        self.plotOutput=self.yamlFile['plotOutput']
        self.restartOutput=self.yamlFile['restartOutput']
        target.write("time.stop_time \t\t\t = -1\n")
        target.write("time.max_step \t\t\t = %g\n"%(self.timeSteps))
        target.write("time.initial_dt \t\t\t = 0.1\n")
        target.write("time.fixed_dt \t\t\t = -1\n")
        target.write("time.cfl \t\t\t = 0.9\n")
        target.write('time.plot_interval \t\t\t = %g\n'%(self.plotOutput))
        target.write("time.checkpoint_interval \t\t\t = %g\n"%(self.restartOutput))

    def createSolverInfo(self,target,terrain=-1,turbine=-1):
        self.caseWindspeedX=self.yamlFile['windX']
        self.caseWindspeedY=self.yamlFile['windY']
        self.caseWindspeedZ=self.yamlFile['windZ']                
        target.write("# incflo \n")
        if(terrain==1 and turbine==1):
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
        target.write('incflo.godunov_type \t\t\t = "weno_z"\n')
        target.write("incflo.diffusion_type \t\t\t = 2\n")

    def createAMRTransport(self,target):
        target.write("# transport equation parameters \n")
        target.write("transport.model \t\t\t = ConstTransport\n")
        target.write("transport.viscosity \t\t\t = 1e-5\n")
        target.write("transport.laminar_prandtl \t\t\t = 0.7\n")
        target.write("transport.turbulent_prandtl \t\t\t = 0.333\n")     

    def createAMRTurbulence(self,target):
        target.write("# turbulence equation parameters \n")
        target.write("turbulence.model \t\t\t = Kosovic\n")
        target.write("Kosovic.refMOL \t\t\t = -1e30\n")

    def createAMRABLData(self,target,iomode=-1,fluctuations=1):
        self.refTemperature=self.yamlFile["refTemperature"]
        self.refRoughness=self.yamlFile["refRoughness"]
        self.refHeatFlux=self.yamlFile["refHeatflux"]
        target.write("# Atmospheric boundary layer\n")
        if(fluctuations==1):
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
        target.write("ABL.kappa \t\t\t = .41\n")
        target.write("ABL.normal_direction \t\t\t = 2\n")
        target.write("ABL.reference_temperature \t\t\t = %g\n"%(self.refTemperature))
        target.write("ABL.stats_output_format \t\t\t = netcdf\n")
        target.write("ABL.surface_roughness_z0 \t\t\t = %g\n"%(self.refRoughness))
        target.write("ABL.temperature_heights \t\t\t = 0  800 900 1900\n")
        target.write("ABL.temperature_values  \t\t\t = 300 300 308 311\n")
        target.write("ABL.wall_shear_stress_type \t\t\t = local\n")
        target.write("ABL.surface_temp_flux \t\t\t = %g\n"%(self.refHeatFlux))
        if(iomode==0):
            target.write('ABL.bndry_file \t\t\t = "bndry_files"\n')
            target.write("ABL.bndry_write_frequency \t\t\t = 100\n")
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
            startTime=0.25*(self.timeSteps/dt)
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
        refLat=self.yamlFile["refLat"]
        refPeriod=self.yamlFile["refPeriod"]
        if(terrain==1 and turbine==1):
            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing ActuatorForcing \n")
        elif(terrain==1):    
            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm DragForcing\n")
        else:
            target.write("ICNS.source_terms \t\t\t = BoussinesqBuoyancy CoriolisForcing GeostrophicForcing RayleighDamping NonLinearSGSTerm\n")
        target.write("BoussinesqBuoyancy.reference_temperature \t\t\t = %g\n"%(self.refTemperature))
        target.write("BoussinesqBuoyancy.thermal_expansion_coeff \t\t\t = %g\n"%(1.0/self.refTemperature))
        target.write("CoriolisForcing.east_vector \t\t\t = 1.0 0.0 0.0 \n")
        target.write("CoriolisForcing.north_vector \t\t\t = 0.0 1.0 0.0 \n")
        target.write("CoriolisForcing.latitude \t\t\t = %g \n"%(refLat))
        target.write("CoriolisForcing.rotational_time_period \t\t\t = %g \n"%(refPeriod))
        target.write("GeostrophicForcing.geostrophic_wind \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        target.write("RayleighDamping.reference_velocity \t\t\t = %g %g %g\n"%(self.caseWindspeedX,self.caseWindspeedY,self.caseWindspeedZ))
        target.write("RayleighDamping.length_sloped_damping \t\t\t = 400\n")
        target.write("RayleighDamping.length_complete_damping \t\t\t = 200\n")
        target.write("RayleighDamping.time_scale \t\t\t = 5.0\n")     

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
        if(modify==1):
            target.write("# Tolerance \n")
            target.write("nodal_proj.mg_rtol \t\t\t = 1.0e-3\n")
            target.write("nodal_proj.maxiter \t\t\t = 360\n")
            target.write("nodal_proj.mg_atol = 1.0e-6\n")
            target.write("mac_proj.mg_rtol = 1.0e-3\n")
            target.write("mac_proj.mg_atol = 1.0e-6\n") 
            target.write("mac_proj.maxiter \t\t\t = 360\n")
        else:
            target.write("# Tolerance \n")
            target.write("nodal_proj.mg_rtol \t\t\t = 1.0e-4\n")
            target.write("nodal_proj.mg_atol = 1.0e-6\n")
            target.write("mac_proj.mg_rtol = 1.0e-4\n")
            target.write("mac_proj.mg_atol = 1.0e-6\n")    
    
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

    def closeAMRFiles(self):
        self.amrPrecursorFile.close()
        try:
            self.amrTerrainFile.close()
        except:
            pass

    def writeTerrainData(self,folder):
        target=Path(self.caseParent,self.caseName,folder,"terrain.amrwind").open("w")
        x1=self.terrainX1.flatten(order='F')
        x2=self.terrainX2.flatten(order='F')
        x3=self.terrainX3.flatten(order='F')
        xmin=min(x1)
        xmax=max(x1)
        ymin=min(x2)
        ymax=max(x2)
        xleft=xmin+1200
        xright=xmax-1200
        yleft=ymin+1200
        yright=ymax-1200
        zmean=np.mean(x3)
        # for i in range(0,len(x1)):
        #     if(x1[i]>xleft and x1[i]<xright and x2[i]>yleft and x2[i]<yright):
        #         pass            
        #     else:
        #         x3[i]=zmean    
        data=np.column_stack([x1,x2,x3])
        import pyvista as pv
        mesh=pv.PolyData(data)
        mesh['elevation']=data[:,2]
        surf = mesh.delaunay_2d()
        totalPoints=100
        decimateFactor=0.95
        while(totalPoints<100000):
            smoothData=surf.decimate(decimateFactor)
            totalPoints=len(smoothData.points[:,0])
            if(totalPoints<100000):
                decimateFactor=decimateFactor-0.05
            print(decimateFactor,totalPoints)
        #smoothData=surf.smooth(1008)  
        print("Reduced points from ",len(x1)," to ",len(smoothData.points[:,0]))
        surf.save(Path(self.caseParent,self.caseName,folder,"terrain.vtk").as_posix())
        target=Path(self.caseParent,self.caseName,folder,"terrain.amrwind").open("w")
        self.smoothTerrainX1=smoothData.points[:,0]
        self.smoothTerrainX2=smoothData.points[:,1]
        self.smoothTerrainX3=smoothData.points[:,2]
        for i in range(0,len(smoothData.points[:,0])):
             target.write("%g %g %g\n"%(smoothData.points[i,0],smoothData.points[i,1],smoothData.points[i,2]))
        target.close()
        smoothData.save(Path(self.caseParent,self.caseName,folder,"smoothterrain.vtk").as_posix())
        # target=Path(self.caseParent,self.caseName,"terrain","terrain.amrwind").open("w")
        # for i in range(0,len(x1)):
        #      target.write("%g %g %g\n"%(smoothData.points[i,0],smoothData.points[i,1],smoothData.points[i,2]))
        # target.close()
        # smoothData.save(Path(self.caseParent,self.caseName,"terrain","smoothterrain.vtk").as_posix())
        # target=Path(self.caseParent,self.caseName,"terrain","terrain.amrwind").open("w")
        # for i in range(0,len(x1)):
        #     if(x1[i]>xleft and x1[i]<xright and x2[i]>yleft and x2[i]<yright):
        #         smoothData.points[i,0]=x1[i]
        #         smoothData.points[i,1]=x2[i]
        #         smoothData.points[i,2]=x3[i]
        #         target.write("%g %g %g\n"%(x1[i],x2[i],x3[i])) 
        #     else:
        #      target.write("%g %g %g\n"%(smoothData.points[i,0],smoothData.points[i,1],smoothData.points[i,2]))
        # target.close()
        # smoothData.save(Path(self.caseParent,self.caseName,"terrain","smoothterrain.vtk").as_posix())

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
        self.turbineMarkType=self.yamlFile['turbineMarkType']
        if(self.turbineMarkType=='database'):
            xmin=np.amin(self.caseLatList)+0.05
            ymin=np.amin(self.caseLonList)+0.05
            xmax=np.amax(self.caseLatList)-0.05
            ymax=np.amax(self.caseLonList)-0.05
            import pandas as pd
            df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1',usecols=[1,2])   
            data=df.to_numpy()
            self.caseTurbineLat=[]
            self.caseTurbineLon=[]
            lon=data[:,0]
            lat=data[:,1]
            df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1')
            data=df.to_numpy()
            # data=np.genfromtxt("turbine.csv",dtype=str,delimiter=",",invalid_raise = False,skip_header=1)
            location=data[:,0]
            index=0
            for i in range(0,len(lat)):
                if(lat[i]>xmin and lat[i]<xmax and lon[i]>ymin and lon[i]<ymax):
                    index=index+1
                    self.caseTurbineLat.append(lat[i])
                    self.caseTurbineLon.append(lon[i]) 
                    #print(index,xmin,xmax,ymin,ymax,location[i],lat[i],lon[i])
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
            for i in range(0,len(self.caseTurbineLon)):
                target=Path(self.caseParent,self.caseName,"terrainTurbine","turbine"+str(i+1)+"Details.info").open("w")
                xturb,yturb=self.srtm.to_xy(self.caseTurbineLat[i],self.caseTurbineLon[i])
                residual=1000000
                for k in range(0,len(xTerrain)):
                    error=np.sqrt((xTerrain[k]-xturb)**2+(yTerrain[k]-yturb)**2)
                    if(error<residual):
                        residual=error
                        zturb=zTerrain[k]
                        kloc=k
                zTurbineLoc.append(zturb)
                self.turbineX1.append(xturb)
                self.turbineX2.append(yturb)
                self.turbineX3.append(zturb)
                index=index+1
                #print(index,xturb,yturb,zturb,xTerrain[kloc],yTerrain[kloc],zTerrain[kloc],residual)
                self.createDefaultTurbine(xturb,yturb,zturb,index,target)
                target.close()
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
        target.write("tagging.labels =")
        index=0
        for i in range(0,len(self.caseTurbineLon)):
                target.write(" g%g gg%g ggg%g"%(i,i,i))
        target.write("\n")
        target.close()
        # Write Refinement Files
        tempTurbineRadius=max(50,0.4*globalMindistance)
        for i in range(0,len(self.caseTurbineLon)):
            target=Path(self.caseParent,self.caseName,"terrainTurbine","turbineRefinement"+str(i+1)+"Details.info").open("w")
            xstart=self.turbineX1[i]-4*tempTurbineRadius
            ystart=self.turbineX2[i]-4*tempTurbineRadius
            zstart=self.turbineX3[i]
            target.write("tagging.g%g.type \t\t\t = GeometryRefinement\n"%(i))
            target.write("tagging.g%g.shapes \t\t\t = b%g\n"%(i,i))
            target.write("tagging.g%g.level \t\t\t = 0\n"%(i))
            target.write("tagging.g%g.b%g.type \t\t\t = box\n"%(i,i))
            target.write("tagging.g%g.b%g.origin = %g %g %g \n"%(i,i,xstart,ystart,zstart))
            target.write("tagging.g%g.b%g.xaxis =  %g %g %g\n"%(i,i,4*tempTurbineRadius,0,0))
            target.write("tagging.g%g.b%g.yaxis =  %g %g %g\n"%(i,i,0,4*tempTurbineRadius,0))
            target.write("tagging.g%g.b%g.zaxis = %g %g %g\n"%(i,i,0,0,400))
            xstart=self.turbineX1[i]-2*tempTurbineRadius
            ystart=self.turbineX2[i]-2*tempTurbineRadius
            zstart=self.turbineX3[i]
            target.write("tagging.gg%g.type \t\t\t = GeometryRefinement\n"%(i))
            target.write("tagging.gg%g.shapes \t\t\t = bb%g\n"%(i,i))
            target.write("tagging.gg%g.min_level \t\t\t = 0\n"%(i))
            target.write("tagging.gg%g.max_level \t\t\t = 1\n"%(i))
            target.write("tagging.gg%g.bb%g.type \t\t\t = box\n"%(i,i))
            target.write("tagging.gg%g.bb%g.origin = %g %g %g \n"%(i,i,xstart,ystart,zstart))
            target.write("tagging.gg%g.bb%g.xaxis =  %g %g %g\n"%(i,i,2*tempTurbineRadius,0,0))
            target.write("tagging.gg%g.bb%g.yaxis =  %g %g %g\n"%(i,i,0,2*tempTurbineRadius,0))
            target.write("tagging.gg%g.bb%g.zaxis = %g %g %g\n"%(i,i,0,0,200))
            xstart=self.turbineX1[i]-tempTurbineRadius
            ystart=self.turbineX2[i]-tempTurbineRadius
            zstart=self.turbineX3[i]
            target.write("tagging.ggg%g.type \t\t\t = GeometryRefinement\n"%(i))
            target.write("tagging.ggg%g.shapes \t\t\t = bbb%g\n"%(i,i))
            target.write("tagging.ggg%g.min_level \t\t\t = 0\n"%(i))
            target.write("tagging.ggg%g.max_level \t\t\t = 2\n"%(i))
            target.write("tagging.ggg%g.bbb%g.type \t\t\t = box\n"%(i,i))
            target.write("tagging.ggg%g.bbb%g.origin = %g %g %g \n"%(i,i,xstart,ystart,zstart))
            target.write("tagging.ggg%g.bbb%g.xaxis =  %g %g %g\n"%(i,i,tempTurbineRadius,0,0))
            target.write("tagging.ggg%g.bbb%g.yaxis =  %g %g %g\n"%(i,i,0,tempTurbineRadius,0))
            target.write("tagging.ggg%g.bbb%g.zaxis = %g %g %g\n"%(i,i,0,0,160))
            target.close()
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


    def createDefaultTurbine(self,xi,yi,zi,index,target):
        string="Actuator.Turb"+str(index)
        target.write("# Turbine %g\n"%(index))
        target.write(string+".type           = JoukowskyDisk\n")
        target.write(string+".base_position  = %g %g %g\n"%(xi,yi,zi))
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
        target.write(string+".rpm            = 8 8 8 9.172052607 10.17611854 10.99428938 11.62116715 12.05261594 12.28578923 12.31914861 12.48582322 12.85143216 13.413563 ")
        target.write("14.16850788 15.11128509 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026965 15.16026 ")
        target.write(" 15.16026965 15.16026965 15.16026965 15.16026965\n")
        target.write(string+".num_points_r   = 5\n")
        target.write(string+".num_points_t   = 5\n")
        target.write(string+".num_blades     = 3\n")
        target.write(string+".vortex_core_size = 13.0\n")
        target.write(string+".use_tip_correction = true\n")
        target.write(string+".use_root_correction = true\n")

    def plotTurbines(self):
        import pyvista as pv 
        import pandas as pd
        print("Plotting Turbines")
        # df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1',usecols=[1,2])
        # data=df.to_numpy()        
        # strings=['Biglow','PaTu','Hay Canyon']
        # xlat=[]
        # xlon=[]
        # lon=data[:,0]
        # lat=data[:,1]
        # df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1')
        # data=df.to_numpy()  
        # location=data[:,0]
        # for turbineLocation in strings:
        #     for i in range(0,len(location)):
        #         if(turbineLocation in location[i]):
        #             print(location[i])
        #             xlat.append(lat[i])
        #             xlon.append(lon[i]) 
        # print(len(xlat))
        #data=np.genfromtxt("turbine.csv",skip_header=1,delimiter=",")
        x1=self.terrainX1.flatten(order='F')
        x2=self.terrainX2.flatten(order='F')
        x3=self.terrainX3.flatten(order='F')
        data=np.column_stack([x1,x2,x3])
        pl = pv.Plotter()
        mesh2=pv.PolyData(data)
        mesh2['elevation']=data[:,2]
        surf = mesh2.delaunay_2d()
        pl.add_mesh(surf)
        # xlat=np.asarray(xlat)
        # xlon=np.asarray(xlon)
        # xlat=np.asarray(self.caseTurbineLat)
        # xlon=np.asarray(self.caseTurbineLon)
        # zloc=0*xlat
        # x1=self.caseLatList.flatten(order='F')
        # x2=self.caseLonList.flatten(order='F')
        # #print(len(xlat),len(xlon),len(zloc),len(x1),len(x2))
        # for i in range(0,len(xlat)):
        #     residual=10000000
        #     error=0
        #     for j in range(0,len(x1)):
        #         error=np.sqrt((x1[j]-xlat[i])**2+(x2[j]-xlon[i])**2)
        #         if(residual>error):
        #             residual=error
        #             zloc[i]=x3[j]+80
            #print(i,residual,zloc[i])
        # for i in range(0,len(xlat)):
        #         print("Plotting turbine",i+1)
        #         xlat[i],xlon[i]=self.srtm.to_xy(xlat[i],xlon[i])
        #         newDisk=pv.Disc(center=(xlat[i],xlon[i],zloc[i]),inner=8,outer=50,normal=(1.0, 0.0,0.0), r_res=1, c_res=24)
        #         pl.add_mesh(newDisk)
        #         if(i==0):
        #             globalBox=newDisk
        #         else:
        #             localBox=newDisk
        #             tempbox=globalBox.merge([localBox])
        #             globalBox=tempbox
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
