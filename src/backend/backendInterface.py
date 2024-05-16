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
        self.srtm = SRTM(bounds,fpath=srtm_output.as_posix(),product=product)
        self.srtm.download()
        x1,x2,x3 = self.srtm.to_terrain(dx,dy)
        cornerCut=12
        self.caseX1=x1[cornerCut:x1.shape[0]-cornerCut,cornerCut:x1.shape[1]-cornerCut]
        self.caseX2=x2[cornerCut:x2.shape[0]-cornerCut,cornerCut:x2.shape[1]-cornerCut]
        self.caseX3=x3[cornerCut:x3.shape[0]-cornerCut,cornerCut:x3.shape[1]-cornerCut]
        self.caseLatList=0*self.caseX1
        self.caseLonList=0*self.caseX2
        for i in range(0,self.caseX1.shape[0]):
            for j in range(0,self.caseX1.shape[1]):
                self.caseLatList[i,j],self.caseLonList[i,j]=self.srtm.to_latlon(self.caseX1[i,j],self.caseX2[i,j])
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
        self.createAMRGeometry(self.amrPrecursorFile,1)
        self.createAMRGeometry(self.amrTerrainFile,-1)
        self.createAMRGrid(self.amrPrecursorFile)
        self.createAMRGrid(self.amrTerrainFile)
        self.createAMRTime(self.amrPrecursorFile)
        self.createAMRTime(self.amrTerrainFile)
        self.createSolverInfo(self.amrPrecursorFile)
        self.createSolverInfo(self.amrTerrainFile,1)
        self.createAMRTransport(self.amrPrecursorFile)
        self.createAMRTransport(self.amrTerrainFile)
        self.createAMRTurbulence(self.amrPrecursorFile)
        self.createAMRTurbulence(self.amrTerrainFile)
        self.createAMRABLData(self.amrPrecursorFile,0,1)
        self.createAMRABLData(self.amrTerrainFile,1,0)  
        self.createAMRSourceTerm(self.amrPrecursorFile)
        self.createAMRSourceTerm(self.amrTerrainFile,1)      
        self.createAMRBC(self.amrPrecursorFile)
        self.createAMRBC(self.amrTerrainFile,1)
        self.createAMRTolerance(self.amrPrecursorFile)
        self.createAMRTolerance(self.amrTerrainFile)
        self.plotTurbines()

    def createAMRGeometry(self,target,periodic=-1):
        target.write("# Geometry\n")
        minX=np.amin(self.caseX1)
        minY=np.amin(self.caseX2)
        minZ=np.amin(self.caseX3)
        maxX=np.amax(self.caseX1)
        maxY=np.amax(self.caseX2)
        maxZ=max(2*np.amax(self.caseX3),1600)
        target.write("geometry.prob_lo \t\t\t = %g %g %g \n"%(minX,minY,minZ))
        target.write("geometry.prob_hi \t\t\t = %g %g %g \n"%(maxX,maxY,maxZ))
        if(periodic==1):
            target.write("geometry.is_periodic \t\t\t = 1 1 0\n")
        else:
            target.write("geometry.is_periodic \t\t\t = 0 0 0\n")

    def createAMRGrid(self,target):
        self.caseCellSize=self.yamlFile['cellSize']
        nx=int((np.amax(self.caseX1)-np.amin(self.caseX1))/self.caseCellSize)
        while (nx%8 !=0):
            nx=nx+1
        print(nx)
        ny=int((np.amax(self.caseX2)-np.amin(self.caseX2))/self.caseCellSize)
        while (ny%8 !=0):
            ny=ny+1
        print(ny)
        nz=int((max(2*np.amax(self.caseX3),1600)-np.amin(self.caseX3))/self.caseCellSize)
        while (nz%8 !=0):
            nz=nz+1
        print(nz)
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
            UPeriod=int((np.amax(self.caseX1)-np.amin(self.caseX1))/200)
            VPeriod=int((np.amax(self.caseX1)-np.amin(self.caseX1))/200)
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
        target.write("ABL.stats_output_format \t\t\t = native\n")
        target.write("ABL.surface_roughness_z0 \t\t\t = %g\n"%(self.refRoughness))
        target.write("ABL.temperature_heights \t\t\t = 0  800 900 1900\n")
        target.write("ABL.temperature_values  \t\t\t = 300 300 308 311\n")
        target.write("ABL.wall_shear_stress_type \t\t\t = local\n")
        target.write("ABL.surface_temp_flux \t\t\t = %g\n"%(self.refHeatFlux))
        if(iomode==0):
            target.write("ABL.bndry_io_mode \t\t\t = 0\n")
            if(self.caseWindspeedX>0 and self.caseWindspeedY>0):
                target.write("ABL.bndry_planes \t\t\t = xlo ylo \n")     
            elif(self.caseWindspeedX>0 and self.caseWindspeedY<0):
                target.write("ABL.bndry_planes \t\t\t = xlo yhi \n")      
            if(self.caseWindspeedX<0 and self.caseWindspeedY>0):
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
        elif(iomode==2):
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
            if(self.caseWindspeedX>0):
                target.write('xlo.type \t\t\t = "mass_inflow"\n')
                target.write("xlo.density \t\t\t = 1.225\n")
                target.write("xlo.temperature \t\t\t = 300\n")
                target.write('xhi.type \t\t\t = "pressure_outflow"\n')
            else:
                target.write('xhi.type \t\t\t = "mass_inflow"\n')
                target.write("xhi.density \t\t\t = 1.225\n")
                target.write("xhi.temperature \t\t\t = 300\n")
                target.write('xlo.type \t\t\t = "pressure_outflow"\n')   
            if(self.caseWindspeedY>0):
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

    def createAMRTolerance(self,target):
        target.write("# Tolerance \n")
        target.write("nodal_proj.mg_rtol \t\t\t = 1.0e-4\n")
        target.write("nodal_proj.mg_atol = 1.0e-6\n")
        target.write("mac_proj.mg_rtol = 1.0e-4\n")
        target.write("mac_proj.mg_atol = 1.0e-6\n")        

    def closeAMRFiles(self):
        self.amrPrecursorFile.close()
        self.amrTerrainFile.close()

    def plotTurbines(self):
        import pyvista as pv 
        import pandas as pd
        df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1',usecols=[1,2])
        data=df.to_numpy()        
        strings=['Biglow','PaTu','Hay Canyon']
        xlat=[]
        xlon=[]
        lon=data[:,0]
        lat=data[:,1]
        df=pd.read_csv('turbine.csv', sep=',',encoding='latin-1')
        data=df.to_numpy()  
        location=data[:,0]
        for turbineLocation in strings:
            for i in range(0,len(location)):
                if(turbineLocation in location[i]):
                    print(location[i])
                    xlat.append(lat[i])
                    xlon.append(lon[i]) 
        print(len(xlat))
        #data=np.genfromtxt("turbine.csv",skip_header=1,delimiter=",")
        x1=self.caseX1.flatten(order='F')
        x2=self.caseX2.flatten(order='F')
        x3=self.caseX3.flatten(order='F')
        data=np.column_stack([x1,x2,x3])
        pl = pv.Plotter()
        mesh2=pv.PolyData(data)
        mesh2['elevation']=data[:,2]
        surf = mesh2.delaunay_2d()
        pl.add_mesh(surf)
        xlat=np.asarray(xlat)
        xlon=np.asarray(xlon)
        zloc=0*xlat
        x1=self.caseLatList.flatten(order='F')
        x2=self.caseLonList.flatten(order='F')
        for i in range(0,len(xlat)):
            residual=10000000
            error=0
            for j in range(0,len(x1)):
                error=np.sqrt((x1[j]-xlat[i])**2+(x2[j]-xlon[i])**2)
                if(residual>error):
                    residual=error
                    zloc[i]=x3[j]
            print(i,residual,zloc[i])
        for i in range(0,len(xlat)):
                xlat[i],xlon[i]=self.srtm.to_xy(xlat[i],xlon[i])
        data=np.column_stack([xlat,xlon,zloc])
        mesh1=pv.PolyData(data)
        pl.add_mesh(mesh1,render_points_as_spheres=True)        
        pl.view_xy()
        pl.show_axes()
        pl.show()


from sys import argv 
amrRef=amrBackend(argv[1])
amrRef.createDomain()
amrRef.createAMRFiles()
