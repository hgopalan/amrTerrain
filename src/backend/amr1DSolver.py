'''



'''

import numpy as np 
import time 
from joblib import Parallel, delayed

class amr1dSolver:

    def __init__(self,nz,zupper,z0,terrain,pathToWrite):
        self.setupgrid(nz,zupper)
        self.initialize_variables()
        self.allocate_memory()
        self.allocate_mean_memory()
        self.z0=z0
        self.initialize_remaining_variables(terrain,pathToWrite)

    def setupgrid(self,nz,zupper):
        self.zupper=zupper
        self.nz=nz
    
    def initialize_variables(self):
        self.ug=10
        self.vg=0
        self.t_ref=300
        self.delta_inversion=8
        self.inversion_height=468
        self.inversion_width=83
        self.lapse_rate=0.003
    
    def allocate_memory(self):
        self.ux=np.zeros(self.nz)
        self.uy=np.zeros(self.nz)
        self.temperature=np.zeros(self.nz)
        self.tke=np.zeros(self.nz)
        self.nut=np.zeros(self.nz)
        self.nutPrime=np.zeros(self.nz)
        self.sigmaT=np.zeros(self.nz)
        self.Rt=np.zeros(self.nz)

    def allocate_mean_memory(self):
        self.umean=np.zeros(self.nz)
        self.vmean=np.zeros(self.nz)
        self.tmean=np.zeros(self.nz)
        self.tkemean=np.zeros(self.nz)
        self.nutmean=np.zeros(self.nz)

    def initialize_remaining_variables(self,terrain,pathToWrite):
        self.pblh=1000
        self.Qh=0
        self.Qb=0
        self.ref_mol=-1e+30
        self.lscale=np.zeros(self.nz)
        self.coriolis=0
        self.ustar=0.41
        self.thetastar=0.0
        self.phiM=1
        self.start_time=0
        self.end_time=20000
        self.tke_init=0.4
        self.heat_flux_mode=1
        self.mode_value=0.0 
        self.z=np.linspace(0,self.zupper,self.nz)
        self.dz=self.z[1]-self.z[0]
        self.terrain_height=terrain
        self.lower=np.abs(self.z - terrain).argmin()
        self.writePath=pathToWrite

    def initialize_physics(self,ug,vg,t_ref,tke_init,ref_ustar,pblh):
        self.ug=ug
        self.vg=vg
        self.t_ref=t_ref
        self.tke_init=tke_init
        self.ustar=ref_ustar
        self.pblh=pblh

    def return_heights(self):
        return self.z

    def return_windspeed(self):
        return self.ux,self.uy

    def reinitialize_windspeed(self,ux,uy):
        self.ux=ux
        self.uy=uy

    def initialize_coriolis(self,lat):
        omega=7.292115e-5
        self.coriolis=2*omega*np.sin(lat*np.pi/180)

    def heat_flux_model(self,mode,value):
        self.heat_flux_mode=mode 
        self.mode_value=float(value) 

    def temperature_inversion(self,height,width,strength,rate):
        self.delta_inversion=strength
        self.inversion_height=height
        self.inversion_width=width
        self.lapse_rate=rate 

    def run_simulation(self,end_time,convergence=1e-4,forceMetMastWind=False,metMastHeight=100,metMastX=10,metmastY=0):
        self.counter=0
        self.start_time=0
        self.err=convergence 
        self.converge=False 
        self.initialize_simulation()
        self.counter=0
        self.end_time=end_time
        self.zloc=[]
        if(forceMetMastWind):
            self.zlist=[metMastHeight]
            for i in range(0,len(self.zlist)):
                self.zloc.append(np.abs(self.z - self.zlist[i]).argmin())
                print("Index:",self.zloc)
            self.uxlist=[metMastX]
            self.vxlist=[metmastY]
        self.dt=0.8*(self.z[1]-self.z[0])/np.sqrt(self.ug**2+self.vg**2)
        self.errvelx=np.sum(self.ux)
        self.errvely=np.sum(self.uy)
        self.errtke=np.sum(self.tke)
        self.errnut=np.sum(self.nut)
        self.errtemp=np.sum(self.temperature)
        while(self.start_time<=self.end_time):
            self.start_time+=self.dt
            self.compute_similarity()
            for i in range(self.lower+1,self.nz-1):
                self.update_windspeed_x(i,self.dt) 
                self.update_windspeed_y(i,self.dt) 
                self.update_temperature(i,self.dt)
                self.update_turbulence(i,self.dt) 
            self.dt=0.8*(self.z[1]-self.z[0])/np.sqrt(max(self.ux)**2+max(self.uy)**2)
            self.counter+=1
            self.compute_average()
            self.compute_error()
            if(self.converge):
                print("Converged:")
                break
        # Normalize Average 
        self.umean=self.umean/self.counter
        self.vmean=self.vmean/self.counter
        self.tmean=self.tmean/self.counter
        self.tkemean=self.tkemean/self.counter
        self.nutmean=self.nutmean/self.counter
        self.compute_zi()
        target=open(self.writePath,"w")
        for i in range(0,len(self.z)):
            target.write("%g %g %g %g %g %g %g %g \n"%(self.z[i]-self.terrain_height,self.ux[i],self.uy[i],0.0,self.temperature[i],self.tke[i],self.nut[i],self.lscale[i]))
        target.close()

    def initialize_simulation(self):
        for i in range(0,self.nz):
            self.ux[i]=self.ug
            self.uy[i]=self.vg
            if(self.z[i]<=self.inversion_height):
                self.temperature[i]=self.t_ref
            elif(self.z[i]>self.inversion_height and self.z[i]<=(self.inversion_height+self.inversion_width)):
                self.temperature[i]=self.t_ref+(self.z[i]-self.inversion_height)*0.08
            else:
                if(self.z[i]<=self.z[self.lower]):
                    self.temperature[i]=self.t_ref
                else:
                    self.temperature[i]=self.temperature[i-1]+self.lapse_rate*max(self.z[i]-self.z[i-1],0.0)
            self.nut[i]=0.41**2*self.z[i]*(1-min(self.z[i],790)/800)**2 #1e-5
            self.tke[i]=0.41**2*self.z[i]*(1-min(self.z[i],790)/800)**2 #0.1
            self.lscale[i]=self.z[i]*(1-min(self.z[i],790)/800)**2 #0.1
            self.nutPrime[i]=self.nut[i]
            self.sigmaT[i]=1.0
        if(self.lower==0):
            pass
        else:
            self.ux[0:self.lower]*=0.0
            self.uy[0:self.lower]*=0.0
            self.temperature[0:self.lower]=self.temperature[self.lower]+0*self.temperature[0:self.lower]
            self.nut[0:self.lower]*=0.0
            self.tke[0:self.lower]*=0.0
            self.lscale[0:self.lower]*=0.0

    def compute_similarity(self):
        M1=np.sqrt(self.ux[self.lower+1]**2+self.uy[self.lower+1]**2)
        self.ustar = 0.41*M1/np.log(self.z[1]/self.z0)
        iter=0
        psi_h=0
        psi_m=0
        iter=0
        error=25
        self.mo_length=self.ref_mol
        # Qh Specified and compulte wall temperature and MOL 
        if(self.heat_flux_mode==1):
            self.Qh=self.mode_value
            while (iter <= 25 and error>1e-5):
                utau_iter = self.ustar
                self.temperature[self.lower] = self.Qh*(np.log(self.z[1]/self.z0) - psi_h)/(self.ustar*0.41)+self.temperature[self.lower+1]
                if(abs(self.Qh)>1e-5):
                    self.mo_length = -self.ustar**3*self.temperature[self.lower+1]/(0.41*9.81*self.Qh)
                    zeta = self.z[1] / self.mo_length
                else: 
                    self.mo_length=-1e30
                    zeta = 0.0;
                if(zeta>=0):
                    psi_m=-5*zeta
                    psi_h=-5*zeta
                else:
                    x = np.sqrt(1 - 16 * zeta)
                    psi_h=2.0 *np.log(0.5 * (1 + x))
                    x = np.sqrt(np.sqrt(1 - 16 * zeta))
                    psi_m= np.log(0.5*(1+x**2)*0.25*(1+x)**2) - \
                            2.0 * np.arctan(x) + 0.5*np.pi
                self.ustar = 0.41 * M1 / (np.log(self.z[1] / self.z0) - psi_m)
                iter+=1
                error=abs(self.ustar-utau_iter)
        # Surface temperature specified. Compute sensible heat flux and MOL 
        elif(self.heat_flux_mode==2):
            self.temperature[0]=self.mode_value
            if(self.temperature[self.lower]==self.t_ref):
                self.Qh=(self.temperature[self.lower]-self.temperature[self.lower+1])*(self.ustar*0.41)/(np.log(self.z[1]/self.z0) - psi_h)
                self.mo_length=-1e30
            else:
                while (iter <= 25 and error>1e-5):
                    utau_iter = self.ustar
                    self.Qh=(self.temperature[self.lower]-self.temperature[self.lower+1])*(self.ustar*0.41)/(np.log(self.z[1]/self.z0) - psi_h)
                    if(abs(self.Qh)>1e-5):
                        self.mo_length = -self.ustar**3*self.temperature[self.lower+1]/(0.41*9.81*self.Qh)
                        zeta = self.z[1] / self.mo_length
                    else: 
                        self.mo_length=-1e30
                        zeta = 0.0;
                    if(zeta>=0):
                        psi_m=-5*zeta
                        psi_h=-5*zeta
                    else:
                        x = np.sqrt(1 - 16 * zeta)
                        psi_h=2.0 *np.log(0.5 * (1 + x))
                        x = np.sqrt(np.sqrt(1 - 16 * zeta))
                        psi_m= np.log(0.5*(1+x**2)*0.25*(1+x)**2) - \
                                2.0 * np.arctan(x) + 0.5*np.pi
                    self.ustar = 0.41 * M1 / (np.log(self.z[1] / self.z0) - psi_m)
                    iter+=1
                    error=abs(self.ustar-utau_iter)
        # Heating/Cooling Rate Specified 
        elif(self.heat_flux_mode==3):
            self.temperature[self.lower]=self.t_ref+self.mode_value*self.start_time/3600
            if(self.temperature[self.lower]==self.t_ref):
                self.Qh=(self.temperature[self.lower]-self.temperature[self.lower+1])*(self.ustar*0.41)/(np.log(self.z[1]/self.z0) - psi_h)
                self.mo_length=-1e30
            else:
                while (iter <= 25 and error>1e-5):
                    utau_iter = self.ustar
                    self.Qh=(self.temperature[self.lower]-self.temperature[self.lower+1])*(self.ustar*0.41)/(np.log(self.z[1]/self.z0) - psi_h)
                    if(abs(self.Qh)>1e-5):
                        self.mo_length = -self.ustar**3*self.temperature[self.lower+1]/(0.41*9.81*self.Qh)
                        zeta = self.z[1] / self.mo_length
                    else: 
                        self.mo_length=-1e30
                        zeta = 0.0;
                    if(zeta>=0):
                        psi_m=-5*zeta
                        psi_h=-5*zeta
                    else:
                        x = np.sqrt(1 - 16 * zeta)
                        psi_h=2.0 *np.log(0.5 * (1 + x))
                        x = np.sqrt(np.sqrt(1 - 16 * zeta))
                        psi_m= np.log(0.5*(1+x**2)*0.25*(1+x)**2) - \
                                2.0 * np.arctan(x) + 0.5*np.pi
                    self.ustar = 0.41 * M1 / (np.log(self.z[1] / self.z0) - psi_m)
                    iter+=1
                    error=abs(self.ustar-utau_iter)
        # MOL Specified. Compute surface flux and wall temperature 
        # Used for Operational Runs 
        else:
            self.mo_length=self.mode_value
            zeta=self.z[1]/self.mo_length
            error=100
            # ustar 
            if(zeta>0):
                psi_m=-5*zeta
                psi_h=-5*zeta
            else:
                x = np.sqrt(1 - 16 * zeta)
                psi_h=2.0 *np.log(0.5 * (1 + x))
                x = np.sqrt(np.sqrt(1 - 16 * zeta))
                psi_m= np.log(0.5*(1+x**2)*0.25*(1+x)**2) - \
                        2.0 * np.arctan(x) + 0.5*np.pi
            self.ustar=M1*0.41/(np.log(self.z[1]/self.z0)-psi_m)
            self.thetastar=self.ustar**2*self.temperature[self.lower+1]/(0.41*9.81*self.mo_length)
            self.temperature[self.lower]=self.temperature[self.lower+1]-self.thetastar/0.41*(np.log(self.z[1]/self.z0)-psi_h)
            self.Qh=-self.ustar*self.thetastar
        if(self.mo_length<0):
            self.phi_m=(1-16*self.z0/self.mo_length)**(-0.25)
        else:
            self.phi_m=1+5*self.z0/self.mo_length
        self.nut[self.lower]=self.ustar*0.41*self.z0/self.phi_m
        M0=M1-self.ustar/0.41*self.phi_m
        self.ux[self.lower]=M0*self.ux[self.lower+1]/M1
        self.uy[self.lower]=M0*self.uy[self.lower+1]/M1
        self.Qb=9.81/self.t_ref*self.Qh   
        self.tke[self.lower]=self.ustar**2/0.556**2+(max(self.Qb,0)*0.41*self.z[1]/0.556**3)**(2.0/3.0)
        self.lscale[self.lower]=0
        self.lscale[self.nz-1]=self.lscale[self.nz-2]
        self.tke[self.nz-1]=self.tke[self.nz-2]
        self.nut[self.nz-1]=self.nut[self.nz-2]
        self.ux[self.nz-1]=self.ux[self.nz-2]
        self.uy[self.nz-1]=self.uy[self.nz-2]
        self.temperature[self.nz-1]=self.temperature[self.nz-2]
        if(self.lower==0):
            pass
        else:
            self.ux[0:self.lower]*=0.0
            self.uy[0:self.lower]*=0.0
            self.temperature[0:self.lower]=self.temperature[self.lower]+0*self.temperature[0:self.lower]
            self.nut[0:self.lower]*=0.0
            self.tke[0:self.lower]*=0.0

    def update_windspeed_x(self,i,dt):
        term1=self.nut[i]*(self.ux[i+1]-2*self.ux[i]+self.ux[i-1])/self.dz**2
        if(i==1):
            dudz=1/self.dz*(self.ux[i+1]-self.ux[i])
        else:
            dudz=1/self.dz*(self.ux[i+1]-self.ux[i])
        term2=0.5/self.dz*(self.nut[i+1]-self.nut[i-1])*dudz
        coriolis=self.coriolis*(self.uy[i])
        geostrophic=-self.coriolis*self.vg
        # Forcing 
        forcing=0
        for j in range(0,len(self.zloc)):
            if(i==self.zloc[j]):
                forcing=-(self.ux[i]-self.uxlist[j])/self.dt  
        self.ux[i]=self.ux[i]+dt*(term1+term2+coriolis+geostrophic+forcing)
        return 1 
    
    def update_windspeed_y(self,i,dt):
        term1=self.nut[i]*(self.uy[i+1]-2*self.uy[i]+self.uy[i-1])/self.dz**2
        if(i==1):
            dvdz=1/self.dz*(self.uy[i+1]-self.uy[i])
        else:
            dvdz=1/self.dz*(self.uy[i+1]-self.uy[i])
        term2=0.5/self.dz*(self.nut[i+1]-self.nut[i-1])*dvdz
        coriolis=-self.coriolis*(self.ux[i])
        geostrophic=self.coriolis*self.ug
        # Forcing 
        forcing=0
        for j in range(0,len(self.zloc)):
            if(i==self.zloc[j]):
                forcing=-(self.uy[i]-self.vxlist[j])/self.dt
        self.uy[i]=self.uy[i]+dt*(term1+term2+coriolis+geostrophic+forcing)
        return 1 

    def update_temperature(self,i,dt):
        term1=self.nut[i]/self.sigmaT[i]*(self.temperature[i+1]-2*self.temperature[i]+self.temperature[i-1])/self.dz**2
        term2=1/self.dz*(self.nut[i]/self.sigmaT[i]-self.nut[i-1]/self.sigmaT[i-1])*1/self.dz*(self.temperature[i]-self.temperature[i-1])
        self.temperature[i]=self.temperature[i]+dt*(term1+term2)
        return 1 

    def update_turbulence(self,i,dt):
        term1=self.nut[i]*(self.tke[i+1]-2*self.tke[i]+self.tke[i-1])/self.dz**2
        term2=1/self.dz*(self.nut[i]-self.nut[i-1])*1/self.dz*(self.tke[i]-self.tke[i-1])
        production=self.nut[i]*(1/self.dz)**2*((self.ux[i]-self.ux[i-1])**2+(self.uy[i]-self.uy[i-1])**2)
        # Original RANS Model 
        lturb=0.41*(self.z[i]-self.z[self.lower])
        lmax=min(0.00027*np.sqrt(self.ug**2+self.vg**2)/abs(self.coriolis),100)
        invLshear=1.0/lturb**2+1.0/lmax**2
        lshear=np.sqrt(1.0/invLshear)
        stratification=9.81*1/(self.dz*self.t_ref)*(self.temperature[i]-self.temperature[i-1])
        dissipation=0.556**3*self.tke[i]**1.5/self.lscale[i]
        Rt=(self.tke[i]/dissipation)**2*stratification
        if(Rt<-1):
            Rt=max(Rt,Rt-(1+Rt)**2/(Rt-1))
        self.Rt[i]=Rt
        buoyancy=-self.nutPrime[i]*stratification
        if(Rt>0):
            lbuoyancy=0.25*np.sqrt(self.tke[i])/np.sqrt(max(stratification,1e-15))
            invLscale=1.0/lshear**2+1.0/lbuoyancy**2
            lscale=np.sqrt(1.0/invLscale)
        else:
            lscale=lshear*np.sqrt(1-0.556**6/0.35**2*Rt)
        diffusion=term1+term2
        self.tke[i]=self.tke[i]+dt*(production+buoyancy-dissipation+diffusion)
        self.tke[i]=max(self.tke[i],1e-15)
        cmu=(0.556+0.108*Rt)/(1+0.308*Rt+0.00837*Rt**2)
        self.nut[i]=cmu*np.sqrt(self.tke[i])*lscale
        cmuprime=0.556/(1+0.277*Rt)
        self.sigmaT[i]=(1+0.193*Rt)/(1+0.0302*Rt)
        self.nutPrime[i]=cmuprime*np.sqrt(self.tke[i])*lscale
        self.lscale[i]=lscale
        return 1 
    
    def compute_zi(self):
        M=self.ux**2+self.uy**2
        Rib=9.81*self.z/self.temperature[self.lower]*(self.temperature-self.temperature[self.lower])/(M+1e-5)
        Ric=0.5
        zi=max(self.z)
        for i in range(self.lower,len(Rib)):
            if(Rib[i]>Ric):
                zi=self.z[i]
                break

    def compute_average(self):
        self.umean+=self.ux
        self.vmean+=self.uy
        self.tmean+=self.temperature
        self.tkemean+=self.tke
        self.nutmean+=self.nut

    def compute_error(self):
        errvelx=abs(np.sum(self.ux)-self.errvelx)
        errvely=abs(np.sum(self.uy)-self.errvely)
        errtke=abs(np.sum(self.tke)-self.errtke)
        errnut=abs(np.sum(self.nut)-self.errnut)
        errtemp=abs(np.sum(self.temperature)-self.errtemp)
        if(self.counter%5000==0):
            print("Residual:%g %g %g %g %g %g"%(self.start_time,errvelx,errvely,errtemp,errtke,errnut))
        if(errvelx<self.err and errvely<self.err):
            self.converge=True
        self.errvelx=np.sum(self.ux)
        self.errvely=np.sum(self.uy)
        self.errtke=np.sum(self.tke)
        self.errnut=np.sum(self.nut)
        self.errtemp=np.sum(self.temperature)


    def print_mean_results(self):
        import matplotlib.pylab as plt 
        #plt.plot(self.umean,self.z)
        plt.plot(self.vmean,self.z)
        #M=np.sqrt(self.umean**2+self.vmean**2)
        #plt.plot(M,self.z)
        plt.savefig("1DWindSpeed_mean.pdf")
        plt.figure(2)
        plt.plot(self.nutmean,self.z)
        plt.savefig("1DNut_mean.pdf")
        plt.figure(3)
        plt.plot(self.tkemean,self.z)
        plt.savefig("1DTKE_mean.pdf")
        plt.figure(4)
        plt.plot(self.lscale,self.z)
        plt.savefig("1DLt.pdf")
        plt.figure(5)
        plt.plot(self.tmean,self.z )
        plt.savefig("1DTemp_mean.pdf")
        plt.figure(6)
        zeta = self.z / self.mo_length
        psi_m=0+0*zeta 
        for i in range(0,self.nz):
            if(zeta[i]>0):
                psi_m=-5*zeta[i]
            else:
                x = np.sqrt(np.sqrt(1 - 16 * zeta[i]))
                psi_m[i]= np.log(0.5*(1+x**2)*0.25*(1+x)**2) - \
                        2.0 * np.arctan(x) + 0.5*np.pi
        plt.semilogx(self.z,np.sqrt(self.umean**2+self.vmean**2))
        plt.semilogx(self.z,self.ustar/0.41*(np.log((self.z+self.z0)/self.z0)-psi_m),'ro')
        plt.xlim(5,200)
        plt.savefig("1DLogLaw_mean.pdf")
        plt.close('all')    
        #plt.show()
    
    def print_instantaneous_results(self):
        import matplotlib.pylab as plt 
        #plt.plot(self.ux,self.z)
        plt.plot(self.uy,self.z)
        #M=np.sqrt(self.ux**2+self.uy**2)
        #plt.plot(M,self.z)
        plt.savefig("1DWindSpeedY.pdf")
        plt.figure(6)
        plt.plot(self.ux,self.z)
        plt.plot(self.uy,self.z)
        M=np.sqrt(self.ux**2+self.uy**2)
        plt.plot(M,self.z)
        plt.savefig("1DWindSpeedX.pdf")
        plt.figure(2)
        plt.plot(self.nut,self.z)
        plt.savefig("1DNut.pdf")
        plt.figure(3)
        plt.plot(self.tke,self.z)
        plt.savefig("1DTKE.pdf")
        plt.figure(5)
        plt.plot(self.temperature,self.z )
        plt.savefig("1DTemp.pdf")
        plt.figure(7)
        plt.plot(np.arctan(self.uy/self.ux)*180/np.pi,self.z)
        plt.savefig("1DWindDir.pdf")
        plt.close('all')   


