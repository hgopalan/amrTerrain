# Handcooked AEP calculator 
#!/usr/bin/env python3
# This ia an experimental tool to write WRG file from AMR-Wind particle format output
# Currently it writes only for one height level 
# Run: python windmap.py parentcasedir 
# It is assumed that parent case dir has folders terrain_0, terrain_20, terrain_40, ...
import numpy as np
import matplotlib.pyplot as plt
import sys


from amrex_particle import AmrexParticleFile

loc_dir="."
pfile = AmrexParticleFile(loc_dir)
folders=[0,20,40,60,80,100,120,140,160,180,200]
heights=[10,100]
times="01200"
backup="00600"
try:
    data=np.genfromtxt(sys.argv[1]+"/utm.info",delimiter=",",dtype=float)
    x1=data[0,1]
    x2=data[0,2]
    refx3=data[0,3]
    xref=data[1,1]
    yref=data[1,2]
    print(x1,x2,refx3,xref,yref)
except:
    print("utm.info file not available so assuming met mast is at 0,0,10 and utm reference is 0,0")
    x1=0
    x2=0
    refx3=10
    xref=0
    yref=0
folders[0]="0_new"
for dir in folders:
    x3=refx3+heights[0]
    try:
        file=sys.argv[1]+"/terrain_"+str(dir)+"/post_processing/terrain"+str(heights[0])+times
        pt_reference=pfile.special_load(file)
    except:
        file=sys.argv[1]+"/terrain_"+str(dir)+"/post_processing/terrain"+str(heights[0])+backup
        pt_reference=pfile.special_load(file)
    print(file)
    pt_reference.parse_header()
    pt_reference.load_binary_data()
    data = pt_reference.df
    x_oo = data.xco
    y_oo = data.yco
    z_oo = data.zco
    u_oo = data.velocityx
    v_oo = data.velocityy
    w_oo = data.velocityz
    # Reorder arrays
    ind_0 = np.argsort(x_oo) 
    x = np.zeros(np.shape(x_oo))
    y = np.zeros(np.shape(x_oo))
    z = np.zeros(np.shape(x_oo))
    u = np.zeros(np.shape(x_oo))
    v = np.zeros(np.shape(x_oo))
    error=10000
    for i in range(0, len(ind_0)): 
        x[i]= x_oo[ind_0[i]]
        y[i]= y_oo[ind_0[i]]
        z[i]= z_oo[ind_0[i]]
        u[i] = u_oo[ind_0[i]]
        v[i] = v_oo[ind_0[i]]
        residual=np.sqrt((x[i]-x1)**2+(y[i]-x2)**2+(z[i]-x3)**2)
        if(residual<error):
            error=residual
            xloc=x[i]
            yloc=y[i]
            zloc=z[i]
            uref=np.sqrt(u[i]**2+v[i]**2)
            u1ref=u[i]
            u2ref=v[i]
            angle=np.arctan2(v[i],u[i])*180/np.pi
            if angle < 0:
                angle+= 360
    try:
        file=sys.argv[1]+"/terrain_"+str(dir)+"/post_processing/terrain"+str(heights[1])+times
        pt = pfile.special_load(file)
    except:
        file=sys.argv[1]+"/terrain_"+str(dir)+"/post_processing/terrain"+str(heights[1])+backup
        pt = pfile.special_load(file)
    print(file)
    x3=refx3+heights[1]
    pt.parse_header()
    pt.load_binary_data()
    data = pt.df
    x_oo = data.xco
    y_oo = data.yco
    z_oo = data.zco
    u_oo = data.velocityx
    v_oo = data.velocityy
    w_oo = data.velocityz
    # Reorder arrays
    ind_0 = np.argsort(x_oo) 
    x = np.zeros(np.shape(x_oo))
    y = np.zeros(np.shape(x_oo))
    z = np.zeros(np.shape(x_oo))
    u = np.zeros(np.shape(x_oo))
    v = np.zeros(np.shape(x_oo))
    error=10000
    for i in range(0, len(ind_0)): 
        x[i]= x_oo[ind_0[i]]
        y[i]= y_oo[ind_0[i]]
        z[i]= z_oo[ind_0[i]]
        u[i] = u_oo[ind_0[i]]
        v[i] = v_oo[ind_0[i]]
    data=np.column_stack([x,y,z])
    import pyvista as pv 
    mesh=pv.PolyData(data)
    mesh["speed-up"]=np.sqrt(u[:]**2+v[:]**2)/uref
    modified_angle=(np.arctan2(v[:],u[:]))*180/np.pi
    counter=0
    for i in range(0,len(modified_angle)):
        if(modified_angle[i]<0):
            modified_angle[i]+=360
        if(abs(modified_angle[i]-angle)>60):
             #print(modified_angle[i],angle)
             counter=counter+1
             modified_angle[i]=0.0
             continue 
        modified_angle[i]=modified_angle[i]-angle        
    print("Fixed:",counter)
    mesh["Direction"]=modified_angle
    mesh.save("file_"+str(dir)+".vtk")
    print(np.min(modified_angle),np.max(modified_angle))

speed_up=np.zeros((len(folders),len(ind_0)))
global_dir=np.zeros((len(folders),len(ind_0)))
x=np.zeros(len(ind_0))
y=np.zeros(len(ind_0))
z=np.zeros(len(ind_0))
for i in range(0,len(folders)):
    data=pv.read("file_"+str(folders[i])+".vtk")
    if(i==0):
        x=data.points[:,0]
        y=data.points[:,1]
        z=data.points[:,2]
    speed_up[i,:]=data["speed-up"]
    global_dir[i,:]=data["Direction"]
data=np.genfromtxt(str(sys.argv[2]),dtype=float,skip_header=1,delimiter=",",encoding='utf-8')
windspeed=data[:,1]
winddir=data[:,21]
import scipy.stats as stats
shape, loc, scale = stats.weibull_min.fit(windspeed, floc=0)
print(scale,shape)
tempwindspeed=windspeed
# A - scale 
# k - shape 
global_scale=np.zeros(len(ind_0))
global_shape=np.zeros(len(ind_0))
local_scale=np.zeros((18,len(ind_0)))
local_shape=np.zeros((18,len(ind_0)))
# https://www.ecampmany.com/vortex/wrgfile.htm
# Line 1 Nx Ny Xmin Ymin CellSize 
# Line 2 Pt 1 
# Line 3 Pt 2 ...
# Pt 1 Details 
# Width Item
# 10    String 
# 10    X-coordinate (easting) of the site [m]
# 10    Y-coordinate (northing) of the site [m]
# 8     Z-coordinate (elevation) of the site [m]
# 5     Height above ground level [m a.g.l.]
# 5     Weibull A-parameter for the total distribution [ms-1]
# 6     Weibull k-parameter for the total distribution
# 15    Power density [Wm-2] or power production [Why-1]
# 3     Number of sectors
# 4     Frequency of occurrence for sector #1 [%·10]
# 4     Weibull A-parameter for sector #1 [ms-1·10]
# 5     Weibull k-parameter for sector #1 [·100]
target=open("test.wrg","w")

target.write("%g %g %g %g %d\n"%(len(np.unique(x)),len(np.unique(y)),np.min(x)+xref,np.min(y)+yref,(np.max(x)-np.min(x))/len(np.unique(x))))
for j in range(0,len(ind_0)):
    tempwindspeed=[]
    localsector1=[] # [0-20]
    localsector2=[] # [20-40]
    localsector3=[] # [40-60]
    localsector4=[] # [60-80]
    localsector5=[] # [80-100]
    localsector6=[] # [100-120]
    localsector7=[] # [120-140]
    localsector8=[] # [140-160]
    localsector9=[] # [160-180]
    localsector10=[] # [180-200]
    localsector11=[] # [200-220]
    localsector12=[] # [220-240]
    localsector13=[] # [240-260]
    localsector14=[] # [260-280]
    localsector15= [] # [280-300]
    localsector16= [] #[300-320]
    localsector17= [] # [320-340]
    localsector18= [] # [340-360]
    pt_speeds=speed_up[:,j]
    pt_angle=global_dir[:,j]
    for i in range(0,len(windspeed)):
        # pt_angle=np.abs(folders+global_dir[:,j])
        # idx = np.argmin(np.abs(pt_angle - winddir[i]))
        speed_up_coefficient=np.interp(winddir[i],folders,pt_speeds)
        angle_correction=np.interp(winddir[i],folders,pt_angle)
        angle_correction+=winddir[i]
        tempwindspeed.append(windspeed[i]*speed_up_coefficient)
        if(angle_correction>=0 and angle_correction<=20):
            localsector1.append(windspeed[i]*speed_up_coefficient)
        elif(angle_correction>20 and angle_correction<=40):
            localsector2.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>40 and angle_correction<=60):
            localsector3.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>60 and angle_correction<=80):
            localsector4.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>80 and angle_correction<=100):
            localsector5.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>100 and angle_correction<=120):
            localsector6.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>120 and angle_correction<=140):
            localsector7.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>140 and angle_correction<=160):
            localsector8.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>160 and angle_correction<=180):
            localsector9.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>180 and angle_correction<=200):
            localsector10.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>200 and angle_correction<=220):
            localsector11.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>220 and angle_correction<=240):
            localsector12.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>240 and angle_correction<=260):
            localsector13.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>260 and angle_correction<=280):
            localsector14.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>280 and angle_correction<=300):
            localsector15.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>300 and angle_correction<=320):
            localsector16.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>320 and angle_correction<=340):
            localsector17.append(windspeed[i]*speed_up_coefficient) 
        elif(angle_correction>340 and angle_correction<=360 or angle_correction<0):
            localsector18.append(windspeed[i]*speed_up_coefficient) 
    shape, loc, scale = stats.weibull_min.fit(tempwindspeed, floc=0)
    global_scale[j]=round(scale,1)
    global_shape[j]=round(shape,1)
    # Sector Weibull 
    # Sector 1 
    if(len(localsector1)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector1,floc=0)
        local_scale[0,j]=int(scale*10)
        local_shape[0,j]=int(shape*100)
    # Sector 2 
    if(len(localsector2)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector2, floc=0)
        local_scale[1,j]=int(scale*10)
        local_shape[1,j]=int(shape*100)
    # Sector 3
    if(len(localsector3)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector3, floc=0)
        local_scale[2,j]=int(scale*10)
        local_shape[2,j]=int(shape*100)
    # Sector 4
    if(len(localsector4)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector4, floc=0)
        local_scale[3,j]=int(scale*10)
        local_shape[3,j]=int(shape*100)
    # Sector 5
    if(len(localsector5)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector5,floc=0)
        local_scale[4,j]=int(scale*10)
        local_shape[4,j]=int(shape*100)
    # Sector 6 
    if(len(localsector6)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector6, floc=0)
        local_scale[5,j]=int(scale*10)
        local_shape[5,j]=int(shape*100)
    # Sector 7
    if(len(localsector7)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector7, floc=0)
        local_scale[6,j]=int(scale*10)
        local_shape[6,j]=int(shape*100)
    # Sector 8
    if(len(localsector8)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector8, floc=0)
        local_scale[7,j]=int(scale*10)
        local_shape[7,j]=int(shape*100)
    # Sector 9
    if(len(localsector9)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector9, floc=0)
        local_scale[8,j]=int(scale*10)
        local_shape[8,j]=int(shape*100)
    # Sector 10
    if(len(localsector10)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector10, floc=0)
        local_scale[9,j]=int(scale*10)
        local_shape[9,j]=int(shape*100)
    # Sector 11
    if(len(localsector11)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector11,floc=0)
        local_scale[10,j]=int(scale*10)
        local_shape[10,j]=int(shape*100)
    # Sector 12 
    if(len(localsector12)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector12, floc=0)
        local_scale[11,j]=int(scale*10)
        local_shape[11,j]=int(shape*100)
    # Sector 13
    if(len(localsector13)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector13, floc=0)
        local_scale[12,j]=int(scale*10)
        local_shape[12,j]=int(shape*100)
    # Sector 14
    if(len(localsector14)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector14, floc=0)
        local_scale[13,j]=int(scale*10)
        local_shape[13,j]=int(shape*100)
    # Sector 15
    if(len(localsector15)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector15,floc=0)
        local_scale[14,j]=int(scale*10)
        local_shape[14,j]=int(shape*100)
    # Sector 16 
    if(len(localsector16)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector16, floc=0)
        local_scale[15,j]=int(scale*10)
        local_shape[15,j]=int(shape*100)
    # Sector 17
    if(len(localsector17)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector17, floc=0)
        local_scale[16,j]=int(scale*10)
        local_shape[16,j]=int(shape*100)
    # Sector 18
    if(len(localsector18)>0):
        shape, loc, scale = stats.weibull_min.fit(localsector18, floc=0)
        local_scale[17,j]=int(scale*10)
        local_shape[17,j]=int(shape*100)
    frequency1=float(len(localsector1)/len(windspeed)*10)
    frequency2=float((len(localsector2)/len(windspeed)*10))
    frequency3=float(len(localsector3)/len(windspeed)*10)
    frequency4=float((len(localsector4)/len(windspeed)*10))
    frequency5=float(len(localsector5)/len(windspeed)*10)
    frequency6=float((len(localsector6)/len(windspeed)*10))
    frequency7=float(len(localsector7)/len(windspeed)*10)
    frequency8=float((len(localsector8)/len(windspeed)*10))
    frequency9=float(len(localsector9)/len(windspeed)*10)
    frequency10=float((len(localsector10)/len(windspeed)*10))
    frequency11=float(len(localsector11)/len(windspeed)*10)
    frequency12=float((len(localsector12)/len(windspeed)*10))
    frequency13=float(len(localsector13)/len(windspeed)*10)
    frequency14=float((len(localsector14)/len(windspeed)*10))
    frequency15=float(len(localsector15)/len(windspeed)*10)
    frequency16=float((len(localsector16)/len(windspeed)*10))
    frequency17=float(len(localsector17)/len(windspeed)*10)
    frequency18=float((len(localsector18)/len(windspeed)*10))
    power_density=int(0.5*1.225*global_scale[j]**3)
    wrgstring="%10s%10s%10s"
    wrgstring+="%8s%5s%5s"
    wrgstring+="%6s%15s%3s"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d"
    wrgstring+="%4.1g%4d%5d\n"
    target.write(wrgstring%("GridPt",round(x[j]+xref,2),round(y[j]+yref,2), \
                                round(z[j]-heights[1],2),round(heights[1],1),global_scale[j], \
                                global_shape[j],power_density,18, \
                                frequency1,local_scale[0,j],local_shape[0,j], \
                                frequency2,local_scale[1,j],local_shape[1,j], \
                                frequency3,local_scale[2,j],local_shape[2,j], \
                                frequency4,local_scale[3,j],local_shape[3,j], \
                                frequency5,local_scale[4,j],local_shape[4,j], \
                                frequency6,local_scale[5,j],local_shape[5,j], \
                                frequency7,local_scale[6,j],local_shape[6,j], \
                                frequency8,local_scale[7,j],local_shape[7,j], \
                                frequency9,local_scale[8,j],local_shape[8,j], \
                                frequency10,local_scale[9,j],local_shape[9,j], \
                                frequency11,local_scale[10,j],local_shape[10,j], \
                                frequency12,local_scale[11,j],local_shape[11,j], \
                                frequency13,local_scale[12,j],local_shape[12,j], \
                                frequency14,local_scale[13,j],local_shape[13,j], \
                                frequency15,local_scale[14,j],local_shape[14,j], \
                                frequency16,local_scale[15,j],local_shape[15,j], \
                                frequency17,local_scale[16,j],local_shape[16,j], \
                                frequency18,local_scale[17,j],local_shape[17,j]))

    # target.write(wrgstring%("GridPt",round(x[j],2),round(y[j],2), \
    #                             round(z[j]-heights[1],2),round(heights[1],1),global_scale[i],\
    #                             global_shape[j],power_density,2, \
    #                             frequency0,local_scale[0,j],local_shape[0,j],\
    #                             frequency1,local_scale[1,j],local_shape[1,j]),\
    #                             frequency2,local_scale[2,j],local_shape[2,j],\
    #                             frequency3,local_scale[3,j],local_shape[3,j])
    print(j,len(ind_0))
target.close()




            
    # Weibull 
    # data=np.genfromtxt(str(sys.argv[4]),dtype=float,skip_header=1,delimiter=",",encoding='utf-8')
    # windspeed=data[:,1]

    # import scipy.stats as stats
    # shape, loc, scale = stats.weibull_min.fit(windspeed, floc=0)
    # print(shape,loc,scale)

    # for i in range(0,len(ind_0)):
    #     tempspeed=np.sqrt(u[i]**2+v[i]**2)/uref
    #     tempspeed=windspeed*tempspeed
    #     shape, loc, scale = stats.weibull_min.fit(tempspeed, floc=0)
    #     print(len(ind_0),i,shape,loc,scale)



