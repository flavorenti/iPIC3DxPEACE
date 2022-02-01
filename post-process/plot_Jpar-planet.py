import fibo as fb
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from vectorfield import *

cycle = 6000
compute_Flag = False
imaging_Flag = True

######################################################
#----path-to-run-------------------
data_address = sys.argv[1]

#----create-your-objects-----------
fibo_name = 'ipic3d'
alldata = fb.fibo(fibo_name)
from_VTK = fb.from_VTK(data_address)

#----load-metadata------------------
from_VTK.get_meta(silent=False)
alldata.meta = from_VTK.meta
x = np.linspace(0,from_VTK.meta['xl'],from_VTK.meta['nx']) 
y = np.linspace(0,from_VTK.meta['yl'],from_VTK.meta['ny']) 
z = np.linspace(0,from_VTK.meta['zl'],from_VTK.meta['nz'])
xc = from_VTK.meta['xc']
yc = from_VTK.meta['yc']
zc = from_VTK.meta['zc']+from_VTK.meta['Doff']
x = (-x + xc )/from_VTK.meta['R']
y = (-y + yc )/from_VTK.meta['R']
z = ( z - zc )/from_VTK.meta['R']
dx = from_VTK.meta['dx']/from_VTK.meta['R']
dy = from_VTK.meta['dy']/from_VTK.meta['R']
dz = from_VTK.meta['dz']/from_VTK.meta['R']
i_good = []
lot = []
lat = []
sph = []
Jpar_good = []

#----miscellaneous------------------
def scalr(a):
    return np.sqrt(a[0]**2+a[1]**2+a[2]**2)

#----load-Jpar------------------
if compute_Flag:
    print('LOAD')
    conv_J = -(4.*np.pi*30.)*(400./0.027)*(0.16) #nA/m^2
    Jparx,Jpary,Jparz = conv_J*from_VTK.get_vect(alldata.meta['name']+'_Jpargf3_'+str(cycle),cycle,fibo_obj=None,tar_var='JparB',silent=True)

#----loop-on-the-box------------
    print('COMPUTE')
    for ix in range(len(x)):
        for iy in range(len(y)):
            for iz in range(len(z)):
                rad = np.sqrt(x[ix]**2+y[iy]**2+z[iz]**2)
                rho = np.sqrt(x[ix]**2+y[iy]**2)
                senphi = y[iy]/rho
                cosphi = x[ix]/rho
                costhe = z[iz]/rad
                if 1 < rad < (1+2.*dx) :
                    i_good.append([ix,iy,iz])
                    sph = np.append(sph,senphi)
                    lot = np.append(lot,cosphi)
                    lat = np.append(lat,costhe)
                    Jpar_good = np.append(Jpar_good,(Jparx[ix,iy,iz]*x[ix]+Jpary[ix,iy,iz]*y[iy]+Jparz[ix,iy,iz]*z[iz])/rad)

    #----print-to-text-----
    print('PRINT')
    out = open('output_plot_Jpargf3.txt','w')
    out.write('# sen(phi)\t cos(phi)\t cos(theta)\t Jpar[nA/m2]')
    for i in range(len(lot)):
        out.write('%.4f\t%.4f\t%.4f\t%.4f\n'%(sph[i],lot[i],lat[i],Jpar_good[i]))
    out.close()


#----sample-on-image---------
Nx=50
Ny=50
grid_x = np.linspace(-1,1,Nx)
grid_y = np.linspace(0,1,Ny)
maap = np.zeros([Nx,Ny])
if imaging_Flag:
    senphi, cosphi,costhe,Jpar = np.loadtxt('output_plot_Jpargf3.txt',unpack=True)
    for ix in range(0,Nx-1):
        for iy in range(0,Ny-1):
            igood = (costhe<grid_y[iy+1])*(costhe>grid_y[iy])*(cosphi<grid_x[ix+1])*(cosphi>grid_x[ix])*(senphi>0)
            maap[ix,iy] = np.mean(Jpar[igood])
    vectorfield(maap,grid_x[:Nx-1],grid_y[:Ny-1],'cos(phi)','cos(theta)','Jparallel map planet DAWN','./images/Jpar_planet_dawn_map_N1',filter_sigma=0,clim_min=-3000,clim_max=3000,Log=False)

