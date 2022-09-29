import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import numpy as np
import glob
import os
import sys
import h5py as h5
import fibo as fb

###################################################################
######################## DEFINE PARAMETERS ########################
###################################################################
# path to simu pcls data
data_path = sys.argv[1]+'dataP/'
cycle = 5999
nrun = 2
# path to trajectory files and text with list
traj_path = './trajectory_files/'+sys.argv[2]+'-flyby-1st_MSO-trajectory1.txt'
good_path = sys.argv[1]+'texts/output_dataP-new_'+sys.argv[2]+'.txt'
output_file = open(sys.argv[1]+'texts/output_good-pcls_%i_%i_'%(nrun,cycle)+sys.argv[2]+'.txt','w')
# width sphere of integration [R]
dr = 0.03
###################################################################
######################## END DEFINITION ###########################
###################################################################

# header output file
output_file.write('# s/c trajectory from file '+traj_path+'\n')
output_file.write('# list good procs from file '+good_path+'\n')
output_file.write('# procs found in path '+data_path+'nrun=%i cycle=%i\n'%(nrun,cycle))
output_file.write('# sphere of radius dr='+str(dr)+'\n')
output_file.write('# coordinates of the particles in MSO ref. system\n\n')
output_file.write('# time[min]\t xp[R]\t yp[R]\t zp[R]\t log10(E/eV) \n')

# conversion factors from code units to SI values
conv_T = 21.5/(0.031**2) #eV 

# miscellaneous
from_VTK = fb.from_VTK(data_path)
from_VTK.get_meta(silent=False)
x = np.linspace(0,from_VTK.meta['xl'],from_VTK.meta['nx'])
y = np.linspace(0,from_VTK.meta['yl'],from_VTK.meta['ny'])
z = np.linspace(0,from_VTK.meta['zl'],from_VTK.meta['nz'])
xc = from_VTK.meta['xc']
yc = from_VTK.meta['yc']
zc = from_VTK.meta['zc']+from_VTK.meta['Doff']
x = (-x + xc )/from_VTK.meta['R']
y = (-y + yc )/from_VTK.meta['R']
z = ( z - zc )/from_VTK.meta['R']

# open the trajectory file
tt,xx,yy,zz = np.loadtxt(traj_path,unpack=True)
dd = np.sqrt(xx**2+yy**2+zz**2)
i_CA = np.where(dd==min(dd))
tt = (tt-tt[i_CA])/60.

# loop on the trajectory file (get in and out limits)
enter=0
for ii in range(0,len(tt)):
    if((xx[ii]>min(x))*(yy[ii]>min(y))*(zz[ii]>min(z))*(xx[ii]<max(x))*(yy[ii]<max(y))*(zz[ii]<max(z))):
        if(enter==0):
            imin=ii
        enter=1
    elif(enter==1):
        imax=ii
        break

# cut trajectory to inside the box
tt=tt[imin:imax]
xx=xx[imin:imax]
yy=yy[imin:imax]
zz=zz[imin:imax]
dd=dd[imin:imax]

# open file listing good data
good_data_files = np.loadtxt(good_path,unpack=True,dtype='int')

# module function
def scalr(a):
    return np.sqrt(a[0]**2+a[1]**2+a[2]**2)

# open particle files
for good_id in good_data_files:

    # file
    d = h5.File(data_path+'restart%i_%i.hdf'%(good_id,nrun),'r')

    # positions
    x_pcls = d['particles/species_0/x/cycle_'+str(cycle)][:]
    y_pcls = d['particles/species_0/y/cycle_'+str(cycle)][:]
    z_pcls = d['particles/species_0/z/cycle_'+str(cycle)][:]

    # velocity
    u_pcls = d['particles/species_0/u/cycle_'+str(cycle)][:]
    v_pcls = d['particles/species_0/v/cycle_'+str(cycle)][:]
    w_pcls = d['particles/species_0/w/cycle_'+str(cycle)][:]

    # energies
    e_pcls = conv_T*(u_pcls**2+v_pcls**2+w_pcls**2)

    # print
    print('###\n###\nUnpack file %i #pcls loaded = %i\n###\n###'%(good_id,len(x_pcls)))

    # delete array not needed
    del d, u_pcls, v_pcls, w_pcls

    # change units x,y,z pcls to MSO
    x_pcls = (-x_pcls + xc )/from_VTK.meta['R']
    y_pcls = (-y_pcls + yc )/from_VTK.meta['R']
    z_pcls = ( z_pcls - zc )/from_VTK.meta['R']

    # loop on the trajectory file - compute good pcls
    for ii in range(0,len(tt)-1):

        versor = [xx[ii+1]-xx[ii],yy[ii+1]-yy[ii],zz[ii+1]-zz[ii]]
        versor = versor/scalr(versor)

        distance = np.sqrt( (1.-versor[0]**2)*(xx[ii]-x_pcls)**2 + (1.-versor[1]**2)*(yy[ii]-y_pcls)**2 + (1.-versor[2]**2)*(zz[ii]-z_pcls)**2 )

        idx_sphere = np.where(distance<dr)

        if len(idx_sphere[0])>0 :
            print(xx[ii],yy[ii],zz[ii])
            print('G0D ii=%i #pcls inside sphere=%i'%(ii,len(idx_sphere[0])))
            for idx in idx_sphere[0]: 
                output_file.write('%.2f\t %.6f\t %.6f\t %.6f\t %.8f\n'%(tt[ii],x_pcls[idx],y_pcls[idx],z_pcls[idx],np.log10(e_pcls[idx])))

        elif ii%100==0:
            print(xx[ii],yy[ii],zz[ii])
            print('BAD ii=%i #pcls inside sphere=%i'%(ii,len(idx_sphere[0])))

    del x_pcls, y_pcls, z_pcls, e_pcls

# close the file
output_file.close()
