import fibo as fb
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

###############################
### define input parameters ###
###############################
data_address = sys.argv[1]
save_name = 'PR1-minusBz-e-'
# theta [0,180] --> [sud,nord]
# phi [-180,180] --> [left,right]
theta_lim = [90,95]
phi_lim = [-135,-130]
species = 'e'
###############################

# read metadata
from_VTK = fb.from_VTK(data_address)
from_VTK.get_meta(silent=True)

# read file with particles
cycle, x, y, z, u, v ,w ,q = np.loadtxt(data_address+'RemovedParticles3.txt',unpack=True,skiprows=2)
print('STEP0 len(cycle)=',len(cycle))

# cut file in time
clim = cycle[0]+35
cycle = cycle[np.where(cycle<clim)]
x = x[np.where(cycle<clim)]
y = y[np.where(cycle<clim)]
z = z[np.where(cycle<clim)]
u = u[np.where(cycle<clim)]
v = v[np.where(cycle<clim)]
w = w[np.where(cycle<clim)]
q = q[np.where(cycle<clim)]
print('STEP1 t<tmax len(cycle)=',len(cycle))

# remove supra-luminous pcls
slow = np.where(((u**2+v**2+w**2)<0.3))#*((u**2+v**2+w**2)>0.18)) #add 4keV filter
cycle = cycle[slow]
x = x[slow]
y = y[slow]
z = z[slow]
u = u[slow]
v = v[slow]
w = w[slow]
q = q[slow]
print('STEP2 v<c/3 len(cycle)=',len(cycle))

# choose one species
if species=='e':
    index_q = np.where(q<0)
else:
    if species=='i':
        index_q = np.where(q>0)
    else:
        print('ERROR wrong species specification -->'+species)
cycle = cycle[index_q]
x = x[index_q]
y = y[index_q]
z = z[index_q]
u = u[index_q]
v = v[index_q]
w = w[index_q]
q = q[index_q]
print('STEP3 q<0 len(cycle)=',len(cycle))


# latitude-longitude pcls
xc = from_VTK.meta['xc']
yc = from_VTK.meta['yc']
zc = from_VTK.meta['zc']+from_VTK.meta['Doff']
R = from_VTK.meta['R']
x_mso = (-x + xc )/R
y_mso = (-y + yc )/R
z_mso = ( z - zc )/R
theta = np.arctan2(np.sqrt(x_mso**2+y_mso**2),-z_mso)*180./3.14
phi   = np.arctan2(y_mso,x_mso)*180./3.14

# apply filter region
index_f = np.where((theta>theta_lim[0])*(theta<theta_lim[1])*(phi>phi_lim[0])*(phi<phi_lim[1]))
cycle = cycle[index_f]
x = x[index_f]
y = y[index_f]
z = z[index_f]
u = u[index_f]
v = v[index_f]
w = w[index_f]
q = q[index_f]
print('STEP4 (x,y,z) len(cycle)=',len(cycle))

# Open .txt files for output
file_name = './pcls/'+save_name+'_%s-%s_%s-%s.txt'%(theta_lim[0],theta_lim[1],phi_lim[0],phi_lim[1])
ff = open(file_name,'w')
for i in range(0,len(q)):
    ff.write('%i\t'%cycle[i])
    ff.write('%.5f\t'%x[i])
    ff.write('%.5f\t'%y[i])
    ff.write('%.5f\t'%z[i])
    ff.write('%.8f\t'%u[i])
    ff.write('%.8f\t'%v[i])
    ff.write('%.8f\t'%w[i])
    ff.write('%.5e\n'%q[i])
ff.close()
