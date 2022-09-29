import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import sys
import fibo as fb

######################################################
################ cross sections EII ##################
######################################################
Name = ['H','Li','Na','K','Rb','Cs','He','Ne','Ar','Kr','Xe','Ti','Mn','Fe','Ni','Cu','Pd','Ag','W','Au','Be','Mg','Al','Si','Hg','Pb','U','Ca']
I = [13.6,5.4,5.1,4.3,4.2,3.9,24.6,21.6,15.8,14.0,12.1,6.8,7.4,7.9,7.7,7.7,8.3,7.6,8.0,9.2,9.3,7.6,6.0,8.2,10.4,7.4,5.7,6.1]
alpha = [0.83,5.72,9.56,6.54,4.83,3.87,0.37,0.37,2.92,3.51,4.30,19.1,8.39,14.80,6.04,6.86,3.09,7.65,7.12,16.5,3.22,13.7,11.6,9.97,1.00,12.80,5.04,11.45]
beta = [0.35,0.50,0.52,0.36,0.21,0.13,0.29,0.14,0.29,0.27,0.26,0.65,0.41,1.15,0.41,0.65,0.15,0.57,0.38,0.27,0.34,0.71,0.34,0.50,0.22,0.59,0.33,0.41]
gamma = [1.91,1.67,1.90,1.57,1.82,1.81,1.91,2.00,1.86,1.80,1.76,1.85,1.62,1.44,1.86,1.52,1.89,1.46,1.62,1.86,2.20,1.87,1.80,1.61,1.74,1.52,1.73,2.22]
def cross_section_EII(eps,sp):
    if sp in Name:
        i = np.where(np.asarray(Name) == sp)[0][0]
        return 1e-16*alpha[i]*((eps/I[i])-1.)/np.power((1.+(beta[i]*((eps/I[i])-1.))),gamma[i])
    elif sp=='O':
        en_arr, sigma_arr = np.loadtxt('Sigma_EII-O.txt',unpack=True,usecols=(0,1))
        f = interpolate.interp1d(en_arr, sigma_arr, kind='linear',bounds_error=False,fill_value=0.)
        return 1e-16*f(eps)
    else:
        print('EII Cross section for '+sp+' NOT FOUND!!!')
        return None
######################################################
######################################################
######################################################

######################################################
################ cross sections XRF ##################
######################################################
path = '/home/flavorenti/Downloads/NIST164/'
name_file = 'NIST164_XRF_'
def cross_section_XRF(eps,sp):
    en_arr,sigma_arr = np.loadtxt(path+name_file+sp+'.txt',unpack=True,usecols=(0,1))
    f = interpolate.interp1d(en_arr, sigma_arr, kind='linear',bounds_error=False,fill_value=0.)
    return f(eps)
######################################################
######################################################
######################################################

### define input parameters ###
data_address = sys.argv[1]
species = 'Na'
save_name = 'map_PR1-plusBz_XRF-'+species
Nbinx = 30 
Nbiny = 15
Nener = 400
ppc = 64.
n_sw = 30.
T_isw = 21.5
T_esw = 21.5
conv_di_to_cm = 23.7e5
conv_wpi_to_s = 2.8e-3
###############################

# read metadata
from_VTK = fb.from_VTK(data_address)
from_VTK.get_meta(silent=False)

# read file with particles
cycle, x, y, z, u, v ,w ,q = np.loadtxt(data_address+'RemovedParticles3.txt',unpack=True,skiprows=2)

# remove supra-luminous pcls
slow = np.where((u**2+v**2+w**2)<0.3)
cycle = cycle[slow]
x = x[slow]
y = y[slow]
z = z[slow]
u = u[slow]
v = v[slow]
w = w[slow]
q = q[slow]

# select only electrons
x = x[q<0]
y = y[q<0]
z = z[q<0]
u = u[q<0]
v = v[q<0]
w = w[q<0]
q = q[q<0]

# define useful metadata
dx = from_VTK.meta['dx']
dy = from_VTK.meta['dy']
dz = from_VTK.meta['dz']
dt = from_VTK.meta['dt']
xc = from_VTK.meta['xc']
yc = from_VTK.meta['yc']
zc = from_VTK.meta['zc']+from_VTK.meta['Doff']
R = from_VTK.meta['R']
vth_e = from_VTK.meta['vths'][0]
vth_i = from_VTK.meta['vths'][1]

# convert pcls position to MSO coord.
x = (-x + xc )/R
y = (-y + yc )/R
z = ( z - zc )/R

# latitude-longitude pcls
theta = np.arctan2(np.sqrt(x**2+y**2),-z)*180./3.14
phi   = np.arctan2(y,x)*180./3.14

# array of cycle print pcls
cycles = np.unique(cycle) 
print('MIN Cycles=',min(cycles))
print('MAX Cycles=',max(cycles))
print('LEN Cycles=',len(cycles))
print('LEN Electrons=',len(cycle[np.where(q<0)]))
print('\n')

# change units speed [mic^2]-->[eV]
u = u*np.sqrt(T_esw)/vth_e
v = v*np.sqrt(T_esw)/vth_e
w = w*np.sqrt(T_esw)/vth_e

# kinetic energy pcls [eV]
Ekin = (u**2+v**2+w**2)
print('MIN pcls Energy=',min(Ekin))
print('MAX pcls Energy=',max(Ekin))
energy_arr = np.linspace(0.,max(Ekin),Nener)
dE = (max(energy_arr)-0.)/Nener

# compute cross_section
if 'XRF' in save_name:
    sigma_arr = cross_section_XRF(energy_arr,species)
if 'EII' in save_name:
    sigma_arr = cross_section_EII(energy_arr,species)

# initializ void histogram
histo = np.zeros((Nbinx,Nbiny))
binx = np.linspace(-180,180,Nbinx+1)
biny = np.linspace(0,180,Nbiny+1)

# conversion factor flux (from macropcls/solid_angle to pcls/cm2/s)
def convF(lat):
    conv_pcls = n_sw*dx*dy*dz/ppc*(conv_di_to_cm**3)
    deltaT = len(cycles)*dt*(conv_wpi_to_s)
    deltaS = 2.*(np.pi*R*conv_di_to_cm)**2/Nbinx/Nbiny*abs(np.sin(lat*3.14/180.))
    return conv_pcls/deltaT/deltaS

# compute Flux*sigma (histogram)
for ix in range(0,Nbinx):
    for iy in range(0,Nbiny):

        # particles in given spatial bin
        index = np.where((phi>=binx[ix])*(phi<binx[ix+1])*(theta>=biny[iy])*(theta<biny[iy+1]))
        Ekin_1 = Ekin[index]
        hist_ener_1, temp = np.histogram(Ekin_1,bins=np.linspace(0,max(Ekin),Nener+1))
        
        # plot histo energy for check
        '''
        plt.plot(energy_arr,hist_ener_1)
        plt.plot(energy_arr,sigma_arr/np.max(sigma_arr))
        plt.show()
        '''

        # compute rate as integral
        rate = np.sum(hist_ener_1*sigma_arr)*dE*convF((biny[iy]+biny[iy+1])/2.)

        # fill histogram
        histo[ix,iy] = rate


# center of the bins size Nx, Ny
binx_cen = np.zeros(int(Nbinx))
for iix in range(int(Nbinx)):
    binx_cen[iix] = (binx[iix]+binx[iix+1])/2.
biny_cen = np.zeros(int(Nbiny))
for iiy in range(int(Nbiny)):
    biny_cen[iiy] = (biny[iiy]+biny[iiy+1])/2.

# local time - x-axis [0,24] hours
lot = (binx_cen+180.)*24./360.
# latitude - y-axis [-90,90] deg
lat = biny_cen-90.

# Open .txt files for output
file_name='./texts_rates/'+save_name+'.txt'

ff = open(file_name,'w')

# Header file:
ff.write('#\n')
if 'EII' in save_name:
    ff.write('# Electrons precipitation process EII for species '+species+'\n')
if 'XRF' in save_name:
    ff.write('# Electrons precipitation process XRF for species '+species+'\n')
if 'plusBz' in save_name:
    ff.write('# Simulation setup: northward IMF, same Solar wind as Aizawa et al. (2021)\n')
if 'minusBz' in save_name:
    ff.write('# Simulation setup: southward IMF, same Solar wind as Aizawa et al. (2021)\n')
ff.write('# Simulation code: iPIC3D\n')
ff.write('# Simulation path: '+data_address+'\n')
ff.write('# Signature: Federico Lavorenti\n')
ff.write('#\n')
ff.write('nan\t')

# def separator
sep = '\t'

# print first row 
for ix,llot in enumerate(lot):
    ff.write('%.3f'%llot+sep)
ff.write('\n')

# print histogram
for iy,llat in enumerate(lat):
    ff.write('%.3f'%llat+sep)
    for ix,llot in enumerate(lot):
        ff.write('%.3e'%(histo[ix,iy])+sep)
    ff.write('\n')

# close file
ff.close()
