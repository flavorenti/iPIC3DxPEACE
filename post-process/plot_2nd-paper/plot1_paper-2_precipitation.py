import fibo as fb
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

### define input parameters ###
data_address = sys.argv[1]
save_address = './images/'
save_name = 'Histo-PR0-minusBz-1D'
cusp_file = 'output_cusps_PR0-plusBz_5800.txt'
Nbinx = 50.
Nbiny = 26.
species = ['e','i']
RemovedParticlesCycle = 2
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
cycle, x, y, z, u, v ,w ,q = np.loadtxt(data_address+'RemovedParticles.txt',unpack=True,skiprows=2)

# cut file 
clim = 5810
cycle = cycle[np.where(cycle<clim)]
x = x[np.where(cycle<clim)]
y = y[np.where(cycle<clim)]
z = z[np.where(cycle<clim)]
u = u[np.where(cycle<clim)]
v = v[np.where(cycle<clim)]
w = w[np.where(cycle<clim)]
q = q[np.where(cycle<clim)]
'''
# each two cycles
cycle = cycle[::2]
x = x[::2]
y = y[::2]
z = z[::2]
u = u[::2]
v = v[::2]
w = w[::2]
q = q[::2]
'''

# remove supra-luminous pcls
slow = np.where((u**2+v**2+w**2)<1)
cycle = cycle[slow]
x = x[slow]
y = y[slow]
z = z[slow]
u = u[slow]
v = v[slow]
w = w[slow]
q = q[slow]

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

# remove pcls with not-too-much radial speed
rad_lim=0.01
index_radiality = np.abs(x*u+y*v+z*w)/np.sqrt(x**2+y**2+z**2)/np.sqrt(u**2+v**2+w**2)
cycle = cycle[np.where(index_radiality>rad_lim)]
x = x[np.where(index_radiality>rad_lim)]
y = y[np.where(index_radiality>rad_lim)]
z = z[np.where(index_radiality>rad_lim)]
u = u[np.where(index_radiality>rad_lim)]
v = v[np.where(index_radiality>rad_lim)]
w = w[np.where(index_radiality>rad_lim)]
q = q[np.where(index_radiality>rad_lim)]

# latitude-longitude pcls
theta = np.arctan2(np.sqrt(x**2+y**2),-z)*180./3.14
phi   = np.arctan2(y,x)*180./3.14

# array of cycle print pcls
maxc = int(np.max(cycle))
minc = int(np.min(cycle))
cycles = range(minc,maxc,RemovedParticlesCycle)
print(minc,maxc,len(cycles))

# kinetic energy pcls [mi*c^2]
Ekin = (u**2+v**2+w**2)

# kinetic energy pcls [eV]
Ekin[q>0] = Ekin[q>0]*T_isw/(vth_i**2)
Ekin[q<0] = Ekin[q<0]*T_esw/(vth_e**2)

# conversion factor flux (from macropcls/solid_angle to pcls/cm2/s)
def convF(lat):
    conv_pcls = n_sw*dx*dy*dz/ppc*(conv_di_to_cm**3)
    deltaT = len(cycles)*dt*(conv_wpi_to_s)
    deltaS = 2.*(np.pi*R*conv_di_to_cm)**2/Nbinx/Nbiny*abs(np.sin(lat*3.14/180.))
    return conv_pcls/deltaT/deltaS

# compute Flux (histogram)
for i,sp in enumerate(species):
    if i==0:
        index = np.where((q<0))
    if i==1:
        index = np.where((q>0))
    histo = plt.hist2d(phi[index],theta[index],bins=[Nbinx,Nbiny],range=[[-180,180],[0,180]])
    Flux = histo[0]
    Temp = np.zeros(np.shape(Flux))

    # limits of the bins, size Nx+1, Ny+1
    binx_lim = histo[1]
    biny_lim = histo[2]

    # center of the bins size Nx, Ny
    binx_cen = np.zeros(int(Nbinx))
    for iix in range(int(Nbinx)):
        binx_cen[iix] = (binx_lim[iix]+binx_lim[iix+1])/2.
    biny_cen = np.zeros(int(Nbiny))
    for iiy in range(int(Nbiny)):
        biny_cen[iiy] = (biny_lim[iiy]+biny_lim[iiy+1])/2.

    # convert units flux [pcls/cm2/s]
    for iy,lat in enumerate(biny_cen):
        print('latitude=',lat,'conversion Flux=',convF(lat))
        Flux[:,iy] *= convF(lat)

        # compute mean energy per bin
        for ix,lon in enumerate(binx_cen):
            if i==0:
                index = np.where((q<0)*(phi<binx_lim[ix+1])*(phi>=binx_lim[ix])*(theta<biny_lim[iy+1])*(theta>=biny_lim[iy]))
            if i==1:
                index = np.where((q>0)*(phi<binx_lim[ix+1])*(phi>=binx_lim[ix])*(theta<biny_lim[iy+1])*(theta>=biny_lim[iy]))
            if Flux[ix,iy]>1 :
                Temp[ix,iy]=np.mean(Ekin[index])
            else:
                Flux[ix,iy]=1.
                Temp[ix,iy]=1e-3


    ### PRINT INTEGRATED QUANTITIES ###
    sum_Flux = 0.
    for iy,lat in enumerate(biny_cen):
        deltaS = 2.*(np.pi*R*conv_di_to_cm)**2/Nbinx/Nbiny*abs(np.sin(lat*3.14/180.))
        sum_Flux += np.sum(Flux[:,iy])*Nbinx*deltaS 
    print('GOOD Integrated Flux_'+sp+'=%.4e #pcls/s'%sum_Flux)

    print('BAD Integrated Flux_'+sp+'=%.4e #pcls/s'%np.sum(Flux))
    print('BAD Integrated Energy_'+sp+'=%.4e #eV'%np.sum(Temp))
    print('BAD Integrated EnergyFlux_'+sp+'=%.4e #eV/s'%np.sum(Flux*Temp))

    ### PLOT THE DATA ###

    # load the cusps
    phi_c,theta_c,comment = np.loadtxt(cusp_file,unpack=True,dtype='str')
    lot_cusp = np.array(phi_c,dtype='float')*24./(2*np.pi)
    lat_cusp = 90.-(np.array(theta_c,dtype='float')*180./np.pi)
    x_cusp=[]
    y_cusp=[]
    for i in range(0,len(phi_c)):
        if lot_cusp[i] not in x_cusp:
            if len(lat_cusp[(phi_c==phi_c[i])*(lat_cusp>0)])>0:
                x_cusp.append(lot_cusp[i])
                y_cusp.append(np.min(lat_cusp[(phi_c==phi_c[i])*(lat_cusp>0)]))
            if len(lat_cusp[(phi_c==phi_c[i])*(lat_cusp<0)])>0:
                x_cusp.append(lot_cusp[i])
                y_cusp.append(np.max(lat_cusp[(phi_c==phi_c[i])*(lat_cusp<0)]))
    x_cusp = np.array(x_cusp)
    y_cusp = np.array(y_cusp)
    x_cuspN = x_cusp[np.where(y_cusp>0)]
    y_cuspN = y_cusp[np.where(y_cusp>0)]
    x_cuspS = x_cusp[np.where(y_cusp<0)]
    y_cuspS = y_cusp[np.where(y_cusp<0)]

    # local time - x-axis [0,24] hours
    lot = (binx_cen+180.)*24./360.
    # latitude - y-axis [-90,90] deg
    lat = biny_cen-90.
    # mesh grid 
    Lot,Lat = np.meshgrid(lot,lat)

    fig = plt.figure(figsize=(7,6))
    plt.axis('off')
    ax1 = fig.add_subplot(211)#, projection='mollweide')
    if (sp=='i'):
        plt.title('Precipitating protons Flux',fontsize=14)
    if (sp=='e'):
        plt.title('Precipitating electrons Flux',fontsize=14)
    #ax1.pcolormesh(Lot,Lat,np.random.rand(Nbiny, Nbinx),cmap=cm.Reds)
    plt.imshow(np.transpose(Flux),cmap=cm.Reds,norm=matplotlib.colors.LogNorm(),extent=(0,24,-90,90),origin='lower',aspect=0.07)
    cb=plt.colorbar()
    plt.clim(1e6,1e9)
    cb.set_label(r'$\frac{pcls}{cm^2 s}$',rotation=0,fontsize=16,labelpad=18)
    cb.ax.tick_params(labelsize=14)
    plt.ylabel(r'$Latitude$ (deg.)',fontsize=14)
    plt.plot(x_cuspN,y_cuspN,linestyle='--',marker='',color='black')
    plt.plot(x_cuspS,y_cuspS,linestyle='--',marker='',color='black')
    plt.xlim(0,24)
    plt.ylim(-90,90)
    plt.xticks([0,6,12,18,24])
    plt.yticks([-90,-60,-30,0,30,60,90])
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
 
    ax2 = fig.add_subplot(212)#, projection='mollweide')
    if (sp=='i'):
        plt.title('Precipitating protons Energy',fontsize=14)
    if (sp=='e'):
        plt.title('Precipitating electrons Energy',fontsize=14)
    #ax2.pcolormesh(Lot,Lat,np.transpose(Temp),cmap=cm.coolwarm)
    plt.imshow(np.transpose(Temp),cmap=cm.coolwarm,norm=matplotlib.colors.LogNorm(),extent=(0,24,-90,90),origin='lower',aspect=0.07)
    cb=plt.colorbar()
    plt.clim(2e1,4e4)
    cb.ax.tick_params(labelsize=14)
    cb.set_label(r'$eV$',rotation=0,fontsize=16,labelpad=18)
    plt.ylabel(r'$Latitude$ (deg.)',fontsize=14)
    plt.xlabel('Local Time',fontsize=14)
    plt.plot(x_cuspN,y_cuspN,linestyle='--',marker='',color='black')
    plt.plot(x_cuspS,y_cuspS,linestyle='--',marker='',color='black')
    plt.xlim(0,24)
    plt.ylim(-90,90)
    plt.xticks([0,6,12,18,24])
    plt.yticks([-90,-60,-30,0,30,60,90])
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
 

    plt.tight_layout()
    if sp=='i':
        plt.savefig(save_address+save_name+'-new_protons_0-inf.png',format='png')
    if sp=='e':
        plt.savefig(save_address+save_name+'-new_electrons_0-inf.png',format='png')
    plt.close()
