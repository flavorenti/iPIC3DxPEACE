import fibo as fb
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

### define input parameters ###
data_address = sys.argv[1]
save_name = 'map_PR1-minusBz_0-inf_xLiz'
Nbinx = 50  #200
Nbiny = 25  #101
species = ['e','i']
filter_region = ''#CN,CS,DR,DayNorth,DaySouth,Dawn,Dusk,Right,Left,Up,Down
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

# cut file 
clim = cycle[0]+1000000#35
cycle = cycle[np.where(cycle<clim)]
x = x[np.where(cycle<clim)]
y = y[np.where(cycle<clim)]
z = z[np.where(cycle<clim)]
u = u[np.where(cycle<clim)]
v = v[np.where(cycle<clim)]
w = w[np.where(cycle<clim)]
q = q[np.where(cycle<clim)]

# remove supra-luminous pcls
slow = np.where((u**2+v**2+w**2)<0.3)
#slow = np.where(((q<0)*((u**2+v**2+w**2)<0.0045))+((q>0)*((u**2+v**2+w**2)<0.000225))) #add 0-100eV+0-500eV filter
#slow = np.where(((q<0)*((u**2+v**2+w**2)>0.0045)*((u**2+v**2+w**2)<0.18))+((q>0)*((u**2+v**2+w**2)>0.000225)*((u**2+v**2+w**2)<0.000675))) #add 100-4000eV+500-1500eV filter
#slow = np.where(((u**2+v**2+w**2)<0.3)*(((q<0)*((u**2+v**2+w**2)>0.18))+((q>0)*((u**2+v**2+w**2)>0.000675)))) # add 4keV-1.5keV filter 
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

# remove pcls with not-too-much radial speed (0 remove nothing, 1 remove everything)
rad_lim=0.0
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

# apply filter region (if needed)
if filter_region is not '':
    if filter_region=='CN':
        index_f = np.where((theta>160.)*(theta<170.)*(phi>-5.)*(phi<5.))
        #index_f = np.where((theta>150.)*(phi>-45.)*(phi<45.))
    if filter_region=='CS':
        index_f = np.where((theta>10.)*(theta<20.)*(phi>-5.)*(phi<5.))
        #index_f = np.where((theta<60.)*(phi>-90.)*(phi<90.))
    if filter_region=='DR':
        index_f = np.where((theta>85.)*(theta<95.)*(phi>-55.)*(phi<-45.))
        #index_f = np.where((theta>60.)*(theta<150.)*(phi<-45.))
    if filter_region=='DayNorth':
        index_f = np.where((theta>90.)*(phi>-90.)*(phi<90.))
    if filter_region=='DaySouth':
        index_f = np.where((theta<90.)*(phi>-90.)*(phi<90.))
    if filter_region=='Dawn':
        index_f = np.where((phi<-90.))
    if filter_region=='Dusk':
        index_f = np.where((phi>90.))
    if filter_region=='Right':
        index_f = np.where((phi>0.))
    if filter_region=='Left':
        index_f = np.where((phi<0.))
    if filter_region=='Up':
        index_f = np.where((theta>90.))
    if filter_region=='Down':
        index_f = np.where((theta<90.))
    theta = theta[index_f]
    phi = phi[index_f]
    cycle = cycle[index_f]
    x = x[index_f]
    y = y[index_f]
    z = z[index_f]
    u = u[index_f]
    v = v[index_f]
    w = w[index_f]
    q = q[index_f]


# array of cycle print pcls
cycles = np.unique(cycle) 
print('MIN Cycles=',min(cycles))
print('MAX Cycles=',max(cycles))
print('LEN Cycles=',len(cycles))
print('LEN Protons=',len(cycle[np.where(q>0)]))
print('LEN Electrons=',len(cycle[np.where(q<0)]))
print('\n\n')

# change units speed [mic^2]-->[eV]
u[q>0] = u[q>0]*np.sqrt(T_isw)/vth_i
v[q>0] = v[q>0]*np.sqrt(T_isw)/vth_i
w[q>0] = w[q>0]*np.sqrt(T_isw)/vth_i
u[q<0] = u[q<0]*np.sqrt(T_esw)/vth_e
v[q<0] = v[q<0]*np.sqrt(T_esw)/vth_e
w[q<0] = w[q<0]*np.sqrt(T_esw)/vth_e

# kinetic energy pcls [eV]
Ekin = (u**2+v**2+w**2)
Ekin_rad = (u*x+v*y+w*z)**2/(x**2+y**2+z**2)

# conversion factor flux (from macropcls/solid_angle to pcls/cm2/s)
def convF(lat):
    conv_pcls = n_sw*dx*dy*dz/ppc*(conv_di_to_cm**3)
    deltaT = len(cycles)*dt*(conv_wpi_to_s)
    deltaS = 2.*(np.pi*R*conv_di_to_cm)**2/Nbinx/Nbiny*abs(np.sin(lat*3.14/180.))
    return conv_pcls/deltaT/deltaS

# integrate flux over surface
def integral(Flux,binx_cent,biny_cent):
    sum_Flux = 0.
    for iy,lat in enumerate(biny_cen):
        deltaS = 2.*(np.pi*R*conv_di_to_cm)**2/Nbinx/Nbiny*abs(np.sin(lat*3.14/180.))
        for ix,lon in enumerate(binx_cent):
            sum_Flux += np.sum(Flux[ix,iy])*Nbinx*deltaS
    return sum_Flux 

# compute Flux (histogram)
for i,sp in enumerate(species):
    if i==0:
        index_q = np.where((q<0))
    if i==1:
        index_q = np.where((q>0))
    histo = plt.hist2d(phi[index_q],theta[index_q],bins=[Nbinx,Nbiny],range=[[-180,180],[0,180]])
    Flux = histo[0]
    
    # remove low #pcls cells
    Flux[np.where(Flux<64)]==np.nan

    # init empty arrays
    Temp = np.zeros(np.shape(Flux))
    Pram = np.zeros(np.shape(Flux))
    Ener = np.zeros(np.shape(Flux))
    ERad = np.zeros(np.shape(Flux))

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
        if iy%2==0:
            print('latitude=',lat,'conversion Flux=',convF(lat))
        Flux[:,iy] *= convF(lat)

        # compute mean energy per bin
        for ix,lon in enumerate(binx_cen):
            if i==0:
                index = np.where((q<0)*(phi<binx_lim[ix+1])*(phi>=binx_lim[ix])*(theta<biny_lim[iy+1])*(theta>=biny_lim[iy]))
            if i==1:
                index = np.where((q>0)*(phi<binx_lim[ix+1])*(phi>=binx_lim[ix])*(theta<biny_lim[iy+1])*(theta>=biny_lim[iy]))
            if Flux[ix,iy]>0 :
                Vavg_x=np.mean(u[index])
                Vavg_y=np.mean(v[index])
                Vavg_z=np.mean(w[index])
                Ekin_trasl = ((u[index]-Vavg_x)**2+(v[index]-Vavg_y)**2+(w[index]-Vavg_z)**2)

                Pram[ix,iy]=(Vavg_x**2+Vavg_y**2+Vavg_z**2)
                Temp[ix,iy]=np.mean(Ekin_trasl)
                Ener[ix,iy]=np.mean(Ekin[index])
                ERad[ix,iy]=np.mean(Ekin_rad[index])
            else:
                Pram[ix,iy]=np.nan
                Temp[ix,iy]=np.nan
                Ener[ix,iy]=np.nan
                ERad[ix,iy]=np.nan

    # local time - x-axis [0,24] hours
    if 'xLiz' in save_name:
        lot = (binx_cen+180.)
    else:
        lot = (binx_cen+180.)*24./360.
    # latitude - y-axis [-90,90] deg
    lat = biny_cen-90.
    # print lat lot
    print(lot)
    print(lat)

    # Open .txt files for output
    files_names=[]
    for string in ['Flux','Ener']:#'Temp','Pram','Ener','ERad']:
        if 'xLiz' in save_name:
            files_names.append('./'+save_name+'_'+string+'-'+sp+filter_region+'.csv')
        else:
            files_names.append('./texts/'+save_name+'_'+string+'-'+sp+filter_region+'.txt')

    # loop over files
    for ff_name in files_names:

        ff = open(ff_name,'w')

        # Header file:
        ff.write('#\n')
        if sp=='i':
            ff.write('# Protons precipitation map\n')
        if sp=='e':
            ff.write('# Electrons precipitation map\n')
        if 'plusBz' in save_name:
            ff.write('# Simulation setup: northward IMF, same Solar wind as Aizawa et al. (2021)\n')
        if 'minusBz' in save_name:
            ff.write('# Simulation setup: southward IMF, same Solar wind as Aizawa et al. (2021)\n')
        ff.write('# Simulation code: iPIC3D\n')
        ff.write('# Simulation path: '+data_address+'\n')
        '''
        if 'Flux' in ff_name:
            ff.write('# Output: Flux, units 1/cm2/s\n')
            ff.write('# Output minimum: 1e-6\n')
        else:
            if 'Ener' in ff_name:
                ff.write('# Output: Mean energy (temp+ram) pcls per bin, units eV\n')
            elif 'Temp' in ff_name:
                ff.write('# Output: Mean energy (Temperature) pcls per bin, units eV\n')
            elif 'Pram' in ff_name:
                ff.write('# Output: Mean energy (Ram pressure) pcls per bin, units eV\n')
            elif 'ERad' in ff_name:
                ff.write('# Output: Mean energy (radial kinetic) pcls per bin, units eV\n')
            ff.write('# Output minimum: 1e-3\n')
        '''
        if 'xLiz' in ff_name:
            ff.write('# Grid X-Axis (first row): local time [deg] from 0 to 360 with %i equispaced bins\n'%Nbinx)
        else:
            ff.write('# Grid X-Axis (first row): local time [hours] from 0 to 24 with %i equispaced bins\n'%Nbinx)
        ff.write('# Grid Y-Axis (first column): latitude [deg.] from -90 to 90 with %i equispaced bins\n'%Nbiny)
        ff.write('# Signature: Federico Lavorenti\n')
        ff.write('#\n')
        ff.write('nan\t')

        # define separator 
        if 'xLiz' in save_name:
            sep = ','
        else:
            sep = '\t'

        # print first row 
        for ix,llot in enumerate(lot):
            ff.write('%.3f'%llot+sep)
        ff.write('\n')

        # print histogram
        for iy,llat in enumerate(lat):
            ff.write('%.3f'%llat+sep)
            for ix,llot in enumerate(lot):
                if 'Flux' in ff_name:
                    if 'xLizzzzz' in ff_name:
                        ex = 2.*(np.pi*2.4e8)**2/Nbinx/Nbiny*abs(np.sin(llat*3.14/180.))
                    else:
                        ex = 1.
                    ff.write('%.3e'%(Flux[ix,iy]*ex)+sep)
                elif 'Ener' in ff_name:
                    ff.write('%.3e'%(Ener[ix,iy]*ex)+sep)
                elif 'Temp' in ff_name:
                    ff.write('%.3e'%(Temp[ix,iy]*ex)+sep)
                elif 'Pram' in ff_name:
                    ff.write('%.3e'%(Pram[ix,iy]*ex)+sep)
                elif 'ERad' in ff_name:
                    ff.write('%.3e'%(ERad[ix,iy]*ex)+sep)
            ff.write('\n')

        # close file
        ff.close()

    # close histo2d
    plt.close()

    ### PRINT INTEGRAL TOTAL ###
    Integral_Flux = integral(Flux,binx_cen,biny_cen)
    print('Flux integrated species '+sp+' '+filter_region+' = %.3e 1/s'%Integral_Flux)
    ### PRINT MEAN ENERGY ###
    Mean_Ener = np.nanmean(Ener)
    Mean_Ener1 = np.nanmean(Ekin[index_q])
    Mean_Temp = np.nanmean(Temp)
    Mean_Pram = np.nanmean(Pram)
    Mean_ERad = np.nanmean(ERad)
    print('Mean Energy species '+sp+' '+filter_region+' = %.3e eV'%Mean_Ener)
    print('Mean Energy1 species '+sp+' '+filter_region+' = %.3e eV'%Mean_Ener1)
    print('Mean Temperature species '+sp+' '+filter_region+' = %.3e eV'%Mean_Temp)
    print('Mean RAM Pressure species '+sp+' '+filter_region+' = %.3e eV'%Mean_Pram)
    print('Mean Radial Energy species '+sp+' '+filter_region+' = %.3e eV'%Mean_ERad)

    ### PLOT D.F. ENERGY WHOLE REGION ###
    x_axis = np.linspace(0,6e3,200)
    plt.title('Norm E.D.F. '+filter_region+' '+sp+' Q=%.2e 1/s'%Integral_Flux,fontsize=15)
    histo = plt.hist(Ekin[index_q],bins=x_axis,range=[np.min(x_axis),np.max(x_axis)], weights=(np.ones(len(Ekin[index_q]))/len(Ekin[index_q])) )
    histo_df = histo[0]
    id_max = np.argmax(histo_df)
    x_max = x_axis[id_max]
    y_max = histo_df[id_max]
    y_fwhm = y_max/2.
    #id_fwhm1 = np.argmin((histo_df[:id_max]-y_fwhm)**2)
    #id_fwhm2 = np.argmin((histo_df[id_max:]-y_fwhm)**2)
    #x_fwhm1 = x_axis[:id_max][id_fwhm1]
    #x_fwhm2 = x_axis[id_max:][id_fwhm2]
    # plot Planck distribution
    if sp =='e':
        x_array = histo[1]
        histo_df[np.where(histo_df<1e-8)]=1e-8
        def Planck(x,T,n):
            return n*(x/T)**3/(np.exp(x/T)-1.)*15./(np.pi**4)
        def Maxw(x,T,n):
            return n*np.sqrt(x/T)*np.exp(-x/T)*2./np.sqrt(np.pi)
        def doubleMaxw(x,T1,n1,T2,n2):
            return (n1*np.sqrt(x/T1)*np.exp(-x/T1)+n2*np.sqrt(x/T2)*np.exp(-x/T2))*2./np.sqrt(np.pi)
        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(doubleMaxw, x_array[1:], histo_df, p0=[400.,1e-3,70.,1e-2])
        chi2 = np.sum( ((histo_df-doubleMaxw(x_array[1:],popt[0],popt[1],popt[2],popt[3]))**2)/histo_df )
        plt.plot(x_array,doubleMaxw(x_array,popt[0],popt[1],popt[2],popt[3]),label='biMaxw T=%.0f,%.0f nI/nII=%.2e\n chi2=%.2e'%(popt[0],popt[2],popt[1]/popt[3],chi2))
        popt, pcov = curve_fit(Maxw, x_array[1:], histo_df, p0=[60.,1e-2])
        chi2 = np.sum( ((histo_df-Maxw(x_array[1:],popt[0],popt[1]))**2)/histo_df )
        plt.plot(x_array,Maxw(x_array,popt[0],popt[1]),label='Maxw T=%.0f\n chi2=%.2e'%(popt[0],chi2))
        plt.legend(loc=0)
    # done Planck
    plt.vlines(x_max,1e-5,y_max,color='red',linestyle='--')
    #plt.hlines(y_fwhm,x_fwhm1,x_fwhm2,color='orange',linestyle='--')
    #plt.vlines([x_fwhm1,x_fwhm2],1e-5,y_fwhm,color='orange',linestyle='--')
    plt.vlines(Mean_Ener,1e-5,y_max*3.,color='green')
    plt.vlines(Mean_Ener1,1e-5,y_max*3.,color='green',linestyle='--')
    plt.yscale('Log')
    plt.ylim(1e-3,3.*y_max)
    plt.xlabel('Energy [eV]',fontsize=16)
    plt.ylabel('#pcls/bin/Q',fontsize=16)
    plt.savefig('./images/'+save_name+'_NDF'+sp+'_'+filter_region+'_biMaxw1.png',format='png')
    plt.show()
    plt.close()
