import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.ndimage.filters import gaussian_filter as gf
import sys


# either Mariner or Bepi
mission = sys.argv[1]

# path to file with pcls
if mission=='Mariner':
    runN_path_file = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/Normal-newbox/texts/output_good-pcls_2_5999_Mariner.txt'
    runS_path_file = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/minusBz-newbox/texts/output_good-pcls_2_5999_Mariner.txt'
if mission=='Bepi':
    runN_path_file = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/Normal-newbox/texts/output_good-pcls_2_5999_Bepi.txt'
    runS_path_file = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/minusBz-newbox/texts/output_good-pcls_2_5999_Bepi1.txt'

# Mariner10 MFB1 data
if mission=='Mariner':
    #~~ PLS ~~
    emin = 13.4
    emax = 687.
    Ne = 100
    # time bins
    Nt = 250
    tmin = -23.
    tmax = +19.
    #~~ M10-MFB1 --
    tBS1 = -19.5
    dBS1 = +1.
    tMP1 = -9.7
    dMP1 = 0.3
    tMP2 = +8.
    dMP2 = 0.6
    tBS2 = +10.5
    dBS2 = +2.
    # string plot
    label = 'Mariner10'
    # remove duplicates?
    apply_rm_dupl=False

# BepiColombo MFB1 data
if mission=='Bepi':
    #~~ MEA-SW ~~
    emin = 3.
    emax = 3000.
    Ne = 126
    #~~ MEA-MS ~~
    #emin = 3.
    #emax = 26010.
    #Ne = 126
    # time bins
    Nt = 420
    tmin = -40.
    tmax = +29.
    #~~ BC-MFB1 ~~
    tBS1 = -29
    dBS1 = +10.
    tMP1 = -10.
    dMP1 = +1.
    tMP2 = +7.
    dMP2 = +1.
    tBS2 = +11.
    dBS2 = +1.
    # string plot
    label = 'BepiColombo'
    # FOV~~MEA1~~
    theta = [15,120]
    phi   = [80,100]
    apply_fov = False
    # remove duplicates?
    apply_rm_dupl=True

# load pcls
def load_pcls(path):

    # load file
    timeX_, xp_, yp_, zp_, ep_ = np.loadtxt(path,unpack=True)

    # remove duplicates
    if apply_rm_dupl:
        print('Apply Remove Duplicates filter!')
        values, idx = np.unique([timeX_,xp_,yp_,zp_,ep_],axis=1,return_index=True)
        timeX1 = timeX_[idx]
        timeX= []
        xp = xp_[idx]
        yp = yp_[idx]
        zp = zp_[idx]
        ep1 = ep_[idx]
        ep = []
    else:
        timeX1 = timeX_
        timeX= []
        xp = xp_
        yp = yp_
        zp = zp_
        ep1 = ep_
        ep = []

    # instrument FOV (only for Bepi)
    if mission=='Bepi' and apply_fov:
        fov = []
        print('Apply FOV filter!')
        for t in np.unique(timeX1):
            it     = (timeX1==t)
            t_good = timeX1[it]
            x_good = xp[it]-np.mean(xp[it])
            y_good = yp[it]-np.mean(yp[it])
            z_good = zp[it]-np.mean(zp[it])
            e_good = ep1[it]
            dist = np.sqrt(x_good**2+y_good**2+z_good**2)
            thetap = np.arctan(z_good/dist)*180./np.pi
            phip   = np.arctan(y_good/x_good)*180./np.pi
            timeX  = np.append(timeX, t_good[(thetap>theta[0])*(thetap<theta[1])*(phip>phi[0])*(phip<phi[1])])
            ep     = np.append(ep,    e_good[(thetap>theta[0])*(thetap<theta[1])*(phip>phi[0])*(phip<phi[1])])            
    else:
        timeX = timeX1
        ep = ep1

    # histogram 2d
    hist, xbins, ybins = np.histogram2d(timeX, ep, bins=[Nt,Ne], range=[[tmin,tmax],[np.log10(emin), np.log10(emax)]])

    # remove zero values
    hist[np.where(hist==0)]=1.

    # compute density
    dens = np.mean(hist,axis=1)

    # remove holes
    for ix in range(1,Nt-1):
        imax = np.argmax([dens[ix-1],dens[ix],dens[ix+1]])
        hist[ix,:] = hist[ix-1+imax,:]
    hist = hist[::2,:]

    # re-compute density
    dens = np.mean(hist,axis=1)
    dens = dens/dens[-1]*30. 

    return xbins, ybins, hist, dens


# call two runs loader
xbinsN, ybinsN, histN, densN = load_pcls(runN_path_file)
xbinsS, ybinsS, histS, densS = load_pcls(runS_path_file)


# figure with subplots
fig = plt.figure(figsize=(16,7.5))
plt.suptitle(label+' MFB#1',fontsize=27)

# plot RunN hist
ax1 = plt.subplot(221)
plt.title('Northward IMF (RunN)',fontsize=23)
im1 = ax1.imshow(np.transpose(histN), cmap='jet', interpolation='none', origin='lower', norm=colors.LogNorm(vmin=np.min(histN), vmax=np.max(histN)), extent=[xbinsN[0], xbinsN[-1], ybinsN[0], ybinsN[-1]], aspect='auto')
ax1.set_ylabel(r'log10(E/eV)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax1.vlines([0.],ybinsN[0],ybinsN[-1],linestyle='--',color='red')
ax1.vlines([tBS1+.5*dBS1,tBS2+.5*dBS2],ybinsN[0],ybinsN[-1],alpha=0.5,linestyle='--',color='black')
ax1.vlines([tMP1+.5*dMP1,tMP2+.5*dMP2],ybinsN[0],ybinsN[-1],alpha=0.5,linestyle='--',color='black')
plt.fill_between([tMP1,tMP1+dMP1],ybinsN[0],ybinsN[-1],alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tMP2,tMP2+dMP2],ybinsN[0],ybinsN[-1],alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS1,tBS1+dBS1],ybinsN[0],ybinsN[-1],alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS2,tBS2+dBS2],ybinsN[0],ybinsN[-1],alpha=0.3,color='grey',linestyle='--')

# plot RunS hist
ax1 = plt.subplot(222)
plt.title('Southward IMF (RunS)',fontsize=23)
im1 = ax1.imshow(np.transpose(histS), cmap='jet', interpolation='none', origin='lower', norm=colors.LogNorm(vmin=np.min(histS), vmax=np.max(histS)), extent=[xbinsS[0], xbinsS[-1], ybinsS[0], ybinsS[-1]], aspect='auto')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax1.vlines([0.],ybinsS[0],ybinsS[-1],linestyle='--',color='red')
ax1.vlines([tMP1+.5*dMP1,tMP2+.5*dMP2],ybinsS[0],ybinsS[-1],alpha=0.5,linestyle='--',color='black')
ax1.vlines([tBS1+.5*dBS1,tBS2+.5*dBS2],ybinsS[0],ybinsS[-1],alpha=0.5,linestyle='--',color='black')
plt.fill_between([tMP1,tMP1+dMP1],ybinsS[0],ybinsS[-1],alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tMP2,tMP2+dMP2],ybinsS[0],ybinsS[-1],alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS1,tBS1+dBS1],ybinsS[0],ybinsS[-1],alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS2,tBS2+dBS2],ybinsS[0],ybinsS[-1],alpha=0.3,color='grey',linestyle='--')

# plot RunN dens
ax2 = plt.subplot(223)
ax2.plot(np.linspace(tmin,tmax,len(densN)), densN)
ax2.set_ylabel(r'$n_e$ [$cm^{-3}$]',fontsize=16)
ax2.set_xlabel('time [CA-min]',fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(xbinsN[0], xbinsN[-1])
plt.ylim(0,np.max(densS)*1.2)
ax2.vlines([0.],min(densS),2*max(densS),linestyle='--',color='red')
ax2.vlines([tMP1+.5*dMP1,tMP2+.5*dMP2],min(densS),2*max(densS),alpha=0.5,linestyle='--',color='black')
ax2.vlines([tBS1+.5*dBS1,tBS2+.5*dBS2],min(densS),2*max(densS),alpha=0.5,linestyle='--',color='black')
plt.fill_between([tMP1,tMP1+dMP1],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tMP2,tMP2+dMP2],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS1,tBS1+dBS1],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS2,tBS2+dBS2],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.text(tBS1+.5*dBS1-3.,np.max(densS)*1.1,'BS',color='black',fontsize=11,weight='bold')
plt.text(tMP1-3,np.max(densS)*1.1,'MP',color='black',fontsize=11,weight='bold')
plt.text(-3,np.max(densS)*1.1,'CA',color='red',fontsize=11,weight='bold')
plt.text(tMP2-3,np.max(densS)*1.1,'MP',color='black',fontsize=11,weight='bold')
plt.text(tBS2-2.2,np.max(densS)*1.1,'BS',color='black',fontsize=11,weight='bold')

# plot RunS dens
ax2 = plt.subplot(224)
ax2.plot(np.linspace(tmin,tmax,len(densS)), densS)
ax2.set_xlabel('time [CA-min]',fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(xbinsS[0], xbinsS[-1])
plt.ylim(0,np.max(densS)*1.2)
ax2.vlines([0.],min(densS),2*max(densS),linestyle='--',color='red')
ax2.vlines([tMP1+.5*dMP1,tMP2+.5*dMP2],min(densS),2*max(densS),alpha=0.5,linestyle='--',color='black')
ax2.vlines([tBS1+.5*dBS1,tBS2+.5*dBS2],min(densS),2*max(densS),alpha=0.5,linestyle='--',color='black')
plt.fill_between([tMP1,tMP1+dMP1],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tMP2,tMP2+dMP2],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS1,tBS1+dBS1],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.fill_between([tBS2,tBS2+dBS2],min(densS),2*max(densS),alpha=0.3,color='grey',linestyle='--')
plt.text(tBS1+.5*dBS1-3,np.max(densS)*1.1,'BS',color='black',fontsize=11,weight='bold')
plt.text(tMP1-3,np.max(densS)*1.1,'MP',color='black',fontsize=11,weight='bold')
plt.text(-3,np.max(densS)*1.1,'CA',color='red',fontsize=11,weight='bold')
plt.text(tMP2-3,np.max(densS)*1.1,'MP',color='black',fontsize=11,weight='bold')
plt.text(tBS2-2.2,np.max(densS)*1.1,'BS',color='black',fontsize=11,weight='bold')


# colorbar
cax = fig.add_axes([0.92, 0.53, 0.01, 0.35])
cbar = plt.colorbar(im1,cax=cax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('counts [arb. units]',fontsize=16)

# show fig
plt.savefig('./images/plot5_dfe_'+mission+'_new1.png',format='png')

