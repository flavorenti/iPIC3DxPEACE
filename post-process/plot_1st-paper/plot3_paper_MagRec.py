import fibo as fb
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge
from scipy.ndimage.filters import gaussian_filter as gf

# define path to simulations
data_address1 = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/Normal-newbox/data1/'
data_address2 = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/minusBz-newbox/data1/'

# define path to simu and cycle to plot
cycle_array = [2200,4300,6500]
xlim=[4.1,-5.1]
ylim=[3.1,-3.1]
zlim=[-3.1,3.1]

# conversion factors from code units to SI values
conv_T = 21.5/(0.0031**2) #eV
conv_B = 20./0.00562 #nT
conv_J = -(4.*np.pi*30.)*(400./0.027)*(0.16) #nA/m^2

# figure to use
fig = plt.figure(figsize=(15,8.4))
plt.text(0.3,1.03,'Electron temperature RunN',fontsize=25)
plt.text(0.3,0.47,'Electron temperature RunS',fontsize=25)
plt.axis('off')
i=1

# loop on two simus
for path in [data_address1,data_address2]:

    # load simu
    from_VTK = fb.from_VTK(path)
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
    ix = [np.argmin(abs(x-xlim[0])), np.argmin(abs(x-xlim[1]))]
    iy = [np.argmin(abs(y-ylim[0])), np.argmin(abs(y-ylim[1]))]
    iz = [np.argmin(abs(z-zlim[0])), np.argmin(abs(z-zlim[1]))]
    x = x[ix[0]:ix[1]]
    y = y[iy[0]:iy[1]]
    z = z[iz[0]:iz[1]]

    # init empty variable
    time = []

    # cycle over time
    for seg in from_VTK.meta['segcycles']:
        if seg in cycle_array :

            # load files 2D
            print('LOAD-->%i'%seg)
            Tpar = conv_T*from_VTK.get_scal(from_VTK.meta['name']+'_Tpare0dp_%i'%seg,seg,fibo_obj=None,tar_var='Tpar')[ix[0]:ix[1],0,iz[0]:iz[1]]
            Tper = conv_T*from_VTK.get_scal(from_VTK.meta['name']+'_Tpere0dp_%i'%seg,seg,fibo_obj=None,tar_var='Tper')[ix[0]:ix[1],0,iz[0]:iz[1]]/2.
            rhoedp = from_VTK.get_scal(from_VTK.meta['name']+'_rhoe0dp_%i'%seg,seg,fibo_obj=None,tar_var='rhoe')[ix[0]:ix[1],0,iz[0]:iz[1]]
            Jedp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%seg,seg,fibo_obj=None,tar_var='Je')[1][ix[0]:ix[1],0,iz[0]:iz[1]]
            Jidp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%seg,seg,fibo_obj=None,tar_var='Ji')[1][ix[0]:ix[1],0,iz[0]:iz[1]]
            Bxdp  = -conv_B*from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%seg,seg,fibo_obj=None,tar_var='B')[0][ix[0]:ix[1],0,iz[0]:iz[1]]
            Bzdp  = conv_B*from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%seg,seg,fibo_obj=None,tar_var='B')[2][ix[0]:ix[1],0,iz[0]:iz[1]]

            # remove cells with low #pcls
            Tpar[np.where(rhoedp>-1/128.)]=np.nan
            Tper[np.where(rhoedp>-1/128.)]=np.nan

            # Compute Tiso
            Tiso = (Tpar+2.*Tper)/3.
            
            # Compute total current
            Jdp = gf(Jidp+Jedp,1)

            # check values Tpar e Tperp
            print(Tpar[10:15,10])
            print(Tper[10:15,10])

            # cumpute time
            time.append(seg*from_VTK.meta['dt']*from_VTK.meta['Vx0']/from_VTK.meta['R'])

            # create plot
            ax = fig.add_subplot(2,3,i)
            plt.text(0.97*max(x),1.07*max(z),r'$t$=%.1fR/$V_x$'%time[-1],fontsize=14)
            im1 = plt.imshow(np.transpose(gf(Tiso,1)), cmap='jet', origin='lower', norm=colors.LogNorm(vmin=1, vmax=3000), extent=(max(x),min(x),min(z),max(z)), aspect=1, alpha=1, zorder=0)
            X,Y=np.meshgrid(x,z)
            #cont = plt.contour(X,Y, np.transpose(Jdp), [800], color='purple')
            plt.streamplot(x, z, np.transpose(Bxdp), np.transpose(Bzdp), zorder=1, color='grey')#, start_points=seeds, linewidth=1, arrowstyle='->', minlength=0.5, maxlength=50, color='black', zorder=1)
            plt.xlim(xlim[0],xlim[1])

            w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
            w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
            ax.add_artist(w1)
            ax.add_artist(w2)

            if i<4:
                plt.xticks(fontsize=14)
            if i==1 or i==4:
                plt.yticks(fontsize=14)
                ax.set_ylabel(r'$z_{_{MSO}}$[R]',fontsize=16)
            else:        
                plt.yticks([])
            if i>3:
                plt.xticks(fontsize=14)
                ax.set_xlabel(r'$x_{_{MSO}}$[R]',fontsize=16)
            if i==3 or i==7 :
                cb1 = plt.colorbar(im1, cax=cax1, orientation='vertical')
                cb1.ax.tick_params(labelsize=14)
                cb1.set_label('$T_e$[eV]',fontsize=16)
            plt.clim(1,3000)

            i+=1


# plot figure
plt.subplots_adjust(wspace=0.1, hspace=0.2)
plt.savefig('./images/plot3_paper_MagRec_Tiso_new8.png',format='png')
plt.show()
plt.close()
cb1.remove()
