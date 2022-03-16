import fibo as fb
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge
from scipy.ndimage.filters import gaussian_filter as gf

# define path to simu and cycle to plot
data_address = sys.argv[1]
cycle_min = 500
cycle_max = 10000
xlim=[3.1,-5.6]
ylim=[3.1,-3.1]
zlim=[-3.1,3.1]

# conversion factors from code units to SI values
conv_T = 21.5/(0.0031**2) #eV
conv_B = 20./0.00562 #nT
conv_J = -(4.*np.pi*30.)*(400./0.027)*(0.16) #nA/m^2

# miscellaneous
from_VTK = fb.from_VTK(data_address)
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
time = []
i=1
def scalr(a):
    return np.sqrt(a[0]**2+a[1]**2+a[2]**2)
fig = plt.figure(figsize=(16,8))
plt.axis('off')
plt.title('Electron temperature RunS',fontsize=25)

# cycle over time
for seg in from_VTK.meta['segcycles'][::2]:
  if cycle_min<=seg<=cycle_max :
    print('LOAD-->%i'%seg)
    Tpar = conv_T*from_VTK.get_scal(from_VTK.meta['name']+'_Tpare0eq_%i'%seg,seg,fibo_obj=None,tar_var='Tpar')[ix[0]:ix[1],iy[0]:iy[1],0]
    Tper = conv_T*from_VTK.get_scal(from_VTK.meta['name']+'_Tpere0eq_%i'%seg,seg,fibo_obj=None,tar_var='Tper')[ix[0]:ix[1],iy[0]:iy[1],0]/2.
    Bdp = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%seg,seg,fibo_obj=None,tar_var='B'))[ix[0]:ix[1],0,iz[0]:iz[1]]
    Jidp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%seg,seg,fibo_obj=None,tar_var='Ji')[1]
    Jedp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%seg,seg,fibo_obj=None,tar_var='Je')[1]
    Jdp = gf(Jedp,1)[ix[0]:ix[1],0,iz[0]:iz[1]]

    Tpar[abs(Tpar)>1e6]=1000
    Tper[abs(Tper)>1e6]=1000

    Tiso = (Tpar+2.*Tper)/3.

    time.append(seg*from_VTK.meta['dt']*from_VTK.meta['Vx0']/from_VTK.meta['R'])

    # plot the figure
    if i%7==0 : 
        plt.subplots_adjust(wspace=0.1, hspace=0.01)
        plt.savefig('./images/plot3_paper_MagRec_Tiso_4_%i.png'%seg,format='png')
        plt.close()
        cb1.remove
        i=1
        fig = plt.figure(figsize=(16,8))
        plt.axis('off')
        plt.title('Electron temperature RunS',fontsize=25)
       
    ax = fig.add_subplot(2,3,i)
    plt.text(0.95*max(x),0.75*max(z),r'$t$=%.1fR/$V_x$'%time[-1],fontsize=16)
    im1 = plt.imshow(np.transpose(gf(Tiso,1)), cmap='jet', origin='lower', norm=colors.LogNorm(vmin=1, vmax=3000), extent=(max(x),min(x),min(z),max(z)), aspect=1, alpha=1)
    X,Y=np.meshgrid(x,y)
    cont = plt.contour(X,Y[::-1], np.transpose(Jdp), [300])

    w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
    w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
    ax.add_artist(w1)
    ax.add_artist(w2)
    if i<4:
        plt.xticks([])
    if i==1 or i==4:
        plt.yticks(fontsize=14)
        ax.set_ylabel(r'$z_{_{MSO}}$[R]',fontsize=16)
    else:        
        plt.yticks([])
    if i>3:
        plt.xticks(fontsize=14)
        ax.set_xlabel(r'$x_{_{MSO}}$[R]',fontsize=16)
    if i==6 :
        cax1 = fig.add_axes([0.93, 0.3, 0.01, 0.4])
        cb1 = plt.colorbar(im1, cax=cax1, orientation='vertical')
        cb1.ax.tick_params(labelsize=14)
    plt.clim(1,3000)

    i+=1
