import fibo as fb
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge

# define path to simu and cycle to plot
data_address = sys.argv[1]
cycle = 6000
xlim=[3.1,0.8]

# conversion factors from code units to SI values
conv_B = 20./0.00562 #nT
conv_J = -(4.*np.pi*30.)*(400./0.027)*(0.16) #nA/m^2
conv_N = 30. #cm-3

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
iyc = len(y)/2
izc = len(z)/2
x = x[ix[0]:ix[1]]
def scalr(a):
    return np.sqrt(a[0]**2+a[1]**2+a[2]**2)

# load fields already in 2D
Bxdp,Bydp,Bzdp = conv_B*from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%cycle,cycle,fibo_obj=None,tar_var='B')
Bdp = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%cycle,cycle,fibo_obj=None,tar_var='B'))
Jedp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[1]
Jidp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[1]
Nidp = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoi1dp_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoi1')
Nedp = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoe0dp_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoe0')

# plot the figure
fig = plt.figure(figsize=(16,8))
plt.axis('off')
plt.title('Dayside cut RunN',fontsize=25)
# Density
ax = fig.add_subplot(311)
plt.plot(x,Nidp[ix[0]:ix[1],0,izc],label='$n_i$[cm-3]',color='orange')
plt.plot(x,-Nedp[ix[0]:ix[1],0,izc],label='$n_e$[cm-3]',color='blue')
plt.legend(loc=2,fontsize=16)
plt.xlim(max(x),min(x))
plt.ylim(0,150)
plt.xticks([])
plt.fill_between([2.03,2.21],y1=1e10,y2=-1e10,alpha=0.1,color='red')
plt.fill_between([1.33,1.47],y1=1e10,y2=-1e10,alpha=0.1,color='red')
## Magnetic field
ax = fig.add_subplot(312)
plt.plot(x,Bdp[ix[0]:ix[1],0,izc],label='$|B|$[nT]',color='black')
plt.plot(x,Bxdp[ix[0]:ix[1],0,izc],label='$B_x$',color='red')
plt.plot(x,Bydp[ix[0]:ix[1],0,izc],label='$B_y$',color='green')
plt.plot(x,Bzdp[ix[0]:ix[1],0,izc],label='$B_z$',color='blue')
plt.hlines(0,0,10,color='grey',linestyle='--')
plt.legend(loc=2,fontsize=16)
plt.xlim(max(x),min(x))
plt.ylim(-200,200)
plt.xticks([])
plt.fill_between([2.03,2.21],y1=1e10,y2=-1e10,alpha=0.1,color='red')
plt.fill_between([1.33,1.47],y1=1e10,y2=-1e10,alpha=0.1,color='red')
## Currents
ax = fig.add_subplot(313)
plt.plot(x,(Jidp+Jedp)[ix[0]:ix[1],0,izc],label='$J_y$[nA/$m^2$]',color='black')
plt.plot(x,Jidp[ix[0]:ix[1],0,izc],label='$J_{y,i}$',color='orange')
plt.plot(x,Jedp[ix[0]:ix[1],0,izc],label='$J_{y,e}$',color='blue')
plt.hlines(0,0,10,color='grey',linestyle='--')
plt.legend(loc=2,fontsize=16)
plt.xlabel(r'$x_{_{MSO}}$[R]',fontsize=16)
plt.xlim(max(x),min(x))
plt.ylim(-2500,6500)
plt.fill_between([2.03,2.21],y1=1e10,y2=-1e10,alpha=0.1,color='red')
plt.fill_between([1.33,1.47],y1=1e10,y2=-1e10,alpha=0.1,color='red')
'''
im1 = plt.imshow(np.transpose(Ndp[ix[0]:ix[1],0,iz[0]:iz[1]]), cmap='jet', origin='lower', norm=colors.LogNorm(vmin=1., vmax=120.), extent=(max(x),min(x),min(z),max(z)), aspect=1)
plt.plot(xx_bs,zz_bs,linestyle='--',color='black')
plt.plot(xx_mp,zz_mp,linestyle='-',color='black')
w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
ax.add_artist(w1)
ax.add_artist(w2)
plt.clim(1.,120.)
plt.xticks([])
plt.yticks(fontsize=14)
plt.text(0.8*max(x),0.75*max(z),r'$n_i$[cm-3]',fontsize=14)
plt.ylabel(r'$z_{_{MSO}}$[R]',fontsize=14)
ax = fig.add_subplot(234)
plt.imshow(np.transpose(Neq[ix[0]:ix[1],iy[0]:iy[1],0]), cmap='jet', origin='lower', norm=colors.LogNorm(vmin=1., vmax=120.), extent=(max(x),min(x),max(y),min(y)), aspect=1)
plt.plot(xx_bs,yy_bs,linestyle='--',color='black')
plt.plot(xx_mp,yy_mp,linestyle='-',color='black')
w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
ax.add_artist(w1)
ax.add_artist(w2)
plt.clim(1.,120.)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.text(0.8*max(x),0.75*min(y),r'$n_i$[cm-3]',fontsize=14)
plt.ylabel(r'$y_{_{MSO}}$[R]',fontsize=14)
#magnetic field
ax = fig.add_subplot(232)
im2 = plt.imshow(np.transpose(Bdp[ix[0]:ix[1],0,iz[0]:iz[1]]), cmap='seismic', origin='lower', norm=colors.LogNorm(vmin=1e-1, vmax=1e3), extent=(max(x),min(x),min(z),max(z)), aspect=1)
plt.plot(xx_bs,zz_bs,linestyle='--',color='black')
plt.plot(xx_mp,zz_mp,linestyle='-',color='black')
w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
ax.add_artist(w1)
ax.add_artist(w2)
plt.clim(1e-1,1e3)
plt.xticks([])
plt.yticks([])
plt.text(0.8*max(x),0.75*max(z),r'$|B|$[nT]',fontsize=14)
ax = fig.add_subplot(235)
plt.imshow(np.transpose(Beq[ix[0]:ix[1],iy[0]:iy[1],0]), cmap='seismic', origin='lower', norm=colors.LogNorm(vmin=1e-1, vmax=1e3), extent=(max(x),min(x),max(y),min(y)), aspect=1)
plt.plot(xx_bs,yy_bs,linestyle='--',color='black')
plt.plot(xx_mp,yy_mp,linestyle='-',color='black')
w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
ax.add_artist(w1)
ax.add_artist(w2)
plt.clim(1e-1,1e3)
plt.xticks(fontsize=14)
plt.yticks([])
plt.text(0.8*max(x),0.75*min(y),r'$|B|$[nT]',fontsize=14)
#current
ax = fig.add_subplot(233)
im3 = plt.imshow(np.transpose(Jdp[ix[0]:ix[1],0,iz[0]:iz[1]]), cmap='seismic', origin='lower', extent=(max(x),min(x),min(z),max(z)), aspect=1)
plt.plot(xx_bs,zz_bs,linestyle='--',color='black')
plt.plot(xx_mp,zz_mp,linestyle='-',color='black')
w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
ax.add_artist(w1)
ax.add_artist(w2)
plt.clim(-2e3,2e3)
plt.xticks([])
plt.yticks([])
plt.text(0.8*max(x),0.75*max(z),r'$J_y$[nA/$m^2$]',fontsize=14)
ax = fig.add_subplot(236)
plt.imshow(np.transpose(Jeq[ix[0]:ix[1],iy[0]:iy[1],0]), cmap='seismic', origin='lower', extent=(max(x),min(x),max(y),min(y)), aspect=1)
plt.plot(xx_bs,yy_bs,linestyle='--',color='black')
plt.plot(xx_mp,yy_mp,linestyle='-',color='black')
w1 = Wedge((0.,0.), 1., 90., 270., fc='black')
w2 = Wedge((0.,0.), 1., 270., 90., fc='grey')
ax.add_artist(w1)
ax.add_artist(w2)
plt.clim(-2e3,2e3)
plt.xticks(fontsize=14)
plt.yticks([])
plt.text(0.8*max(x),0.75*min(y),r'$J_y$[nA/$m^2$]',fontsize=14)

#colorbar in right places
cax1 = fig.add_axes([0.128, 0.05, 0.237, 0.015])
cbar = plt.colorbar(im1, cax=cax1, orientation='horizontal')
cbar.ax.tick_params(labelsize=12)
cax2 = fig.add_axes([0.395, 0.05, 0.237, 0.015])
cbar = plt.colorbar(im2, cax=cax2, orientation='horizontal')
cbar.ax.tick_params(labelsize=12)
cax3 = fig.add_axes([0.668, 0.05, 0.237, 0.015])
cbar = plt.colorbar(im3, cax=cax3, orientation='horizontal')
cbar.ax.tick_params(labelsize=12)
'''
#reduce spacing subplots
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('./images/plot2_paper_nBJ_2N.png',format='png')

