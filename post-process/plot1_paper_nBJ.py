import fibo as fb
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge

# define path to simu and cycle to plot
data_address = sys.argv[1]
cycle = 6000
xlim=[4.1,-6.1]
ylim=[4.1,-4.1]
zlim=[-4.1,4.1]

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
iy = [np.argmin(abs(y-ylim[0])), np.argmin(abs(y-ylim[1]))]
iz = [np.argmin(abs(z-zlim[0])), np.argmin(abs(z-zlim[1]))]
print(min(x),max(x),min(y),max(y),min(z),max(z))
x = x[ix[0]:ix[1]]
y = y[iy[0]:iy[1]]
z = z[iz[0]:iz[1]]
def scalr(a):
    return np.sqrt(a[0]**2+a[1]**2+a[2]**2)
theta = np.linspace(-np.pi/1.1,np.pi/1.1,10000)
Rmp   = 1.45*(2./(1.+np.cos(theta)))**(0.5)
Rbs   = 2.75*1.04/(1+1.04*np.cos(theta))
x_mp  =  Rmp*np.cos(theta)
y_mp  =  Rmp*np.sin(theta)
z_mp  =  0.2+Rmp*np.sin(theta)
x_bs  =  0.5+Rbs*np.cos(theta)
y_bs  =  Rbs*np.sin(theta)
z_bs  =  0.2+Rbs*np.sin(theta)
inside_mp = (x_mp<max(x))*(x_mp>min(x))*(y_mp<max(y))*(y_mp>min(y))*(z_mp<max(z))*(z_mp>min(z))
inside_bs = (x_bs<max(x))*(x_bs>min(x))*(y_bs<max(y))*(y_bs>min(y))*(z_bs<max(z))*(z_bs>min(z))
xx_mp = x_mp[inside_mp]
yy_mp = y_mp[inside_mp]
zz_mp = z_mp[inside_mp]
xx_bs = x_bs[inside_bs]
yy_bs = y_bs[inside_bs]
zz_bs = z_bs[inside_bs]


# load fields already in 2D
Bdp = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%cycle,cycle,fibo_obj=None,tar_var='B'))
Beq = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Beq_%i'%cycle,cycle,fibo_obj=None,tar_var='B'))
Jdp  = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[1]
Jdp += conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[1]
Jeq  = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jieq_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[1]
Jeq += conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jeeq_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[1]
Ndp = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoi1dp_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoi1')
Neq = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoi1eq_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoi1')

# plot the figure
fig = plt.figure(figsize=(18,10))
plt.axis('off')
plt.title('Northward IMF (RunN)',fontsize=25)
#density
ax = fig.add_subplot(231)
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

#reduce spacing subplots
plt.subplots_adjust(wspace=0.1, hspace=0.01)

plt.savefig('./images/plot1_paper_nBJ_6N.png',format='png')

