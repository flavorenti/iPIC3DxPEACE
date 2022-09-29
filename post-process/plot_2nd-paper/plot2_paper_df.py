import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 

#################################################
#file_name = 'PR1-minusBz-e-_90-95_-135--130.txt'
#theta = 92.5*np.pi/180.
#phi = -132.5*np.pi/180.
################################################
#file_name = 'PR1-minusBz-e-_165-170_-10--5.txt'
#theta = 167.5*np.pi/180.
#phi = -7.5*np.pi/180.
################################################
#file_name = 'PR1-plusBz-e-_165-170_-10--5.txt'
#theta = 167.5*np.pi/180.
#phi = -7.5*np.pi/180.
################################################
vlim=15.

# versore radiale
rx = np.sin(theta)*np.cos(phi)
ry = np.sin(theta)*np.sin(-phi)
rz = np.cos(theta)

# load pcls file
c,x,y,z,u,v,w,q = np.loadtxt('./pcls/'+file_name,unpack=True)

# change units in vthe
u /= 0.031
v /= 0.031
w /= 0.031

plt.hist2d(u,v,norm=mpl.colors.LogNorm(),range=[[-vlim,vlim],[-vlim,vlim]],bins=[20,20])
plt.plot([-vlim*rx,0,vlim*rx],[-vlim*ry,0,vlim*ry],color='red',alpha=0.5)
plt.vlines(0,-vlim,vlim,color='grey',linestyle='--')
plt.hlines(0,-vlim,vlim,color='grey',linestyle='--')
plt.xlabel(r'$v_x/v_{the}$',fontsize=20)
plt.ylabel(r'$v_y/v_{the}$',fontsize=20)
plt.title(file_name,fontsize=16)
plt.colorbar()
plt.savefig('./images/pcls-VxVy_'+file_name.replace('txt','png'),format='png')
plt.show()

plt.hist2d(u,w,norm=mpl.colors.LogNorm(),range=[[-vlim,vlim],[-vlim,vlim]],bins=[20,20])
plt.plot([-vlim*rx,0,vlim*rx],[-vlim*rz,0,vlim*rz],color='red',alpha=0.5)
plt.vlines(0,-vlim,vlim,color='grey',linestyle='--')
plt.hlines(0,-vlim,vlim,color='grey',linestyle='--')
plt.xlabel(r'$v_x/v_{the}$',fontsize=20)
plt.ylabel(r'$v_z/v_{the}$',fontsize=20)
plt.title(file_name,fontsize=16)
plt.colorbar()
plt.savefig('./images/pcls-VxVz_'+file_name.replace('txt','png'),format='png')
plt.show()

plt.hist2d(v,w,norm=mpl.colors.LogNorm(),range=[[-vlim,vlim],[-vlim,vlim]],bins=[20,20])
plt.plot([-vlim*ry,0,vlim*ry],[-vlim*rz,0,vlim*rz],color='red',alpha=0.5)
plt.vlines(0,-vlim,vlim,color='grey',linestyle='--')
plt.hlines(0,-vlim,vlim,color='grey',linestyle='--')
plt.xlabel(r'$v_y/v_{the}$',fontsize=20)
plt.ylabel(r'$v_z/v_{the}$',fontsize=20)
plt.title(file_name,fontsize=16)
plt.colorbar()
plt.savefig('./images/pcls-VyVz_'+file_name.replace('txt','png'),format='png')
plt.show()
