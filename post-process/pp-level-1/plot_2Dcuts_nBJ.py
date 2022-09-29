import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fibo as fb
import sys
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge
from scipy.ndimage.filters import gaussian_filter as gf

# define path to simu and cycle to plot
data_address = sys.argv[1]
output_path = sys.argv[2]

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
# limits colorbar
lim_min = []
lim_max = []
 
# loop on cycles
for cycle in from_VTK.meta['segcycles']:

    # load fields already in 2D
    Bdp = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%cycle,cycle,fibo_obj=None,tar_var='B'))[:,0,:]
    Beq = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Beq_%i'%cycle,cycle,fibo_obj=None,tar_var='B'))[:,:,0]
    Jidp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[1][:,0,:]
    Jedp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[1][:,0,:]
    Jdp = Jidp + Jedp
    Jieq = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jieq_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[1][:,:,0]
    Jeeq = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jeeq_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[1][:,:,0]
    Jeq = Jieq + Jeeq
    Ndp = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoi1dp_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoi1')[:,0,:]
    Neq = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoi1eq_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoi1')[:,:,0]

    # miscellanous
    fields_to_plot = [Ndp,Bdp,Jdp,Neq,Beq,Jeq]

    # plot the figure
    fig = plt.figure(figsize=(17,10))
    for i,field in enumerate(fields_to_plot):
        ax = fig.add_subplot(2,3,i+1)
        plt.imshow(np.transpose(field), origin='lower', extent=(max(x),min(x),min(z),max(z)), aspect=1)
        plt.plot(xx_bs,zz_bs,linestyle='--',color='black')
        plt.plot(xx_mp,zz_mp,linestyle='-',color='black')
        plt.colorbar()

        if cycle==0 :
            plt.clim(np.min(field),np.max(field))
            lim_min.append(np.min(field))
            lim_max.append(np.max(field))
        else:
            plt.clim(lim_min[i],lim_max[i])

    #reduce spacing subplots
    plt.subplots_adjust(wspace=0.1, hspace=0.01)
    plt.savefig(output_path+'/plot_2Dcuts_nBJ_%i.png'%cycle,format='png')
    plt.close()
