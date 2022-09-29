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
# profiles MP, BS
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

# loop on species
for sp in from_VTK.meta['species']:
    # limits colorabar
    lim_min = []
    lim_max = []
    # loop on cycles
    for cycle in from_VTK.meta['segcycles']:
        # load fields already in 2D
        TXXdp = from_VTK.get_scal(from_VTK.meta['name']+'_TXX'+sp+'dp_%i'%cycle,cycle,fibo_obj=None,tar_var='TXX')[:,0,:]
        TXXeq = from_VTK.get_scal(from_VTK.meta['name']+'_TXX'+sp+'eq_%i'%cycle,cycle,fibo_obj=None,tar_var='TXX')[:,:,0]
        TYYdp = from_VTK.get_scal(from_VTK.meta['name']+'_TYY'+sp+'dp_%i'%cycle,cycle,fibo_obj=None,tar_var='TYY')[:,0,:]
        TYYeq = from_VTK.get_scal(from_VTK.meta['name']+'_TYY'+sp+'eq_%i'%cycle,cycle,fibo_obj=None,tar_var='TYY')[:,:,0]
        TZZdp = from_VTK.get_scal(from_VTK.meta['name']+'_TZZ'+sp+'dp_%i'%cycle,cycle,fibo_obj=None,tar_var='TZZ')[:,0,:]
        TZZeq = from_VTK.get_scal(from_VTK.meta['name']+'_TZZ'+sp+'eq_%i'%cycle,cycle,fibo_obj=None,tar_var='TZZ')[:,:,0]

        # miscellanous
        fields_to_plot = [TXXdp,TYYdp,TZZdp,TXXeq,TYYeq,TZZeq]

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
        plt.savefig(output_path+'/plot_2Dcuts_T'+sp+'_%i.png'%cycle,format='png')
        plt.close()
