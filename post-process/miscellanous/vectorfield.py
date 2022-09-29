import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.style.use('classic')

import matplotlib.cm as cm
from matplotlib import colors

import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf


#2D plot using the imshow function change the xlim, ylim parameter if needed
def vectorfield(a,x,y,xlabel,ylabel,title,name,filter_sigma=0,clim_min=None,clim_max=None,Log=False):

    figa = plt.figure()
    a = np.transpose(a)
    ax  = figa.add_subplot(111)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_title(title, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_xlim(min(x),max(x))
    ax.set_ylim(min(y),max(y))

    if ((clim_min==None or clim_max==None) and Log==True):
        print('ERROR: pls pass arguments clim_min,clim_max not None\n\n')
	return 

    if (filter_sigma>0):
        a=gf(a,filter_sigma)

    if(Log):
        plt.imshow(a, cmap=cm.seismic, origin='lower', extent=(min(x),max(x),min(y),max(y)), aspect='auto', norm=colors.LogNorm(vmin=clim_min/2., vmax=clim_max*2.))
    else:
        plt.imshow(a, cmap=cm.seismic, origin='lower', extent=(min(x),max(x),min(y),max(y)), aspect='auto')

    cb = plt.colorbar()
    plt.clim(clim_min,clim_max)

    plt.savefig(name+'.png', format='png', dpi=200) 
    cb.remove()
    plt.close() 
