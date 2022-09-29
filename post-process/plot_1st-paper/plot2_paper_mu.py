import fibo as fb
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge

# define path to simu and cycle to plot
data_address = sys.argv[1]
cycle_min = 0
cycle_max = 10000

# conversion factors from code units to SI values
conv_mu = np.sqrt(0.00562)/0.0031 #mu_SW

# miscellaneous
from_VTK = fb.from_VTK(data_address)
from_VTK.get_meta(silent=False)
time = []

# cycle over time
for seg in from_VTK.meta['segcycles']:
  if cycle_min<=seg<=cycle_max :
    print('LOAD-->%i'%seg)
    mu_e = conv_mu*from_VTK.get_scal(from_VTK.meta['name']+'_mue0_%i'%seg,seg,fibo_obj=None,tar_var='mue0'))
    mu_i = conv_mu*from_VTK.get_scal(from_VTK.meta['name']+'_mui1_%i'%seg,seg,fibo_obj=None,tar_var='mui1'))
    time.append(seg*from_VTK.meta['dt'])

# plot the figure
fig = plt.figure(figsize=(10,4))
ax = fig.add_subplot(111)
ax.plot(t,mu_e,color='blue',marker='o',linestyle='-')
ax.plot(t,mu_i,color='red',marker='o',linestyle='-')
plt.xlabel('$t\omega_{pi}$',fontsize=18)
plt.text(0.8*max(t),0.75*max(mu_i),r'$\mu_i$',fontsize=18)

plt.savefig('./images/plot2_paper_mu_1.png',format='png')

