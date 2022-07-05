import fibo as fb
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.patches import Wedge

# define path to simu and cycle to plot
cycle_min = 6500
cycle_max = 6600
d_cycle = 100
xlim=[3.02,0.98]

# conversion factors from code units to SI values
conv_B = 20./0.00562 #nT
conv_J = -(4.*np.pi*30.)*(400./0.027)*(0.16) #nA/m^2
conv_N = 30. #cm-3
conv_Pe  = 100.*(4.*np.pi)/(3.8e-4)  #nPa
conv_Pi  = 1.*(4.*np.pi)/(3.8e-4)  #nPa

# do figure for plot
fig = plt.figure(figsize=(15,8))
plt.axis('off')
 
##########
## MAIN ##
##########
def main(i,data_address):

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

    # loop on cycles
    for cycle in range(cycle_min,cycle_max,d_cycle):

        # load fields already in 2D
        print('LOAD')
        Bxdp,Bydp,Bzdp = conv_B*from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%cycle,cycle,fibo_obj=None,tar_var='B')
        Bdp = conv_B*scalr(from_VTK.get_vect(from_VTK.meta['name']+'_Bdp_%i'%cycle,cycle,fibo_obj=None,tar_var='B'))
        Jedp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[1]
        Jidp = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[1]
        Jedpx = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jedp_%i'%cycle,cycle,fibo_obj=None,tar_var='Je')[0]
        Jidpx = conv_J*from_VTK.get_vect(from_VTK.meta['name']+'_Jidp_%i'%cycle,cycle,fibo_obj=None,tar_var='Ji')[0]
        Nidp = conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoi1dp_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoi1')
        Nedp = -conv_N*from_VTK.get_scal(from_VTK.meta['name']+'_rhoe0dp_%i'%cycle,cycle,fibo_obj=None,tar_var='rhoe0')
        PeXXdp = -conv_Pe*from_VTK.get_scal(from_VTK.meta['name']+'_PXXe0dp_%i'%cycle,cycle,fibo_obj=None,tar_var='PXXe0')
        PiXXdp = -conv_Pi*from_VTK.get_scal(from_VTK.meta['name']+'_PXXi1dp_%i'%cycle,cycle,fibo_obj=None,tar_var='PXXi1')

        # compute pressures [nPa]
        mu0 = 4.*np.pi*1e-7
        mi = 1.7e-27
        me = 1.7e-28
        q = 1.6e-19
        Nidp[Nidp==0.]=1e-3
        Nedp[Nedp==0.]=1e-3
        Pmagdp = 1e-9*0.5*Bdp**2/mu0
        Pramidp = 1e-15*mi*Jidpx**2/Nidp/(q**2) 
        Pramedp = 1e-15*me*Jedpx**2/Nedp/(q**2) 
        Pkinedp = PeXXdp-Pramedp
        Pkinidp = PiXXdp-Pramidp

        # bow shock
        print('FIND BS')
        if i==1:
            i_BS = np.argmax((Jidp+Jedp)[ix[0]:ix[1],0,izc])
        if i==2:
            i_BS = np.argmin((Jidp+Jedp)[ix[0]:ix[1],0,izc])
        # magnetopause
        print('FIND MP')
        i_MP = np.argmin((Pmagdp[ix[0]:ix[1],0,izc]-8.16)**2)

        # plot
        # Density
        imag = i
        ax = fig.add_subplot(4,2,imag)
        plt.plot(x,Nidp[ix[0]:ix[1],0,izc],label='$n_i$[cm$^{-3}$]',color='orange')
        plt.plot(x,Nedp[ix[0]:ix[1],0,izc],label='$n_e$[cm$^{-3}$]',color='blue')
        if imag==1 :
            plt.title('Dayside cut RunN',fontsize=25)
        if imag==2 :
            plt.title('Dayside cut RunS',fontsize=25)
        plt.legend(loc=2,fontsize=14)
        plt.xlim(max(x),min(x))
        plt.ylim(0,150)
        plt.xticks([])
        plt.yticks(fontsize=14)
        if imag%2==0 :
            plt.yticks([])
        plt.vlines([x[i_BS],x[i_MP]],0,150,color='red',linestyle='-',alpha=1.)
        plt.vlines([1.9,1.45],0,150,color='red',linestyle='--',alpha=1.)
        #plt.fill_between([2.03,2.21],y1=1e10,y2=-1e10,alpha=0.1,color='red')
        #plt.fill_between([1.33,1.47],y1=1e10,y2=-1e10,alpha=0.1,color='red')
        ## Magnetic field
        imag += 2
        ax = fig.add_subplot(4,2,imag)
        plt.plot(x,Bdp[ix[0]:ix[1],0,izc],label='$|B|$[nT]',color='black')
        plt.plot(x,Bxdp[ix[0]:ix[1],0,izc],label='$B_x$',color='red')
        plt.plot(x,Bydp[ix[0]:ix[1],0,izc],label='$B_y$',color='green')
        plt.plot(x,Bzdp[ix[0]:ix[1],0,izc],label='$B_z$',color='blue')
        plt.hlines(0,0,10,color='grey',linestyle='--')
        plt.legend(loc=2,fontsize=14)
        plt.xlim(max(x),min(x))
        plt.ylim(-200,200)
        plt.xticks([])
        plt.yticks(fontsize=14)
        if imag%2==0 :
            plt.yticks([])
        plt.vlines([x[i_BS],x[i_MP]],-200,200,color='red',linestyle='-',alpha=1.)
        plt.vlines([1.9,1.45],-200,200,color='red',linestyle='--',alpha=1.)
        ## Currents
        imag += 2
        ax = fig.add_subplot(4,2,imag)
        plt.plot(x,(Jidp+Jedp)[ix[0]:ix[1],0,izc],label='$J_y$[nA/m$^2$]',color='black')
        plt.plot(x,Jidp[ix[0]:ix[1],0,izc],label='$J_{y,i}$',color='orange')
        plt.plot(x,Jedp[ix[0]:ix[1],0,izc],label='$J_{y,e}$',color='blue')
        plt.hlines(0,0,10,color='grey',linestyle='--')
        plt.legend(loc=2,fontsize=14)
        plt.xlim(max(x),min(x))
        plt.ylim(-2500,5000)
        plt.xticks([])
        plt.yticks(fontsize=14)
        if imag%2==0 :
            plt.yticks([])
        plt.vlines([1.9,1.45],-2500,5000,color='red',linestyle='--',alpha=1.)
        plt.vlines([x[i_BS],x[i_MP]],-2500,5000,color='red',linestyle='-',alpha=1.)
        ## Pressures
        imag += 2
        ax = fig.add_subplot(4,2,imag)
        plt.plot(x,(Pmagdp+Pramidp+Pramedp)[ix[0]:ix[1],0,izc],label='$P_{mag}+P_{ram}$[nPa]',color='black')
        plt.plot(x,Pmagdp[ix[0]:ix[1],0,izc],label='$P_{mag}$',color='grey')
        plt.plot(x,(Pramidp+Pramedp)[ix[0]:ix[1],0,izc],label='$P_{ram}$',color='orange')
        #plt.plot(x,(Pkinidp+Pkinedp)[ix[0]:ix[1],0,izc],label='$P_{kin}$',color='red')
        plt.hlines(0,0,10,color='grey',linestyle='--')
        plt.legend(loc=2,fontsize=14)
        plt.xlabel(r'$x_{_{MSO}}$[R]',fontsize=16)
        plt.xlim(max(x),min(x))
        plt.ylim(0,20)
        plt.vlines([1.9,1.45],0,20,color='red',linestyle='--',alpha=1.)
        plt.vlines([x[i_BS],x[i_MP]],0,20,color='red',linestyle='-',alpha=1.)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        if imag%2==0 :
            plt.yticks([])


# call mains
main(1,'/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/Normal-newbox/data1')
main(2,'/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/minusBz-newbox/data1')
        
#reduce spacing subplots
plt.subplots_adjust(wspace=0.1, hspace=0.1)
# save fig
plt.savefig('./images/plot2_paper_nBJ_new1.png',format='png')
plt.close()
