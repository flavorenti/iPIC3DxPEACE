from vectorfield import *
from scipy import ndimage
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
from vtk_routines import *
import time
import matplotlib.pyplot as plt


# main of the code: PR_number  = number production run (0,1,2,0-bigbox)
#                   fact       = scaling factor (PR0,1,2 = 1,2,4) for bigbox this is 2
#                   run_number = number of the run we want to plot
def main(PR_number,fact,run_number,Ncycle):

    # define pars simulation as global variables
    global x,y,z,dt,ypl,zpl,qome

    # path to the simulation data
    data_path = '../data_simu_iPIC3D/Mercury_SaeInit/PR'+PR_number+'/run'+str(run_number)+'/data/'

    # box parameters, code units di,wci,c
    x   = np.linspace(0,10.,256*fact)
    y   = np.linspace(0,8.,208*fact)
    z   = np.linspace(0,8.,208*fact)

    # double box for the bigbox run
    if(PR_number=='0-bigbox'):
        x=x*2.
        y=y*2.
        z=z*2.

    # timestep wci-1
    dt = 0.5
    # position planet center
    xpl = 0.4*max(x)
    ypl = 0.5*max(y)
    zpl = 0.5*max(z)#-0.2
    # q/m electrons
    qome = -100.

    # start time
    start = time.time()

    # open the vtk files
    if(PR_number!='2'):
        Pxx_vtk = open_vtk_scal(data_path+'Dipole3D_PXXe0_'+str(Ncycle)+'.vtk')
        Pyy_vtk = open_vtk_scal(data_path+'Dipole3D_PYYe0_'+str(Ncycle)+'.vtk')
        Pzz_vtk = open_vtk_scal(data_path+'Dipole3D_PZZe0_'+str(Ncycle)+'.vtk')
        rhoi_vtk = open_vtk_scal(data_path+'Dipole3D_rhoi1_'+str(Ncycle)+'.vtk')
        rhoe_vtk = open_vtk_scal(data_path+'Dipole3D_rhoe0_'+str(Ncycle)+'.vtk')
        ni =  convert_vtk_scal(rhoi_vtk)[:,:,:]
        ne = -convert_vtk_scal(rhoe_vtk)[:,:,:]
        Px  =  convert_vtk_scal(Pxx_vtk)*4.*np.pi/qome
        Py  =  convert_vtk_scal(Pyy_vtk)*4.*np.pi/qome
        Pz  =  convert_vtk_scal(Pzz_vtk)*4.*np.pi/qome
    Ji_vtk  = open_vtk_vect(data_path+'Dipole3D_Ji_'+str(Ncycle)+'.vtk')
    Je_vtk  = open_vtk_vect(data_path+'Dipole3D_Je_'+str(Ncycle)+'.vtk')
    print('open B')
    B_vtk   = open_vtk_vect(data_path+'Dipole3D_B_'+str(Ncycle)+'.vtk')
    Jx = convert_vtk_vect(Je_vtk,'Je')[0][:,:,:]
    Jy = convert_vtk_vect(Je_vtk,'Je')[1][:,:,:] 
    Jz = convert_vtk_vect(Je_vtk,'Je')[2][:,:,:]
    modJ = np.sqrt(Jx**2+Jy**2+Jz**2)
    Bx = convert_vtk_vect(B_vtk,'B')[0][:,:,:]
    print('bx done')
    By = convert_vtk_vect(B_vtk,'B')[1][:,:,:]
    print('by done')
    Bz = convert_vtk_vect(B_vtk,'B')[2][:,:,:]
    print('bz done')
    modB = np.sqrt(Bx**2+By**2+Bz**2)
    print('3D fields are open')

    # create maxwell d.f. electrons using n,J,T from simu
    v_plot = np.linspace(-0.1,0.1,2000) #c
    # retrieve thermal speed
    Px    = Px + (4.*np.pi*Jx)**2/gf(ne,2)/qome #subtract the mean flow to get real P
    Py    = Py + (4.*np.pi*Jy)**2/gf(ne,2)/qome 
    Pz    = Pz + (4.*np.pi*Jz)**2/gf(ne,2)/qome 
    Tx    = abs(Px/gf(ne,2)*qome)                    # divide by density to get T
    Ty    = abs(Py/gf(ne,2)*qome)      
    Tz    = abs(Pz/gf(ne,2)*qome)      
    vthx  = np.sqrt(Tx)                         # sqrt to get thermal speed
    vthy  = np.sqrt(Ty)
    vthz  = np.sqrt(Tz)
    # retrieve drift J
    Vx = -4.*np.pi*Jx/gf(ne,2) #c 
    Vy = -4.*np.pi*Jy/gf(ne,2) #c 
    Vz = -4.*np.pi*Jz/gf(ne,2) #c 

    def maxw_v(v,n,J,vth):
        return n*np.exp(-((v-J)/vth)**2)

    # distribution function in energy E=m*v^2 T=me*vthe^2
    e_plot = np.linspace(0.2e-5,2.e-3,2000)  #me*c^2
    e_plot = e_plot*500.e3 #eV
    def maxw_E(E,n,vth):
        Temp = 21.5*(vth/0.031)**2  #eV
        #return n*np.exp(-E*abs(qome)/vth**2)*np.sqrt(E)
        return n*np.exp(-E/Temp)*np.sqrt(E/Temp)

    #testx = maxw_v(v_plot,ne[24,104,104],Vx[24,104,104],vthx[24,104,104]) 
    #testy = maxw_v(v_plot,ne[24,104,104],Vy[24,104,104],vthy[24,104,104]) 
    #testz = maxw_v(v_plot,ne[24,104,104],Vz[24,104,104],vthz[24,104,104]) 

    #plt.subplot(311)
    #plt.plot(v_plot, testx)
    #plt.subplot(312)
    #plt.plot(v_plot, testy)
    #plt.subplot(313)
    #plt.plot(v_plot, testz)
    #plt.show()

    # filter the 3D files
    #ni = gf(ni,2)
    #ne = gf(ne,2)
    #Bx = gf(Bx,2)
    #By = gf(By,2)
    #Bz = gf(Bz,2)

    # open the trajectory file
    tt,bbX,bxX,byX,bzX,xx,yy,zz,n1X,n2X,n3X,T1X,T2X,T3X = np.loadtxt('Mariner-flyby-1st_full-dataset.txt',unpack=True)
    d = np.sqrt(xx**2+yy**2+zz**2)
    i_CA = np.where(d==min(d))
    tt = (tt-tt[i_CA])/60.
    modbX = np.sqrt(bxX**2+byX**2+bzX**2)
    # pass from MSO -> box coordinate system
    xx = -xx + xpl
    yy = -yy + ypl
    zz = +zz + zpl
    bxX = -bxX
    byX = -byX
    # define trajectory arrays to fill in later
    NN = len(xx)
    bxbx = []
    byby = []
    bzbz = []
    bb   = []
    nini = []
    nene = []
    clock= []
    cone = []
    lenn = 330
    if(PR_number=='0-bigbox'):
        lenn=532
    DFxx = np.zeros(2000*lenn).reshape(2000,lenn) #330 for PR0-PR1 -- 532 for PR0-bigbox
    DFyy = np.zeros(2000*lenn).reshape(2000,lenn)
    DFzz = np.zeros(2000*lenn).reshape(2000,lenn)
    vxvx = (np.roll(xx,-1)-xx)*2440./(np.roll(tt,-1)-tt)/60. #km/s 
    vyvy = (np.roll(yy,-1)-yy)*2440./(np.roll(tt,-1)-tt)/60. #km/s 
    vzvz = (np.roll(zz,-1)-zz)*2440./(np.roll(tt,-1)-tt)/60. #km/s 
    vv   = np.sqrt(vxvx**2+vyvy**2+vzvz**2)

    # angles from marinerX data
    cloX = 180./np.pi*np.arctan(byX/bzX)
    conX = 180./np.pi*np.arccos(bxX/modbX)

    enter=0
    ip=0
    # loop on the trajectory file downloaded from amda
    for ii in range(NN):
        if((xx[ii]>0)*(yy[ii]>0)*(zz[ii]>0)*(xx[ii]<max(x))*(yy[ii]<max(y))*(zz[ii]<max(z))):
            ix = np.where(abs(x-xx[ii])==min(abs(x-xx[ii])))
            iy = np.where(abs(y-yy[ii])==min(abs(y-yy[ii])))
            iz = np.where(abs(z-zz[ii])==min(abs(z-zz[ii])))
            bxbx = np.append(bxbx,Bx[ix,iy,iz]/0.0056*20.)
            byby = np.append(byby,By[ix,iy,iz]/0.0056*20.)
            bzbz = np.append(bzbz,Bz[ix,iy,iz]/0.0056*20.)
            bb   = np.append(bb,modB[ix,iy,iz]/0.0056*20.)

            vth_para = np.sqrt(Tx[ix,iy,iz]*Bx[ix,iy,iz]**2 + Ty[ix,iy,iz]*By[ix,iy,iz]**2 + Tz[ix,iy,iz]*Bz[ix,iy,iz]**2)/modB[ix,iy,iz]
            vth_perp = np.sqrt(Tx[ix,iy,iz]*(1.-(Bx[ix,iy,iz]/modB[ix,iy,iz])**2) + Ty[ix,iy,iz]*(1.-(By[ix,iy,iz]/modB[ix,iy,iz])**2) + Tz[ix,iy,iz]*(1.-(Bz[ix,iy,iz]/modB[ix,iy,iz])**2))/np.sqrt(2.)

            #DFx = maxw_v(v_plot,ne[ix,iy,iz],Vx[ix,iy,iz],vthx[ix,iy,iz])
            #DFy = maxw_v(v_plot,ne[ix,iy,iz],Vy[ix,iy,iz],vthy[ix,iy,iz])
            #DFz = maxw_v(v_plot,ne[ix,iy,iz],Vz[ix,iy,iz],vthz[ix,iy,iz])
            DFx = maxw_E(e_plot,ne[ix,iy,iz],vth_perp)
            DFy = maxw_E(e_plot,ne[ix,iy,iz],vth_para)
            DFxx[:,ip] = DFx#/np.max(DFx)
            DFyy[:,ip] = DFy#/np.max(DFy)
            #DFzz[:,ip] = DFz/np.max(DFz)
            ip=ip+1
            print(ip)

            clock= np.append(clock,180./np.pi*np.arctan(By[ix,iy,iz]/Bz[ix,iy,iz])) 
            cone = np.append(cone ,180./np.pi*np.arccos(Bx[ix,iy,iz]/modB[ix,iy,iz])) 

            nini = np.append(nini,ni[ix,iy,iz]*30.) 
            nene = np.append(nene,ne[ix,iy,iz]*30.) 
            if(enter==0):
                imin=ii
            enter=1
        elif(enter==1):
            imax=ii
            break

    # end time and print
    end = time.time()
    print('\n### COMPUTATIONAL TIME PR'+str(PR_number)+'-run'+str(run_number)+' CYCLE = %.2f ###\n'%(end-start))

    # bow shock and MP crossings
    tBS1 = -19.5
    tMP1 = -10.0
    tMP2 = + 7.5
    tBS2 = +11.5

    # plot magnetic fields
    fig = plt.figure(1)
    plt.subplot(411)
    plt.plot(tt[imin:imax],bb,color='black',label='|B|simu')
    plt.plot(tt[imin:imax],bbX[imin:imax],color='black',label='|B|obs',alpha=0.5)
    plt.legend(loc=0,fontsize=10)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],0,140,linestyle='--',color='grey')
    plt.xlim(tt[imin],tt[imax])
    plt.xticks([])
    plt.subplot(412)
    plt.plot(tt[imin:imax],bxX[imin:imax],color='red',label=r'$Bx_{s/c}$',alpha=0.5)
    plt.plot(tt[imin:imax],byX[imin:imax],color='blue',label=r'$By_{s/c}$',alpha=0.5)
    plt.plot(tt[imin:imax],bzX[imin:imax],color='green',label=r'$Bz_{s/c}$',alpha=0.5)

    #plt.plot(tt[imin:imax],bxX[imin:imax],color='red',alpha=0.5)
    plt.legend(loc=0,ncol=3,fontsize=10)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],-40,100,linestyle='--',color='grey')
    plt.xlim(tt[imin],tt[imax])
    plt.xticks([])
    '''
    plt.subplot(513)
    plt.plot(tt[imin:imax],gf(byby,1),color='green',label='By')
    plt.plot(tt[imin:imax],byX[imin:imax],color='green',alpha=0.5)
    plt.legend(loc=0,fontsize=10)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],-40,60,linestyle='--',color='grey')
    plt.xlim(tt[imin],tt[imax])
    plt.xticks([])
    plt.subplot(514)
    plt.plot(tt[imin:imax],gf(bzbz,1),color='blue',label='Bz')
    plt.plot(tt[imin:imax],bzX[imin:imax],color='blue',alpha=0.5)
    plt.legend(loc=0,fontsize=10)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],-40,140,linestyle='--',color='grey')
    plt.xlim(tt[imin],tt[imax])
    plt.xticks([])
    '''

    # plot distribution functions
    ax1 = plt.subplot(413)
    im1=ax1.imshow(DFxx, cmap=cm.jet,  origin='lower',norm=colors.LogNorm(vmin=5e-5, vmax=5e0), extent=(tt[imin],tt[imax],10,800),aspect='auto')
    plt.ylim(10,800)
    plt.yscale('Log')
    plt.ylabel(r'$E/eV$,$T_{\bot}$',fontsize=16)
    plt.grid()
    plt.xticks([])

    # plot distribution functions
    ax2 = plt.subplot(414)
    im2=ax2.imshow(DFyy, cmap=cm.jet,  origin='lower',norm=colors.LogNorm(vmin=5e-5, vmax=5e0), extent=(tt[imin],tt[imax],10,800),aspect='auto')
    plt.ylim(10,800)
    plt.yscale('Log')
    plt.ylabel(r'$E/eV$,$T_{\parallel}$',fontsize=16)
    plt.grid()
    #plt.xticks([])

    # plot distribution functions
    '''
    ax = plt.subplot(515)
    im=ax.imshow(DFzz, cmap=cm.jet,  origin='lower',norm=colors.LogNorm(vmin=1e-5, vmax=5e0), extent=(tt[imin],tt[imax],10,800),aspect='auto')
    plt.ylim(10,800)
    plt.yscale('Log')
    plt.ylabel(r'$E/eV$',fontsize=16)
    '''
    plt.xlabel('time [CA-min]',fontsize=16)
    plt.grid()
    cax = fig.add_axes([0.915, 0.25, 0.01, 0.45])
    cbar = plt.colorbar(im2,cax=cax,format='%.0e')
    cbar.ax.tick_params(labelsize=8)
    '''
    # plot densities
    plt.subplot(515)
    if(PR_number!='2'):
        plt.plot(tt[imin:imax],gf(nene,1),color='blue',label='ne simu')
    plt.plot(tt[imin:imax],n1X[imin:imax],color='blue',label='ne obs',alpha=0.5)
    plt.legend(loc=0,fontsize=10)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],0,90,linestyle='--',color='grey')
    plt.xlim(tt[imin],tt[imax])
    '''
    plt.savefig('images/plot_MarinerX_DF-9-2_PR'+str(PR_number)+'_'+str(Ncycle)+'.png',format='png')
    plt.close()
    '''
    # plot angles
    plt.subplot(211)
    plt.plot(tt[imin:imax],clock,color='grey',label='simu')
    plt.plot(tt[imin:imax],cloX[imin:imax],color='grey',label='obs',alpha=0.5)
    plt.ylabel('clock angle')
    plt.legend(loc=0,fontsize=10)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],-160,160,linestyle='--',color='grey')
    plt.ylim(-180,180)
    plt.xlim(tt[imin],tt[imax])
    plt.subplot(212)
    plt.plot(tt[imin:imax],cone,color='orange')
    plt.plot(tt[imin:imax],conX[imin:imax],color='orange',alpha=0.5)
    plt.ylabel('cone angle')
    plt.ylim(0,180)
    plt.vlines([tBS1,tMP1,tBS2,tMP2],-160,160,linestyle='--',color='grey')
    plt.xlim(tt[imin],tt[imax])
    plt.xlabel('time [CA-min]',fontsize=16)

    plt.savefig('images/plot_MarinerX_9-angles_PR'+str(PR_number)+'_'+str(Ncycle)+'.png',format='png')
    plt.close()

    print('i0=267\nimin=%d\nimax=%d'%(imin,imax))

    plt.hist2d(cloX[imax:imax+300],conX[imax:imax+300],bins=100,range=[[-360.,360.],[-360.,360.]])
    plt.colorbar()
    plt.xlabel('clock angle',fontsize=16)
    plt.ylabel('cone angle',fontsize=16)
    plt.show()
    '''
    return np.array([xx[imin:imax],yy[imin:imax],zz[imin:imax]])





# routine to show which patch the trajectory passes trough --> optimize download restart
def pcl_patch(PR_number,fact,run_number,Ncycle,traj):

    XLEN = 16
    YLEN = 16
    ZLEN = 16

    if(PR_number=='2'):
        XLEN=32

    Nx_ = int(len(x)/XLEN)
    Ny_ = int(len(y)/YLEN)
    Nz_ = int(len(z)/ZLEN)

    # open file for output
    output = open('good_patches_MarinerX_PR'+str(PR_number)+'.txt','w')

    # array to plot
    patch_list_x = []
    patch_list_y = []
    patch_list_z = []
    color_list = []

    # loop on the processes
    ip=0
    ic=0
    for ipx in range(0,XLEN):
        ix1 = ipx*(Nx_-1)
        ix2 = ix1+Nx_+1
        for ipy in range(0,YLEN):
            iy1 = ipy*(Ny_-1)
            iy2 = iy1+Ny_+1
            for ipz in range(0,ZLEN):
                iz1 = ipz*(Nz_-1)
                iz2 = iz1+Nz_+1

                # look if any of the trajectory points are in this patch
                igoodx = (traj[0]>x[ix1])*(traj[0]<x[ix2])
                igoody = (traj[1]>y[iy1])*(traj[1]<y[iy2])
                igoodz = (traj[2]>z[iz1])*(traj[2]<z[iz2])
                booll  = igoodx*igoody*igoodz

                # write output
                if(booll.any()):
                    output.write('%d\n'%ip)
                    ic=ic+1
                    color_list.append('red')
                else:
                    color_list.append('blue')

                # increase patch index
                ip=ip+1

                # plot a sort of grid
                patch_list_x.append(x[ix1])
                patch_list_y.append(y[iy1])
                patch_list_z.append(z[iz1])



    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    plt.title('MPI processes + flyby Ncross/Nproc=%d/%d'%(ic,ip))
    for i in range(len(patch_list_x)):
        if(i%2==0):
            ax.scatter(patch_list_x[i],patch_list_y[i],patch_list_z[i],marker='o',color=color_list[i],alpha=0.4)

    ax.set_xlabel('x[Rm]',fontsize=16)
    ax.set_ylabel('y[Rm]',fontsize=16)
    ax.set_zlabel('z[Rm]',fontsize=16)

    plt.show()
    plt.savefig('images/3d_plot_patches_'+str(PR_number)+'_'+str(Ncycle)+'.png',format='png')
    plt.close()
                
                



#traj=main('1',2,1,7500)
#pcl_patch('1',2,1,7500,traj)

main('0',1,1,5500)
#main('0-bigbox',2,1,9500)
main('1',2,1,7500)
#main('2',4,28,11200)



