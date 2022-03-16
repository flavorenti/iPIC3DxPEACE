from vectorfield import *
from scipy import ndimage
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
from vtk_routines import *
import time


# main of the code: PR_number  = number production run (0,1,2,0-bigbox)
#                   fact       = scaling factor (PR0,1,2 = 1,2,4) for bigbox this is 2
#                   run_number = number of the run we want to plot
#                   Nstart     = start cycle for the plots
#                   Nend       = end cycle for the plots
#                   step       = number of cycles between two plots
#                   fig_min, fig_max = define the code of th pictures to save
def main(args):

    # extract arguments from list args
    run_number = args[2]
    Nstart     = args[3]
    Nend       = args[4]
    step       = args[5]
    fig_min    = args[6]
    fig_max    = args[7]

    # define pars simulation as global variables
    global x,y,z,dt,xpl,ypl,zpl,qome,iyc,izc

    # pass also fact and PR_number in global
    global PR_number, fact
    PR_number = str(args[0])
    fact = args[1]

    # path to simulation data, s/c trajectory file, where to save plot
    global data_path, bepi_path, save_path
    data_path = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR'+PR_number+'/run'+str(run_number)+'/data/'
    bepi_path = './Bepi-flyby-1st_MSO-trajectory1.txt'
    save_path = '/home/flavorenti/Bureau/my_scripts_only/images/'

    # box parameters, code units di,wci,c
    x   = np.linspace(0,10.,256*fact)
    y   = np.linspace(0,8.,208*fact)
    z   = np.linspace(0,8.,208*fact)
    # bigbox has double size
    if(PR_number=='0-bigbox'):
        x=x*2.
        y=y*2.
        z=z*2.
    # timestep wci-1
    dt = 0.5
    # position planet center
    xpl = 0.5*max(x)*4./5.
    ypl = 0.5*max(y)
    zpl = 0.5*max(z)-0.2
    # q/m electrons
    qome = -100.
    # cut indexes
    iyc = len(y)/2 # *2 for bigbox-cut
    izc = len(z)/2 # *2 for bigbox-cut

    # loop over the cycles
    for Ncycle in range(Nstart,Nend,step):

        # start time
        start = time.time()

        if(fig_min<=1 and fig_max>=1):
            # open the currents files (only once/cycle!)
            print('Open Vectors [Ji,Je] for cycle=%d/%d'%(Ncycle,Nend))
            Ji_vtk = open_vtk_vect(data_path+'Dipole3D_Ji_'+str(Ncycle)+'.vtk')
            Je_vtk = open_vtk_vect(data_path+'Dipole3D_Je_'+str(Ncycle)+'.vtk')
            print('Save Figure1 for cycle=%d/%d'%(Ncycle,Nend))
            Ji = Figure1_Ji(Ji_vtk,Ncycle)
            Je = Figure1_Je(Je_vtk,Ncycle)

        if(fig_min<=2 and fig_max>=2):
            # open the densities files (only once/cycle!)
            print('Open Scalars [rhoi,rhoe] for cycle=%d/%d'%(Ncycle,Nend))
            rhoi_vtk = open_vtk_scal(data_path+'Dipole3D_rhoi1_'+str(Ncycle)+'.vtk')
            rhoe_vtk = open_vtk_scal(data_path+'Dipole3D_rhoe0_'+str(Ncycle)+'.vtk')
            print('Save Figure2 for cycle=%d/%d'%(Ncycle,Nend))
            n  = Figure2(rhoi_vtk,rhoe_vtk,Ncycle)

        if(fig_min<=3 and fig_max>=3):
            # open pressure files
            print('Open Scalars [rhoi,rhoe] for cycle=%d/%d'%(Ncycle,Nend))
            Pixx_vtk = open_vtk_scal(data_path+'Dipole3D_PXXi1_'+str(Ncycle)+'.vtk')
            Piyy_vtk = open_vtk_scal(data_path+'Dipole3D_PYYi1_'+str(Ncycle)+'.vtk')
            Pizz_vtk = open_vtk_scal(data_path+'Dipole3D_PZZi1_'+str(Ncycle)+'.vtk')
            Pixy_vtk = open_vtk_scal(data_path+'Dipole3D_PXYi1_'+str(Ncycle)+'.vtk')
            Pixz_vtk = open_vtk_scal(data_path+'Dipole3D_PXZi1_'+str(Ncycle)+'.vtk')
            Piyz_vtk = open_vtk_scal(data_path+'Dipole3D_PYZi1_'+str(Ncycle)+'.vtk')
            Pexx_vtk = open_vtk_scal(data_path+'Dipole3D_PXXe0_'+str(Ncycle)+'.vtk')
            Peyy_vtk = open_vtk_scal(data_path+'Dipole3D_PYYe0_'+str(Ncycle)+'.vtk')
            Pezz_vtk = open_vtk_scal(data_path+'Dipole3D_PZZe0_'+str(Ncycle)+'.vtk')
            Pexy_vtk = open_vtk_scal(data_path+'Dipole3D_PXYe0_'+str(Ncycle)+'.vtk')
            Pexz_vtk = open_vtk_scal(data_path+'Dipole3D_PXZe0_'+str(Ncycle)+'.vtk')
            Peyz_vtk = open_vtk_scal(data_path+'Dipole3D_PYZe0_'+str(Ncycle)+'.vtk')
            print('Save Figure3 for cycle=%d/%d'%(Ncycle,Nend))
            T  = Figure3(Pixx_vtk,Piyy_vtk,Pizz_vtk,Pixy_vtk,Pixz_vtk,Piyz_vtk,Pexx_vtk,Peyy_vtk,Pezz_vtk,Pexy_vtk,Pexz_vtk,Peyz_vtk,Ji[0],Ji[1],Ji[2],Ji[3],Ji[4],Ji[5],Je[0],Je[1],Je[2],Je[3],Je[4],Je[5],n[0],n[1],n[2],n[3],Ncycle)

        if(fig_min<=4 and fig_max>=4):
            # open electric field file
            print('Open Vectors [E] for cycle=%d/%d'%(Ncycle,Nend))
            E_vtk  = open_vtk_vect(data_path+'Dipole3D_E_'+str(Ncycle)+'.vtk')
            print('Save Figure4 for cycle=%d/%d'%(Ncycle,Nend))
            E = Figure4(E_vtk,Ncycle)

        if(fig_min<=5 and fig_max>=5):
            # open magnetic field file
            print('Open Vectors [B] for cycle=%d/%d'%(Ncycle,Nend))
            B_vtk  = open_vtk_vect(data_path+'Dipole3D_B_'+str(Ncycle)+'.vtk')
            print('Save Figure5 for cycle=%d/%d'%(Ncycle,Nend))
            B = Figure5(B_vtk,Ncycle)
        
        if(fig_min<=6 and fig_max>=6):
            # compute velocity from current and density
            print('Save Figure6 for cycle=%d/%d'%(Ncycle,Nend))
            Figure6_Ui(Ji,n,Ncycle)
            Figure6_Ue(Je,n,Ncycle)

        if(fig_min<=7 and fig_max>=7):
            # compute length scales di, rhoi, debye length
            print('Save Figure7 for cycle=%d/%d'%(Ncycle,Nend))
            Figure7(n[0],n[1],n[2],n[3],B[0],B[1],T[0],T[1],T[4],T[5],Ncycle)

        if(fig_min<=8 and fig_max>=8):
            # compute temperature anisotropy
            print('Save Figure8 for cycle=%d/%d'%(Ncycle,Nend))
            Figure8(T[2],T[3],T[6],T[7],B[2],B[3],Ncycle)

        # end time and print
        end = time.time()
        print('\n### COMPUTATIONAL TIME PR'+str(PR_number)+'-run'+str(run_number)+' CYCLE%d/%d = %.2f ###\n'%(Ncycle,Nend,end-start))

        # delete arrays
        if(fig_min<=1 and fig_max>=1):
            del Ji_vtk,Je_vtk,Ji,Je
        if(fig_min<=2 and fig_max>=2):
            del rhoi_vtk,rhoe_vtk,n
        if(fig_min<=3 and fig_max>=3):
            del Pixx_vtk,Piyy_vtk,Pizz_vtk,Pixz_vtk,Pixy_vtk,Piyz_vtk,Pexx_vtk,Peyy_vtk,Pezz_vtk,Pexz_vtk,Pexy_vtk,Peyz_vtk,T
        if(fig_min<=1 and fig_max>=1):
            del E_vtk,E
        if(fig_min<=5 and fig_max>=5):
            del B_vtk,B



# unpack a 3d vtk file to a more handy 2d cut: data_file = the vtk file opened above
#                                              iy        =  slice to take in y, if None means cut in z
#                                              iz        =  slice to take in z, if None means cut in y
#                                              string    =  name of the array to unpack (B,E,Ji,Je,ni,ne,Pixx,Pexx,etc.)
def unpack(data_file,iy,iz,string):

    # XY cut of the 3d array
    if(iz==None):
        # vectors
        if(string=='Bx'):
            return convert_vtk_vect(data_file,'B')[0][::fact,iy,::fact]
        if(string=='By'):
            return convert_vtk_vect(data_file,'B')[1][::fact,iy,::fact]
        if(string=='Bz'):
            return convert_vtk_vect(data_file,'B')[2][::fact,iy,::fact]
        if(string=='Ex'):
            return convert_vtk_vect(data_file,'E')[0][::fact,iy,::fact]
        if(string=='Ey'):
            return convert_vtk_vect(data_file,'E')[1][::fact,iy,::fact]
        if(string=='Ez'):
            return convert_vtk_vect(data_file,'E')[2][::fact,iy,::fact]
        if(string=='Jxi'):
            return convert_vtk_vect(data_file,'Ji')[0][::fact,iy,::fact]
        if(string=='Jyi'):
            return convert_vtk_vect(data_file,'Ji')[1][::fact,iy,::fact]
        if(string=='Jzi'):
            return convert_vtk_vect(data_file,'Ji')[2][::fact,iy,::fact]
        if(string=='Jxe'):
            return convert_vtk_vect(data_file,'Je')[0][::fact,iy,::fact]
        if(string=='Jye'):
            return convert_vtk_vect(data_file,'Je')[1][::fact,iy,::fact]
        if(string=='Jze'):
            return convert_vtk_vect(data_file,'Je')[2][::fact,iy,::fact]
        # scalars
        scal = convert_vtk_scal(data_file)
        return scal[::fact,iy,::fact]
        
    # XZ cut of the 3d array
    if(iy==None):
        # vectors
        if(string=='Bx'):
            return convert_vtk_vect(data_file,'B')[0][::fact,::fact,iz]
        if(string=='By'):
            return convert_vtk_vect(data_file,'B')[1][::fact,::fact,iz]
        if(string=='Bz'):
            return convert_vtk_vect(data_file,'B')[2][::fact,::fact,iz]
        if(string=='Ex'):
            return convert_vtk_vect(data_file,'E')[0][::fact,::fact,iz]
        if(string=='Ey'):
            return convert_vtk_vect(data_file,'E')[1][::fact,::fact,iz]
        if(string=='Ez'):
            return convert_vtk_vect(data_file,'E')[2][::fact,::fact,iz]
        if(string=='Jxi'):
            return convert_vtk_vect(data_file,'Ji')[0][::fact,::fact,iz]
        if(string=='Jyi'):
            return convert_vtk_vect(data_file,'Ji')[1][::fact,::fact,iz]
        if(string=='Jzi'):
            return convert_vtk_vect(data_file,'Ji')[2][::fact,::fact,iz]
        if(string=='Jxe'):
            return convert_vtk_vect(data_file,'Je')[0][::fact,::fact,iz]
        if(string=='Jye'):
            return convert_vtk_vect(data_file,'Je')[1][::fact,::fact,iz]
        if(string=='Jze'):
            return convert_vtk_vect(data_file,'Je')[2][::fact,::fact,iz]
        # scalars
        scal = convert_vtk_scal(data_file)
        return scal[::fact,::fact,iz]


# add plot to figure, with (1) labels axis (bool option)
#                          (2) text
#                          (3) colorbar with specified colormap and limits (bool option)
#                          (4) planet disk (offset, radius=1)
#                          (5) logaritmic colorscale (bool option)
#                          (6) BShock and MPause profiles (bool option)
#                          (7) flyby trajectory (bool option) 
def add_plot(fig,position,field,x,y,xlabel,ylabel,text,ccmap,cmin,cmax,offset=0.,xlab=False,ylab=False,cmap=False,log=False,add_profiles=True,flyby_file=None):
    fig.add_subplot(position)
    a = field#[102:358,104:312 this only for bigbox-cut]
    a = np.transpose(a)
    #a = gf(a,fact)
    if(log==False):
        im = plt.imshow(a, cmap=ccmap, origin='lower', extent=(min(x),max(x),min(y),max(y)), aspect=1)
    else:
       im = plt.imshow(a, cmap=ccmap, origin='lower', norm=colors.LogNorm(vmin=cmin, vmax=cmax), extent=(min(x),max(x),min(y),max(y)), aspect=1)
    planet = plt.Circle((xpl, ypl-offset), 1., color='grey', fill=True, linewidth=2.)
    ax = fig.gca()
    ax.add_patch(planet)
    plt.text(0.1*max(x),0.75*max(y),text,fontsize=16)
    plt.clim(cmin,cmax)
    if(xlab):
        plt.xlabel(xlabel,fontsize=14)
    if(ylab):
        plt.ylabel(ylabel,fontsize=14)
    if(add_profiles==True):
        theta = np.linspace(-3.,3.,1000)
        Rmp     = 1.45*(2./(1.+np.cos(theta)))**(0.5)
        Rbs     = 2.75*1.04/(1+1.04*np.cos(theta))
        xx     =  xpl-Rmp*np.cos(theta)
        zz     =  ypl+Rmp*np.sin(theta)
        plt.plot(xx,zz,linestyle='--',color='black',marker='',linewidth=2)
        xx     =  xpl-0.5-Rbs*np.cos(theta)
        zz     =  ypl+Rbs*np.sin(theta)
        plt.plot(xx,zz,linestyle='--',color='grey',marker='',linewidth=2)
        plt.ylim(0,np.max(y))
        plt.xlim(0,np.max(x))
    if(flyby_file!=None):
        tt,xx,yy,zz = np.loadtxt(flyby_file,unpack=True)
        ica = np.argmin(xx**2+yy**2+zz**2)
        tt = tt-tt[ica]
        xx=-xx
        xx+=xpl
        yy+=ypl
        zz+=zpl
        for it in range(-4,5):
            tgood=tt[ica+it*300]
            if(ylabel=='y[Rm]'):
                plt.plot(xx[ica+it*300],yy[ica+it*300],marker='o',label='time=%.1f m'%(tgood/60.))
            if(ylabel=='z[Rm]'):
                plt.plot(xx[ica+it*300],zz[ica+it*300],marker='o',label='time=%.1f m'%(tgood/60.))
        if(ylabel=='y[Rm]'):    
            plt.plot(xx,yy,color='yellow')
        if(ylabel=='z[Rm]'):   
            plt.plot(xx,zz,color='yellow')
        plt.ylim(0,np.max(y))
        plt.xlim(0,np.max(x))
        if(cmap==True):
            plt.legend(loc=0,fontsize=10)
    if(cmap==True and log==False):
        cax = fig.add_axes([0.915, 0.2, 0.015, 0.6])
        cbar = plt.colorbar(im, cax=cax, format='%.0e')
        cbar.ax.tick_params(labelsize=10)
    if(cmap==True and log==True):
        cax = fig.add_axes([0.915, 0.2, 0.015, 0.6])
        cbar = plt.colorbar(im, cax=cax)
    if(cmap==5):
        cax = fig.add_axes([0.673, 0.5, 0.22, 0.01])
        cbar = plt.colorbar(im, cax=cax, format='%.0e', orientation='horizontal')
        cbar.ax.tick_params(labelsize=10)


# plot and save Fig1-1 showing cuts of 3 components of Ji
def Figure1_Ji(Ji_vtk,Ncycle):
        ### 1st figure Ji ####
        fig1 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Ion Currents t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Jxi_xz = 4.*np.pi*unpack(Ji_vtk,iyc,None,'Jxi')
        Jyi_xz = 4.*np.pi*unpack(Ji_vtk,iyc,None,'Jyi')
        Jzi_xz = 4.*np.pi*unpack(Ji_vtk,iyc,None,'Jzi')
        Jxi_xy = 4.*np.pi*unpack(Ji_vtk,None,izc,'Jxi')
        Jyi_xy = 4.*np.pi*unpack(Ji_vtk,None,izc,'Jyi')
        Jzi_xy = 4.*np.pi*unpack(Ji_vtk,None,izc,'Jzi')
        add_plot(fig1,231,Jxi_xz,x,z,'x[Rm]','z[Rm]','Jxi','jet',-0.04,0.04,offset=0.2,ylab=True)
        add_plot(fig1,232,Jyi_xz,x,z,'x[Rm]','z[Rm]','Jyi','jet',-0.04,0.04,offset=0.2,cmap=True)
        add_plot(fig1,233,Jzi_xz,x,z,'x[Rm]','z[Rm]','Jzi','jet',-0.04,0.04,offset=0.2)
        add_plot(fig1,234,Jxi_xy,x,y,'x[Rm]','y[Rm]','Jxi','jet',-0.04,0.04,xlab=True,ylab=True)
        add_plot(fig1,235,Jyi_xy,x,y,'x[Rm]','y[Rm]','Jyi','jet',-0.04,0.04,xlab=True)
        add_plot(fig1,236,Jzi_xy,x,y,'x[Rm]','y[Rm]','Jzi','jet',-0.04,0.04,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig1-1_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        return np.array([Jxi_xz,Jxi_xy,Jyi_xz,Jyi_xy,Jzi_xz,Jzi_xy])

# plot and save Fig1-2 showing cuts of 3 components of Je
def Figure1_Je(Je_vtk,Ncycle):
        ### 1st figure Je ####
        fig1 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Electron Currents t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Jxe_xz = 4.*np.pi*unpack(Je_vtk,iyc,None,'Jxe')
        Jye_xz = 4.*np.pi*unpack(Je_vtk,iyc,None,'Jye')
        Jze_xz = 4.*np.pi*unpack(Je_vtk,iyc,None,'Jze')
        Jxe_xy = 4.*np.pi*unpack(Je_vtk,None,izc,'Jxe')
        Jye_xy = 4.*np.pi*unpack(Je_vtk,None,izc,'Jye')
        Jze_xy = 4.*np.pi*unpack(Je_vtk,None,izc,'Jze')
        add_plot(fig1,231,Jxe_xz,x,z,'x[Rm]','z[Rm]','Jxe','jet_r',-0.04,0.04,offset=0.2,ylab=True)
        add_plot(fig1,232,Jye_xz,x,z,'x[Rm]','z[Rm]','Jye','jet_r',-0.04,0.04,offset=0.2,cmap=True)
        add_plot(fig1,233,Jze_xz,x,z,'x[Rm]','z[Rm]','Jze','jet_r',-0.04,0.04,offset=0.2)
        add_plot(fig1,234,Jxe_xy,x,y,'x[Rm]','y[Rm]','Jxe','jet_r',-0.04,0.04,xlab=True,ylab=True)
        add_plot(fig1,235,Jye_xy,x,y,'x[Rm]','y[Rm]','Jye','jet_r',-0.04,0.04,xlab=True)
        add_plot(fig1,236,Jze_xy,x,y,'x[Rm]','y[Rm]','Jze','jet_r',-0.04,0.04,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig1-2_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        return np.array([Jxe_xz,Jxe_xy,Jye_xz,Jye_xy,Jze_xz,Jze_xy])


# plot and save Fig2 showing cuts of ni,ne
def Figure2(rhoi_vtk,rhoe_vtk,Ncycle):
        ### 2nd figure n ####
        fig2 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Densities t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        ni_xz =  unpack(rhoi_vtk,iyc,None,'ni')
        ne_xz = -unpack(rhoe_vtk,iyc,None,'ne')
        ni_xy =  unpack(rhoi_vtk,None,izc,'ni')
        ne_xy = -unpack(rhoe_vtk,None,izc,'ne')
        add_plot(fig2,231,30.*ni_xz,x,z,'x[Rm]','z[Rm]','ni[cm-3]','jet',1.,120.,offset=0.2,ylab=True,log=True,flyby_file=bepi_path)
        add_plot(fig2,232,30.*ne_xz,x,z,'x[Rm]','z[Rm]','ne[cm-3]','jet',1.,120.,offset=0.2,cmap=True,log=True,flyby_file=bepi_path)
        add_plot(fig2,234,30.*ni_xy,x,y,'x[Rm]','y[Rm]','ni','jet',1.,120.,xlab=True,ylab=True,log=True,flyby_file=bepi_path)
        add_plot(fig2,235,30.*ne_xy,x,y,'x[Rm]','y[Rm]','ne','jet',1.,120.,xlab=True,log=True,flyby_file=bepi_path)
        add_plot(fig2,233,gf((ni_xz-ne_xz),fact)/gf(ni_xz,fact),x,z,'x[Rm]','z[Rm]',r'$\Delta$n/n','BrBG',-1e-1,1e-1,offset=0.2,cmap=5)
        add_plot(fig2,236,gf((ni_xy-ne_xy),fact)/gf(ni_xy,fact),x,y,'x[Rm]','y[Rm]',r'$\Delta$n/n','BrBG',-1e-1,1e-1,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig2_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        ni_xz[np.where(ni_xz==0)]=1e-6
        ni_xy[np.where(ni_xy==0)]=1e-6
        ne_xz[np.where(ne_xz==0)]=1e-6
        ne_xy[np.where(ne_xy==0)]=1e-6
        return np.array([gf(ni_xz,fact),gf(ni_xy,fact),gf(ne_xz,fact),gf(ne_xy,fact)])


# plot and save Fig3 showing cuts of P = (Pxx+Pyy+Pzz)/3
# also plot figure with T = P/n
# also compute full tensor of cut Temperatures (ions and electrons) and return it
def Figure3(Pixx_vtk,Piyy_vtk,Pizz_vtk,Pixy_vtk,Pixz_vtk,Piyz_vtk,Pexx_vtk,Peyy_vtk,Pezz_vtk,Pexy_vtk,Pexz_vtk,Peyz_vtk,Jxi_xz,Jxi_xy,Jyi_xz,Jyi_xy,Jzi_xz,Jzi_xy,Jxe_xz,Jxe_xy,Jye_xz,Jye_xy,Jze_xz,Jze_xy,ni_xz,ni_xy,ne_xz,ne_xy,Ncycle): 
        ### 3rd figure P ####
        fig3 = plt.figure()
        plt.title('Mercury SHOTS-Init Pressures t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Pixx_xz = unpack(Pixx_vtk,iyc,None,'Pixx')*4.*np.pi
        Pexx_xz = unpack(Pexx_vtk,iyc,None,'Pexx')*4.*np.pi/qome
        Pixx_xy = unpack(Pixx_vtk,None,izc,'Pixx')*4.*np.pi
        Pexx_xy = unpack(Pexx_vtk,None,izc,'Pexx')*4.*np.pi/qome
        Piyy_xz = unpack(Piyy_vtk,iyc,None,'Piyy')*4.*np.pi
        Peyy_xz = unpack(Peyy_vtk,iyc,None,'Peyy')*4.*np.pi/qome
        Piyy_xy = unpack(Piyy_vtk,None,izc,'Piyy')*4.*np.pi
        Peyy_xy = unpack(Peyy_vtk,None,izc,'Peyy')*4.*np.pi/qome
        Pizz_xz = unpack(Pizz_vtk,iyc,None,'Pizz')*4.*np.pi
        Pezz_xz = unpack(Pezz_vtk,iyc,None,'Pezz')*4.*np.pi/qome
        Pizz_xy = unpack(Pizz_vtk,None,izc,'Pizz')*4.*np.pi
        Pezz_xy = unpack(Pezz_vtk,None,izc,'Pezz')*4.*np.pi/qome
        Pixx_xz = (Pixx_xz-(Jxi_xz**2/gf(ni_xz,fact)))
        Pexx_xz = (Pexx_xz+(Jxe_xz**2/gf(ne_xz,fact))/qome)
        Pixx_xy = (Pixx_xy-(Jxi_xy**2/gf(ni_xy,fact)))
        Pexx_xy = (Pexx_xy+(Jxe_xy**2/gf(ne_xy,fact))/qome)
        Piyy_xz = (Piyy_xz-(Jyi_xz**2/gf(ni_xz,fact)))
        Peyy_xz = (Peyy_xz+(Jye_xz**2/gf(ne_xz,fact))/qome)
        Piyy_xy = (Piyy_xy-(Jyi_xy**2/gf(ni_xy,fact)))
        Peyy_xy = (Peyy_xy+(Jye_xy**2/gf(ne_xy,fact))/qome)
        Pizz_xz = (Pizz_xz-(Jzi_xz**2/gf(ni_xz,fact)))
        Pezz_xz = (Pezz_xz+(Jze_xz**2/gf(ne_xz,fact))/qome)
        Pizz_xy = (Pizz_xy-(Jzi_xy**2/gf(ni_xy,fact)))
        Pezz_xy = (Pezz_xy+(Jze_xy**2/gf(ne_xy,fact))/qome)
        Pi_xz = (Pixx_xz+Piyy_xz+Pizz_xz)/3.
        Pe_xz = (Pexx_xz+Peyy_xz+Pezz_xz)/3.
        Pi_xy = (Pixx_xy+Piyy_xy+Pizz_xy)/3.
        Pe_xy = (Pexx_xy+Peyy_xy+Pezz_xy)/3.
        add_plot(fig3,221,Pi_xz,x,z,'x[Rm]','z[Rm]','Pi','copper',0,1e-4,offset=0.2,ylab=True)
        add_plot(fig3,222,abs(Pe_xz),x,z,'x[Rm]','z[Rm]','Pe','copper',0,1e-4,offset=0.2,cmap=True)
        add_plot(fig3,223,Pi_xy,x,y,'x[Rm]','y[Rm]','Pi','copper',0,1e-4,xlab=True,ylab=True)
        add_plot(fig3,224,abs(Pe_xy),x,y,'x[Rm]','y[Rm]','Pe','copper',0,1e-4,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig3-P_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        ### 3th figure T ###
        fig3 = plt.figure()
        plt.title('Mercury SHOTS-Init Temperatures t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Ti_xz = Pi_xz/gf(ni_xz,fact)
        Te_xz = Pe_xz/gf(ne_xz,fact)
        Ti_xy = Pi_xy/gf(ni_xy,fact)
        Te_xy = Pe_xy/gf(ne_xy,fact)
        add_plot(fig3,221,Ti_xz,x,z,'x[Rm]','z[Rm]','Ti','copper',0,1e-4,offset=0.2,ylab=True)
        add_plot(fig3,222,abs(Te_xz),x,z,'x[Rm]','z[Rm]','Te','copper',0,1e-4,offset=0.2,cmap=True)
        add_plot(fig3,223,Ti_xy,x,y,'x[rm]','y[Rm]','Ti','copper',0,1e-4,xlab=True,ylab=True)
        add_plot(fig3,224,abs(Te_xy),x,y,'x[Rm]','y[Rm]','Te','copper',0,1e-4,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig3-T_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        ### make T tensor to return it ###
        Tensi_xz = np.zeros((3,3,len(x)/fact,len(y)/fact))
        Tensi_xy = np.zeros((3,3,len(x)/fact,len(y)/fact))
        Tense_xz = np.zeros((3,3,len(x)/fact,len(y)/fact))
        Tense_xy = np.zeros((3,3,len(x)/fact,len(y)/fact))
        Pixy_xz = unpack(Pixy_vtk,iyc,None,'Pixy')*4.*np.pi
        Pixz_xz = unpack(Pixz_vtk,iyc,None,'Pixz')*4.*np.pi
        Piyz_xz = unpack(Piyz_vtk,iyc,None,'Piyz')*4.*np.pi
        Pixy_xy = unpack(Pixy_vtk,None,izc,'Pixy')*4.*np.pi
        Pixz_xy = unpack(Pixz_vtk,None,izc,'Pixz')*4.*np.pi
        Piyz_xy = unpack(Piyz_vtk,None,izc,'Piyz')*4.*np.pi
        Pexy_xz = unpack(Pexy_vtk,iyc,None,'Pexy')*4.*np.pi/qome
        Pexz_xz = unpack(Pexz_vtk,iyc,None,'Pexz')*4.*np.pi/qome
        Peyz_xz = unpack(Peyz_vtk,iyc,None,'Peyz')*4.*np.pi/qome
        Pexy_xy = unpack(Pexy_vtk,None,izc,'Pexy')*4.*np.pi/qome
        Pexz_xy = unpack(Pexz_vtk,None,izc,'Pexz')*4.*np.pi/qome
        Peyz_xy = unpack(Peyz_vtk,None,izc,'Peyz')*4.*np.pi/qome
        Tensi_xz[0,0,:,:] = Pixx_xz/gf(ni_xz,fact) # Temp ions cut xz 
        Tensi_xz[0,1,:,:] = Pixy_xz/gf(ni_xz,fact)
        Tensi_xz[0,2,:,:] = Pixz_xz/gf(ni_xz,fact)
        Tensi_xz[1,0,:,:] = Tensi_xz[0,1,:,:]
        Tensi_xz[1,1,:,:] = Piyy_xz/gf(ni_xz,fact)
        Tensi_xz[1,2,:,:] = Piyz_xz/gf(ni_xz,fact)
        Tensi_xz[2,0,:,:] = Tensi_xz[0,2,:,:]
        Tensi_xz[2,1,:,:] = Tensi_xz[1,2,:,:]
        Tensi_xz[2,2,:,:] = Pizz_xy/gf(ni_xy,fact)
        Tensi_xy[0,0,:,:] = Pixx_xy/gf(ni_xy,fact) # Temp ions cut xy
        Tensi_xy[0,1,:,:] = Pixy_xy/gf(ni_xy,fact)
        Tensi_xy[0,2,:,:] = Pixz_xy/gf(ni_xy,fact)
        Tensi_xy[1,0,:,:] = Tensi_xy[0,1,:,:]
        Tensi_xy[1,1,:,:] = Piyy_xy/gf(ni_xy,fact)
        Tensi_xy[1,2,:,:] = Piyz_xy/gf(ni_xy,fact)
        Tensi_xy[2,0,:,:] = Tensi_xy[0,2,:,:]
        Tensi_xy[2,1,:,:] = Tensi_xy[1,2,:,:]
        Tensi_xy[2,2,:,:] = Pizz_xy/gf(ni_xy,fact)
        Tense_xz[0,0,:,:] = Pexx_xz/gf(ne_xz,fact) # Temp electrons cut xz 
        Tense_xz[0,1,:,:] = Pexy_xz/gf(ne_xz,fact)
        Tense_xz[0,2,:,:] = Pexz_xz/gf(ne_xz,fact)
        Tense_xz[1,0,:,:] = Tense_xz[0,1,:,:]
        Tense_xz[1,1,:,:] = Peyy_xz/gf(ne_xz,fact)
        Tense_xz[1,2,:,:] = Peyz_xz/gf(ne_xz,fact)
        Tense_xz[2,0,:,:] = Tense_xz[0,2,:,:]
        Tense_xz[2,1,:,:] = Tense_xz[1,2,:,:]
        Tense_xz[2,2,:,:] = Pezz_xy/gf(ne_xy,fact)
        Tense_xy[0,0,:,:] = Pexx_xy/gf(ne_xy,fact) # Temp electrons cut xy
        Tense_xy[0,1,:,:] = Pexy_xy/gf(ne_xy,fact)
        Tense_xy[0,2,:,:] = Pexz_xy/gf(ne_xy,fact)
        Tense_xy[1,0,:,:] = Tense_xy[0,1,:,:]
        Tense_xy[1,1,:,:] = Peyy_xy/gf(ne_xy,fact)
        Tense_xy[1,2,:,:] = Peyz_xy/gf(ne_xy,fact)
        Tense_xy[2,0,:,:] = Tense_xy[0,2,:,:]
        Tense_xy[2,1,:,:] = Tense_xy[1,2,:,:]
        Tense_xy[2,2,:,:] = Pezz_xy/gf(ne_xy,fact)
        return np.array([Ti_xz,Ti_xy,Tensi_xz,Tensi_xy,Te_xz,Te_xy,Tense_xz,Tense_xy])


# plot and save cuts of Electric field
def Figure4(E_vtk,Ncycle):
        ### 4th figure E ####
        fig4 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Electric field t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Ex_xz = unpack(E_vtk,iyc,None,'Ex')
        Ey_xz = unpack(E_vtk,iyc,None,'Ey')
        Ez_xz = unpack(E_vtk,iyc,None,'Ez')
        Ex_xy = unpack(E_vtk,None,izc,'Ex')
        Ey_xy = unpack(E_vtk,None,izc,'Ey')
        Ez_xy = unpack(E_vtk,None,izc,'Ez')
        add_plot(fig4,231,Ex_xz,x,z,'x[Rm]','z[Rm]','Ex','jet',-6e-4,6e-4,offset=0.2,ylab=True)
        add_plot(fig4,232,Ey_xz,x,z,'x[Rm]','z[Rm]','Ey','jet',-6e-4,6e-4,offset=0.2,cmap=True)
        add_plot(fig4,233,Ez_xz,x,z,'x[Rm]','z[Rm]','Ez','jet',-6e-4,6e-4,offset=0.2)
        add_plot(fig4,234,Ex_xy,x,y,'x[Rm]','y[Rm]','Ex','jet',-6e-4,6e-4,xlab=True,ylab=True)
        add_plot(fig4,235,Ey_xy,x,y,'x[Rm]','y[Rm]','Ey','jet',-6e-4,6e-4,xlab=True)
        add_plot(fig4,236,Ez_xy,x,y,'x[Rm]','y[Rm]','Ez','jet',-6e-4,6e-4,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig4_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        return np.array([Ex_xz, Ex_xy, Ey_xz, Ey_xy, Ez_xz, Ez_xy])


# plot and save magnetic field cut
# return magnetic field versor (B/|B|)
def Figure5(B_vtk,Ncycle):        
        ### 5th figure B ####
        fig5 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Magnetic field t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Bx_xz = gf(unpack(B_vtk,iyc,None,'Bx'),2*fact)
        By_xz = gf(unpack(B_vtk,iyc,None,'By'),2*fact)
        Bz_xz = gf(unpack(B_vtk,iyc,None,'Bz'),2*fact)
        Bx_xy = gf(unpack(B_vtk,None,izc,'Bx'),2*fact)
        By_xy = gf(unpack(B_vtk,None,izc,'By'),2*fact)
        Bz_xy = gf(unpack(B_vtk,None,izc,'Bz'),2*fact)
        add_plot(fig5,231,Bx_xz,x,z,'x[Rm]','z[Rm]','Bx','jet',-0.03,0.03,offset=0.2,ylab=True)
        add_plot(fig5,232,By_xz,x,z,'x[Rm]','z[Rm]','By','jet',-0.03,0.03,offset=0.2,cmap=True)
        add_plot(fig5,233,Bz_xz,x,z,'x[Rm]','z[Rm]','Bz','jet',-0.03,0.03,offset=0.2)
        add_plot(fig5,234,Bx_xy,x,y,'x[Rm]','y[Rm]','Bx','jet',-0.03,0.03,xlab=True,ylab=True)
        add_plot(fig5,235,By_xy,x,y,'x[Rm]','y[Rm]','By','jet',-0.03,0.03,xlab=True)
        add_plot(fig5,236,Bz_xy,x,y,'x[Rm]','y[Rm]','Bz','jet',-0.03,0.03,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig5_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        ### compute versor and return it ###
        modB_xz = np.sqrt(Bx_xz**2 + By_xz**2 + Bz_xz**2)
        Bvers_xz = np.array([Bx_xz/modB_xz,By_xz/modB_xz,Bz_xz/modB_xz])
        modB_xy = np.sqrt(Bx_xy**2 + By_xy**2 + Bz_xy**2)
        Bvers_xy = np.array([Bx_xy/modB_xy,By_xy/modB_xy,Bz_xy/modB_xy])
        return np.array([modB_xz,modB_xy,Bvers_xz,Bvers_xy])
       

# plot and save Fig6-1 showing cuts of 3 components of Ui
def Figure6_Ui(Ji,n,Ncycle):
        ### 1st figure Ji ####
        fig1 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Ion Velocities t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Uxi_xz = Ji[0]/n[0]
        Uxi_xy = Ji[1]/n[1]
        Uyi_xz = Ji[2]/n[0]
        Uyi_xy = Ji[3]/n[1]
        Uzi_xz = Ji[4]/n[0]
        Uzi_xy = Ji[5]/n[1]
        add_plot(fig1,231,Uxi_xz,x,z,'x[Rm]','z[Rm]','Uxi','jet',-0.04,0.04,offset=0.2,ylab=True)
        add_plot(fig1,232,Uyi_xz,x,z,'x[Rm]','z[Rm]','Uyi','jet',-0.04,0.04,offset=0.2,cmap=True)
        add_plot(fig1,233,Uzi_xz,x,z,'x[Rm]','z[Rm]','Uzi','jet',-0.04,0.04,offset=0.2)
        add_plot(fig1,234,Uxi_xy,x,y,'x[Rm]','y[Rm]','Uxi','jet',-0.04,0.04,xlab=True,ylab=True)
        add_plot(fig1,235,Uyi_xy,x,y,'x[Rm]','y[Rm]','Uyi','jet',-0.04,0.04,xlab=True)
        add_plot(fig1,236,Uzi_xy,x,y,'x[Rm]','y[Rm]','Uzi','jet',-0.04,0.04,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig6-1_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del Uxi_xz,Uxi_xy,Uyi_xz,Uyi_xy,Uzi_xz,Uzi_xy

# plot and save Fig6-2 showing cuts of 3 components of Ue
def Figure6_Ue(Je,n,Ncycle):
        ### 1st figure Ji ####
        fig1 = plt.figure(figsize=(20,10))
        plt.title('Mercury SHOTS-Init Electron Velocities t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Uxe_xz = Je[0]/n[2]
        Uxe_xy = Je[1]/n[3]
        Uye_xz = Je[2]/n[2]
        Uye_xy = Je[3]/n[3]
        Uze_xz = Je[4]/n[2]
        Uze_xy = Je[5]/n[3]
        add_plot(fig1,231,Uxe_xz,x,z,'x[Rm]','z[Rm]','Uxe','jet',-0.04,0.04,offset=0.2,ylab=True)
        add_plot(fig1,232,Uye_xz,x,z,'x[Rm]','z[Rm]','Uye','jet',-0.04,0.04,offset=0.2,cmap=True)
        add_plot(fig1,233,Uze_xz,x,z,'x[Rm]','z[Rm]','Uze','jet',-0.04,0.04,offset=0.2)
        add_plot(fig1,234,Uxe_xy,x,y,'x[Rm]','y[Rm]','Uxe','jet',-0.04,0.04,xlab=True,ylab=True)
        add_plot(fig1,235,Uye_xy,x,y,'x[Rm]','y[Rm]','Uye','jet',-0.04,0.04,xlab=True)
        add_plot(fig1,236,Uze_xy,x,y,'x[Rm]','y[Rm]','Uze','jet',-0.04,0.04,xlab=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig6-2_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del Uxe_xz,Uxe_xy,Uye_xz,Uye_xy,Uze_xz,Uze_xy

# plot and save Fig7 cartography of spatial scales Debye length, skin depth, larmor radius
def Figure7(ni_xz,ni_xy,ne_xz,ne_xy,modB_xz,modB_xy,Ti_xz,Ti_xy,Te_xz,Te_xy,Ncycle):
        dx = 0.1       #in di
        B0 = 0.00562   #Bsw in code units
        T0 = 0.0031**2 #Tsw in code units
        di_xz = 1./np.sqrt(ni_xz)/dx
        di_xy = 1./np.sqrt(ni_xy)/dx
        ri_xz = np.sqrt(Ti_xz)/(modB_xz)/dx
        ri_xy = np.sqrt(Ti_xy)/(modB_xy)/dx
        ld_xz = np.sqrt(abs(qome*Te_xz/ne_xz))/dx
        ld_xy = np.sqrt(abs(qome*Te_xy/ne_xy))/dx
        fpefce_xz = 20.*np.sqrt(ne_xz)/modB_xz/np.sqrt(abs(qome))
        fpefce_xy = 20.*np.sqrt(ne_xy)/modB_xy/np.sqrt(abs(qome))
        fig1 = plt.figure()
        plt.title('Mercury full-kinetic global simulation northward IMF',fontsize=20)#ales t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        #add_plot(fig1,231,di_xz,x,z,'x[Rm]','z[Rm]','di/dx','seismic_r',1e-1,1e3,offset=0.2,ylab=True,log=True)
        #add_plot(fig1,232,ri_xz,x,z,'x[Rm]','z[Rm]','ri/dx','seismic_r',1e-1,1e3,offset=0.2,cmap=True,log=True,add_profiles=True)
        add_plot(fig1,121,fpefce_xz,x,z,'$x[R_m]$','$z[R_m]$',r'$f_{pe}/f_{ce}$','seismic_r',1e-1,5e3,offset=0.2,xlab=True,ylab=True,log=True)
        #add_plot(fig1,234,di_xy,x,y,'x[Rm]','y[Rm]','di/dx','seismic_r',1e-1,1e3,xlab=True,ylab=True,log=True)
        #add_plot(fig1,235,ri_xy,x,y,'x[Rm]','y[Rm]','ri/dx','seismic_r',1e-1,1e3,xlab=True,log=True)
        add_plot(fig1,122,fpefce_xy,x,y,'$x[R_m]$','$y[R_m]$',r'$f_{pe}/f_{ce}$','seismic_r',1e-1,5e3,xlab=True,ylab=True,log=True,cmap=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig7_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()

# plot and save Fig8 anisotropy for ions and electrons
def Figure8(Tensi_xz,Tensi_xy,Tense_xz,Tense_xy,Bvers_xz,Bvers_xy,Ncycle):
        Tpari_xz = 0.
        Tpari_xy = 0.
        Tperi_xz = 0.
        Tperi_xy = 0.
        Tpare_xz = 0.
        Tpare_xy = 0.
        Tpere_xz = 0.
        Tpere_xy = 0.
        for i in range(0,3):
            for j in range(0,3):
                Tpari_xz += Tensi_xz[i,j,:,:]*Bvers_xz[i]*Bvers_xz[j] 
                Tpari_xy += Tensi_xy[i,j,:,:]*Bvers_xy[i]*Bvers_xy[j] 
                Tperi_xz += gf(Tensi_xz[i,j,:,:],fact)*(int(i==j)-Bvers_xz[i]*Bvers_xz[j])
                Tperi_xy += gf(Tensi_xy[i,j,:,:],fact)*(int(i==j)-Bvers_xy[i]*Bvers_xy[j])
                Tpare_xz += Tense_xz[i,j,:,:]*Bvers_xz[i]*Bvers_xz[j]
                Tpare_xy += Tense_xy[i,j,:,:]*Bvers_xy[i]*Bvers_xy[j] 
                Tpere_xz += gf(Tense_xz[i,j,:,:],fact)*(int(i==j)-Bvers_xz[i]*Bvers_xz[j])
                Tpere_xy += gf(Tense_xy[i,j,:,:],fact)*(int(i==j)-Bvers_xy[i]*Bvers_xy[j])
        fig1 = plt.figure()
        plt.title(r'Mercury SHOTS-Init Anisotropy $T_{\parallel}/T_{\bot}$ t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        add_plot(fig1,221,Tpari_xz/Tperi_xz,x,z,'x[Rm]','z[Rm]',r'','BrBG',1.e-2,1.e2,offset=0.2,ylab=True,log=True)
        add_plot(fig1,222,Tpare_xz/Tpere_xz,x,z,'x[Rm]','z[Rm]',r'','BrBG',1.e-2,1.e2,offset=0.2,cmap=True,add_profiles=True,log=True)
        add_plot(fig1,223,Tpari_xy/Tperi_xy,x,y,'x[Rm]','z[Rm]',r'','BrBG',1.e-2,1.e2,offset=0.2,xlab=True,ylab=True,log=True)
        add_plot(fig1,224,Tpare_xy/Tpere_xy,x,y,'x[Rm]','y[Rm]',r'','BrBG',1.e-2,1.e2,offset=0.2,xlab=True,log=True)
        plt.savefig(save_path+'Mercury_SaeInit_Fig8_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
