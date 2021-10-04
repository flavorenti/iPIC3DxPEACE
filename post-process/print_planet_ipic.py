from vectorfield import *
from scipy import ndimage
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
import glob
import os
from vtk_routines import *



def unpack(data_path,Ncycle,iy,iz,string):

    print('Unpack Ncycle =',Ncycle)
    print('slicing at iy =',iy)
    print('slicing at iz =',iz)

    B_path = data_path+'Dipole3D_B_'+str(Ncycle)+'.vtk'
    E_path = data_path+'Dipole3D_E_'+str(Ncycle)+'.vtk'
    Ji_path = data_path+'Dipole3D_Ji_'+str(Ncycle)+'.vtk'
    Je_path = data_path+'Dipole3D_Je_'+str(Ncycle)+'.vtk'
    rhoi_path = data_path+'Dipole3D_rhoi1_'+str(Ncycle)+'.vtk'
    rhoe_path = data_path+'Dipole3D_rhoe0_'+str(Ncycle)+'.vtk'
    Pxxe_path = data_path+'Dipole3D_PXXe0_'+str(Ncycle)+'.vtk'
    Pyye_path = data_path+'Dipole3D_PYYe0_'+str(Ncycle)+'.vtk'
    Pzze_path = data_path+'Dipole3D_PZZe0_'+str(Ncycle)+'.vtk'
    Pxye_path = data_path+'Dipole3D_PXYe0_'+str(Ncycle)+'.vtk'
    Pyze_path = data_path+'Dipole3D_PYZe0_'+str(Ncycle)+'.vtk'
    Pxze_path = data_path+'Dipole3D_PXZe0_'+str(Ncycle)+'.vtk'
    Pxxi_path = data_path+'Dipole3D_PXXi1_'+str(Ncycle)+'.vtk'
    Pyyi_path = data_path+'Dipole3D_PYYi1_'+str(Ncycle)+'.vtk'
    Pzzi_path = data_path+'Dipole3D_PZZi1_'+str(Ncycle)+'.vtk'
    Pxyi_path = data_path+'Dipole3D_PXYi1_'+str(Ncycle)+'.vtk'
    Pyzi_path = data_path+'Dipole3D_PYZi1_'+str(Ncycle)+'.vtk'
    Pxzi_path = data_path+'Dipole3D_PXZi1_'+str(Ncycle)+'.vtk'

    if(iz==None):
        if(string=='Bx'):
            return convert_vtk_vect(B_path,'B')[0][:,iy,:]
        if(string=='By'):
            return convert_vtk_vect(B_path,'B')[1][:,iy,:]
        if(string=='Bz'):
            return convert_vtk_vect(B_path,'B')[2][:,iy,:]
        if(string=='Ex'):
            return convert_vtk_vect(E_path,'E')[0][:,iy,:]
        if(string=='Ey'):
            return convert_vtk_vect(E_path,'E')[1][:,iy,:]
        if(string=='Ez'):
            return convert_vtk_vect(E_path,'E')[2][:,iy,:]
        if(string=='Jxi'):
            return convert_vtk_vect(Ji_path,'Ji')[0][:,iy,:]
        if(string=='Jyi'):
            return convert_vtk_vect(Ji_path,'Ji')[1][:,iy,:]
        if(string=='Jzi'):
            return convert_vtk_vect(Ji_path,'Ji')[2][:,iy,:]
        if(string=='Jxe'):
            return convert_vtk_vect(Je_path,'Je')[0][:,iy,:]
        if(string=='Jye'):
            return convert_vtk_vect(Je_path,'Je')[1][:,iy,:]
        if(string=='Jze'):
            return convert_vtk_vect(Je_path,'Je')[2][:,iy,:]
        if(string=='Uxi'):
            return ( convert_vtk_vect(Ji_path,'Ji')[0]/(convert_vtk_scal(rhoi_path)+1e-6) )[:,iy,:]
        if(string=='Uyi'):
            return ( convert_vtk_vect(Ji_path,'Ji')[1]/(convert_vtk_scal(rhoi_path)+1e-6) )[:,iy,:]
        if(string=='Uzi'):
            return ( convert_vtk_vect(Ji_path,'Ji')[2]/(convert_vtk_scal(rhoi_path)+1e-6) )[:,iy,:]
        if(string=='Uxe'):
            return ( convert_vtk_vect(Je_path,'Je')[0]/(convert_vtk_scal(rhoe_path)+1e-6) )[:,iy,:]
        if(string=='Uye'):
            return ( convert_vtk_vect(Je_path,'Je')[1]/(convert_vtk_scal(rhoe_path)+1e-6) )[:,iy,:]
        if(string=='Uze'):
            return ( convert_vtk_vect(Je_path,'Je')[2]/(convert_vtk_scal(rhoe_path)+1e-6) )[:,iy,:]
        if(string=='ni'):
            return convert_vtk_scal(rhoi_path)[:,iy,:]
        if(string=='ne'):
            return -convert_vtk_scal(rhoe_path)[:,iy,:]
        if(string=='Pi'):
            return (1./3.)*( convert_vtk_scal(Pxxi_path)+convert_vtk_scal(Pyyi_path)+convert_vtk_scal(Pzzi_path) )[:,iy,:]
        if(string=='Pe'):
            return (1./3.)*( convert_vtk_scal(Pxxe_path)+convert_vtk_scal(Pyye_path)+convert_vtk_scal(Pzze_path) )[:,iy,:]

    if(iy==None):
        if(string=='Bx'):
            return convert_vtk_vect(B_path,'B')[0][:,:,iz]
        if(string=='By'):
            return convert_vtk_vect(B_path,'B')[1][:,:,iz]
        if(string=='Bz'):
            return convert_vtk_vect(B_path,'B')[2][:,:,iz]
        if(string=='Ex'):
            return convert_vtk_vect(E_path,'E')[0][:,:,iz]
        if(string=='Ey'):
            return convert_vtk_vect(E_path,'E')[1][:,:,iz]
        if(string=='Ez'):
            return convert_vtk_vect(E_path,'E')[2][:,:,iz]
        if(string=='Jxi'):
            return convert_vtk_vect(Ji_path,'Ji')[0][:,:,iz]
        if(string=='Jyi'):
            return convert_vtk_vect(Ji_path,'Ji')[1][:,:,iz]
        if(string=='Jzi'):
            return convert_vtk_vect(Ji_path,'Ji')[2][:,:,iz]
        if(string=='Jxe'):
            return convert_vtk_vect(Je_path,'Je')[0][:,:,iz]
        if(string=='Jye'):
            return convert_vtk_vect(Je_path,'Je')[1][:,:,iz]
        if(string=='Jze'):
            return convert_vtk_vect(Je_path,'Je')[2][:,:,iz]
        if(string=='Uxi'):
            return ( convert_vtk_vect(Ji_path,'Ji')[0]/(convert_vtk_scal(rhoi_path)+1e-6) )[:,:,iz]
        if(string=='Uyi'):
            return ( convert_vtk_vect(Ji_path,'Ji')[1]/(convert_vtk_scal(rhoi_path)+1e-6) )[:,:,iz]
        if(string=='Uzi'):
            return ( convert_vtk_vect(Ji_path,'Ji')[2]/(convert_vtk_scal(rhoi_path)+1e-6) )[:,:,iz]
        if(string=='Uxe'):
            return ( convert_vtk_vect(Je_path,'Je')[0]/(convert_vtk_scal(rhoe_path)+1e-6) )[:,:,iz]
        if(string=='Uye'):
            return ( convert_vtk_vect(Je_path,'Je')[1]/(convert_vtk_scal(rhoe_path)+1e-6) )[:,:,iz]
        if(string=='Uze'):
            return ( convert_vtk_vect(Je_path,'Je')[2]/(convert_vtk_scal(rhoe_path)+1e-6) )[:,:,iz]
        if(string=='ni'):
            return convert_vtk_scal(rhoi_path)[:,:,iz]
        if(string=='ne'):
            return -convert_vtk_scal(rhoe_path)[:,:,iz]
        if(string=='Pi'):
            return (1./3.)*( convert_vtk_scal(Pxxi_path)+convert_vtk_scal(Pyyi_path)+convert_vtk_scal(Pzzi_path) )[:,:,iz]
        if(string=='Pe'):
            return (1./3.)*( convert_vtk_scal(Pxxe_path)+convert_vtk_scal(Pyye_path)+convert_vtk_scal(Pzze_path) )[:,:,iz]




def add_plot(fig,position,field,x,y,yplanet,xlabel,ylabel,text,colors,cmin,cmax,xlab=False,ylab=False,cmap=False):
    fig.add_subplot(position)
    a = np.transpose(field)
    im = plt.imshow(a, cmap=colors, origin='lower', extent=(min(x),max(x),min(y),max(y)), aspect=4./5.)
    planet = plt.Circle((4./5.*max(x)/2., yplanet), 1., color='grey')
    ax = fig.gca()
    ax.add_patch(planet)
    plt.text(0.1*max(x),0.75*max(y),text,fontsize=16)
    plt.clim(cmin,cmax)
    if(cmap):
        cax = fig.add_axes([0.915, 0.2, 0.015, 0.6])
        cbar = plt.colorbar(im, cax=cax, format='%.0e')
        cbar.ax.tick_params(labelsize=11)
    if(xlab):
        plt.xlabel(xlabel)
    if(ylab):
        plt.ylabel(ylabel)



def print_run(PR_number,run_number,fact,Nstart,Nend,step):

    data_path = '/ccc/scratch/cont005/gen12622/lavorenf/Mercury_SaeInit/PR'+str(PR_number)+'/run'+str(run_number)+'/data/'
    x   = np.linspace(0,10.,256*fact)
    y   = np.linspace(0,8.,208*fact)
    z   = np.linspace(0,8.,208*fact)
    dt  = 0.5
    ypl = 0.5*max(y)
    zpl = 0.5*max(z)-0.2

    for Ncycle in range(Nstart,Nend,step):

        ### 1st figure Ui ####
        fig1 = plt.figure()
        plt.title('Mercury SHOTS-Init Ion Velocities t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Jx_xz = unpack(data_path,Ncycle,len(y)/2,None,'Jxi')
        Jy_xz = unpack(data_path,Ncycle,len(y)/2,None,'Jyi')
        Jz_xz = unpack(data_path,Ncycle,len(y)/2,None,'Jzi')
        Jx_xy = unpack(data_path,Ncycle,None,len(z)/2,'Jxi')
        Jy_xy = unpack(data_path,Ncycle,None,len(z)/2,'Jyi')
        Jz_xy = unpack(data_path,Ncycle,None,len(z)/2,'Jzi')
        add_plot(fig1,231,Jx_xz,x,z,zpl,'x[di]','z[di]','Jxi','jet',-0.003,0.003,ylab=True,cmap=True)
        add_plot(fig1,232,Jy_xz,x,z,zpl,'x[di]','z[di]','Jyi','jet',-0.003,0.003)
        add_plot(fig1,233,Jz_xz,x,z,zpl,'x[di]','z[di]','Jzi','jet',-0.003,0.003)
        add_plot(fig1,234,Jx_xy,x,y,ypl,'x[di]','y[di]','Jxi','jet',-0.003,0.003,xlab=True,ylab=True)
        add_plot(fig1,235,Jy_xy,x,y,ypl,'x[di]','y[di]','Jyi','jet',-0.003,0.003,xlab=True)
        add_plot(fig1,236,Jz_xy,x,y,ypl,'x[di]','y[di]','Jzi','jet',-0.003,0.003,xlab=True)
        plt.savefig('images/Mercury_SaeInit_Fig1-1_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del Jx_xz, Jy_xz, Jz_xz, Jx_xy, Jy_xy, Jz_xy

        ### 1st figure Ue ####
        fig1 = plt.figure()
        plt.title('Mercury SHOTS-Init Electron Velocities t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Jx_xz = unpack(data_path,Ncycle,len(y)/2,None,'Jxe')
        Jy_xz = unpack(data_path,Ncycle,len(y)/2,None,'Jye')
        Jz_xz = unpack(data_path,Ncycle,len(y)/2,None,'Jze')
        Jx_xy = unpack(data_path,Ncycle,None,len(z)/2,'Jxe')
        Jy_xy = unpack(data_path,Ncycle,None,len(z)/2,'Jye')
        Jz_xy = unpack(data_path,Ncycle,None,len(z)/2,'Jze')
        add_plot(fig1,231,Jx_xz,x,z,zpl,'x[di]','z[di]','Jxe','jet',-0.003,0.003,ylab=True,cmap=True)
        add_plot(fig1,232,Jy_xz,x,z,zpl,'x[di]','z[di]','Jye','jet',-0.003,0.003)
        add_plot(fig1,233,Jz_xz,x,z,zpl,'x[di]','z[di]','Jze','jet',-0.003,0.003)
        add_plot(fig1,234,Jx_xy,x,y,ypl,'x[di]','y[di]','Jxe','jet',-0.003,0.003,xlab=True,ylab=True)
        add_plot(fig1,235,Jy_xy,x,y,ypl,'x[di]','y[di]','Jye','jet',-0.003,0.003,xlab=True)
        add_plot(fig1,236,Jz_xy,x,y,ypl,'x[di]','y[di]','Jze','jet',-0.003,0.003,xlab=True)
        plt.savefig('images/Mercury_SaeInit_Fig1-2_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del Jx_xz, Jy_xz, Jz_xz, Jx_xy, Jy_xy, Jz_xy

        ### 2nd figure n ####
        fig2 = plt.figure()
        plt.title('Mercury SHOTS-Init Densities t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        ni_xz = unpack(data_path,Ncycle,len(y)/2,None,'ni')
        ne_xz = unpack(data_path,Ncycle,len(y)/2,None,'ne')
        ni_xy = unpack(data_path,Ncycle,None,len(z)/2,'ni')
        ne_xy = unpack(data_path,Ncycle,None,len(z)/2,'ne')
        add_plot(fig2,221,ni_xz,x,z,zpl,'x[di]','z[di]','ni','seismic',0,3,ylab=True,cmap=True)
        add_plot(fig2,222,ne_xz,x,z,zpl,'x[di]','z[di]','ne','seismic',0,3)
        add_plot(fig2,223,ni_xy,x,y,ypl,'x[di]','y[di]','ni','seismic',0,3,xlab=True,ylab=True)
        add_plot(fig2,224,ne_xy,x,y,ypl,'x[di]','y[di]','ne','seismic',0,3,xlab=True)
        plt.savefig('images/Mercury_SaeInit_Fig2_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()

        ### 3rd figure T ####
        fig3 = plt.figure()
        plt.title('Mercury SHOTS-Init Temperatures t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Ti_xz = unpack(data_path,Ncycle,len(y)/2,None,'Pi')/gf(ni_xz,fact*2)
        Te_xz = -unpack(data_path,Ncycle,len(y)/2,None,'Pe')/gf(ne_xz,fact*2)
        Ti_xy = unpack(data_path,Ncycle,None,len(z)/2,'Pi')/gf(ni_xy,fact*2)
        Te_xy = -unpack(data_path,Ncycle,None,len(z)/2,'Pe')/gf(ne_xy,fact*2)
        add_plot(fig3,221,Ti_xz,x,z,zpl,'x[di]','z[di]','Ti','copper',0,5e-5,ylab=True,cmap=True)
        add_plot(fig3,222,Te_xz,x,z,zpl,'x[di]','z[di]','Te','copper',0,5e-4)
        add_plot(fig3,223,Ti_xy,x,y,ypl,'x[di]','y[di]','Ti','copper',0,5e-5,xlab=True,ylab=True)
        add_plot(fig3,224,Te_xy,x,y,ypl,'x[di]','y[di]','Te','copper',0,5e-4,xlab=True)
        plt.savefig('images/Mercury_SaeInit_Fig3_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del ni_xz, ne_xz, ni_xy, ne_xy, Ti_xz, Te_xz, Ti_xy, Te_xy

        ### 4th figure E ####
        fig4 = plt.figure()
        plt.title('Mercury SHOTS-Init Electric field t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Ex_xz = unpack(data_path,Ncycle,len(y)/2,None,'Ex')
        Ey_xz = unpack(data_path,Ncycle,len(y)/2,None,'Ey')
        Ez_xz = unpack(data_path,Ncycle,len(y)/2,None,'Ez')
        Ex_xy = unpack(data_path,Ncycle,None,len(z)/2,'Ex')
        Ey_xy = unpack(data_path,Ncycle,None,len(z)/2,'Ey')
        Ez_xy = unpack(data_path,Ncycle,None,len(z)/2,'Ez')
        add_plot(fig4,231,Ex_xz,x,z,zpl,'x[di]','z[di]','Ex','jet',-0.001,0.001,ylab=True,cmap=True)
        add_plot(fig4,232,Ey_xz,x,z,zpl,'x[di]','z[di]','Ey','jet',-0.001,0.001)
        add_plot(fig4,233,Ez_xz,x,z,zpl,'x[di]','z[di]','Ez','jet',-0.001,0.001)
        add_plot(fig4,234,Ex_xy,x,y,ypl,'x[di]','y[di]','Ex','jet',-0.001,0.001,xlab=True,ylab=True)
        add_plot(fig4,235,Ey_xy,x,y,ypl,'x[di]','y[di]','Ey','jet',-0.001,0.001,xlab=True)
        add_plot(fig4,236,Ez_xy,x,y,ypl,'x[di]','y[di]','Ez','jet',-0.001,0.001,xlab=True)
        plt.savefig('images/Mercury_SaeInit_Fig4_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del Ex_xz, Ey_xz, Ez_xz, Ex_xy, Ey_xy, Ez_xy

        ### 5th figure B ####
        fig5 = plt.figure()
        plt.title('Mercury SHOTS-Init Magnetic field t=%d wpi'%(Ncycle*dt))
        plt.axis('off')
        Bx_xz = unpack(data_path,Ncycle,len(y)/2,None,'Bx')
        By_xz = unpack(data_path,Ncycle,len(y)/2,None,'By')
        Bz_xz = unpack(data_path,Ncycle,len(y)/2,None,'Bz')
        Bx_xy = unpack(data_path,Ncycle,None,len(z)/2,'Bx')
        By_xy = unpack(data_path,Ncycle,None,len(z)/2,'By')
        Bz_xy = unpack(data_path,Ncycle,None,len(z)/2,'Bz')
        add_plot(fig5,231,Bx_xz,x,z,zpl,'x[di]','z[di]','Bx','jet',-0.01,0.01,ylab=True,cmap=True)
        add_plot(fig5,232,By_xz,x,z,zpl,'x[di]','z[di]','By','jet',-0.01,0.01)
        add_plot(fig5,233,Bz_xz,x,z,zpl,'x[di]','z[di]','Bz','jet',-0.01,0.01)
        add_plot(fig5,234,Bx_xy,x,y,ypl,'x[di]','y[di]','Bx','jet',-0.01,0.01,xlab=True,ylab=True)
        add_plot(fig5,235,By_xy,x,y,ypl,'x[di]','y[di]','By','jet',-0.01,0.01,xlab=True)
        add_plot(fig5,236,Bz_xy,x,y,ypl,'x[di]','y[di]','Bz','jet',-0.01,0.01,xlab=True)
        plt.savefig('images/Mercury_SaeInit_Fig5_PR'+str(PR_number)+'_'+str(Ncycle*dt)+'.png',format='png')
        plt.close()
        del Bx_xz, By_xz, Bz_xz, Bx_xy, By_xy, Bz_xy


'''
# PR0 run0;1
for irun in range(0,2):
    print('@@@ PR0 - run%d @@@\n'%irun)
    print_run('0',irun,1,3800*irun,3700*(irun+1),100)
'''

# PR1 run0;1
for irun in range(0,2):
    print('@@@ PR1 - run%d @@@\n'%irun)
    print_run('1',irun,2,400+5000*irun,400+5000*(irun+1),100)
                                                                                         
