import fibo as fb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
import sys

file_path = './'
file_name = 'Mariner-flyby-1st_full-dataset.txt'
data_path = sys.argv[1]
Plot_patch= True
XLEN = 16
YLEN = 16
ZLEN = 16


# return the trajectory inside the simulation box
def trajectory_inside(from_VTK):

    # arrays box simulation
    x = from_VTK.meta['x']
    y = from_VTK.meta['y']
    z = from_VTK.meta['z']

    # open the trajectory file (MSO coord.)
    tt,bbX,bxX,byX,bzX,xx,yy,zz,n1X,n2X,n3X,T1X,T2X,T3X = np.loadtxt(file_path+file_name,unpack=True)

    # distance from the planet center
    d = np.sqrt(xx**2+yy**2+zz**2)
    i_CA = np.where(d==min(d))
    # time c.a. in minutes
    tt = (tt-tt[i_CA])/60.

    # pass from MSO -> box coordinate system
    xx = xx*from_VTK.meta['R']
    yy = yy*from_VTK.meta['R']
    zz = zz*from_VTK.meta['R']
    xx = -xx + from_VTK.meta['xc']
    yy = -yy + from_VTK.meta['yc']
    zz = +zz + from_VTK.meta['zc'] + from_VTK.meta['Doff']
    bxX = -bxX
    byX = -byX
 
    enter=0
    ip=0
     # loop on the trajectory file downloaded from amda
    for ii in range(len(xx)):
        if((xx[ii]>0)*(yy[ii]>0)*(zz[ii]>0)*(xx[ii]<max(x))*(yy[ii]<max(y))*(zz[ii]<max(z))):
            ix = np.where(abs(x-xx[ii])==min(abs(x-xx[ii])))
            iy = np.where(abs(y-yy[ii])==min(abs(y-yy[ii])))
            iz = np.where(abs(z-zz[ii])==min(abs(z-zz[ii])))
            if(enter==0):
                imin=ii
            enter=1
        elif(enter==1):
            imax=ii
            break

    return np.array([xx[imin:imax],yy[imin:imax],zz[imin:imax]])


# return and print the good patche in .txt and plot them 
def pcl_patch(from_VTK,traj):

    # arrays box simulation
    x = from_VTK.meta['x']
    y = from_VTK.meta['y']
    z = from_VTK.meta['z']

    Nx_ = int(from_VTK.meta['nx']/XLEN)
    Ny_ = int(from_VTK.meta['ny']/YLEN)
    Nz_ = int(from_VTK.meta['nz']/ZLEN)

    # open file for output
    output = open('./output_dataP.txt','w')
    output.write('# '+file_path+file_name+'\n# '+data_path+'\n')

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

    # 3d plot
    if Plot_patch :
        fig = plt.figure()
        ax = Axes3D(fig)
        plt.title('MPI processes + flyby Ncross/Nproc=%d/%d'%(ic,ip))
        for i in range(len(patch_list_x)):
            if(i%2==0):
                ax.scatter(patch_list_x[i],patch_list_y[i],patch_list_z[i],marker='o',color=color_list[i],alpha=0.4)

        ax.set_xlabel('x[Rm]',fontsize=16)
        ax.set_ylabel('y[Rm]',fontsize=16)
        ax.set_zlabel('z[Rm]',fontsize=16)

        plt.show()
        #plt.savefig('./images/3d_plot_patches_1.png',format='png')
        plt.close()
                

# main of the code
def main():
    from_VTK = fb.from_VTK(data_path)
    from_VTK.get_meta(silent=False)
    traj=trajectory_inside(from_VTK)
    pcl_patch(from_VTK,traj)

# call main
main()
