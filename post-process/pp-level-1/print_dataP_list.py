import matplotlib
#matplotlib.use('agg')
import fibo as fb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
import sys

data_path = sys.argv[1]
output_path = sys.argv[2]
Plot_patch= True
iplot=1
if sys.argv[3]=='Mariner' :
    file_path = '../trajectory_files/Mariner-flyby-1st_MSO-trajectory1.txt'
elif sys.argv[3]=='Bepi' :
    file_path = '../trajectory_files/Bepi-flyby-1st_MSO-trajectory1.txt'
elif sys.argv[3]=='nose' :
    file_path='../trajectory_files/noseMP_trajectory1.txt'
elif sys.argv[3]=='tail' :
    file_path='../trajectory_files/tailXline_trajectory1.txt'


# return the trajectory inside the simulation box
def trajectory_inside(from_VTK):

    # arrays box simulation
    x = from_VTK.meta['x']
    y = from_VTK.meta['y']
    z = from_VTK.meta['z']

    # open the trajectory file (MSO coord.)
    tt,xx,yy,zz = np.loadtxt(file_path,unpack=True)

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
    
    '''
    output = open('Bepi_box-coord_PR1-newbox.txt','w')
    output.write('# time[min]\t x[di]\t y[di]\t z[di]\n')
    print('OUTPUT START')
    for ix in range(0,len(xx),10):
        output.write('%f\t %f\t %f\t %f\n'%(tt[ix],xx[ix],yy[ix],zz[ix]))
    output.close()
    print('OUTPUT DONE')
    '''

    enter=0
    ip=0
     # loop on the trajectory file downloaded from amda
    for ii in range(len(xx)):
        if((xx[ii]>0)*(yy[ii]>0)*(zz[ii]>0)*(xx[ii]<max(x))*(yy[ii]<max(y))*(zz[ii]<max(z))):
            if(enter==0):
                imin=ii
            enter=1
        elif(enter==1):
            imax=ii
            break

    if ((sys.argv[3]=='Mariner') or (sys.argv[3]=='Bepi')):
        return np.array([xx[imin:imax],yy[imin:imax],zz[imin:imax]])
    else:
        return np.array([xx,yy,zz])

# return and print the good patche in .txt and plot them 
def pcl_patch(from_VTK,traj):

    # arrays box simulation
    x = from_VTK.meta['x']
    y = from_VTK.meta['y']
    z = from_VTK.meta['z']

    Nx_ = int(from_VTK.meta['nx']/from_VTK.meta['XLEN'])
    Ny_ = int(from_VTK.meta['ny']/from_VTK.meta['YLEN'])
    Nz_ = int(from_VTK.meta['nz']/from_VTK.meta['ZLEN'])

    # need some spacing to be sure IMPORTANT!
    dr = np.max( [max(x)/from_VTK.meta['XLEN'],max(y)/from_VTK.meta['YLEN'],max(z)/from_VTK.meta['ZLEN']] )/2.

    # open file for output
    output = open(output_path+'/texts/output_dataP-new_'+sys.argv[3]+'.txt','w')

    # array to plot
    patch_list_x = []
    patch_list_y = []
    patch_list_z = []
    color_list = []

    # loop on the processes
    ip=0
    ic=0
    for ipx in range(0,from_VTK.meta['XLEN']):
        ix1 = ipx*(Nx_-1)
        ix2 = ix1+Nx_
        for ipy in range(0,from_VTK.meta['YLEN']):
            iy1 = ipy*(Ny_-1)
            iy2 = iy1+Ny_
            for ipz in range(0,from_VTK.meta['ZLEN']):
                iz1 = ipz*(Nz_-1)
                iz2 = iz1+Nz_
                
                # look if any of the trajectory points are in this patch
                igoodx = ((traj[0]-dr)>x[ix1])*((traj[0]-dr)<x[ix2]) + ((traj[0]+dr)>x[ix1])*((traj[0]+dr)<x[ix2])
                igoody = ((traj[1]-dr)>y[iy1])*((traj[1]-dr)<y[iy2]) + ((traj[1]+dr)>y[iy1])*((traj[1]+dr)<y[iy2])
                igoodz = ((traj[2]-dr)>z[iz1])*((traj[2]-dr)<z[iz2]) + ((traj[2]+dr)>z[iz1])*((traj[2]+dr)<z[iz2])
                booll  = igoodx*igoody*igoodz

                # write output
                if(booll.any()):
                    output.write('%i\n'%ip)
                    ic=ic+1
                    color_list.append('red')
                else:
                    color_list.append('blue')

                # increase patch index
                ip=ip+1

                # plot a sort of grid
                patch_list_x.append((x[ix1]+x[ix2])/2.)
                patch_list_y.append((y[iy1]+y[iy2])/2.)
                patch_list_z.append((z[iz1]+z[iz2])/2.)

    # 3d plot
    if Plot_patch :
        fig = plt.figure()
        ax = Axes3D(fig)
        plt.title('MPI processes + flyby Ncross/Nproc=%d/%d'%(ic,ip))
        for i in range(len(patch_list_x)):
            if(i%iplot==0):
                if color_list[i]=='red':
                    ax.scatter(patch_list_x[i],patch_list_y[i],patch_list_z[i],marker='s',s=150,color=color_list[i],alpha=0.4)
        
        ax.plot(traj[0],traj[1],traj[2],color='grey')
        ax.set_xlabel('x[di]',fontsize=16)
        ax.set_ylabel('y[di]',fontsize=16)
        ax.set_zlabel('z[di]',fontsize=16)

        plt.xlim(min(x),max(x))
        plt.ylim(min(y),max(y))
        ax.set_zlim(min(z),max(z))
        plt.show()
        plt.savefig(output_path+'/images/3d_plot_patches-new_'+sys.argv[3]+'.png',format='png')
        plt.close()

    # close txt file
    output.close()
                

# main of the code
def main():
    from_VTK = fb.from_VTK(data_path)
    from_VTK.get_meta(silent=False)
    traj=trajectory_inside(from_VTK)
    pcl_patch(from_VTK,traj)

# call main
main()
