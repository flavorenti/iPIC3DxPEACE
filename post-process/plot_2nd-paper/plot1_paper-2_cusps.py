import matplotlib
#matplotlib.use('agg')
import fibo as fb
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
import sys
from matplotlib import colors
from matplotlib.patches import Wedge
from scipy.ndimage.filters import gaussian_filter as gf
from scipy.interpolate import RegularGridInterpolator

# define path to simu and cycle to plot
data_address = sys.argv[1]
cycle = 6500
radius = 5.501
Nbinx = 100
Nbiny = 50
ds = 0.15 #di
Nmax = 1500

# miscellaneous
from_VTK = fb.from_VTK(data_address)
from_VTK.get_meta(silent=False)
x = np.linspace(0,from_VTK.meta['xl'],from_VTK.meta['nx'])
y = np.linspace(0,from_VTK.meta['yl'],from_VTK.meta['ny'])
z = np.linspace(0,from_VTK.meta['zl'],from_VTK.meta['nz'])
xc = from_VTK.meta['xc']
yc = from_VTK.meta['yc']
zc = from_VTK.meta['zc']+from_VTK.meta['Doff']
dx = from_VTK.meta['dx'] 
dy = from_VTK.meta['dy'] 
dz = from_VTK.meta['dz']
# load fields in 3D
print('LOAD')
Bx,By,Bz = from_VTK.get_vect(from_VTK.meta['name']+'_B_%i'%cycle,cycle,fibo_obj=None,tar_var='B')

# NEW VERSION = interpolation
fx  = RegularGridInterpolator((x,y,z), Bx)
fy  = RegularGridInterpolator((x,y,z), By)
fz  = RegularGridInterpolator((x,y,z), Bz)

# main routine of the code
def integrator(r,ds):

    # void fields to be filled
    Bline = np.ones((Nbinx,Nbiny,Nmax,2))
    traj = np.ones((3,Nbinx,Nbiny,Nmax,2))
    ibreak = np.ones((Nbinx,Nbiny,2),dtype=int)*Nmax

    # loop on points (along phi)
    for ipx in range(0,Nbinx):
        # loop latitude (theta)
        for ipy in range(0,Nbiny):

            phi=float(ipx)/float(Nbinx)*2.*np.pi
            print('point phi=%.3f'%phi)
            theta=float(ipy)/float(Nbiny)*np.pi
            print('point theta=%.3f'%theta)
            
            rho = radius*np.sin(theta)
            x00 = xc+rho*np.cos(phi)
            y00 = yc+rho*np.sin(phi)
            z00 = zc+radius*np.cos(theta)

            # two loops forward and backward
            for il in range(2):
                # integrator field-line
                for n in range(0,Nmax):
                
                    if n==0:
                        x0=x00
                        y0=y00
                        z0=z00

                    point = np.array([[x0,y0,z0],[0,0,0]])

                    Bxval = fx(point)[0]
                    Byval = fy(point)[0]
                    Bzval = fz(point)[0]
                    Bval = np.sqrt(Bxval**2+Byval**2+Bzval**2)

                    Bline[ipx,ipy,n,il]=Bval
                    traj[0][ipx,ipy,n,il]=x0
                    traj[1][ipx,ipy,n,il]=y0
                    traj[2][ipx,ipy,n,il]=z0

                    x0 += Bxval/Bval*ds 
                    y0 += Byval/Bval*ds 
                    z0 += Bzval/Bval*ds

                    if np.sqrt((x0-xc)**2+(y0-yc)**2+(z0-zc)**2)<from_VTK.meta['R']:
                        ibreak[ipx,ipy,il]=n
                        break

                    elif (x0>np.max(x) or x0<0 or y0>np.max(y) or y0<0 or z0>np.max(z) or z0<0):
                        ibreak[ipx,ipy,il]=-n
                        break
                ds=-ds

    return Bline, traj, ibreak



Bline,traj,ibreak = integrator(radius,ds)

# print in txt file
out = open('output_cusps_PR1-minusBz_%i.txt'%cycle,'w')
out.write('#phi\t theta\t comment\n')

# plot data
fig = plt.figure(1)
ax = Axes3D(fig)
for ipx in range(Nbinx):
    for ipy in range(Nbiny):
        phi=float(ipx)/float(Nbinx)*2.*np.pi
        theta=float(ipy)/float(Nbiny)*np.pi
        rho = radius*np.sin(theta)
        x00 = xc+rho*np.cos(phi)
        y00 = yc+rho*np.sin(phi)
        z00 = zc+radius*np.cos(theta)

        if ibreak[ipx,ipy,0]<0:
            ax.scatter(x00,y00,z00,color='orange',marker='o')
            out.write('%.4f\t %.4f\t forward-outbox\n'%(phi,theta))
        elif ibreak[ipx,ipy,0]==Nmax:
            ax.scatter(x00,y00,z00,color='red',marker='o')
            out.write('%.4f\t %.4f\t forward-eol\n'%(phi,theta))
        elif ibreak[ipx,ipy,1]<0:
            ax.scatter(x00,y00,z00,color='purple',marker='o')
            out.write('%.4f\t %.4f\t backward-outbox\n'%(phi,theta))
        elif ibreak[ipx,ipy,1]==Nmax:
            ax.scatter(x00,y00,z00,color='blue',marker='o')
            out.write('%.4f\t %.4f\t backward-eol\n'%(phi,theta))

            #ax.plot(traj[0][ipx,ipy,:abs(ibreak[ipx,ipy,0]),0],traj[1][ipx,ipy,:abs(ibreak[ipx,ipy,0]),0],traj[2][ipx,ipy,:abs(ibreak[ipx,ipy,0]),0],color='red')
            #ax.plot(traj[0][ipx,ipy,:abs(ibreak[ipx,ipy,1]),1],traj[1][ipx,ipy,:abs(ibreak[ipx,ipy,1]),1],traj[2][ipx,ipy,:abs(ibreak[ipx,ipy,1]),1],color='blue')

ax.set_xlabel('x[di]',fontsize=16)
ax.set_ylabel('y[di]',fontsize=16)
ax.set_zlabel('z[di]',fontsize=16)
plt.xlim(xc-10,xc+10)
plt.ylim(yc-10,yc+10)
ax.set_zlim(zc-10,zc+10)
plt.savefig('./images/plot1_paper2_3D-Blines_PR1-minusBz_%i.png'%cycle,format='png')
plt.close()
out.close()

