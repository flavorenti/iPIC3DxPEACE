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
radius = 5.5*1.1
Npoints=6
ds = 0.02
Nmax = 200

# conversion factors from code units to SI values
conv_B = 20./0.00562 #nT

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
    Bline = np.ones((Npoints,Nmax,2))
    traj = np.ones((3,Npoints,Nmax,2))
    ibreak = np.ones((Npoints,2),dtype=int)*Nmax

    # loop on points (along phi)
    for ip in range(0,Npoints):

        phi=float(ip)/float(Npoints)*2.*np.pi
        print('point phi=%.3f'%phi)
        x00 = xc+r*np.cos(phi)
        y00 = yc+r*np.sin(phi)
        z00 = zc-from_VTK.meta['Doff']

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

                Bline[ip,n,il]=(conv_B*Bval)
                traj[0][ip,n,il]=x0
                traj[1][ip,n,il]=y0
                traj[2][ip,n,il]=z0

                x0 += Bxval/Bval*ds 
                y0 += Byval/Bval*ds 
                z0 += Bzval/Bval*ds

                if np.sqrt((x0-xc)**2+(y0-yc)**2+(z0-zc)**2)<from_VTK.meta['R']:
                    ibreak[ip,il]=n
                    break

            ds=-ds

    return Bline, traj, ibreak



Bline1,traj1,ibreak1 = integrator(radius,ds)
Bline2,traj2,ibreak2 = integrator(radius+0.2,ds)
Bline3,traj3,ibreak3 = integrator(radius+0.4,ds)

'''
fig = plt.figure(1)
ax = Axes3D(fig)
for ip in range(Npoints):
    ax.plot(traj[0][ip,:ibreak[ip,0],0],traj[1][ip,:ibreak[ip,0],0],traj[2][ip,:ibreak[ip,0],0],color='red')
    ax.plot(traj[0][ip,:ibreak[ip,1],1],traj[1][ip,:ibreak[ip,1],1],traj[2][ip,:ibreak[ip,1],1],color='blue')
ax.set_xlabel('x[di]',fontsize=16)
ax.set_ylabel('y[di]',fontsize=16)
ax.set_zlabel('z[di]',fontsize=16)
plt.xlim(xc-10,xc+10)
plt.ylim(yc-10,yc+10)
ax.set_zlim(zc-10,zc+10)
plt.savefig('./images/plot4_paper_3D-Blines_1_%i.png'%cycle,format='png')
plt.close()
'''

pl.figure(figsize=(15,6))
colors = pl.cm.jet(np.linspace(0,1,Npoints))

pl.subplot(131)
line1 = []
line2 = []
pl.title('Radius=1.1R',fontsize=18)
for ip in range(Npoints):
    theta_lim = np.arcsin( np.sqrt( Bline1[ip,0,0]/200. ) )*180./3.14
    l1, =pl.plot(np.arange(ibreak1[ip,0])*ds/5.5,Bline1[ip,:ibreak1[ip,0],0],label=r'$|\theta_0|$<%.0f'%theta_lim,color=colors[ip])
    line1.append(l1)

leg1 = pl.legend(handles=line1,loc='lower left',fontsize=12)
pl.gca().add_artist(leg1)

for ip in range(Npoints):
    l2, = pl.plot(-np.arange(ibreak1[ip,1])*ds/5.5,Bline1[ip,:ibreak1[ip,1],1],label=r'LT=%i'%(ip*24/Npoints),color=colors[ip])
    line2.append(l2)

leg2 = pl.legend(handles=line2,loc='lower right',fontsize=12)

pl.ylabel('|B(s)| [nT]',fontsize=16)
pl.xlabel('s [R]',fontsize=16)
plt.xlim(-0.65,0.65)
plt.ylim(80,250)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

pl.subplot(132)
line1 = []
line2 = []
pl.title('Radius=1.3R',fontsize=18)
for ip in range(Npoints):
    theta_lim = np.arcsin( np.sqrt( Bline2[ip,0,0]/200. ) )*180./3.14
    l1, = pl.plot(np.arange(ibreak2[ip,0])*ds/5.5,Bline2[ip,:ibreak2[ip,0],0],label=r'$|\theta_0|$<%.0f'%theta_lim,color=colors[ip])
    l2, = pl.plot(-np.arange(ibreak2[ip,1])*ds/5.5,Bline2[ip,:ibreak2[ip,1],1],color=colors[ip])
    line1.append(l1)
    line2.append(l2)

pl.xlabel('s [R]',fontsize=16)
pl.legend(handles=line1,loc='lower left',fontsize=12)
plt.xlim(-0.65,0.65)
plt.ylim(80,250)
plt.xticks(fontsize=14)
plt.yticks([])

pl.subplot(133)
line1 = []
line2 = []
pl.title('Radius=1.5R',fontsize=18)
for ip in range(Npoints):
    theta_lim = np.arcsin( np.sqrt( Bline3[ip,0,0]/200. ) )*180./3.14
    l1, = pl.plot(np.arange(ibreak3[ip,0])*ds/5.5,Bline3[ip,:ibreak3[ip,0],0],label=r'$|\theta_0|$<%.0f'%theta_lim,color=colors[ip])
    l2, = pl.plot(-np.arange(ibreak3[ip,1])*ds/5.5,Bline3[ip,:ibreak3[ip,1],1],color=colors[ip])
    line1.append(l1)
    line2.append(l2)

pl.xlabel('s [R]',fontsize=16)
pl.legend(handles=line1,loc='lower left',fontsize=12)
plt.xlim(-0.65,0.65)
plt.ylim(80,250)
plt.xticks(fontsize=14)
plt.yticks([])

'''
pl.subplot(312)
gradB = (np.roll(Bline,-1,axis=1)-Bline)/ds
grad2B = (np.roll(gradB,-1,axis=1)-gradB)/ds
for ip in range(Npoints):
    pl.plot(np.arange(ibreak[ip,0])*ds,gradB[ip,:ibreak[ip,0],0],color=colors[ip])
    pl.plot(-np.arange(ibreak[ip,1])*ds,gradB[ip,:ibreak[ip,1],1],color=colors[ip])
pl.ylabel('dB/ds [nT/di]')
pl.xlabel('B-line length [di]')
pl.subplot(313)
for ip in range(Npoints):
    pl.plot(np.arange(ibreak[ip,0])*ds,grad2B[ip,:ibreak[ip,0],0],color=colors[ip])
    pl.plot(-np.arange(ibreak[ip,1])*ds,grad2B[ip,:ibreak[ip,1],1],color=colors[ip])
pl.ylabel(r'$d^2B/ds^2$ [nT/$d_i^2$]')
pl.xlabel('B-line length [di]')
'''
pl.subplots_adjust(wspace=0.08)
pl.savefig('./images/plot4_paper_Blines_6_%i.png'%cycle,format='png')
plt.close()
