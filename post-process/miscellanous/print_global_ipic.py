from vectorfield import *
from vtk_routines import *


# where is your data
file_path = '/ccc/scratch/cont005/gen12622/lavorenf/Mercury_SaeInit/PR0/run0/data/'
# what kind of simulation you did
simu = 'Dipole3D'
# what are the box dimensions
x = np.linspace(0,10.,256)
y = np.linspace(0,8.,208)
z = y
# cycles you want to plot (start end step)
start = 1000
end   = 3000
step  = 200


def divergence(vect,resx,resy,resz):
	nx,ny,nz = np.shape(vect[0][...])
	div = 0.5*resx*(np.roll(vect[0][...],-1,axis=0)-np.roll(vect[0][...],+1,axis=0))+0.5*resy*(np.roll(vect[1][...],-1,axis=1)-np.roll(vect[1][...],+1,axis=1))+0.5*resz*(np.roll(vect[2][...],-1,axis=2)-np.roll(vect[2][...],+1,axis=2))
	div[0,:,:]=0.
	div[:,0,:]=0.
	div[:,:,0]=0.
	div[nx-1,:,:]=0.
	div[:,ny-1,:]=0.
	div[:,:,nz-1]=0.
	return div

def rotore(vect,resx,resy,resz):
	nx,ny,nz = np.shape(vect[0][...])
	rot_x = 0.5*resy*(np.roll(vect[2][...],-1,axis=1)-np.roll(vect[2][...],+1,axis=1)) - 0.5*resz*(np.roll(vect[1][...],-1,axis=2)-np.roll(vect[1][...],+1,axis=2)) 
	rot_y =-0.5*resx*(np.roll(vect[2][...],-1,axis=0)-np.roll(vect[2][...],+1,axis=0)) + 0.5*resz*(np.roll(vect[0][...],-1,axis=2)-np.roll(vect[0][...],+1,axis=2)) 
	rot_z = 0.5*resx*(np.roll(vect[1][...],-1,axis=0)-np.roll(vect[1][...],+1,axis=0)) - 0.5*resy*(np.roll(vect[0][...],-1,axis=1)-np.roll(vect[0][...],+1,axis=1)) 
	rot_x[0,:,:]=0.
	rot_x[:,0,:]=0.
	rot_x[:,:,0]=0.
	rot_x[nx-1,:,:]=0.
	rot_x[:,ny-1,:]=0.
	rot_x[:,:,nz-1]=0.
	rot_y[0,:,:]=0.
	rot_y[:,0,:]=0.
	rot_y[:,:,0]=0.
	rot_y[nx-1,:,:]=0.
	rot_y[:,ny-1,:]=0.
	rot_y[:,:,nz-1]=0.
	rot_z[0,:,:]=0.
	rot_z[:,0,:]=0.
	rot_z[:,:,0]=0.
	rot_z[nx-1,:,:]=0.
	rot_z[:,ny-1,:]=0.
	rot_z[:,:,nz-1]=0.
	return np.array([rot_x,rot_y,rot_z])



# C'est parti!!!
cycle=start
for i in range(start,end,step):

	print('** print step %d/%d **'%(cycle,end))

	rhoe_vtk = open_vtk_scal(file_path+simu+'_rhoe0_'+str(cycle)+'.vtk')
	rhoi_vtk = open_vtk_scal(file_path+simu+'_rhoi1_'+str(cycle)+'.vtk')
	Je_vtk = open_vtk_vect(file_path+simu+'_Je_'+str(cycle)+'.vtk')
	Ji_vtk = open_vtk_vect(file_path+simu+'_Ji_'+str(cycle)+'.vtk') 
	B_vtk = open_vtk_vect(file_path+simu+'_B_'+str(cycle)+'.vtk')
	E_vtk = open_vtk_vect(file_path+simu+'_E_'+str(cycle)+'.vtk')

	rhoe = convert_vtk_scal(rhoe_vtk)
	rhoi = convert_vtk_scal(rhoi_vtk)
	Je = convert_vtk_vect(Je_vtk,'Je')
	Ji = convert_vtk_vect(Ji_vtk,'Ji')
	B = convert_vtk_vect(B_vtk,'B')
	E = convert_vtk_vect(E_vtk,'E')
	
	# compute total charge
	TOTrho = rhoi+rhoe
	print_vtk_scal(TOTrho, file_path+simu+'_TOTrho_'+str(cycle)+'.vtk','TOTrho', np.array([x[1],y[1],z[1]]))

	# compute total current
	TOTJ = Ji+Je
	print_vtk_vect(TOTJ, file_path+simu+'_TOTJ_'+str(cycle)+'.vtk', 'TOTJ', np.array([x[1],y[1],z[1]]))

	# compute DIV(B)=0
	divB = divergence(B,1./x[1],1./y[1],1./z[1])
	print_vtk_scal(divB, file_path+simu+'_divB_'+str(cycle)+'.vtk','divB', np.array([x[1],y[1],z[1]]))

	# compute ROT(B)-J=0
	rotB = rotore(B,1./x[1],1./y[1],1./z[1])
	print_vtk_vect(rotB-TOTJ, file_path+simu+'_rotB-J_'+str(cycle)+'.vtk','rotB-J', np.array([x[1],y[1],z[1]]))

	# compute ROT(E)=0
	rotE = rotore(E,1./x[1],1./y[1],1./z[1])
	print_vtk_vect(rotE, file_path+simu+'_rotE_'+str(cycle)+'.vtk','rotE', np.array([x[1],y[1],z[1]]))

	# compute DIV(E)=TOTrho
	divE = divergence(E,1./x[1],1./y[1],1./z[1])
	print_vtk_scal(divE-TOTrho, file_path+simu+'_Gauss_'+str(cycle)+'.vtk', 'Gauss', np.array([x[1],y[1],z[1]]))

	del TOTJ,divB,Je,Ji,rotB,B,divE,rotE,E,rhoe,rhoi,TOTrho

	cycle += step


