import numpy as np
import os
import vtk
from vtk.util import numpy_support as VN
import time

# open vtk scalar file and prepare it for conversion to numpy array
def open_vtk_scal(file_path):

    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName(file_path)
    reader.ReadAllScalarsOn()
    reader.Update()

    return reader.GetOutput()


# open vtk vector file and prepare it for conversion to numpy array
def open_vtk_vect(file_path):

    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName(file_path)
    reader.ReadAllVectorsOn()
    reader.Update()

    return reader.GetOutput()


# conversion to numpy array
def convert_vtk_scal(data_file):

	data = VN.vtk_to_numpy(data_file.GetPointData().GetScalars())
	dims = data_file.GetDimensions()
	data = data.reshape(dims[2],dims[1],dims[0]) # need to regive the right shape to np.array because VN givs it flatten
	data = data.transpose(2,1,0)

	return np.array(data)


def convert_vtk_vect(data_file,     # vtk file prepared above (3 components)
                           string):	# string with the name of the array we want to unpack, e.g. B,E,Ji,Je

    data = VN.vtk_to_numpy(data_file.GetPointData().GetArray(string))
    dims = data_file.GetDimensions()

    data_x = data[:,0].reshape(dims[2],dims[1],dims[0]) 
    data_y = data[:,1].reshape(dims[2],dims[1],dims[0]) 
    data_z = data[:,2].reshape(dims[2],dims[1],dims[0]) 

    data_x = data_x.transpose(2,1,0)
    data_y = data_y.transpose(2,1,0)
    data_z = data_z.transpose(2,1,0)

    return np.array([data_x, data_y, data_z])



def print_vtk(field,                    # scalar field to turn into vtk file
              output_path,              # address path for the output
              tar_name,                 # name for the printed variable into the VTK file
              dl,                       # data interval (typically (dx,dy,dz) for a field described in space)
              digits=None):             # z format you will use for printing

    if len(dl)!=3:
        print('Either your field is not 3D, either your dl is wrong')
        return 0

    shape = np.shape(field)

    if (len(shape)<3) or (len(shape)>4):
        print('Your input field is incorrect.')
        print('Your field must scale either as a scalar field (shape=[nx,ny,nz]), either as a vector field (shape=[3,nx,ny,nz]).')
        return 0

    elif len(shape)==3:
        if digits is not None:
            print_vtk_scal(field,output_path,tar_name,dl,digits=digits)
        else:
            print_vtk_scal(field,output_path,tar_name,dl)

    elif len(shape)==4:
        if shape[0]!=3:
            print('The field format is wrong. Field must be an array of shape=[3,nx,ny,nz].')
        else:
            if digits is not None:
                print_vtk_vect(field,output_path,tar_name,dl,digits_x=digits,digits_y=digits,digits_z=digits)
            else:
                print_vtk_vect(field,output_path,tar_name,dl)

    else:
        print('If you have this message, it means you found a totally unexpected bug in my routine.')
        print('Congratulation!')
        print("Can you please tell me what you've done?")


    return 0


def print_vtk_scal(field,               # scalar field to turn into vtk file
                   output_path,         # address output file 
                   tar_name,            # name for the printed variable into the VTK file
                   dl,                  # data interval (typically (dx,dy,dz) for a field described in space)
                   digits = '%.9f'):    # z format you will use for printing

    #determine the coefficients
    nx,ny,nz = np.shape(field)
    dx = dl[0]
    dy = dl[1]
    dz = dl[2]
    
    wf = open(output_path, 'w')
    wf.write('# vtk DataFile Version 1.0 \n')
    wf.write(tar_name+'\n')
    wf.write('ASCII'+'\n')
    wf.write('DATASET STRUCTURED_POINTS'+'\n')
    wf.write('DIMENSIONS'+' '+str(nx)+' '+str(ny)+' '+str(nz)+'\n')
    wf.write('ORIGIN  0.0  0.0  0.0'+'\n')
    wf.write('SPACING'+' '+str(dx)+' '+str(dy)+' '+str(dz)+'\n')
    wf.write('                '+'\n')
    wf.write('POINT_DATA'+' '+str(nx*ny*nz)+'\n')
    wf.write('SCALARS    '+tar_name+'    float  1'+'\n')
    wf.write('LOOKUP_TABLE default'+'\n')
    
    for iz in range(0, nz):
        #print the first time in y
        for iy in range(0, ny):
            for ix in range(0, nx):
                to_write = digits %field[ix,iy,iz]
                wf.write(to_write+'\n')

    wf.close()
    print('done with the print!'+output_path)
    
    return 0


def print_vtk_vect(field,               # vector field to turn into vtk file
                   output_path,         # address output file 
                   tar_name,            # name for the printed variable into the VTK file
                   dl,                  # data interval (typically (dx,dy,dz) for a field described in space)
                   digits_x = '%.7f',	# x format you will use for printing
                   digits_y = '%.7f',	# y format you will use for printing
                   digits_z = '%.7f'):  # z format you will use for printing
    
    #determine the coefficients
    nx,ny,nz = np.shape(field[0,...])
    dx = dl[0]
    dy = dl[1]
    dz = dl[2]
                
    wf = open(output_path, 'w')
    wf.write('# vtk DataFile Version 1.0 \n')
    wf.write(tar_name+'\n')
    wf.write('ASCII'+'\n')
    wf.write('DATASET STRUCTURED_POINTS \n')
    wf.write('DIMENSIONS'+' '+str(nx)+' '+str(ny)+' '+str(nz)+'\n')
    wf.write('ORIGIN  0.0  0.0  0.0'+'\n')
    wf.write('SPACING'+' '+str(dx)+' '+str(dy)+' '+str(dz)+'\n')
    wf.write('                '+'\n')
    wf.write('POINT_DATA'+' '+str(nx*ny*nz)+'\n')
    wf.write('VECTORS    '+tar_name+'    float'+'\n')
                
    for iz in range(0, nz):
        #print the first time in y
        for iy in range(0, ny):
            for ix in range(0, nx):
                to_write_x = digits_x %field[0,ix,iy,iz]
                to_write_y = digits_y %field[1,ix,iy,iz]
                to_write_z = digits_z %field[2,ix,iy,iz]
                wf.write(to_write_x+' \t'+to_write_y+' \t'+to_write_z+'\n')

    wf.close()
    print('done with the print!'+output_path)

    return 0







#---------------------------------------------------------------------------------------
class from_vtk (object):
    
    def __init__(self,
                 address):		#where are the data we are referring to
        
        self.address = address
            
        self.nx = 0
        self.ny = 0
        self.nz = 0
    
        self.xl = 0
        self.yl = 0
        self.zl = 0
        
    #------------------------------------------------------------
    def get_scal(self,
                 name):					# The vtk file name
            
            
        #create data vector, fill it!
        data_file = open(os.path.normpath(self.address+'/'+name),'r')
            
        data_file.readline()
        print('Begin reading of '+ data_file.readline())
        data_file.readline()
        data_file.readline()
        self.nx, self.ny, self.nz = list(map(int, data_file.readline().split()[1:4]))
        data_file.readline()
        dx, dy, dz = list(map(float, data_file.readline().split()[1:4]))
        data_file.readline()
        data_file.readline()	#NB here you have the nx*ny*nz preduct
        data_file.readline()
        data_file.readline()
        
        self.xl = self.nx*dx
        self.yl = self.ny*dy
        self.zl = self.nz*dz

        
        scal = np.zeros([self.nx,self.ny,self.nz])
            
        for iz in range(self.nz):
            for iy in range(self.ny):
                for ix in range(self.nx):
                    scal[ix,iy,iz] = float(data_file.readline().split())
    
        data_file.close()
        

        return scal
                
        #if not silent: print('done with the reading!')
                
    #------------------------------------------------------------
    def get_vect(self,
                     name):				# The vtk file name
                    
                    
        #create data vector, fill it!
        data_file = open(os.path.normpath(self.address+'/'+name),'r')
            
        data_file.readline()
        print('Begin reading of '+ data_file.readline())
        data_file.readline()
        data_file.readline()
        self.nx, self.ny, self.nz = list(map(int, data_file.readline().split()[1:4]))
        data_file.readline()
        dx, dy, dz = list(map(float, data_file.readline().split()[1:4]))
        data_file.readline()
        data_file.readline()	#NB here you have the nx*ny*nz preduct
        data_file.readline()
    
            
        self.xl = self.nx*dx
        self.yl = self.ny*dy
        self.zl = self.nz*dz
        
            
        vect = np.empty([3,self.nx,self.ny,self.nz])
            
        for iz in range(self.nz):
            for iy in range(self.ny):
                for ix in range(self.nx):
                    aa, bb, cc = list(map(float, data_file.readline().split()))
                    vect[0,ix,iy,iz] = aa
                    vect[1,ix,iy,iz] = bb
                    vect[2,ix,iy,iz] = cc
                
        data_file.close()

                
        return vect

        #if not silent: print('done with the reading!')
