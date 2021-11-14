
import collections 
import cv2
import numpy as np
import os
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndm
import itertools as itt
import struct
import re
import time

class fibo_print: 

  #------------------------------------------------------------ 
  #--------------------routines-for-printing--------------------
  #------------------------------------------------------------  
  def print_vtk_scal(self,  
      address,
      tar_name,
      tar_var,
      digits = '%.9f',
      double_y = False, 
      aka = None):
    """ 
    Prints .vtk file of your tar_var
    
    Parameters :
      - address           [address] address of printing
      - tar_name          [str] name for the printed variable (inside the file)
      - tar_var           [fibo.data] target variable of the procedure
      - digits = '%.9f'   [format str] format you will use for printing
      - double_y = False  [bool] if you want to print twice the box (two boxes close in y)
      - aka = None        [None or str] if not None you give the full name of the .vtk file
    
    """  

    #determine the coefficients
    nx,ny,nz = self.meta['nnn'] 
    dx,dy,dz = self.meta['ddd'] 

    #clean the tar_name from possible characters that mess up with paraview's calculator
    tar_name = re.sub('[()!@#$.,*^]', '', tar_name)
    if tar_name[0].isdigit() : tar_name = '_'+tar_name
    #this is the way you create names of the printed files, if you did not mess up with my standard ...
    if (aka == None) :
      filen = os.path.join(address,self.fibo_name+'_'+tar_name+'.vtk')
    else : filen = os.path.join(address,aka+'.vtk')

    if double_y : ny = ny*2

    wf = open(filen, 'w')
    wf.write('# vtk DataFile Version 1.0 \n')
    wf.write(tar_var+' from '+self.fibo_name+'\n')
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
          to_write = digits %self.get_data(tar_var)[ix,iy,iz]
          wf.write(to_write+'\n')
          #old way: wf.write(str(vals[2,ix,iy,iz,it])+' ')
      
      if double_y :

        #print the second time in y (to keep the res spacing you can't start from zero!)
        for iy in range(0, ny):
          for ix in range(0, nx):
            to_write = digits %self.get_data(tar_var)[ix,iy,iz]
            wf.write(to_write+'\n')
            #old way: wf.write(str(vals[2,ix,iy,iz,it])+' ')
        #wf.write('\n')

    wf.close()
    print('done with the print!')

  #------------------------------------------------------------
  def print_vtk_vect(self,
      address,
      tar_name,
      tar_var_x,
      tar_var_y,
      tar_var_z,
      digits_x = '%.7f',
      digits_y = '%.7f',
      digits_z = '%.7f',
      double_y = False,  
      aka = None):
    """ 
    Prints .vtk file of your tar_var_x, tar_vr_y, tar_var_z
    
    Parameters :
      - address           [address] address of printing
      - tar_name          [str] name for the printed variable (inside the file)
      - tar_var_x         [fibo.data] x target variable of the procedure
      - tar_var_y         [fibo.data] y target variable of the procedure
      - tar_var_z         [fibo.data] z target variable of the procedure
      - digits_x = '%.7f' [format str] format you will use for printing tar_var_x
      - digits_y = '%.7f' [format str] format you will use for printing tar_var_y
      - digits_z = '%.7f' [format str] format you will use for printing tar_var_z
      - double_y = False  [bool] if you want to print twice the box (two boxes close in y)
      - aka = None        [None OR str] you can give the full name of the .vtk file
    
    """  

    #determine the coefficients
    nx,ny,nz = self.meta['nnn'] #np.shape(self.get_data(tar_var_x))
    dx,dy,dz = self.meta['ddd'] 

    #clean the tar_name from possible characters that mess up with paraview's calculator
    tar_name = re.sub('[()!@#$.,*^]', '', tar_name)
    if tar_name[0].isdigit() : tar_name = '_'+tar_name
    #this is the way you create names of the printed files, if you did not mess up with my standard ... (and here I hope you really have not)
    if (aka == None) :
      filen = os.path.join(address,self.fibo_name+'_'+tar_name+'.vtk')
    else : filen = os.path.join(address,aka+'.vtk')

    if double_y : ny = ny*2

    wf = open(filen, 'w')
    wf.write('# vtk DataFile Version 1.0 \n')
    wf.write('('+tar_var_x+','+tar_var_y+','+tar_var_z+') from '+self.fibo_name+'\n')
    wf.write('ASCII'+'\n')
    wf.write('DATASET STRUCTURED_POINTS'+'\n')
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
          to_write_x = digits_x %self.get_data(tar_var_x)[ix,iy,iz]
          to_write_y = digits_y %self.get_data(tar_var_y)[ix,iy,iz]
          to_write_z = digits_z %self.get_data(tar_var_z)[ix,iy,iz]
          wf.write(to_write_x+' \t'+to_write_y+' \t'+to_write_z+'\n')
          #old way: wf.write(str(vals_MCA[2,ix,iy,iz,it])+' ')

      if double_y :

        #print the second time in y (to keep the res spacing you can't start from zero!)
        for iy in range(0, ny):
          for ix in range(0, nx):
            to_write_x = digits_x %self.get_data(tar_var_x)[ix,iy,iz]
            to_write_y = digits_y %self.get_data(tar_var_y)[ix,iy,iz]
            to_write_z = digits_z %self.get_data(tar_var_z)[ix,iy,iz]
            wf.write(to_write_x+' \t'+to_write_y+' \t'+to_write_z+'\n')
            #old way: wf.write(str(vals_MCA[2,ix,iy,iz,it])+' ')
        #wf.write('\n')

    wf.close()
    print('done with the print!')
