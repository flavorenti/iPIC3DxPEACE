
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
import vtk
from vtk.util import numpy_support as VN

class fibo_print: 

  #------------------------------------------------------------ 
  #--------------------routines-for-printing--------------------
  #------------------------------------------------------------  
  def print_vtk_scal(self,  
      tar_var,
      seg,
      address,
      out_name,
      double_y = False, 
      silent   = True):
    """ 
    Prints .vtk file of your tar_var
    
    Parameters :
      - tar_var           [fibo.data] target variable of the procedure
      - seg               [str] cycle of the simulation
      - address           [str] address of printing
      - out_name          [str] name for the printed variable (don't include '.vtk')
      - double_y = False  [bool] if you want to print twice the box (two boxes close in y)
      - silent = True     [bool] print status at the end?
    
    """  
    #determine the coefficients
    nx,ny,nz = self.meta['nnn'] 
    dx,dy,dz = self.meta['ddd'] 

    # add seg to the name of tar_var (fibo) 
    tar_var = tar_var+'%.8i'%int(seg)

    # in case of periodic in y
    if double_y : ny = ny*2

    # convert numpy array to vtk
    vtk_data = VN.numpy_to_vtk(self.get_data(tar_var).ravel(), deep=True, array_type=vtk.VTK_FLOAT)

    # create vtk variable
    vtk_void = vtk.vtkStructuredPoints()
    vtk_void.SetName(tar_var.split('0')[0])
    vtk_void.SetDimensions(nx,ny,nz)
    #vtk_void.SetOrigin(0.,0.,0.)
    vtk_void.SetSpacing(dx,dy,dz)
    vtk_void.AllocateScalars(vtk.VTK_FLOAT, 1)

    # fill vtk variable
    vtk_void.GetPointData().SetScalars(vtk_data)

    # print vtk to output file
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(os.path.join(address,out_name+'.vtk'))
    writer.SetInputData(vtk_void)
    writer.SetFileTypeToBinary()
    writer.Update()
    writer.Write()

    if not silent:
      print('print_vtk_scal> grid dimensions         :  ', self.meta['nnn'])
      print('print_vtk_scal> grid size               :  ', self.meta['lll'])
      print('print_vtk_scal> grid spacing            :  ', self.meta['ddd'])
      print('print_vtk_scal> created file '+out_name+' from fibo.obj '+tar_var)

  #------------------------------------------------------------
  def print_vtk_vect(self,
      tar_var_x,
      tar_var_y,
      tar_var_z,
      seg,
      address,
      out_name,
      double_y = False,  
      silent   = False):
    """ 
    Prints .vtk file of your tar_var_x, tar_var_y, tar_var_z
    
    Parameters :
      - tar_var_x         [fibo.data] x target variable of the procedure
      - tar_var_y         [fibo.data] y target variable of the procedure
      - tar_var_z         [fibo.data] z target variable of the procedure
      - address           [address] address of printing
      - out_name          [str] name for the printed variable (excluding '.vtk')
      - double_y = False  [bool] if you want to print twice the box (two boxes close in y)
      - silent = True     [bool] print status at the end?
 
    """  

    #determine the coefficients
    nx,ny,nz = self.meta['nnn'] #np.shape(self.get_data(tar_var_x))
    dx,dy,dz = self.meta['ddd'] 

    # add seg to the name of tar_var (fibo)
    tar_var_x = tar_var_x+'%.8i'%int(seg)
    tar_var_y = tar_var_y+'%.8i'%int(seg)
    tar_var_z = tar_var_z+'%.8i'%int(seg)

    # in case of periodic in y
    if double_y : ny = ny*2

    # convert numpy array to vtk
    arr_numpy = np.array([self.get_data(tar_var_x),self.get_data(tar_var_y),self.get_data(tar_var_z)]).transpose()
    vtk_data = VN.numpy_to_vtk(arr_numpy.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_data.SetNumberOfComponents(3)

    # create vtk variable
    vtk_void = vtk.vtkStructuredPoints()
    vtk_data.SetName(tar_var_x.split('_x')[0])
    vtk_void.SetDimensions(nx,ny,nz)
    #vtk_void.SetOrigin(0.,0.,0.)
    vtk_void.SetSpacing(dx,dy,dz)
    vtk_void.AllocateScalars(vtk.VTK_FLOAT, 3)

    # fill vtk variable
    vtk_void.GetPointData().SetScalars(vtk_data)

    # print vtk to output file
    writer = vtk.vtkStructuredPointsWriter()
    writer.SetFileName(os.path.join(address,out_name+'.vtk'))
    writer.SetInputData(vtk_void)
    writer.SetFileTypeToBinary()
    writer.Update()
    writer.Write()

    if not silent:    
      print('print_vtk_vect> grid dimensions         :  ', self.meta['nnn'])
      print('print_vtk_vect> grid size               :  ', self.meta['lll'])
      print('print_vtk_vect> grid spacing            :  ', self.meta['ddd'])
      print('print_vtk_vect> created file '+out_name+' from fibo.obj '+tar_var_x+','+tar_var_y+','+tar_var_z)
