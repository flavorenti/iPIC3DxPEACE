
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

class fibo_extract :

  #------------------------------------------------------------ 
  #-----routines-for-basic-extractions-------------------------
  #------------------------------------------------------------
  def extract_line(self,  
      tar_data,  
      line_pnts,  
      line_dims,  
      center,  
      versor):  
    """
    Extracts array values along a specified line
    
    Parameters :
      - tar_data    [fibo.data] target variable
      - line_pnts   [int] length (in points) of the line you ask to extract
      - line_dims   [float] dimension of the line you ask to extract 
      - center      [float,float,float] coordinates of the first and last point 
      - versor      [float,float,float] projections of the versor over which the line is built
    
    Returns :
      - new_data    [fibo.data] new data
      - new_meta    [fibo.meta] new metadata
    
    """

    nx,ny,nz = self.meta['nnn']
    dx,dy,dz = self.meta['ddd']

    # here I determine how many relevant periodic dimensions there are
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    # then double the box in every relevant periodic dimension
    arr = self.get_data(tar_data)
    
    if dim_x : arr = np.tile(arr,(2,1,1))
    if dim_y : arr = np.tile(arr,(1,2,1))
    if dim_z : arr = np.tile(arr,(1,1,2))

    # now let's move the center into the central part of this extended array
    if dim_x : center[0] += (1 - 2*center[0]//nx)*nx
    if dim_y : center[1] += (1 - 2*center[1]//ny)*ny
    if dim_z : center[2] += (1 - 2*center[2]//nz)*nz
    
    # re-scale the in-plane versors by the dimensions
    vec = np.array(versor) * line_dims

    # then, add and subtract the in-plane vectors from the plane center coordinates 
    rx = np.array([center[0]]*2)
    rx[0] += 0.5* ( -vec[0] ) / dx
    rx[1] += 0.5* ( +vec[0] ) / dx
    ry = np.array([center[1]]*2)
    ry[0] += 0.5* ( -vec[1] ) / dy
    ry[1] += 0.5* ( +vec[1] ) / dy
    rz = np.array([center[2]]*2)
    rz[0] += 0.5* ( -vec[2] ) / dz
    rz[1] += 0.5* ( +vec[2] ) / dz
    
    # so that you can finally calculate the new data
    lin_x = np.linspace(rx[0],rx[1],line_pnts,endpoint=True)
    lin_y = np.linspace(ry[0],ry[1],line_pnts,endpoint=True)
    lin_z = np.linspace(rz[0],rz[1],line_pnts,endpoint=True)

    new_data = ndm.map_coordinates(arr, np.vstack((lin_x,lin_y,lin_z)))
    new_data = new_data.reshape((line_pnts,1,1))

    # and the new meta data 
    new_meta = self.meta.copy()
    new_meta['space_dim'] = '1D'
    
    new_meta['nx'] = line_pnts
    new_meta['ny'] = 1
    new_meta['nz'] = 1
    
    new_meta['xl'] = line_dims
    new_meta['yl'] = 2*np.pi
    new_meta['zl'] = 2*np.pi
    
    new_meta['dx'] = new_meta['xl'] / new_meta['nx']
    new_meta['dy'] = new_meta['yl'] / new_meta['ny']
    new_meta['dz'] = new_meta['zl'] / new_meta['nz']
    
    new_meta['nnn'] = (new_meta['nx'],new_meta['ny'],new_meta['nz'])
    new_meta['lll'] = (new_meta['xl'],new_meta['yl'],new_meta['zl'])
    new_meta['ddd'] = (new_meta['dx'],new_meta['dy'],new_meta['dz'])    
    
    del(new_meta['np_col'])
    del(new_meta['np_pla'])
    del(new_meta['np_row'])

    return new_data, new_meta

  #------------------------------------------------------------
  def extract_grid(self,
      tar_data,  
      grid_pnts,
      grid_dims,
      center,
      versor_one,
      versor_two):
    """
    Extracts array values in a specified grid
    
    Parameters :
      - tar_data    [fibo.data] target variable
      - grid_pnts   [int,int] lengths (in points) of the plane you ask to extract
      - grid_dims   [float,float] dimensions of the plane you ask to extract
      - center      [float,float,float] coordinates of the central point 
      - versor_one  [float,float,float] projections of the first in-plane versor
      - versor_two  [float,float,float] projections of the second in-plane versor
    
    Returns :
      - new_data    [fibo.data] new data
      - new_meta    [fibo.meta] new metadata
    
    """

    nx,ny,nz = self.meta['nnn']
    dx,dy,dz = self.meta['ddd']

    # here I determine how many relevant periodic dimensions there are
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    # then double the box in every relevant periodic dimension
    arr = self.get_data(tar_data)
    
    if dim_x : arr = np.tile(arr,(2,1,1))
    if dim_y : arr = np.tile(arr,(1,2,1))
    if dim_z : arr = np.tile(arr,(1,1,2))

    # now let's move the center into the central part of this extended array
    if dim_x : center[0] += (1 - 2*center[0]//nx)*nx
    if dim_y : center[1] += (1 - 2*center[1]//ny)*ny
    if dim_z : center[2] += (1 - 2*center[2]//nz)*nz
    
    # re-scale the in-plane versors by the dimensions
    vec_one = np.array(versor_one) * grid_dims[0]
    vec_two = np.array(versor_two) * grid_dims[1]

    # then, add and subtract the in-plane vectors from the plane center coordinates 
    rx = np.array([float(center[0])]*3)
    rx[0] += 0.5* ( -vec_one[0] -vec_two[0] ) / dx
    rx[1] += 0.5* ( +vec_one[0] -vec_two[0] ) / dx
    rx[2] += 0.5* ( -vec_one[0] +vec_two[0] ) / dx
    ry = np.array([float(center[1])]*3)
    ry[0] += 0.5* ( -vec_one[1] -vec_two[1] ) / dy
    ry[1] += 0.5* ( +vec_one[1] -vec_two[1] ) / dy
    ry[2] += 0.5* ( -vec_one[1] +vec_two[1] ) / dy
    rz = np.array([float(center[2])]*3)
    rz[0] += 0.5* ( -vec_one[2] -vec_two[2] ) / dz
    rz[1] += 0.5* ( +vec_one[2] -vec_two[2] ) / dz
    rz[2] += 0.5* ( -vec_one[2] +vec_two[2] ) / dz
    
    # so that  you can finally calculate the new data
    sec_x,sec_xx = np.meshgrid(np.linspace(rx[0],rx[1],grid_pnts[0]),np.linspace(rx[0],rx[2],grid_pnts[1]))
    sec_y,sec_yy = np.meshgrid(np.linspace(ry[0],ry[1],grid_pnts[0]),np.linspace(ry[0],ry[2],grid_pnts[1]))
    sec_z,sec_zz = np.meshgrid(np.linspace(rz[0],rz[1],grid_pnts[0]),np.linspace(rz[0],rz[2],grid_pnts[1]))
    
    sec_x = np.ndarray.flatten(sec_x+sec_xx) - rx[0]
    sec_y = np.ndarray.flatten(sec_y+sec_yy) - ry[0]
    sec_z = np.ndarray.flatten(sec_z+sec_zz) - rz[0]
    
    new_data = np.transpose(np.reshape(ndm.map_coordinates(arr, np.vstack((sec_x,sec_y,sec_z))),(grid_pnts[1],grid_pnts[0])))
    new_data = new_data.reshape((grid_pnts[0],grid_pnts[1],1))

    # and the new meta data 
    new_meta = self.meta.copy()
    new_meta['space_dim'] = '2D'
    
    new_meta['nx'] = grid_pnts[0]
    new_meta['ny'] = grid_pnts[1]
    new_meta['nz'] = 1
    
    new_meta['xl'] = grid_dims[0]
    new_meta['yl'] = grid_dims[1]
    new_meta['zl'] = 2*np.pi
    
    new_meta['dx'] = new_meta['xl'] / new_meta['nx']
    new_meta['dy'] = new_meta['yl'] / new_meta['ny']
    new_meta['dz'] = new_meta['zl'] / new_meta['nz']
    
    new_meta['nnn'] = (new_meta['nx'],new_meta['ny'],new_meta['nz'])
    new_meta['lll'] = (new_meta['xl'],new_meta['yl'],new_meta['zl'])
    new_meta['ddd'] = (new_meta['dx'],new_meta['dy'],new_meta['dz'])    
    
    del(new_meta['np_col'])
    del(new_meta['np_pla'])
    del(new_meta['np_row'])

    return new_data , new_meta

  #------------------------------------------------------------  
  def extract_range(self,
      tar_data,
      range_x,
      range_y,
      range_z):
    """
    Extract subarray ... 
    
    Parameters :
      - tar_data          [fibo.data] target variable
      - range_x           [int,int] range of points to be kept in x 
      - range_y           [int,int] range of points to be kept in y 
      - range_z           [int,int] range of points to be kept in z 
    
    Returns :
      - new_data          [np.ndarray(range_x,range_y,range_z)] 
      - new_meta          [dict] new metadata
    
    """  
    
    # let's fetch the data and meta
    new_data = self.get_data(tar_data)  
    new_meta = self.meta.copy()
    
    nx,ny,nz = self.meta['nnn']
    dx,dy,dz = self.meta['ddd']
    
    # here I determine how many relevant periodic dimensions there are
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    # then double the box in every relevant periodic dimension
    if dim_x : new_data = np.tile(new_data,(2,1,1))
    if dim_y : new_data = np.tile(new_data,(1,2,1))
    if dim_z : new_data = np.tile(new_data,(1,1,2))
    
    # now let's correct the coordinates of subarray to be extracted
    if dim_x : range_x[0] = (range_x[0] + nx)%nx
    if dim_y : range_y[0] = (range_y[0] + ny)%ny
    if dim_z : range_z[0] = (range_z[0] + nz)%nz

    if dim_x : range_x[1] = range_x[0]+(range_x[1]-range_x[0])%(nx+1)
    if dim_y : range_y[1] = range_y[0]+(range_y[1]-range_y[0])%(ny+1)
    if dim_z : range_z[1] = range_z[0]+(range_z[1]-range_z[0])%(nz+1)
    
    # and finally proceed with the extraction
    new_data = new_data[range_x[0]:range_x[1],range_y[0]:range_y[1],range_z[0]:range_z[1]]
    
    # then let's have the meta fixed
    if np.any([range_x[0]!=0 , range_x[1]!=nx]) :
     
      new_meta['xl'] = (range_x[1] - range_x[0]) * dx
      new_meta['nx'] = (range_x[1] - range_x[0])
       
      new_meta['lll'] = (new_meta['xl'], new_meta['lll'][1], new_meta['lll'][2])
      new_meta['nnn'] = (new_meta['nx'], new_meta['nnn'][1], new_meta['nnn'][2])
      new_meta['ppp'] = (False,          new_meta['ppp'][1], new_meta['ppp'][2])
       
    if np.any([range_y[0]!=0 , range_y[1]!=ny]) :
      
      new_meta['yl'] = (range_y[1] - range_y[0]) * dy
      new_meta['ny'] = (range_y[1] - range_y[0])
      
      new_meta['lll'] = (new_meta['lll'][0], new_meta['yl'], new_meta['lll'][2])
      new_meta['nnn'] = (new_meta['nnn'][0], new_meta['ny'], new_meta['nnn'][2])
      new_meta['ppp'] = (new_meta['ppp'][0], False,          new_meta['ppp'][2])
    
    if np.any([range_z[0]!=0 , range_z[1]!=nz]) :
      
      new_meta['zl'] = (range_z[1] - range_z[0]) * dz
      new_meta['nz'] = (range_z[1] - range_z[0])
      
      new_meta['lll'] = (new_meta['lll'][0], new_meta['lll'][1], new_meta['zl'])
      new_meta['nnn'] = (new_meta['nnn'][0], new_meta['nnn'][1], new_meta['nz'])
      new_meta['ppp'] = (new_meta['ppp'][0], new_meta['ppp'][1], False         )
      
    new_meta['space_dim'] =  str(np.count_nonzero(np.array(new_meta['nnn']) -1))+'D'
    
    return new_data, new_meta

  #------------------------------------------------------------  
  def extract_every(self,  
      tar_data,
      every_x,
      every_y,
      every_z):
    """
    Lowers resolution of array, generating a new array by taking one point every ...
    
    Parameters :
      - tar_data    [fibo.data] target variable
      - every_x     [int] points to be kept in x (one every ratio_x)
      - every_y     [int] points to be kept in y (one every ratio_y)
      - every_z     [int] points to be kept in z (one every ratio_z)
    
    Returns :
      - new_data    [np.ndarray(nx//ratio_x,ny//ratio_y,nz//ratio_z)] new variable
      - new_meta    [dict] new metadata    
    
    """
    
    #determine the shape of the new array
    arr = self.get_data(tar_data)
    nx,ny,nz = self.meta['nnn']

    nnx = int((nx+every_x-1)/every_x)
    nny = int((ny+every_y-1)/every_y)
    nnz = int((nz+every_z-1)/every_z)
    
    #create and fill in the new array
    new_data = np.zeros([nnx,nny,nnz], dtype=float)
    
    for iz in range(0,nz,every_z):
      for iy in range(0,ny,every_y):
        for ix in range(0,nx,every_x):
          new_data[ix//every_x,iy//every_y,iz//every_z] = arr[ix,iy,iz]

    #print('done with the new array!')

    #create new metadata for the object
    new_meta = {} 
    new_meta['nnn'] = (nnx, nny, nnz)
    new_meta['lll'] = self.meta['lll']
    new_meta['ddd'] = tuple( np.array(self.meta['ddd']) * np.array([every_x,every_y,every_z]) )
    new_meta['ppp'] = (np.array([nx,ny,nz]) % np.array([every_x,every_y,every_z])).astype('bool')
    new_meta['ppp'] = tuple( np.logical_not(new_meta['ppp']) )
    
    new_meta['nx'],new_meta['ny'],new_meta['nz'] = new_meta['nnn']
    new_meta['xl'],new_meta['yl'],new_meta['zl'] = new_meta['lll']
    new_meta['dx'],new_meta['dy'],new_meta['dz'] = new_meta['ddd']
    
    return new_data, new_meta
    
  #------------------------------------------------------------  
  def extract_interpol_fourier(self,  
      tar_data,
      ratio_x,
      ratio_y,
      ratio_z):
    """
    Increases resolution of array by perfect Fourier-interpolation
      (obviously, it acts only over periodic dimensions)
    
    Parameters :
      - tar_data    [fibo.data] target variable
      - ratio_x     [int] ratio of point increment in x
      - ratio_y     [int] ratio of point increment in y 
      - ratio_z     [int] ratio of point increment in z
    
    Returns :
      - new_data    [np.ndarray(nx*ratio_x,ny*ratio_y,nz*ratio_z)] new variable
      - new_meta    [dict] new metadata    
    
    """

    nx,ny,nz = self.meta['nnn']
    nnx,nny,nnz = self.meta['nnn']
    
    # here I determine how many relevant periodic dimensions are there
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    # let's perform the Fourier transform
    arr = np.fft.fftn( self.get_data(tar_data)  )

    # here I widen the Fourier-transformed array, adding zeros
    if dim_x : 
      zer = int(nx*(ratio_x-1.))
      arr = np.concatenate((arr[:nx//2,:,:], np.zeros((zer,nny,nnz)),arr[nx//2:,:,:]), axis=0)
      nnx = nx + zer
    if dim_y : 
      zer = int(ny*(ratio_y-1.))
      arr = np.concatenate((arr[:,:ny//2,:], np.zeros((nnx,zer,nnz)),arr[:,ny//2:,:]), axis=1) 
      nny = ny + zer
    if dim_z : 
      zer = int(nz*(ratio_z-1.))
      arr = np.concatenate((arr[:,:,:nz//2], np.zeros((nnx,nny,zer)),arr[:,:,nz//2:]), axis=2) 
      nnz = nz + zer

    # now find the actual ratios 
    ratio_x = nnx / nx
    ratio_y = nny / ny
    ratio_z = nnz / nz
    
    # now let's anti-transform the array: 
    new_data  = np.fft.ifftn( arr ).real
    if dim_x : new_data *= ratio_x 
    if dim_y : new_data *= ratio_y
    if dim_z : new_data *= ratio_z

    # and let's provide the meta-data for this new object 
    new_meta = {} 
    new_meta['nnn'] = (nnx, nny, nnz)
    new_meta['lll'] = self.meta['lll']
    new_meta['ddd'] = tuple( np.array(self.meta['ddd']) / np.array([ratio_x,ratio_y,ratio_z]) )
    new_meta['ppp'] = np.copy(self.meta['ppp'])

    new_meta['nx'],new_meta['ny'],new_meta['nz'] = new_meta['nnn']
    new_meta['xl'],new_meta['yl'],new_meta['zl'] = new_meta['lll']
    new_meta['dx'],new_meta['dy'],new_meta['dz'] = new_meta['ddd']

    return new_data, new_meta
