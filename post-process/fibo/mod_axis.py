
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

class fibo_axis: 

  #------------------------------------------------------------ 
  #-----routines-for-axis-calculations-------------------------
  #------------------------------------------------------------
  def axis_data(self,  
      range_x = None,
      range_y = None,
      range_z = None,
      offset_x = 0,
      offset_y = 0,
      offset_z = 0,
      coor_like = False,
      squa_like = False,
      mesh_like = False,
      wrap_like = False):
    """
    Calculates axes values for some range around a given origin
    
    Parameters :
      - range_x = None    [int,int] x coordinates of the first and last point 
      - range_y = None    [int,int] y coordinates of the first and last point 
      - range_z = None    [int,int] z coordinates of the first and last point
      - offset_x = 0  [int] x coordinate of the origin
      - offset_y = 0  [int] y coordinate of the origin
      - offset_z = 0  [int] z coordinate of the origin
      - coor_like = False [bool] do you want axes to be coordinate-like scaled? 
      - squa_like = False [bool] do you want axes to measure square of distance from origin?
      - mesh_like = False [bool] do you want axes in mesh mode?
    
    Returns :
      - x_coor [fibo.data OR np.ndarray(nx)]
      - y_coor [fibo.data OR np.ndarray(ny)]
      - z_coor [fibo.data OR np.ndarray(nz)]
    
    """
    
    nx,ny,nz = self.meta['nnn']
    
    if range_x is None : range_x = [0,nx]
    if range_y is None : range_y = [0,ny]
    if range_z is None : range_z = [0,nz]
    
    # here I determine how many relevant periodic dimensions there are
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    # for every relevant periodic dimension, change the range if necessary
    if dim_x : range_x[0] = (range_x[0] + nx)%nx
    if dim_y : range_y[0] = (range_y[0] + ny)%ny
    if dim_z : range_z[0] = (range_z[0] + nz)%nz
    
    if dim_x : range_x[1] = range_x[0]+(range_x[1]-range_x[0])%(nx+1)
    if dim_y : range_y[1] = range_y[0]+(range_y[1]-range_y[0])%(ny+1)
    if dim_z : range_z[1] = range_z[0]+(range_z[1]-range_z[0])%(nz+1)
    
    # now let's produce the first version of our output arrays
    x_coor = np.arange( range_x[0], range_x[1] ) - offset_x
    y_coor = np.arange( range_y[0], range_y[1] ) - offset_y
    z_coor = np.arange( range_z[0], range_z[1] ) - offset_z
    
    # and re-change them over every periodic dimension
    if wrap_like : 
      if dim_x : x_coor = x_coor%nx
      if dim_y : y_coor = y_coor%ny
      if dim_z : z_coor = z_coor%nz
    
    # if you want coor-like: multiply by dx,dy,dz
    if coor_like :
      x_coor = x_coor*self.meta['ddd'][0]
      y_coor = y_coor*self.meta['ddd'][1]
      z_coor = z_coor*self.meta['ddd'][2]
    
    # if you want these squared (btw what? - whatever ...)
    if squa_like :
      x_coor = np.square(x_coor)
      y_coor = np.square(y_coor)
      z_coor = np.square(z_coor)
    
    # and if you want the mesh also (notice the weird order - yeah)
    if mesh_like : 
      y_coor, x_coor, z_coor = np.meshgrid(y_coor,x_coor,z_coor)


    return x_coor, y_coor, z_coor

  #------------------------------------------------------------ 
  def axis_line(self,  
      line_pts,  
      range_x,  
      range_y,  
      range_z):  
    """
    Calculates axis values along a specified line
    
    Parameters :
      - line_pts   [int] length (in points) of the line you ask to extract
      - range_x    [int,int] x coordinates of the first and last point 
      - range_y    [int,int] y coordinates of the first and last point 
      - range_z    [int,int] z coordinates of the first and last point
    
    Returns :
      - line       [np.ndarray(line_pts)]
    
    """
    
    line_len  = (self.meta['dx']*(range_x[1]-range_x[0]))**2
    line_len += (self.meta['dy']*(range_y[1]-range_y[0]))**2
    line_len += (self.meta['dz']*(range_z[1]-range_z[0]))**2
    line_len = np.sqrt(line_len)

    return np.linspace(0, line_len, line_pts, endpoint=True)