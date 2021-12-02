
import collections 

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

class fibo_extra : 

  #------------------------------------------------------------
  #-----extra-routines-----------------------------------------
  #------------------------------------------------------------
  def make_spect(self,  
      tar_data,  
      raw = False):    
      
    """ 
    Calculates fourier transform of some variable (real-valued)
    
    Parameters :
      - tar_data      [fibo.data] target field
      - raw = False   [bool] if False gives calc_slided(fou_data,nx//2,ny//2,nz//2)
    
    Returns :
      - fou_data      [new fibo.data]
    
    """

    #create a raw vector (symmetric part of our variable goes to real, 
    # antisymmetric part goes to imaginary, normalization to get integral ...)
    nx,ny,nz = self.meta['nnn']
    fou_data = np.absolute(np.fft.fftn(self.get_data(tar_data))) * 2. / (nx*ny*nz)

    if not raw : 
      if nx != 1 : fou_data = np.tile(fou_data,(2,1,1))
      if ny != 1 : fou_data = np.tile(fou_data,(1,2,1))
      if nz != 1 : fou_data = np.tile(fou_data,(1,1,2))

      fou_data = fou_data[nx//2:nx+nx//2,ny//2:ny+ny//2,nz//2:nz+nz//2]

    return fou_data

  #------------------------------------------------------------
  def make_plane(self,  
      normal):    
      
    """ 
    Calculates two in-plane vectors given the plane normal
    
    Parameters :
      - normal       [float,float,float] normal direction
    
    Returns :
      - versor_one   [float,float,float] first in-plane versor
      - versor_two   [float,float,float] second in-plane versor
    """ 
    
    # first, make the normal unitary
    normal /= np.sqrt((normal**2).sum())
  
    # second, check for the direction which is closest to the plane normal
    # if this is the same as the normal, take one of the two perpendiculars 
    quasi_normal = np.zeros(3) 
    quasi_normal[np.argmax(np.absolute(normal))] = 1.
    if np.all(quasi_normal == normal) : 
        quasi_normal = np.tile(quasi_normal,(2))[1:4]

    # third, find the two in-plane versors by applying the cross product
    vec_one  = np.cross(quasi_normal,normal)
    vec_one /= np.sqrt((vec_one**2).sum())
    vec_two  = np.cross(normal,vec_one)
    vec_two /= np.sqrt((vec_two**2).sum())

    return vec_one, vec_two

  #------------------------------------------------------------
  def stat_minimal(self,  
      tar_data, 
      tar_pnts = None,
      ave_label = 'ave',
      var_label = 'var',
      p12_label = 'p12',
      p50_label = 'p50',
      p88_label = 'p88'):
    """
    Performs average, variance, and the 12%, 50%, 88% percentiles 
    
    Parameters :
      - tar_data              [fibo.data] 
      - tar_pnts = None       [None OR fibo.pnts]
      - ave_label = 'ave'     [None OR str] name of average
      - var_label = 'var'     [None OR str] name of variance
      - p12_label = 'p12'     [None OR str] name of 12% percentile
      - p50_label = 'p50'     [None OR str] name of median
      - p88_label = 'p88'     [None OR str] name of 88% percentile
    
    Returns :
      - stat_dict             [dict] dictionary with {label : value}
    
    """

    tar_data = self.get_data(tar_data)
    tar_pnts = self.get_pnts(tar_pnts)
    if tar_pnts is None : tar_pnts = self.find_box()
    
    stat_dict = {}
    
    if (ave_label is not None) : stat_dict[ave_label] = np.average(self.get_data(tar_data)[list(tar_pnts)])
    if (var_label is not None) : stat_dict[var_label] = np.var(self.get_pnts(tar_data)[list(tar_pnts)])
    if (p12_label is not None) : stat_dict[p12_label] = np.percentile(self.get_data(tar_data)[list(tar_pnts)], 12)
    if (p50_label is not None) : stat_dict[p50_label] = np.percentile(self.get_data(tar_data)[list(tar_pnts)], 50)
    if (p88_label is not None) : stat_dict[p88_label] = np.percentile(self.get_data(tar_data)[list(tar_pnts)], 88)

    return stat_dict 

  #------------------------------------------------------------    
  def frame_MVA(self,  
      tar_var_x,
      tar_var_y,
      tar_var_z,
      pos_perm = 'LMN',
      pos_z = True):
    """
    Finds MVA frame relative to some dataset (can be any shape)
    
    Parameters :
      - tar_var_x        [fibo.data OR array] x component of the field you consider
      - tar_var_y        [fibo.data OR array] y component of the field you consider
      - tar_var_z        [fibo.data OR array] z component of the field you consider
      - pos_perm = 'LMN' ['LMN' OR 'NML'] positive permutation of the indices
      - pos_z = True     [bool] to impose that vec_L's z component is positive 
    
    Returns :
      - vec_N    [np.ndarray(3)]  #memory tip: NaMeLy
      - vec_M    [np.ndarray(3)]  
      - vec_L    [np.ndarray(3)]
      - vals     [np.ndarray(3)]
    
    """

    arrx = np.ndarray.flatten(self.get_data(tar_var_x))
    arry = np.ndarray.flatten(self.get_data(tar_var_y))
    arrz = np.ndarray.flatten(self.get_data(tar_var_z))

    matr = np.zeros([3,3])
    matr[0,0] = np.mean(arrx * arrx) - np.square(np.mean(arrx)) 
    matr[1,1] = np.mean(arry * arry) - np.square(np.mean(arry))
    matr[2,2] = np.mean(arrz * arrz) - np.square(np.mean(arrz))
    matr[0,1] = matr[1,0] = np.mean(arrx * arry) - np.mean(arrx)*np.mean(arry)
    matr[0,2] = matr[2,0] = np.mean(arrx * arrz) - np.mean(arrx)*np.mean(arrz)
    matr[1,2] = matr[2,1] = np.mean(arry * arrz) - np.mean(arry)*np.mean(arrz)

    vals, vecs = np.linalg.eigh(matr)
    vec_N, vec_M, vec_L = vecs[:,0], vecs[:,1], vecs[:,2]
    vec_N = np.reshape(vec_N, (3,1))
    vec_M = np.reshape(vec_M, (3,1))
    vec_L = np.reshape(vec_L, (3,1))

    #post-processing
    if (pos_z and vec_M[2] < 0.) : vec_M = -vec_M
    if pos_perm is 'LMN' : 
      vec_N[0],vec_N[1],vec_N[2] =  self.calc_cross(vec_L[0],vec_M[0],vec_L[1],vec_M[1],vec_L[2],vec_M[2])
    if pos_perm is 'NML' : 
      vec_N[0],vec_N[1],vec_N[2] = -self.calc_cross(vec_L[0],vec_M[0],vec_L[1],vec_M[1],vec_L[2],vec_M[2])

    return vec_N[:,0], vec_M[:,0], vec_L[:,0], vals
