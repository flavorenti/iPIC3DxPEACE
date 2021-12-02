
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

class fibo_get :

  #------------------------------------------------------------
  #-----routines-for-finding-data-and-meta-data----------------
  #------------------------------------------------------------
  def get_time_exit(self,
      tar_segs,
      tar_time,
      criteria = 'closer',
      silent=True):
    """
    Gives you the segment and exact exit number to look at ...
    
    Parameters :
      - tar_segs                [loader.segs] time segments 
      - tar_time                [float] time you want to be 
      - criteria = 'closer'     ['after','before','closer'] how to choose your time exit!
      - Returns :
      - to_load                 [str,int,float] seg to load, exit number, time difference
    
    """

    to_load = [' ',0,float("inf")]

    for seg in tar_segs.keys() :
      for exit_num in range(len(tar_segs[seg])):
        time_diff = float(tar_segs[seg][exit_num]) - tar_time

        if criteria == 'closer' : time_diff = np.absolute(time_diff)
        if criteria == 'before' : time_diff = -time_diff

        if (time_diff >= 0. and time_diff < to_load[2]):
          to_load[0] = seg
          to_load[1] = exit_num
          to_load[2] = time_diff

    if not silent:
      print('fibo_get_time_exit()>tar_time    :',tar_time)
      print('fibo_get_time_exit()>criteria    :',criteria)
      print('fibo_get_time_exit() return -->   ',to_load)

    return to_load

  #------------------------------------------------------------
  def get_data(self,
      tar_var):
    """
    Gives you the data array, overrunning fibo privileged data storage ...
    
    Parameters :
      - tar_var   [str OR np.ndarray OR None]   name in list(fibo.data) OR (nx,ny,nz) array
    
    Returns :
      - data      [fibo.data]
    
    """

    if   isinstance(tar_var, str) : return self.data[tar_var]
    elif isinstance(tar_var, np.ndarray) : return tar_var
    else : 
      print ('??WTW??')
      return None

  #------------------------------------------------------------
  def get_pnts(self,
      tar_var):
    """
    Gives you the data array, overrunning fibo privileged data storage ...
    
    Parameters :
      - tar_var   [str OR np.ndarray OR None]   name in list(fibo.pnts) OR (3,nn) array
    
    Returns :
      - data      [fibo.pnts]
    
    """

    if   isinstance(tar_var, str) : return self.pnts[tar_var]
    elif isinstance(tar_var, np.ndarray) : return tar_var
    else : 
      print ('??WTW??')
      return None
