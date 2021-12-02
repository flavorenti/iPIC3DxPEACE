
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

class fibo_comp: 

  #------------------------------------------------------------ 
  #-------comp-routines-provide-other-quantities--------------- 
  #------------------------------------------------------------  

  def comp_distances(self,
        tar_pnts,
        all_comp): 
    """
    Find distance matrix for a set of points
    
    Parameters: 
      - tar_pnts     [fibo.pnts] original points 
      - all_comp     [bool] restitute all components of distances?

    Returns :
      - dist_mat     [float array] distance matrix, with or without components specified 
      
    """

    nx,ny,nz = self.meta['nnn']
    dx,dy,dz = self.meta['ddd']
    pnts = self.get_pnts(tar_pnts)
    
    # here I determine how many relevant periodic dimensions are there
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    dist_mat_x  = np.tile(pnts[0,:],(pnts.shape[1],1))
    dist_mat_x -= np.transpose(dist_mat_x)
    if dim_x : 
      dist_mat_x[ dist_mat_x >  nx/2 ] -= nx
      dist_mat_x[ dist_mat_x < -nx/2 ] += nx
    
    dist_mat_y  = np.tile(pnts[1,:],(pnts.shape[1],1))
    dist_mat_y -= np.transpose(dist_mat_y)
    if dim_y : 
      dist_mat_y[ dist_mat_y >  ny/2 ] -= ny
      dist_mat_y[ dist_mat_y < -ny/2 ] += ny
    
    dist_mat_z  = np.tile(pnts[2,:],(pnts.shape[1],1))
    dist_mat_z -= np.transpose(dist_mat_z)
    if dim_z : 
      dist_mat_z[ dist_mat_z >  nz/2 ] -= nz
      dist_mat_z[ dist_mat_z < -nz/2 ] += nz
      #anti[2, dist_mat_z > nz/2-1 ] = True 
    
    if all_comp: 
      dist_mat = [dx*dist_mat_x, dy*dist_mat_y, dz*dist_mat_z]
    else: 
      dist_mat =  np.sqrt( np.square( dx*dist_mat_x ) + np.square( dy*dist_mat_y ) + np.square( dz*dist_mat_z ) )
    
    return dist_mat

  #------------------------------------------------------------
  def comp_clusters(self,  
      tar_pnts,
      rad_cluster):
    """ 
    Finds clusters of points in a set
      (ok version, uses distance matrix)
    
    Parameters :
      - tar_pnts        [fibo.pnts] original points 
      - rad_cluster     [int>0] cluster radius 
    
    Returns :
      - cluster_labels  [int>0 array] labels of clusters
      - cluster_number  [int>0] number of clusters
    
    """

    nx,ny,nz = self.meta['nnn']
    pnts = self.get_pnts(tar_pnts)
    
    # original labels of the points
    cluster_labels = np.array([k for k in range(1,pnts.shape[1]+1)])
    
    # here you find all distances
    dist_mat = self.comp_distances(pnts)
    
    # and let's now collapse the labels 
    for k in range(pnts.shape[1]) : 
      pp = np.argwhere(dist_mat[k,:] < rad_cluster)[:,0]
      cluster_labels[ pp ] = np.min( cluster_labels[ pp ] )
    
    # which are the numbers you've used? how many clusters?
    old_labels = np.unique(cluster_labels)
    cluster_number = len(old_labels)
    
    # here you've got all your new labels - let's re-name them
    for kk in range(1,cluster_number+1) :      
      if kk not in old_labels :
        k = old_labels[0]
        pp = np.argwhere(cluster_labels==k)[:,0]
        cluster_labels[pp] = kk        
      old_labels = old_labels[1:]

    return cluster_labels, cluster_number


  #------------------------------------------------------------


  #------------------------------------------------------------  
  def comp_spot_number(self,
      tar_spots, 
      bin_struc = np.ones((3,3,3))):
    """ 
    Finds the number of connected spots
    
    Parameters :
      - tar_spots                     [fibo_data] binary mask individuating the spots
      - bin_struc = np.ones((3,3,3))  [3x3x3 binary matrix] binary structure for blurring
    
    Returns :
      - spot_num                      [int>0] connected spots in your environment
    
    """  
    
    nx,ny,nz = self.meta['nnn']
    arr = self.get_data(tar_spots)
    
    # here I determine how many relevant periodic dimensions are there
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1  
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 
    
    # let's double the array over every relevant periodic dimension
    if dim_x : arr = np.concatenate( [arr[-1:,:,:],arr,arr[:1,:,:]] ,axis=0 )
    if dim_y : arr = np.concatenate( [arr[:,-1:,:],arr,arr[:,:1,:]] ,axis=1 )
    if dim_z : arr = np.concatenate( [arr[:,:,-1:],arr,arr[:,:,:1]] ,axis=2 )
    
    # individuate the spot mask (this is the core of all)
    spot_mask, spot_number = ndm.measurements.label(arr,structure=bin_struc)
    # and fold the boundaries of the spot mask so that it stacks over the original array
    # then you can check the stacked edges for label-correspondences
    # and merge identical labels
    # at the same time, look for the tiling that allows for the complete view of each spot
    
    relabel = np.arange(1,spot_number+1)
    
    def corr_check(relabel,edg,check):
      # if check then perform re-labeling on the edge before checking correlations
      if check: 
        for k in range(len(relabel)) :
          kk = relabel[k]
          if kk != k+1 : 
            edg[edg == k+1] = kk
      # candidates for new labels
      edg_labels = np.min(edg,axis=0)
      # for every value that the label has in edg (except zero)
      for k in np.flip( np.unique(edg)[1:] ) :
        pos_k = np.where(edg == k)
        kk = np.min( edg_labels[pos_k[1:]] )
        # and - if different - substitute the latter to the former, everywhere
        if k != kk :     
          edg[pos_k] = kk
          edg_labels = np.min(edg,axis=0)
          loc_k = np.where(relabel == k)
          relabel[loc_k] = kk 
          #print('changed mask '+str(k)+' to '+str(kk))
    
    if dim_x: corr_check(relabel,np.array([spot_mask[:2,:,:],spot_mask[-2:,:,:]]),False)
    if dim_y: corr_check(relabel,np.array([spot_mask[:,:2,:],spot_mask[:,-2:,:]]),dim_x)
    if dim_z: corr_check(relabel,np.array([spot_mask[:,:,:2],spot_mask[:,:,-2:]]),dim_x|dim_y)

    # we're at the end: just re-name the labels now!
    spot_number = len(np.unique(relabel))

    return spot_number


  #------------------------------------------------------------  
  def comp_region_number(self,
      tar_data,
      lev_list,
      bin_struc = np.ones((3,3,3)),
      bin_dilat = np.ones((3,3,3))):
    """ 
    Finds the number of connected spots in a 2d toroidal environment
    
    Parameters :
      - tar_data                      [fibo_data] target variable
      - lev_list                      [list of floats] levels at which we have to cut 
      - bin_struc = np.ones((3,3,3))  [3x3x3 binary matrix] binary structure for connectivity
      - bin_dilat = np.ones((3,3,3))  [3x3x3 binary matrix] binary structure for dilation
    
    Returns :
      - spot_number                   [int>0] 
    
    """

    array = self.get_data(tar_data)
    reg_num_above = []
    reg_num_below = []

    for lvl in lev_list:

      #calculate the above bulk selector!
      spot_selector = ndm.binary_dilation(array > lvl,structure=bin_dilat)
      regions = self.comp_spot_number(spot_selector,bin_struc=bin_struc)
      #regions = self.calc_spots_2d(spot_selector,struct_ind=struct_ind)[1] #slow!
      reg_num_above.append(regions)

      #calculate the below bulk selector!
      spot_selector = ndm.binary_dilation(array < lvl,structure=bin_dilat)
      regions = self.comp_spot_number(spot_selector,bin_struc=bin_struc)
      #regions = self.calc_spots_2d(spot_selector,struct_ind=struct_ind)[1] #slow! 
      reg_num_below.append(regions)
    
    return reg_num_above, reg_num_below
  
  #------------------------------------------------------------    
  def comp_shell_pop(self,  
      tar_var,      
      offset_x,
      offset_y,
      offset_z,
      lvl_num,
      lvl_max,
      squared = False, 
      density = True):  
    """
    Calculates population in each spherical shell centered in offx, offy, offz
    
    Parameters :
      - tar_var          [fibo.data] target variable
      - offset_x         [int] offset in x for the center of concentric shells
      - offset_y         [int] offset in y for the center of concentric shells
      - offset_z         [int] offset in z for the center of concentric shells
      - lvl_num          [int] number of levels you want
      - lvl_max          [float] distance you want to go sempling from the center
      - squared = False  [bool] wanna label shells by the square of their distance from center?
      - density = True   [bool] wanna divide values in each shell by the number of points in it?
    
    Returns :
      - shell_sum        [np.ndarray(lvl_num)]  shell integral OR shell density
      - shell_num        [np.ndarray(lvl_num)]  number of points in the shell
    
    """

    #create the levels vector, that determines in which shell each point goes
    dl =  lvl_max/lvl_num
    arrx, arry, arrz = self.axis_data(offset_x=offset_x,offset_y=offset_y,offset_z=offset_z,coor_like=True,squa_like=True,mesh_like=True)
    levels = ((arrx +arry +arrz)/dl**2).astype(int)
    if not squared : levels = (np.sqrt(levels)).astype(int)

    #create the shell_pop array and fill it
    shell_pop = np.zeros([lvl_num,2])  

    #strange stuff: the following loop is faster on small arrays, but slower on big ones ...
    #  for il in range(lvl_num):
    #    shell_pop[il,0] = np.sum(self.get_data(tar_var)[levels == il])
    #    shell_pop[il,1] = np.sum(self.levels == il)
    #... than these three loops nested one into the other ...
    for iz in range(self.meta['nz']):
      for iy in range(self.meta['ny']):
        for ix in range(self.meta['nx']): 

          if dl*levels[ix,iy,iz] < lvl_max: #rad_max :
            shell_pop[levels[ix,iy,iz],0] += self.get_data(tar_var)[ix,iy,iz]
            shell_pop[levels[ix,iy,iz],1] += 1
    #... how possible?

    if density : 
      shell_pop[shell_pop[:,1] == 0,1] = 0.1
      shell_pop[:,0] /= shell_pop[:,1]

    #print(time.clock() - time_start)
    return shell_pop[:,0], shell_pop[:,1]  
