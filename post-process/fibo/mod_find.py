
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

class fibo_find:

  #------------------------------------------------------------
  #--the-find-functions-give-you-fibo.pnts---------------------
  #------------------------------------------------------------
  def find_box(self,
      range_x = None,
      range_y = None,
      range_z = None):
    """
    Gives list of points covering all box specified

    Parameters :
      - range_x = None       [int,int] x coordinates of the first and last point 
      - range_y = None       [int,int] y coordinates of the first and last point 
      - range_z = None       [int,int] z coordinates of the first and last point
    
    Returns :
      - tar_box              [fibo.pnts]
    
    """  

    xx,yy,zz = self.axis_data(range_x,range_y,range_z,mesh_like = True)

    xx = xx.flatten()
    yy = yy.flatten()
    zz = zz.flatten()         
    
    return np.array([xx,yy,zz])

  #------------------------------------------------------------
  def find_slice(self,
      tar_pnts):
    """
    From list of vertices gives all points inside 2d polygon in x,y plane
    
    Parameters :
      - tar_pnts        [fibo.pnts] x,y,z vertex coordinates 
    
    Returns :
      - tar_slice       [fibo.pnts]
    
    """
    
    nx,ny,nz = self.meta['nnn']
    
    cut_z = tar_pnts[2,0]
    plane_vertices = np.int32([ np.transpose( tar_pnts[0:2,:]) ])
    plane_selector = np.zeros( [ny,nx], dtype=np.uint8)
    plane_selector = cv2.fillPoly(plane_selector,plane_vertices,1)

    pnts_xy = np.transpose(np.argwhere(plane_selector==1))
    pnts_z = np.ones(np.shape(pnts_xy)[1]) * cut_z
    
    return (np.vstack([pnts_xy[1,:],pnts_xy[0,:],pnts_z[:]])).astype('int')

  #------------------------------------------------------------  
  def find_common(self,
      tar_pnts_1,
      tar_pnts_2):
    """ 
    Finds common points in two sets of coordinates
    
    Parameters :
      - tar_pnts_1      [fibo.pnts] first coordinate set
      - tar_pnts_2      [fibo.pnts] second coordinate set
    
    Returns :
      - ok_pnts        [fibo.pnts] common points
    
    """  

    # let's put the points into binary masks    
    pnts_1 = self.get_pnts(tar_pnts_1)
    pnts_2 = self.get_pnts(tar_pnts_2)
    
    mask_1 = np.zeros(self.meta['nnn'],dtype='bool')
    mask_2 = np.zeros(self.meta['nnn'],dtype='bool')

    mask_1[tuple(pnts_1.astype('int'))] = True
    mask_2[tuple(pnts_2.astype('int'))] = True

    # then compare the masks
    mask = np.logical_and(mask_1,mask_2)
    
    # and return the common points
    return np.array( np.where(mask) )
    
    # old version - should be 
    #set_a = set(tuple(i) for i in np.transpose(self.get_pnts(tar_pnts_1)))
    #set_b = set(tuple(i) for i in np.transpose(self.get_pnts(tar_pnts_2)))
    #myset = set_a.intersection(set_b)
    #return np.transpose(np.array([list(i) for i in myset]))
 

  #------------------------------------------------------------  
  def find_range(self,
      tar_pnts,
      range_x = None,
      range_y = None,
      range_z = None):
    """ 
    Finds points inside a range of coordinates
     (equivalent to find_common(tar_pnts,find_box(range_x,...)) but much faster)
    
    Parameters :
      - tar_pnts     [fibo.pnts] point set
      - range_x      [int,int] 
      - range_y      [int,int] 
      - range_z      [int,int] 
    
    Returns :
      - ok_pnts      [fibo.pnts] points inside range
    
    """  

    nx,ny,nz = self.meta['nnn']
    
    if range_x is None : range_x = [0,nx]
    if range_y is None : range_y = [0,ny]
    if range_z is None : range_z = [0,nz]
    
    mask_pnt = np.zeros(self.meta['nnn'],dtype='bool')
    mask_ran = np.zeros(self.meta['nnn'],dtype='bool')

    # let's put the points into a binary mask
    mask_pnt[tuple(self.get_pnts(tar_pnts).astype('int'))] = True

    # and create another binary mask on the point range 
    ran_x,ran_y,ran_z = self.axis_data(range_x,range_y,range_z,mesh_like=True,wrap_like=True)
    mask_ran[ran_x,ran_y,ran_z] = True

    # then compare the masks
    mask = np.logical_and(mask_pnt,mask_ran)
    
    # and return the common points
    return np.array( np.where(mask) )



  #------------------------------------------------------------
  def find_critical_squarely_2d(self,  
      tar_var,
      cut_z=0,
      aa=1, 
      bb=0,    
      thr=1e-1):
    """ 
    Finds zeros in 2d array by triangulation on neighbors 
      (x-y-periodic only)
    
    Parameters :
      - tar_var      [fibo.data] target field
      - cut_z = 0    [int] z coordinate of the cut considered
      - aa = 1       [int>0] first parameter to determine neighbors' diamond 
      - bb = 0       [int>0] second parameter to determine neighbors' diamond
      - thr=1e-1     [float>0] threshold value for maxima/minima detection 
    
    Returns :
      - coord_Oa     [fibo.pnts] maxima (O-points above)
      - coord_Ob     [fibo.pnts] minima (O-points below)
      - coord_X      [fibo.pnts] saddle (X-points)
    
    """
    nx,ny,nz = self.meta['nnn']
    arr = np.tile(self.get_data(tar_var)[:,:,cut_z],(2,2))

    #create eight (nx,ny) arrays with all differences between one value and its neighbours
    neigh = np.zeros([8,nx,ny])
  
    neigh[0,:,:] = (arr[aa:nx+aa, bb:ny+bb] - arr[0:nx, 0:ny])
    neigh[2,:,:] = (arr[nx-bb:2*nx-bb, aa:ny+aa] - arr[0:nx, 0:ny])
    neigh[4,:,:] = (arr[nx-aa:2*nx-aa, ny-bb:2*ny-bb] - arr[0:nx, 0:ny])
    neigh[6,:,:] = (arr[bb:nx+bb, ny-aa:2*ny-aa] - arr[0:nx, 0:ny])

    if aa > bb :
      neigh[1,:,:] = (arr[ aa-bb:nx+aa-bb, aa+bb:ny+aa+bb] - arr[0:nx, 0:ny])
      neigh[3,:,:] = (arr[ nx-aa-bb:2*nx-aa-bb, aa-bb:ny+aa-bb] - arr[0:nx, 0:ny])
      neigh[5,:,:] = (arr[ nx-aa+bb:2*nx-aa+bb, ny-aa-bb:2*ny-aa-bb] - arr[0:nx, 0:ny])
      neigh[7,:,:] = (arr[ aa+bb:nx+aa+bb, ny-aa+bb:2*ny-aa+bb] - arr[0:nx, 0:ny])
    else : 
      neigh[1,:,:] = (arr[ nx+aa-bb:2*nx+aa-bb, aa+bb:ny+aa+bb] - arr[0:nx, 0:ny])
      neigh[3,:,:] = (arr[ nx-aa-bb:2*nx-aa-bb, ny+aa-bb:2*ny+aa-bb] - arr[0:nx, 0:ny])
      neigh[5,:,:] = (arr[-aa+bb:nx-aa+bb, ny-aa-bb:2*ny-aa-bb] - arr[0:nx, 0:ny])
      neigh[7,:,:] = (arr[aa+bb:nx+aa+bb, -aa+bb:ny-aa+bb] - arr[0:nx, 0:ny])

    #calculate Delta plus and minus
    pos_neighs = np.maximum(neigh, np.zeros([8,nx,ny]))
    Delta_plus  = np.sum(pos_neighs,axis=0)
    Delta_minus = np.sum(pos_neighs-neigh,axis=0)

    #select positive neighbors and calculate sign changes
    pos_sel = pos_neighs > 0. 
    sign_change = np.sum(np.diff(pos_sel,axis=0),axis=0) + 1
    sign_change -= sign_change % 2

    #calculate coordinates of Oa points (maxima), Ob points (minima) and X points (saddles)
    coord_Oa = np.transpose(np.argwhere(np.logical_and(Delta_plus < thr*Delta_minus, sign_change == 0)))
    coord_Ob = np.transpose(np.argwhere(np.logical_and(Delta_minus < thr*Delta_plus, sign_change == 0)))
    coord_X = np.transpose(np.argwhere(sign_change > 2)) #np.logical_and(Delta_plus-Delta_minus < thr*(Delta_plus+Delta_minus), sign_change > 2)))

    coord_Oa = (np.vstack([coord_Oa[0],coord_Oa[1],np.ones(np.shape(coord_Oa)[1]) * cut_z])).astype('int')
    coord_Ob = (np.vstack([coord_Ob[0],coord_Ob[1],np.ones(np.shape(coord_Ob)[1]) * cut_z])).astype('int')
    coord_X = (np.vstack([coord_X[0],coord_X[1],np.ones(np.shape(coord_X)[1]) * cut_z])).astype('int')

    return coord_Oa, coord_Ob, coord_X

  #------------------------------------------------------------
  def find_critical_delaunay_2d(self, 
      tar_var,   
      cut_z=0,
      aa=1, 
      bb=0,  
      thr=1e-1):
    """ 
    Finds zeros in 2d array by triangulation on delaunay-selected neighbours 
      (x-y-periodic only)
    
    Parameters :
      - tar_var      [fibo_var] target field
      - cut_z = 0    [int] z coordinate of the cut considered
      - aa = 1       [int>0] first parameter to determine neighbors' diamond 
      - bb = 0       [int>0] second parameter to determine neighbors' diamond
      - thr=1e-1     [float>0] threshold value for maxima/minima detection 
    
    Returns :
      - coord_Oa     [fibo.pnts] maxima (O-points above)
      - coord_Ob     [fibo.pnts] minima (O-points below)
      - coord_X      [fibo.pnts] saddle (X-points)
    
    """    

    nx,ny,nz = self.meta['nnn'] 
    arr = np.tile(self.get_data(tar_var)[:,:,cut_z],(2,2))

    #create eight (nx,ny) arrays with all differences between one value and its neighbours
    neigh = np.zeros([8,nx,ny])

    neigh[0,:,:] = (arr[aa:nx+aa, bb:ny+bb] - arr[0:nx, 0:ny])
    neigh[2,:,:] = (arr[nx-bb:2*nx-bb, aa:ny+aa] - arr[0:nx, 0:ny])
    neigh[4,:,:] = (arr[nx-aa:2*nx-aa, ny-bb:2*ny-bb] - arr[0:nx, 0:ny])
    neigh[6,:,:] = (arr[bb:nx+bb, ny-aa:2*ny-aa] - arr[0:nx, 0:ny])

    if aa > bb :
      neigh[1,:,:] = (arr[ aa-bb:nx+aa-bb, aa+bb:ny+aa+bb] - arr[0:nx, 0:ny])
      neigh[3,:,:] = (arr[ nx-aa-bb:2*nx-aa-bb, aa-bb:ny+aa-bb] - arr[0:nx, 0:ny])
      neigh[5,:,:] = (arr[ nx-aa+bb:2*nx-aa+bb, ny-aa-bb:2*ny-aa-bb] - arr[0:nx, 0:ny])
      neigh[7,:,:] = (arr[ aa+bb:nx+aa+bb, ny-aa+bb:2*ny-aa+bb] - arr[0:nx, 0:ny])
    else : 
      neigh[1,:,:] = (arr[ nx+aa-bb:2*nx+aa-bb, aa+bb:ny+aa+bb] - arr[0:nx, 0:ny])
      neigh[3,:,:] = (arr[ nx-aa-bb:2*nx-aa-bb, ny+aa-bb:2*ny+aa-bb] - arr[0:nx, 0:ny])
      neigh[5,:,:] = (arr[-aa+bb:nx-aa+bb, ny-aa-bb:2*ny-aa-bb] - arr[0:nx, 0:ny])
      neigh[7,:,:] = (arr[aa+bb:nx+aa+bb, -aa+bb:ny-aa+bb] - arr[0:nx, 0:ny])

    #calculate which neighbors are to be discarded
    #ne_diags = np.transpose(np.argwhere(np.absolute(neigh[1,:,:]) > np.absolute(neigh[2,:,:] - neigh[0,:,:])))
    #nw_diags = np.transpose(np.argwhere(np.absolute(neigh[3,:,:]) >= np.absolute(neigh[4,:,:] - neigh[2,:,:])))
    #sw_diags = np.transpose(np.argwhere(np.absolute(neigh[5,:,:]) > np.absolute(neigh[6,:,:] - neigh[4,:,:])))
    #se_diags = np.transpose(np.argwhere(np.absolute(neigh[7,:,:]) >= np.absolute(neigh[0,:,:] - neigh[6,:,:])))
    
    plane_params_x = np.zeros([4,nx,ny])
    plane_params_y = np.zeros([4,nx,ny])

    plane_params_x[0,:,:] =  neigh[0,:,:] 
    plane_params_x[1,:,:] =  neigh[1,:,:] -neigh[2,:,:]
    plane_params_x[2,:,:] =  neigh[1,:,:]  
    plane_params_x[3,:,:] = -neigh[3,:,:] +neigh[2,:,:]
    
    plane_params_y[0,:,:] = -neigh[0,:,:] +neigh[1,:,:]
    plane_params_y[1,:,:] =  neigh[2,:,:] 
    plane_params_y[2,:,:] = -neigh[3,:,:] 
    plane_params_y[3,:,:] = -neigh[1,:,:] +neigh[2,:,:]

    plane_params_x /= self.meta['dx'] 
    plane_params_y /= self.meta['dy'] 
    
    cos_index = np.ones([2,nx,ny],dtype=float)
    
    cos_index[0,:,:] += np.multiply(plane_params_x[0,:,:] ,plane_params_x[1,:,:])
    cos_index[0,:,:] += np.multiply(plane_params_y[0,:,:] ,plane_params_y[1,:,:])
    cos_index[0,:,:] = np.divide(cos_index[0,:,:], np.sqrt(np.square(plane_params_x[0,:,:]) + np.square(plane_params_y[0,:,:])))
    cos_index[0,:,:] = np.divide(cos_index[0,:,:], np.sqrt(np.square(plane_params_x[1,:,:]) + np.square(plane_params_y[1,:,:])))
    
    cos_index[1,:,:] += np.multiply(plane_params_x[2,:,:] ,plane_params_x[3,:,:])
    cos_index[1,:,:] += np.multiply(plane_params_y[2,:,:] ,plane_params_y[3,:,:])
    cos_index[1,:,:] = np.divide(cos_index[1,:,:], np.sqrt(np.square(plane_params_x[2,:,:]) + np.square(plane_params_y[2,:,:])))
    cos_index[1,:,:] = np.divide(cos_index[1,:,:], np.sqrt(np.square(plane_params_x[3,:,:]) + np.square(plane_params_y[3,:,:])))
    
    diags_sel = (np.absolute(cos_index[1,:,:]) > np.absolute(cos_index[0,:,:]))
    diags_sel = np.tile(diags_sel,(2,2))
    
    ne_diags = np.transpose(np.argwhere(diags_sel[0:nx,0:ny]))
    nw_diags = np.transpose(np.argwhere(np.invert(diags_sel[nx-1:2*nx-1,0:ny])))
    sw_diags = np.transpose(np.argwhere(diags_sel[nx-1:2*nx-1,ny-1:2*ny-1]))
    se_diags = np.transpose(np.argwhere(np.invert(diags_sel[0:nx,ny-1:2*ny-1])))

    #calculate Delta plus and minus, careful about discarded neighbors
    neigh[1,ne_diags[0],ne_diags[1]] = 0.
    neigh[3,nw_diags[0],nw_diags[1]] = 0.
    neigh[5,sw_diags[0],sw_diags[1]] = 0.
    neigh[7,se_diags[0],se_diags[1]] = 0.

    pos_neighs = np.maximum(neigh, np.zeros([8,nx,ny]))
    Delta_plus  = np.sum(pos_neighs,axis=0)
    Delta_minus = np.sum(pos_neighs-neigh,axis=0)

    #select positive neighbors and calculate sign changes, careful about discarded neighbors
    good_neighs = (neigh > 0)

    good_neighs[1,ne_diags[0],ne_diags[1]] = np.logical_and(good_neighs[0,ne_diags[0],ne_diags[1]],good_neighs[2,ne_diags[0],ne_diags[1]])
    good_neighs[3,nw_diags[0],nw_diags[1]] = np.logical_and(good_neighs[2,nw_diags[0],nw_diags[1]],good_neighs[4,nw_diags[0],nw_diags[1]])
    good_neighs[5,sw_diags[0],sw_diags[1]] = np.logical_and(good_neighs[4,sw_diags[0],sw_diags[1]],good_neighs[6,sw_diags[0],sw_diags[1]])
    good_neighs[7,se_diags[0],se_diags[1]] = np.logical_and(good_neighs[6,se_diags[0],se_diags[1]],good_neighs[0,se_diags[0],se_diags[1]])
    good_neighs = np.diff(good_neighs,axis=0)

    sign_change = np.sum(good_neighs,axis=0) + 1
    sign_change -= sign_change % 2

    #calculate coordinates of Oa points (maxima), Ob points (minima) and X points (saddles)
    coord_Oa = np.transpose(np.argwhere(np.logical_and(Delta_plus < thr*Delta_minus, sign_change == 0)))
    coord_Ob = np.transpose(np.argwhere(np.logical_and(Delta_minus < thr*Delta_plus, sign_change == 0)))
    coord_X = np.transpose(np.argwhere(sign_change > 2)) #np.logical_and(Delta_plus-Delta_minus < thr*(Delta_plus+Delta_minus), sign_change > 2)))

    coord_Oa = (np.vstack([coord_Oa[0],coord_Oa[1],np.ones(np.shape(coord_Oa)[1]) * cut_z])).astype('int')
    coord_Ob = (np.vstack([coord_Ob[0],coord_Ob[1],np.ones(np.shape(coord_Ob)[1]) * cut_z])).astype('int')
    coord_X = (np.vstack([coord_X[0],coord_X[1],np.ones(np.shape(coord_X)[1]) * cut_z])).astype('int')

    return coord_Oa, coord_Ob, coord_X

  #------------------------------------------------------------
  def find_critical_starlike_2d(self,  
      tar_var,  
      cut_z=0,
      thr=1e-1):
    """ 
    Finds zeros in 2d array by triangulation on star of neighbours
      (x-y-periodic only)
    
    Parameters :
      - tar_var      [fibo_var] target field
      - cut_z = 0    [int] z coordinate of the cut considered
      - thr=1e-1     [float>0] threshold value for maxima/minima detection
    
    Returns :
      - coord_Oa     [fibo.pnts] maxima (O-points above)
      - coord_Ob     [fibo.pnts] minima (O-points below)
      - coord_X      [fibo.pnts] saddle (X-points)
    
    """

    nx,ny,nz = self.meta['nnn']
    arr = np.tile(self.get_data(tar_var)[:,:,cut_z],(2,2))
    
    #create sixteen (nx,ny) arrays with all differences between one value and its neighbours
    neigh = np.zeros([16,nx,ny])
    neigh[0,:,:] = (arr[1:nx+1,0:ny] - arr[0:nx,0:ny])
    neigh[2,:,:] = (arr[1:nx+1,1:ny+1] - arr[0:nx,0:ny])
    neigh[4,:,:] = (arr[0:nx,1:ny+1] - arr[0:nx,0:ny])
    neigh[6,:,:] = (arr[nx-1:2*nx-1,1:ny+1] - arr[0:nx,0:ny])
    neigh[8,:,:] = (arr[nx-1:2*nx-1,0:ny] - arr[0:nx,0:ny])
    neigh[10,:,:] = (arr[nx-1:2*nx-1,ny-1:2*ny-1] - arr[0:nx,0:ny])
    neigh[12,:,:] = (arr[0:nx,ny-1:2*ny-1] - arr[0:nx,0:ny])
    neigh[14,:,:] = (arr[1:nx+1,ny-1:2*ny-1] - arr[0:nx,0:ny])
    neigh[1,:,:] = (arr[2:nx+2,1:ny+1] - arr[0:nx,0:ny])
    neigh[3,:,:] = (arr[1:nx+1,2:ny+2] - arr[0:nx,0:ny])
    neigh[5,:,:] = (arr[nx-1:2*nx-1,2:ny+2] - arr[0:nx,0:ny])
    neigh[7,:,:] = (arr[nx-2:2*nx-2,1:ny+1] - arr[0:nx,0:ny])
    neigh[9,:,:] = (arr[nx-2:2*nx-2,ny-1:2*ny-1] - arr[0:nx,0:ny])
    neigh[11,:,:] = (arr[nx-1:2*nx-1,ny-2:2*ny-2] - arr[0:nx,0:ny])
    neigh[13,:,:] = (arr[1:nx+1,ny-2:2*ny-2] - arr[0:nx,0:ny])
    neigh[15,:,:] = (arr[2:nx+2,ny-1:2*ny-1] - arr[0:nx,0:ny])

    #calculate Delta plus and minus
    pos_neighs = np.maximum(neigh, np.zeros([16,nx,ny]))
    Delta_plus  = np.sum(pos_neighs,axis=0)
    Delta_minus = np.sum(pos_neighs-neigh,axis=0)

    #select positive neighbors and calculate sign changes
    pos_sel = pos_neighs > 0. 
    sign_change = np.sum(np.diff(pos_sel,axis=0),axis=0) + 1
    sign_change -= sign_change % 2

    #calculate coordinates of Oa points (maxima), Ob points (minima) and X points (saddles)
    coord_Oa = np.transpose(np.argwhere(np.logical_and(Delta_plus < thr*Delta_minus, sign_change == 0)))
    coord_Ob = np.transpose(np.argwhere(np.logical_and(Delta_minus < thr*Delta_plus, sign_change == 0)))
    coord_X = np.transpose(np.argwhere(sign_change > 2))
    
    coord_Oa = (np.vstack([coord_Oa[0],coord_Oa[1],np.ones(np.shape(coord_Oa)[1]) * cut_z])).astype('int')
    coord_Ob = (np.vstack([coord_Ob[0],coord_Ob[1],np.ones(np.shape(coord_Ob)[1]) * cut_z])).astype('int')
    coord_X = (np.vstack([coord_X[0],coord_X[1],np.ones(np.shape(coord_X)[1]) * cut_z])).astype('int')

    return coord_Oa, coord_Ob, coord_X

  #------------------------------------------------------------
  def find_middle_pnts(self,  
      tar_pnts,
      bar_pnt = True,
      med_pnt = True):
    """ 
    Finds relavant points in a given point-cluster
    
    Parameters :
      - tar_pnts      [fibo.pnts] set of points 
      - bar_pnt       [bool] do you want the barycentre?
      - med_pnt       [bool] do you want the medoid?
    
    Returns :
      - bar_coord     [fibo.pnts] barycentre coordinates
      - med_coord     [fibo.pnts] medoid coordintes
    
    """

    nx,ny,nz = self.meta['nnn']
    dx,dy,dz = self.meta['ddd']
    pnts = self.get_pnts(tar_pnts)

    # here I determine how many relevant periodic dimensions are there
    dim_x = self.meta['ppp'][0] and self.meta['nnn'][0] != 1 
    dim_y = self.meta['ppp'][1] and self.meta['nnn'][1] != 1 
    dim_z = self.meta['ppp'][2] and self.meta['nnn'][2] != 1 

    # create the relevant-point coordinate arrays
    med = np.zeros(3)
    bar = np.zeros(3)

    # let's find the medoid first
    # it's easier given that we are in a usually-periodic environment
    dist_mat = self.comp_distances(pnts)
    med_ind = np.argmin( np.sum(dist_mat, axis=0) ) 
    med = pnts[:,med_ind].astype('int') 
  
    # now we can go for the barycentre: 
    if bar_pnt : 
      # let's centre onto the medoid our coordinate system
      clus_coords = pnts - np.tile( med.reshape(3,1),(1,pnts.shape[1])) 
      clus_coords = clus_coords.astype('float')
      
      if dim_x : 
        clus_coords[ 0, clus_coords[0,:] >  nx/2  ] -= nx
        clus_coords[ 0, clus_coords[0,:] < -nx/2  ] += nx    
      if dim_y : 
        clus_coords[ 1, clus_coords[1,:] >  ny/2  ] -= ny
        clus_coords[ 1, clus_coords[1,:] < -ny/2  ] += ny  
      if dim_z :
        clus_coords[ 2, clus_coords[2,:] >  nz/2  ] -= nz
        clus_coords[ 2, clus_coords[2,:] < -nz/2  ] += nz  
      
      # then find the barycentre and back to original coordinates  
      bar = med + np.round( np.array([np.sum(clus_coords[i,:])/pnts.shape[1] for i in range(3)]) )
      if dim_x : bar[0] = bar[0]%nx
      if dim_y : bar[1] = bar[1]%ny
      if dim_z : bar[2] = bar[2]%nz
  
      # note that this procedure might be problematic if the set radius is too big 
      # so let me caluculate the radius here :
      rad_x = np.max(clus_coords[0]) - np.min(clus_coords[0])
      rad_y = np.max(clus_coords[1]) - np.min(clus_coords[1])
      rad_z = np.max(clus_coords[2]) - np.min(clus_coords[2])
  
      # and set an alert: 
      if np.any( [rad_x==nx-1 and dim_x,rad_y==ny-1 and dim_y,rad_z==nz-1 and dim_z] ) : 
        print('compliments! your point set is pretty big - you might consider looking into it in detail to be sure of what happens tho')
    
    # now return the values
    if bar_pnt and not med_pnt : return bar
    if med_pnt and not bar_pnt : return med
    if med_pnt and bar_pnt : return bar, med