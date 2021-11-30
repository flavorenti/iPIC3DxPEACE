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

class fibo_calc: 

  #------------------------------------------------------------
  #--SID:the-calc-functions-give-you-fibo.data-----------------
  #--JOB:some-calc-functions-give-you-fibo.data----------------
  #--JOB:_scalr,_cross,_par_per-have-optional-fill-mode--------
  #------------------------------------------------------------
  def calc_slided(self,  
      tar_var,
      offset_x,
      offset_y,
      offset_z):
    """
    Slides your array so to see the box with another 'origin' ...
    
    Parameters :
      - tar_var    [fibo_var] target variable
      - offset_x   [float] offset we'll introduce in x (move for offx)
      - offset_y   [float] offset we'll introduce in y (move for offy)
      - offset_z   [float] offset we'll introduce in z (move for offz)
    
    Returns :
      - new_var    [fibo.data] same variable but slided
    
    """

    #correct any offset that is bigger than box size
    offset_x = offset_x % self.meta['nx']
    offset_y = offset_y % self.meta['ny']
    offset_z = offset_z % self.meta['nz']

    arr = self.get_data(tar_var)

    if offset_x != 0 : arr = np.tile(arr,(2,1,1))
    if offset_y != 0 : arr = np.tile(arr,(1,2,1))
    if offset_z != 0 : arr = np.tile(arr,(1,1,2))

    return arr[offset_x:self.meta['nx']+offset_x,offset_y:self.meta['ny']+offset_y,offset_z:self.meta['nz']+offset_z]

  #------------------------------------------------------------
  def calc_scalr(self,  
      tar_var_1x,  
      tar_var_2x,  
      tar_var_1y,  
      tar_var_2y,
      tar_var_1z = None,
      tar_var_2z = None,
      seg = None,
      fill = True):  
    """
    Calculates scalar product between two vectors
    
    Parameters :
      - tar_var_1x          [fibo.data]  first target variable: x component
      - tar_var_2x          [fibo.data]  second target variable: x component
      - tar_var_1y          [fibo.data]  first target variable: y component
      - tar_var_2y          [fibo.data]  second target variable: y component
      - tar_var_1z = None   [None OR fibo.data]  first target variable: z component
      - tar_var_2z = None   [None OR fibo.data]  second target variable: z component
      - seg = None          [None or str] segment to use, if None use simply tar_var
      - fill = True         [bool] do you want to add this new array to the class self?

    Returns :
     - new_var             [fibo.data]  scalar product
    
    """
    if fill : new_tar = tar_var_1x.replace('_x','')+tar_var_2x.replace('_x','')+'%.8i'%int(seg)

    if (seg != None) :
      tar_var_1x = tar_var_1x+'%.8i'%int(seg)
      tar_var_2x = tar_var_2x+'%.8i'%int(seg)
      tar_var_1y = tar_var_1y+'%.8i'%int(seg)
      tar_var_2y = tar_var_2y+'%.8i'%int(seg)
      tar_var_1z = tar_var_1z+'%.8i'%int(seg)
      tar_var_2z = tar_var_2z+'%.8i'%int(seg)

    tar_arr  = np.multiply(self.get_data(tar_var_1x),self.get_data(tar_var_2x))
    tar_arr += np.multiply(self.get_data(tar_var_1y),self.get_data(tar_var_2y))
    if not (tar_var_1z is None) :
      tar_arr += np.multiply(self.get_data(tar_var_1z),self.get_data(tar_var_2z))

    if fill : self.data[new_tar] = tar_arr
    else : return tar_arr

  #------------------------------------------------------------
  def calc_cross(self,  
      tar_var_1x,
      tar_var_2x,
      tar_var_1y,  
      tar_var_2y,
      tar_var_1z = None,  
      tar_var_2z = None,
      seg = None,
      fill = True):  
    """
    Calculates cross product between two vectors
    
    Parameters :
      - tar_var_1x          [fibo.data]  first target variable: x component
      - tar_var_2x          [fibo.data]  second target variable: x component
      - tar_var_1y          [fibo.data]  first target variable: y component
      - tar_var_2y          [fibo.data]  second target variable: y component
      - tar_var_1z = None   [None OR fibo.data]  first target variable: z component
      - tar_var_2z = None   [None OR fibo.data]  second target variable: z component
      - seg = None          [None or str] segment to use, if None use simply tar_var
      - fill = True         [bool] do you want to add this new array to the class self?
    
    Returns :
      - new_var_x    [None OR fibo.data]  cross product: x component
      - new_var_y    [None OR fibo.data]  cross product: y component
      - new_var_z    [fibo.data]  cross product: z component
    
    """
    if fill :
      new_tar_x = tar_var_1x.replace('_x','')+'x'+tar_var_2x.replace('_x','')+'_x%.8i'%int(seg)
      new_tar_y = tar_var_1x.replace('_x','')+'x'+tar_var_2x.replace('_x','')+'_y%.8i'%int(seg)
      new_tar_z = tar_var_1x.replace('_x','')+'x'+tar_var_2x.replace('_x','')+'_z%.8i'%int(seg)

    if (seg != None) :
      tar_var_1x = tar_var_1x+'%.8i'%int(seg)
      tar_var_2x = tar_var_2x+'%.8i'%int(seg)
      tar_var_1y = tar_var_1y+'%.8i'%int(seg)
      tar_var_2y = tar_var_2y+'%.8i'%int(seg)
      tar_var_1z = tar_var_1z+'%.8i'%int(seg)
      tar_var_2z = tar_var_2z+'%.8i'%int(seg)

    if not (tar_var_1z is None) :
      tar_arr_x  = np.multiply(self.get_data(tar_var_1y),self.get_data(tar_var_2z))
      tar_arr_x -= np.multiply(self.get_data(tar_var_1z),self.get_data(tar_var_2y))
      tar_arr_y  = np.multiply(self.get_data(tar_var_1z),self.get_data(tar_var_2x))
      tar_arr_y -= np.multiply(self.get_data(tar_var_1x),self.get_data(tar_var_2z))
    tar_arr_z  = np.multiply(self.get_data(tar_var_1x),self.get_data(tar_var_2y))
    tar_arr_z -= np.multiply(self.get_data(tar_var_1y),self.get_data(tar_var_2x))

    if fill :
      self.data[new_tar_z] = tar_arr_z
      if (tar_var_1z != None) :
        self.data[new_tar_y] = tar_arr_y
        self.data[new_tar_x] = tar_arr_x
    elif (tar_var_1z == None) : 
      return None, None, tar_arr_z
    elif (tar_var_1z != None) :
      return tar_arr_x, tar_arr_y, tar_arr_z

  #------------------------------------------------------------
  def calc_project(self,
      tar_var_x,  
      ref_var_1x,  
      ref_var_2x,  
      ref_var_3x,  
      tar_var_y,  
      ref_var_1y,  
      ref_var_2y,  
      ref_var_3y,  
      tar_var_z = None,  
      ref_var_1z = None,  
      ref_var_2z = None,  
      ref_var_3z = None):  
    """
    Projects a vector field over three vector components
    
    Parameters :
      - tar_var_x            [fibo.data] target variable: x component  
      - ref_var_1x           [fibo.data OR np.float] x component of 1st versor 
      - ref_var_2x           [fibo.data OR np.float] x component of 2nd versor 
      - ref_var_3x           [fibo.data OR np.float] x component of 3rd versor 
      - tar_var_y            [fibo.data] target variable: y component
      - ref_var_1y           [fibo.data OR np.float] y component of 1st versor 
      - ref_var_2y           [fibo.data OR np.float] y component of 2nd versor 
      - ref_var_3y           [fibo.data OR np.float] y component of 3st versor 
      - tar_var_z = None     [None OR fibo.data] target variable: z component
      - ref_var_1z = None    [None OR fibo.data OR np.float] z component of 1st versor 
      - ref_var_2z = None    [None OR fibo.data OR np.float] z component of 2nd versor 
      - ref_var_3z = None    [None OR fibo.data OR np.float] z component of 3rd versor 
    
    Returns :
      - new_var_x      [fibo.data] projected field: x component
      - new_var_y      [fibo.data] projected field: y component
      - new_var_z      [None OR fibo.data] projected field: z component
    
    """

    #if the reference is given as single-float array, make it the same shape as the other!  
    if  np.shape(self.get_data(ref_var_1x)) == (1,) :
      arr_shape = np.shape(self.get_data(tar_var_x))

      ref_1x = np.ones(arr_shape) * self.get_data(ref_var_1x)
      ref_2x = np.ones(arr_shape) * self.get_data(ref_var_2x)
      ref_3x = np.ones(arr_shape) * self.get_data(ref_var_3x)
      ref_1y = np.ones(arr_shape) * self.get_data(ref_var_1y)
      ref_2y = np.ones(arr_shape) * self.get_data(ref_var_2y)
      ref_3y = np.ones(arr_shape) * self.get_data(ref_var_3y)
      if tar_var_z is None :
        ref_1z = ref_2z = ref_3z = None
      else :
        ref_1z = np.ones(arr_shape) * self.get_data(ref_var_1z)
        ref_2z = np.ones(arr_shape) * self.get_data(ref_var_2z)
        ref_3z = np.ones(arr_shape) * self.get_data(ref_var_3z)

    #otherwise, just see how the thing projects point by point 
    else : 
      ref_1x = self.get_data(ref_var_1x)
      ref_2x = self.get_data(ref_var_2x)
      ref_3x = self.get_data(ref_var_3x)
      ref_1y = self.get_data(ref_var_1y)
      ref_2y = self.get_data(ref_var_2y)
      ref_3y = self.get_data(ref_var_3y)
      if tar_var_z is None :
        ref_1z = ref_2z = ref_3z = None
      else :
        ref_1z = self.get_data(ref_var_1z)
        ref_2z = self.get_data(ref_var_2z)
        ref_3z = self.get_data(ref_var_3z)

    pro_1 = self.calc_scalr(tar_var_x,ref_1x, tar_var_y,ref_1y, tar_var_z,ref_1z)
    pro_2 = self.calc_scalr(tar_var_x,ref_2x, tar_var_y,ref_2y, tar_var_z,ref_2z)
    pro_3 = self.calc_scalr(tar_var_x,ref_3x, tar_var_y,ref_3y, tar_var_z,ref_3z)
  
    return pro_1, pro_2, pro_3

  #------------------------------------------------------------
  def calc_par_per(self,  
      tar_var_x,  
      ref_var_x,
      tar_var_y,
      ref_var_y,
      tar_var_z = None,
      ref_var_z = None,
      seg = None,
      fill = True):  
    """
    Calculates parallel and perpendicolar parts of tar_var with respect to ref_var
    
    Parameters :
      - tar_var_x         [fibo.data] tar_var: x component
      - ref_var_x         [fibo.data] ref_var: x component
      - tar_var_y         [fibo.data] tar_var: y component
      - ref_var_y         [fibo.data] ref_var: y component
      - tar_var_z = None  [None OR fibo.data] tar_var: z component
      - ref_var_z = None  [None OR fibo.data] ref_var: z component
    
    Returns :
      - par_var_x      [fibo.data] parallel part: x component
      - par_var_y      [fibo.data] parallel part: y component
      - par_var_z      [None OR fibo.data] parallel part: z component
      - per_var_x      [fibo.data] perpendicular part: x component
      - per_var_y      [fibo.data] perpendicular part: y component
      - per_var_z      [None OR fibo.data] perpendicular part: z component

    """
    if fill :
      par_tar_x = tar_var_x.replace('_x','')+'par'+ref_var_x.replace('_x','')+'_x%.8i'%int(seg)
      par_tar_y = tar_var_x.replace('_x','')+'par'+ref_var_x.replace('_x','')+'_y%.8i'%int(seg)
      par_tar_z = tar_var_x.replace('_x','')+'par'+ref_var_x.replace('_x','')+'_z%.8i'%int(seg)
      per_tar_x = par_tar_x.replace('par','per')
      per_tar_y = par_tar_y.replace('par','per')
      per_tar_z = par_tar_z.replace('par','per')

    if (seg != None) :
      tar_var_x = tar_var_x+'%.8i'%int(seg)
      ref_var_x = ref_var_x+'%.8i'%int(seg)
      tar_var_y = tar_var_y+'%.8i'%int(seg)
      ref_var_y = ref_var_y+'%.8i'%int(seg)
      tar_var_z = tar_var_z+'%.8i'%int(seg)
      ref_var_z = ref_var_z+'%.8i'%int(seg)

    fact  = np.reciprocal(self.calc_scalr(ref_var_x,ref_var_x,ref_var_y,ref_var_y,ref_var_z,ref_var_z,fill=False))
    fact *= self.calc_scalr(tar_var_x,ref_var_x,tar_var_y,ref_var_y,tar_var_z,ref_var_z,fill=False)

    par_arrx = self.get_data(ref_var_x)*fact
    par_arry = self.get_data(ref_var_y)*fact
    per_arrx = self.get_data(tar_var_x)-par_arrx
    per_arry = self.get_data(tar_var_y)-par_arry
    if not (tar_var_z is None) : 
      par_arrz = self.get_data(ref_var_z)*fact
      per_arrz = self.get_data(tar_var_z)-par_arrz

    if fill :
      self.data[par_tar_x] = par_arrx
      self.data[par_tar_y] = par_arry
      self.data[per_tar_x] = per_arrx
      self.data[per_tar_y] = per_arry
      if (tar_var_z != None) :
        self.data[par_tar_z] = par_arrz
        self.data[per_tar_z] = per_arrz
    elif (tar_var_z == None) : 
      return par_arrx, par_arry, None, per_arrx, per_arry, None
    elif (tar_var_z != None) :
      return par_arrx, par_arry, par_arrz, per_arrx, per_arry, per_arrz


  #------------------------------------------------------------
  def calc_irr_sol(self,
      tar_var_x,
      tar_var_y,
      tar_var_z = None):
    """
    Calculates irrotational and solenoidal components of the given vector field
    
    Parameters :
      - tar_var_x           [fibo.data] target variable: x component
      - tar_var_y           [fibo.data] target variable: y component
      - tar_var_z = None    [None OR fibo.data] target variable: z component
    
    Returns :
      - irr_var_x           [fibo.data] irrotational part: x component
      - irr_var_y           [fibo.data] irrotational part: y component
      - irr_var_z           [None OR fibo.data] irrotational part: z component
      - sol_var_x           [fibo.data] solenoidal part: x component
      - sol_var_y           [fibo.data] solenoidal part: y component
      - sol_var_z           [None OR fibo.data] solenoidal part: z component
    
    Credit: L.R.
    """

    #determine the shape of the array, correct any offset that is bigger than box size
    nx,ny,nz = np.shape(self.get_data(tar_var_x))
    nxm = nx//2
    nym = ny//2
    nzm = nz//2

    myaxes = ()
    if nx > 1 : myaxes = myaxes + (0,)
    if ny > 1 : myaxes = myaxes + (1,)
    if nz > 1 : myaxes = myaxes + (2,)

    #if not (tar_var_1z is None) :     
    G_x = np.fft.fftn(self.get_data(tar_var_x),axes=myaxes)
    G_y = np.fft.fftn(self.get_data(tar_var_y),axes=myaxes)
    G_z = np.fft.fftn(self.get_data(tar_var_z),axes=myaxes)

    #print(myaxes)
    GPsi = np.zeros([nx,ny,nz],dtype='complex128')
    for ix in range(1,nxm):
      GPsi[ix,:,:] += 1j*G_x[ix,:,:]*ix / self.meta['xl']
      GPsi[nx-ix,:,:] += -1j*G_x[nx-ix,:,:]*ix / self.meta['xl']
    for iy in range(1,nym):
      GPsi[:,iy,:] += 1j*G_y[:,iy,:]* iy / self.meta['yl']
      GPsi[:,ny-iy,:] += -1j*G_y[:,ny-iy,:]* iy / self.meta['yl']
    for iz in range(1,nzm):
      GPsi[:,:,iz] += 1j*G_z[:,:,iz]* iz / self.meta['zl']
      GPsi[:,:,nz-iz] += -1j*G_z[:,:,nz-iz]* iz / self.meta['zl']

    del(G_x, G_y, G_z)
    mesh_x = np.arange(0,nx) / self.meta['xl']
    mesh_y = np.arange(0,ny) / self.meta['yl']
    mesh_z = np.arange(0,nz) / self.meta['zl']
    mesh_y, mesh_x, mesh_z = np.meshgrid(mesh_y,mesh_x,mesh_z)   #achtung!  (y,x,z) is correct order!  
    mesh = np.square(mesh_x) + np.square(mesh_y) + np.square(mesh_z)
    mesh[0,0,0] = 1.

    GPsi = -1j* np.multiply(GPsi, np.reciprocal(mesh))

    # back to physical space
    arrx = (np.fft.ifftn((GPsi * mesh_x),axes=myaxes)).real
    arrxx = self.get_data(tar_var_x)-arrx
    arry = (np.fft.ifftn((GPsi * mesh_y),axes=myaxes)).real
    arryy = self.get_data(tar_var_y)-arry
    if (tar_var_z is None) : 
      arrz = None    #print(time.clock()-time_start)
      arrzz = None
    else : #if (myaxes == (0, 1, 2)) :
      arrz = (np.fft.ifftn((GPsi * mesh_z),axes=myaxes)).real
      arrzz = self.get_data(tar_var_z)-arrz
    #else :
    #  arrz = np.zeros([nx,ny,nz])
    #  arrzz = np.zeros([nx,ny,nz])

    return arrx, arry, arrz, arrxx, arryy, arrzz


  #------------------------------------------------------------
  def calc_gradx(self,
      tar_var,
      periodic = None,
      der_ord = 1,
      brutal = False):
    """ 
    Calculates x gradient component with(out) periodical boundary conditions
    
    Parameters :
      - tar_var           [fibo.data] target variable of the procedure
      - periodic = None   [None OR bool] periodic boundary conditions? if None will consult meta
      - der_ord = 1       [int] order of derivation (0: no derivative, 1: first derivative, 2: second derivative ...)
      - brutal = False     [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - xgr_var           [fibo.data]
    
    """

    #determine the coefficients
    nx,ny,nz = self.meta['nnn']
    if periodic is None : periodic = self.meta['ppp'][0]
    
    if nx == 1 : 
      return np.zeros([nx,ny,nz])

    else: 
    
      # set the coefficients !
      if brutal : 
        oo_6 = self.meta['dx']**(-der_ord) * np.array([1.0,  0.0, -2.0])      
        aa_6 = self.meta['dx']**(-der_ord) * np.array([0.0,  0.5, 1.0 ]) 
        bb_6 = self.meta['dx']**(-der_ord) * np.array([0.0,  0.0, 0.0]) 
        cc_6 = self.meta['dx']**(-der_ord) * np.array([0.0,  0.0, 0.0]) 
      else : 
        oo_6 = self.meta['dx']**(-der_ord) * np.array([1.0,  0.0, -49./18.])      
        aa_6 = self.meta['dx']**(-der_ord) * np.array([0.0,  9.0/12.0,  3./2. ]) 
        bb_6 = self.meta['dx']**(-der_ord) * np.array([0.0, -3.0/20.0, -3./20.]) 
        cc_6 = self.meta['dx']**(-der_ord) * np.array([0.0,  1.0/60.0,  1./90.]) 
    
      #create the new vector, fill it
      if periodic :
        ff = np.tile(self.get_data(tar_var),(2,1,1))
        dx_f  = oo_6[der_ord] * ff[0:nx,:,:] 
        dx_f += aa_6[der_ord] *(ff[1:1+nx,:,:] - ff[nx-1:2*nx-1,:,:])
        dx_f += bb_6[der_ord] *(ff[2:2+nx,:,:] - ff[nx-2:2*nx-2,:,:])
        dx_f += cc_6[der_ord] *(ff[3:3+nx,:,:] - ff[nx-3:2*nx-3,:,:])
      else :
        f = self.get_data(tar_var)  
        dx_f  = np.zeros([nx,ny,nz])  
        dx_f[3:nx-3,:,:]  = oo_6[der_ord] * f[3:nx-3,:,:]
        dx_f[3:nx-3,:,:] += aa_6[der_ord] *(f[4:nx-2,:,:] - f[2:nx-4,:,:])
        dx_f[3:nx-3,:,:] += bb_6[der_ord] *(f[5:nx-1,:,:] - f[1:nx-5,:,:])
        dx_f[3:nx-3,:,:] += cc_6[der_ord] *(f[6:nx-0,:,:] - f[0:nx-6,:,:])

      #print('done with the new array!')
      return dx_f

  #------------------------------------------------------------  
  def calc_grady(self,    
      tar_var,
      periodic = None,
      der_ord = 1,
      brutal = False): 
    """ 
    Calculates y gradient component with(out) periodical boundary conditions
    
    Parameters :
      - tar_var           [fibo.data] target variable of the procedure
      - periodic = None   [None OR bool] periodic boundary conditions? - if None will consult meta
      - der_ord = 1       [int] order of derivation (0: no derivative, 1: first derivative, 2: second derivative ...)
      - brutal = False     [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - ygr_var           [fibo.data]
    
    """

    #determine the coefficients
    nx,ny,nz = self.meta['nnn']
    if periodic is None : periodic = self.meta['ppp'][1]

    if ny == 1 : 
      return np.zeros([nx,ny,nz])

    else: 
    
      # set the coefficients !
      if brutal : 
        oo_6 = self.meta['dy']**(-der_ord) * np.array([1.0,  0.0, -2.0])      
        aa_6 = self.meta['dy']**(-der_ord) * np.array([0.0,  0.5, 1.0 ]) 
        bb_6 = self.meta['dy']**(-der_ord) * np.array([0.0,  0.0, 0.0]) 
        cc_6 = self.meta['dy']**(-der_ord) * np.array([0.0,  0.0, 0.0]) 
      else : 
        oo_6 = self.meta['dy']**(-der_ord) * np.array([1.0,  0.0, -49./18.])      
        aa_6 = self.meta['dy']**(-der_ord) * np.array([0.0,  9.0/12.0,  3./2. ]) 
        bb_6 = self.meta['dy']**(-der_ord) * np.array([0.0, -3.0/20.0, -3./20.]) 
        cc_6 = self.meta['dy']**(-der_ord) * np.array([0.0,  1.0/60.0,  1./90.]) 
    
      #create the new vector, fill it  (periodic / nonperiodic boundary cases ...)  
      if periodic :
        ff = np.tile(self.get_data(tar_var),(1,2,1))
        dy_f  = oo_6[der_ord] * ff[:,0:ny,:] 
        dy_f += aa_6[der_ord] *(ff[:,1:1+ny,:] - ff[:,ny-1:2*ny-1,:])
        dy_f += bb_6[der_ord] *(ff[:,2:2+ny,:] - ff[:,ny-2:2*ny-2,:])
        dy_f += cc_6[der_ord] *(ff[:,3:3+ny,:] - ff[:,ny-3:2*ny-3,:])
      else :
        f = self.get_data(tar_var)  
        dy_f = np.zeros([nx,ny,nz])  
        dy_f[:,3:ny-3,:]  = oo_6[der_ord] * f[:,3:ny-3,:]
        dy_f[:,3:ny-3,:] += aa_6[der_ord] *(f[:,4:ny-2,:] - f[:,2:ny-4,:])
        dy_f[:,3:ny-3,:] += bb_6[der_ord] *(f[:,5:ny-1,:] - f[:,1:ny-5,:])
        dy_f[:,3:ny-3,:] += cc_6[der_ord] *(f[:,6:ny-0,:] - f[:,0:ny-6,:])

      #print('done with the new array!')
      return dy_f

  #------------------------------------------------------------  
  def calc_gradz(self,
      tar_var,
      periodic = None,
      der_ord = 1,
      brutal = False):
    """ 
    Calculates z gradient component with(out) periodical boundary conditions
    
    Parameters :
      - tar_var            [fibo.data] target variable of the procedure
      - periodic = True    [bool] periodic boundary conditions? if None will consult meta
      - der_ord = 1        [int] order of derivation (0: no derivative, 1: first derivative, 2: second derivative ...)
      - brutal = False     [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - zgr_var            [fibo.data]
    
    """

    #determine the coefficients
    nx,ny,nz = self.meta['nnn']
    if periodic is None : periodic = self.meta['ppp'][2]

    if nz == 1 :
      return np.zeros([nx,ny,nz])

    else :

      # set the coefficients !
      if brutal : 
        oo_6 = self.meta['dz']**(-der_ord) * np.array([1.0,  0.0, -2.0])      
        aa_6 = self.meta['dz']**(-der_ord) * np.array([0.0,  0.5, 1.0 ]) 
        bb_6 = self.meta['dz']**(-der_ord) * np.array([0.0,  0.0, 0.0]) 
        cc_6 = self.meta['dz']**(-der_ord) * np.array([0.0,  0.0, 0.0]) 
      else : 
        oo_6 = self.meta['dz']**(-der_ord) * np.array([1.0,  0.0, -49./18.])      
        aa_6 = self.meta['dz']**(-der_ord) * np.array([0.0,  9.0/12.0,  3./2. ]) 
        bb_6 = self.meta['dz']**(-der_ord) * np.array([0.0, -3.0/20.0, -3./20.]) 
        cc_6 = self.meta['dz']**(-der_ord) * np.array([0.0,  1.0/60.0,  1./90.]) 
      
      #create the new vector, fill it
      if periodic :
        ff = np.tile(self.get_data(tar_var),(1,1,2))
        dz_f  = oo_6[der_ord] * ff[:,:,0:nz] 
        dz_f += aa_6[der_ord] *(ff[:,:,1:1+nz] - ff[:,:,nz-1:2*nz-1])
        dz_f += bb_6[der_ord] *(ff[:,:,2:2+nz] - ff[:,:,nz-2:2*nz-2])
        dz_f += cc_6[der_ord] *(ff[:,:,3:3+nz] - ff[:,:,nz-3:2*nz-3])  
      else:
        f = self.get_data(tar_var)  
        dz_f  = np.zeros([nx,ny,nz])  
        dz_f[:,:,3:nz-3]  = oo_6[der_ord] * f[:,:,3:nz-3]
        dz_f[:,:,3:nz-3] += aa_6[der_ord] *(f[:,:,4:nz-2] - f[:,:,2:nz-4])
        dz_f[:,:,3:nz-3] += bb_6[der_ord] *(f[:,:,5:nz-1] - f[:,:,1:nz-5])
        dz_f[:,:,3:nz-3] += cc_6[der_ord] *(f[:,:,6:nz-0] - f[:,:,0:nz-6])

      #print('done with the new array!')
      return dz_f

  #------------------------------------------------------------
  def calc_divr(self,
      tar_var_x,
      tar_var_y,
      tar_var_z,
      perx = None,
      pery = None,
      perz = None,
      brutal = False):      
    """ 
    Calculates the divergence
    
    Parameters :
      - tar_var_x    [fibo.data] target field, x component
      - tar_var_y    [fibo.data] target field, y component
      - tar_var_z    [fibo.data] target field, z component
      - perx = None  [None OR bool] x periodic? (not considered if nx==1)
      - pery = None  [None OR bool] y periodic? (not considered if ny==1)
      - perz = None  [None OR bool] z periodic? (not considered if nz==1)     
      - brutal = False     [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - divr_var     [fibo.data]
    
    """

    #create a raw vector 
    nx,ny,nz = self.meta['nnn'] 
    divr = np.zeros([nx,ny,nz])

    if nx != 1 : divr += self.calc_gradx(tar_var_x,periodic=perx,brutal=brutal)
    if ny != 1 : divr += self.calc_grady(tar_var_y,periodic=pery,brutal=brutal)
    if nz != 1 : divr += self.calc_gradz(tar_var_z,periodic=perz,brutal=brutal)

    return divr

  #------------------------------------------------------------
  def calc_curl(self,  
      tar_var_x,
      tar_var_y,
      tar_var_z,
      perx = None,
      pery = None,
      perz = None,
      brutal = False):   
    """ 
    Calculates the curl
    
    Parameters :
      - tar_var_x    [fibo.data] target field, x component
      - tar_var_y    [fibo.data] target field, y component
      - tar_var_z    [fibo.data] target field, z component
      - perx = None  [None OR bool] x periodic? (not considered if nx==1)
      - pery = None  [None OR bool] y periodic? (not considered if ny==1)
      - perz = None  [None OR bool] z periodic? (not considered if nz==1)
      - brutal = False     [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - curl_var     [fibo.data]
    
    """

    #create a raw vector field
    nx,ny,nz = self.meta['nnn']
    curl_x = np.zeros([nx,ny,nz])
    curl_y = np.zeros([nx,ny,nz])
    curl_z = np.zeros([nx,ny,nz])

    if nx != 1 : 
      curl_y -= self.calc_gradx(tar_var_z,periodic=perx,brutal=brutal)
      curl_z += self.calc_gradx(tar_var_y,periodic=perx,brutal=brutal)
    if ny != 1 : 
      curl_z -= self.calc_grady(tar_var_x,periodic=pery,brutal=brutal)
      curl_x += self.calc_grady(tar_var_z,periodic=pery,brutal=brutal)
    if nz != 1 : 
      curl_x -= self.calc_gradz(tar_var_y,periodic=perz,brutal=brutal)
      curl_y += self.calc_gradz(tar_var_x,periodic=perz,brutal=brutal)

    return curl_x, curl_y, curl_z

  #------------------------------------------------------------
  def calc_lapl(self,  
      tar_var_x,
      tar_var_y,
      tar_var_z,
      perx = None,
      pery = None,
      perz = None,
      brutal = False):     
    """ 
    Calculates the laplacian
    
    Parameters :
      - tar_var_x    [fibo.data] target field, x component
      - tar_var_y    [fibo.data] target field, y component
      - tar_var_z    [fibo.data] target field, z component
      - perx = None  [None OR bool] x periodic? (not considered if nx==1)
      - pery = None  [None OR bool] y periodic? (not considered if ny==1)
      - perz = None  [None OR bool] z periodic? (not considered if nz==1)
      - brutal = False     [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - lapl_var     [fibo.data]
    
    """
    #create a raw vector
    nx,ny,nz = self.meta['nnn']
    lapl = np.zeros([nx,ny,nz])

    if nx != 1 : lapl += self.calc_gradx(tar_var_x,periodic=perx,der_ord=2,brutal=brutal)
    if ny != 1 : lapl += self.calc_grady(tar_var_y,periodic=pery,der_ord=2,brutal=brutal)
    if nz != 1 : lapl += self.calc_gradz(tar_var_z,periodic=perz,der_ord=2,brutal=brutal)

    return lapl

  #------------------------------------------------------------  
  def calc_len_scal(self,
      tar_data,
      relative = True,
      normal_x = None,
      normal_y = None,
      normal_z = None,
      brutal = False):
    """
    Calculates scale lengths for a scalar field
    
    Parameters :
      - tar_data             [fibo.data] target scalar field
      - relative = True      [bool] relative to the local field intensity?
      - normal_x = None      [None OR fibo.data] of local evaluation plane: x component
      - normal_y = None      [None OR fibo.data] of local evaluation plane: y component
      - normal_z = None      [None OR fibo.data] of local evaluation plane: z component
      - brutal = False       [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - tar_length          [fibo.data]
      - tar_versor_x      [fibo.data]
      - tar_versor_y      [fibo.data]
      - tar_versor_z      [fibo.data]
    
    """
    
    nx,ny,nz = self.meta['nnn']
    
    grad = np.zeros([nx,ny,nz,3])
    grad[:,:,:,0] = self.calc_gradx(tar_data,brutal=brutal)
    grad[:,:,:,1] = self.calc_grady(tar_data,brutal=brutal)
    grad[:,:,:,2] = self.calc_gradz(tar_data,brutal=brutal)

    if not (normal_x is None) : 
      #normalize the plane normal versor
      norm = self.calc_scalr(normal_x,normal_x,normal_y,normal_y,normal_z,normal_z)
      norm = np.sqrt(np.reciprocal( norm ))
      normal_x,normal_y,normal_z = normal_x*norm,normal_y*norm,normal_z*norm

      #create the projecting tensor and projec
      plan = np.transpose(np.array([normal_x,normal_y,normal_z]),(1,2,3,0))
      for ix in range(nx):
        for iy in range(ny):
          for iz in range(nz): 
            proj = np.eye(3) - np.outer(plan[ix,iy,iz,:],plan[ix,iy,iz,:])
            grad[ix,iy,iz,:] = np.tensordot(proj,grad[ix,iy,iz,:])

    #calculate the scale length
    length = self.calc_scalr(grad[:,:,:,0],grad[:,:,:,0],grad[:,:,:,1],grad[:,:,:,1],grad[:,:,:,2],grad[:,:,:,2])
    length = np.sqrt(np.reciprocal( length ))
    
    #now the versors 
    grad[:,:,:,0],grad[:,:,:,1],grad[:,:,:,2] = grad[:,:,:,0]*length,grad[:,:,:,1]*length,grad[:,:,:,2]*length

    if relative : length *= self.data[tar_data]

    return length, grad[:,:,:,0],grad[:,:,:,1],grad[:,:,:,2]

  #------------------------------------------------------------  
  def calc_len_vect(self,
      tar_data_x,
      tar_data_y,
      tar_data_z,
      relative = True,
      normal_x = None,
      normal_y = None,
      normal_z = None,
      brutal = False):
    """
    Calculates scale lengths for a vector field
    
    Parameters :
      - tar_data_x           [fibo.data] target vector field: x component
      - tar_data_y           [fibo.data] target vector field: y component
      - tar_data_z           [fibo.data] target vector field: z component
      - relative = True      [bool] relative to the local field intensity?
      - normal_x = None      [None OR fibo.data] of local evaluation plane: x component
      - normal_y = None      [None OR fibo.data] of local evaluation plane: y component
      - normal_z = None      [None OR fibo.data] of local evaluation plane: z component
      - brutal = False       [bool] do you want to revert to a brutal calculation?
    
    Returns :
      - tar_length_A           [fibo.data]
      - tar_length_B           [fibo.data]
      - tar_length_C           [fibo.data]
      - tar_versor_A_x         [fibo.data]
      - tar_versor_A_y         [fibo.data]
      - tar_versor_A_z         [fibo.data]
      - tar_versor_B_x         [fibo.data]
      - tar_versor_B_y         [fibo.data]
      - tar_versor_B_z         [fibo.data]
      - tar_versor_C_x         [fibo.data]
      - tar_versor_C_y         [fibo.data]
      - tar_versor_C_z         [fibo.data]
    
    """

    nx,ny,nz = self.meta['nnn']
    
    grad = np.zeros([nx,ny,nz,3,3])
    grad[:,:,:,0,0] = self.calc_gradx(tar_data_x,brutal=brutal)
    grad[:,:,:,1,0] = self.calc_grady(tar_data_x,brutal=brutal)
    grad[:,:,:,2,0] = self.calc_gradz(tar_data_x,brutal=brutal)
    print('done the gradient of x component!')
    grad[:,:,:,0,1] = self.calc_gradx(tar_data_y,brutal=brutal)
    grad[:,:,:,1,1] = self.calc_grady(tar_data_y,brutal=brutal)
    grad[:,:,:,2,1] = self.calc_gradz(tar_data_y,brutal=brutal)
    print('done the gradient of y component!')
    grad[:,:,:,0,2] = self.calc_gradx(tar_data_z,brutal=brutal)
    grad[:,:,:,1,2] = self.calc_grady(tar_data_z,brutal=brutal)
    grad[:,:,:,2,2] = self.calc_gradz(tar_data_z,brutal=brutal)
    print('done the gradient of z component!')
    
    if len(list({normal_x, normal_y, normal_z})) != 1 : 
      #normalize the plane normal versor
      norm = self.calc_scalr(normal_x,normal_x,normal_y,normal_y,normal_z,normal_z)
      norm = np.sqrt(np.reciprocal( norm ))
      normal_x,normal_y,normal_z = normal_x*norm,normal_y*norm,normal_z*norm

      #create the projecting tensor and projec
      plan = np.transpose(np.array([normal_x,normal_y,normal_z]),(1,2,3,0))
      for ix in range(nx):
        for iy in range(ny):
          for iz in range(nz): 
            proj = np.eye(3) - np.outer(plan[ix,iy,iz,:],plan[ix,iy,iz,:])
            grad[ix,iy,iz,:,:] = np.dot(proj,grad[ix,iy,iz,:,:])
    
    lengthA, lengthB, lengthC = np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz])
    vecA_x, vecA_y, vecA_z = np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz])
    vecB_x, vecB_y, vecB_z = np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz])
    vecC_x, vecC_y, vecC_z = np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz]), np.zeros([nx,ny,nz])

    vals = np.zeros([3])
    vecs = np.zeros([3,3])
    
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          vals,vecs = np.linalg.eigh(np.dot(grad[ix,iy,iz,:,:],np.transpose(grad[ix,iy,iz,:,:])))
          #vals = np.reciprocal(np.sqrt(vals))
          lengthC[ix,iy,iz],lengthB[ix,iy,iz],lengthA[ix,iy,iz] = vals[0],vals[1],vals[2]
          vecC_x[ix,iy,iz],vecB_x[ix,iy,iz],vecA_x[ix,iy,iz] = vecs[0,0],vecs[0,1],vecs[0,2]
          vecC_y[ix,iy,iz],vecB_y[ix,iy,iz],vecA_y[ix,iy,iz] = vecs[1,0],vecs[1,1],vecs[1,2]
          vecC_z[ix,iy,iz],vecB_z[ix,iy,iz],vecA_z[ix,iy,iz] = vecs[2,0],vecs[2,1],vecs[2,2]
    lengthA = np.reciprocal(np.sqrt(lengthA))
    lengthB = np.reciprocal(np.sqrt(lengthB))
    lengthC = np.reciprocal(np.sqrt(lengthC))
    print('done the values and vectors!')
    
    if relative : 
      tar_data_norm = np.sqrt(self.calc_scalr(tar_data_x,tar_data_x,tar_data_y,tar_data_y,tar_data_z,tar_data_z))
      lengthA, lengthB, lengthC = tar_data_norm*lengthA, tar_data_norm*lengthB, tar_data_norm*lengthC

    return lengthA,lengthB,lengthC, vecA_x,vecA_y,vecA_z, vecB_x,vecB_y,vecB_z, vecC_x,vecC_y,vecC_z

  #------------------------------------------------------------  
  def calc_spots(self,
      tar_spots,
      bin_struc = np.ones((3,3,3))):
    """
    Finds connected spots (also with periodic boundaries)
    
    Parameters :
      - tar_spots    [fibo.data] binary mask individuating the spots
      - bin_struc    [bool array] binary structure to characterize the connectivity
    
    Returns :
      - spot_mask    [fibo.data] positive integers individuating the connected spots
      - spot_number  [int>0] nuber of different labels 
      - spot_tiling  [(3,spot_number) int>0 array] minimum tiling to see each connected spot whole
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
    spot_tiling = np.zeros((3,spot_number+1),'int')
    
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
    
    if dim_x: 
      relabel = np.arange(1,spot_number+1)
      # check label correspondences via the other periodic dimensions 
      if dim_y: corr_check(relabel,np.array([spot_mask[:,:2,:],spot_mask[:,-2:,:]]),False)
      if dim_z: corr_check(relabel,np.array([spot_mask[:,:,:2],spot_mask[:,:,-2:]]),dim_y)
      # extract the edges and collapse corresponding labels
      edg = np.array([spot_mask[:2,:,:],spot_mask[-2:,:,:]])
      for k in range(len(relabel)) :
        kk = relabel[k]
        if kk != k+1 : 
          edg[edg == k+1] = kk
      # find out which labels correspond to coast-to-coast crossing spots
      ctc_labels = np.array(list( set(edg[0].flatten()).intersection(set(edg[1].flatten())) ))
      inf_labels = np.unique(edg[0,edg[0]==edg[1]])
      spot_tiling[0,ctc_labels] = 1
      spot_tiling[0,inf_labels] = -spot_number
      # now you may complete the re-labelling :)
      corr_check(relabel,edg,False)
    
    if dim_y: 
      relabel = np.arange(1,spot_number+1)
      # check label correspondences via the other periodic dimensions 
      if dim_z: corr_check(relabel,np.array([spot_mask[:,:,:2],spot_mask[:,:,-2:]]),False)
      if dim_x: corr_check(relabel,np.array([spot_mask[:2,:,:],spot_mask[-2:,:,:]]),dim_z)
      # extract the edges and collapse corresponding labels
      edg = np.array([spot_mask[:,:2,:],spot_mask[:,-2:,:]])
      for k in range(len(relabel)) :
        kk = relabel[k]
        if kk != k+1 : 
          edg[edg == k+1] = kk
      # find out which labels correspond to coast-to-coast crossing spots
      ctc_labels = np.array(list( set(edg[0].flatten()).intersection(set(edg[1].flatten())) ))
      inf_labels = np.unique(edg[0,edg[0]==edg[1]])
      spot_tiling[1,ctc_labels] = 1
      spot_tiling[1,inf_labels] = -spot_number
      # now you may complete the re-labelling :)
      corr_check(relabel,edg,False)
    
    if dim_z: 
      relabel = np.arange(1,spot_number+1)
      # check label correspondences via the other periodic dimensions 
      if dim_x: corr_check(relabel,np.array([spot_mask[:2,:,:],spot_mask[-2:,:,:]]),False)
      if dim_y: corr_check(relabel,np.array([spot_mask[:,:2,:],spot_mask[:,-2:,:]]),dim_x)
      # extract the edges and collapse corresponding labels
      edg = np.array([spot_mask[:,:,:2],spot_mask[:,:,-2:]])
      for k in range(len(relabel)) :
        kk = relabel[k]
        if kk != k+1 : 
          edg[edg == k+1] = kk
      # find out which labels correspond to coast-to-coast crossing spots
      ctc_labels = np.array(list( set(edg[0].flatten()).intersection(set(edg[1].flatten())) ))
      inf_labels = np.unique(edg[0,edg[0]==edg[1]])
      spot_tiling[2,ctc_labels] = 1
      spot_tiling[2,inf_labels] = -spot_number
      # now you may complete the re-labelling :)
      corr_check(relabel,edg,False)
    
    # finally you can trim down your spot_mask
    if dim_x : spot_mask = spot_mask[1:-1,:,:]
    if dim_y : spot_mask = spot_mask[:,1:-1,:]
    if dim_z : spot_mask = spot_mask[:,:,1:-1]
    
    # we're at the end: just re-name the labels now!
    old_labels = np.unique(relabel)
    spot_number = len(old_labels)
    
    for kk in range(1,spot_number+1) :
      if kk not in old_labels : 
        k = old_labels[0]
        loc_k = np.where(relabel == k)
        relabel[loc_k] = kk    
        #print('changed label '+str(k)+' to '+str(kk))
      old_labels = old_labels[1:]
    
    # and change the labels everywhere! 
    for k in range(len(relabel)) :
      kk = relabel[k]
      if kk != k+1 : 
        spot_mask[spot_mask == k+1] = kk
        spot_tiling[:,kk] += spot_tiling[:,k+1]
        spot_tiling[:,k+1] = np.zeros((3),'int')
        #print('relabel '+str(k+1)+' to '+str(kk))
      
    # trim the tiling indices matrix and add the null tile factor 
    spot_tiling  = spot_tiling[:,1:spot_number+1] 
    spot_tiling += np.ones((3,spot_number),'int')
    spot_tiling  = np.maximum(spot_tiling, np.zeros((3,spot_number),'int'))

    return spot_mask, spot_number, spot_tiling


  #------------------------------------------------------------  
  #--JOB:function-added-to-change-units-of-output-iPIC3D-------  
  #------------------------------------------------------------  
  def calc_units(self,
      units,
      fibo_obj = None,
      silent = True):
    """
    Change units for all data stored in fibo_obj
    
    Parameters :
      - units            [str] specify the new units you want to implement (assume to start from code output)
      - fibo_obj = None  [None or fibo] fibo object you want to fill, else return normalization factors
    
    Returns :
      - norm_factors     [float,float,float,float,float] factors to change units of B, E, J, P, rho
    """  
    FourPi = 4.*np.pi
    me     = np.absolute(self.meta['sQOM'][0])
    mi     = np.absolute(self.meta['sQOM'][1])
    qe     = np.sign(self.meta['sQOM'][0])
    qi     = np.sign(self.meta['sQOM'][1])
    Bsw    = np.sqrt(self.meta['Bx0']**2+self.meta['By0']**2+self.meta['Bz0']**2)
    Vsw    = np.sqrt(self.meta['Vx0']**2+self.meta['Vy0']**2+self.meta['Vz0']**2) 
    Esw    = Bsw*Vsw      #NB: this is not really the most general case...

    if mi!=1. or qi!=1 : print('ERROR: qi or mi are not 1 !!!')

    if units == 'iPIC'  : norm_fact = np.array([1.,     1.,     1.,         1.,            1./FourPi]) 
    elif units == 'SW'  : norm_fact = np.array([1./Bsw, 1./Esw, FourPi/Vsw, FourPi/Vsw**2, 1.])
    #elif units == 'SI' : norm_fact = np.array()
    else :
      print('ERROR: arg units ='+units+' not recognized (choose iPIC,SW,SI)')

    if (fibo_obj != None) :
      for key in fibo_obj.data.keys() :
          if 'B'        in key : fibo_obj.data[key] = fibo_obj.data[key]*norm_fact[0]
          elif 'E'   in key : fibo_obj.data[key] = fibo_obj.data[key]*norm_fact[1]
          elif 'J'   in key : fibo_obj.data[key] = fibo_obj.data[key]*norm_fact[2]
          elif 'P'   in key : fibo_obj.data[key] = fibo_obj.data[key]*norm_fact[3]
          elif 'rho' in key : fibo_obj.data[key] = fibo_obj.data[key]*norm_fact[4]
          else: print('ERROR: there is a key='+key+' in data.keys() that does not match!!!')

    if not silent:
      print('calc_units>units         :',units)
      print('calc_units>B_factor      :',norm_fact[0])
      print('calc_units>E_factor      :',norm_fact[1])
      print('calc_units>J_factor      :',norm_fact[2])
      print('calc_units>P_factor      :',norm_fact[3])
      print('calc_units>rho_factor    :',norm_fact[4])

    return norm_fact

  #------------------------------------------------------------  
  #--JOB:func-to-divide-vector-by-scalar-----------------------  
  #------------------------------------------------------------  
  def calc_divid(self,
      tar_var_vect,
      tar_var_scal,
      new_tar,
      seg = None,
      fill = True):

    """
    compute vector / scalar and creates a new fibo.data

    Parameters :
      - tar_var_vect            [str] vector we want to divide (e.g. 'Ji')
      - tar_var_scal            [str] scalar field that divides vect
      - new_tar                 [str] name of new fibo.data created
      - seg = None              [None or str] segment to use, if None use simply tar_var_vect, tar_var_scal
      - fill = True             [bool] do you want to add this new array to the class self?

    Returns :
      - divided numpy array
    """
    tar_var_vect_x = tar_var_vect+'_x'
    tar_var_vect_y = tar_var_vect+'_y'
    tar_var_vect_z = tar_var_vect+'_z'

    if (seg != None) :
      tar_var_vect_x = tar_var_vect_x+'%.8i'%int(seg)
      tar_var_vect_y = tar_var_vect_y+'%.8i'%int(seg)
      tar_var_vect_z = tar_var_vect_z+'%.8i'%int(seg)
      tar_var_scal   = tar_var_scal  +'%.8i'%int(seg)

    tar_arr_x  = np.divide(self.get_data(tar_var_vect_x),self.get_data(tar_var_scal))
    tar_arr_y  = np.divide(self.get_data(tar_var_vect_y),self.get_data(tar_var_scal))
    tar_arr_z  = np.divide(self.get_data(tar_var_vect_z),self.get_data(tar_var_scal))

    if fill: 
      self.data[new_tar+'_x'+'%.8i'%int(seg)] = tar_arr_x
      self.data[new_tar+'_y'+'%.8i'%int(seg)] = tar_arr_y
      self.data[new_tar+'_z'+'%.8i'%int(seg)] = tar_arr_z
    else: 
      return np.array([tar_arr_x,tar_arr_y,tar_arr_z])


