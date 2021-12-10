import collections 
import codecs

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

#---------------------------------------------------------------------------------------
#------fill-fibo-objects-from-various-sources-------------------------------------------
#-------or-simply-get-your-data---------------------------------------------------------
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
class from_VTK (object):

  def __init__(self, 
      address):
    """ 
    Creates the object to retrieve data from VTK files
    
    Parameters :
      - address      [address] where your data are (folder with segs inside)
    
    """
    self.address = address
    self.meta = {}

  #------------------------------------------------------------
  def get_meta(self,  #counts lines in file and calls the appropriate function to get the meta data
      extra_address = '',
      silent = True):
    """
    ------------------------------------------------------------------------------------
      fills the metadata list 
    ------------------------------------------------------------------------------------
    extra_address = ''      [address] to reach any subfolder where your meta-data are
    silent        = True    [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    
    with open(os.path.join(self.address,extra_address,'SimulationData.txt'),'r') as foo:
      line_number = len(foo.readlines())

    if line_number==35 : old_vers=True
    elif line_number>35 : old_vers=False
    
    self.get_meta_A(old_vers,extra_address)

    #---------get-dimensions-of-your-simulation----------

    self.meta['dx'] = self.meta['xl']/self.meta['nx']
    self.meta['dy'] = self.meta['yl']/self.meta['ny']
    self.meta['dz'] = self.meta['zl']/self.meta['nz']

    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz']) 

    self.meta['ppp'] = (False, False, False)       # HARD-CODED !!! 

    self.meta['x'] = np.arange(0.,self.meta['xl'],self.meta['dx'])
    try:
        self.meta['y'] = np.arange(0.,self.meta['yl'],self.meta['dy'])
    except:
        self.meta['y'] = np.array([0.])
    try:
        self.meta['z'] = np.arange(0.,self.meta['zl'],self.meta['dz'])
    except:
        self.meta['z'] = np.array([0.])


    #----------get-time-infos-from all-vtk-files----------------- 

    segments = [f for f in os.listdir(self.address) if f.split('.')[-1]=='vtk']
    for i in range(len(segments)):
        if i == 0:
          self.meta['name'] = segments[i].split('_')[0]
        segments[i] = segments[i].split('_')[-1].split('.')[0]
    segments = set(segments)
    segments = map(str, sorted(map(int, segments)))
    self.meta['segcycles']=[]
    self.meta['segtimes']=[]
    for seg in segments:
      self.meta['segcycles'].append(int(seg))
      self.meta['segtimes'].append(float(seg)*self.meta['dt'])
      

    #----------add-informations-on-species-----------------

    species  = []
    for isp in range(0,self.meta['nss']):
      if self.meta['sQOM'][isp]<0  : species.append('e'+str(isp))
      elif self.meta['sQOM'][isp]>0  : species.append('i'+str(isp))

    self.meta['species']  = species

    if self.meta['ny'] == 1:
      self.meta['space_dim'] = '1D'
    elif self.meta['nz'] == 1:
      self.meta['space_dim'] = '2D'
    else:
      self.meta['space_dim'] = '3D'

    #----------print-summary-----------------

    if not silent : 
      print('iPIC3D> cell number               :  ', self.meta['nnn'])
      print('iPIC3D> domain size               :  ', self.meta['lll'])
      print('iPIC3D> mesh spacing              :  ', self.meta['ddd'])
      print('iPIC3D> periodicity               :  ', self.meta['ppp'])
      print('iPIC3D> time step                 :  ', self.meta['dt'])
      print('iPIC3D> species                   :  ', self.meta['species'])
      for i in range(self.meta['nss']):
        print('          '+species[i]+' charge-over-mass           :  ', self.meta['sQOM'][i])

  #------------------------------------------------------------
  def get_meta_A(self,
      old_vers=False,
      extra_address = ''):
    """
    ------------------------------------------------------------------------------------
      extra routine, reads meta data from SimulationData.txt
    ------------------------------------------------------------------------------------
      old_vers      = False   [bool] is your simu older than 11/2021? means that I have changed the simulationdata.txt in ipic3d
      extra_address = ''      [address] to reach any subfolder where your meta-data are
    ------------------------------------------------------------------------------------
    """

    #get mesh infos from SimulationData.txt
    infos = open(os.path.join(self.address,extra_address,'SimulationData.txt'),'r')

    infos.readline() #---------------------------
    infos.readline() #-  Simulation Parameters  -
    infos.readline() #---------------------------
    self.meta['nss'] = int(infos.readline().split('=')[-1]) #number of species
    stag=[]
    sQOM=[]
    for i in range(self.meta['nss']):
      sQOM.append(float(infos.readline().split('=')[-1]))
    self.meta['sQOM'] = sQOM
    infos.readline() #---------------------------
    self.meta['xl'] = float(infos.readline().split('=')[-1]) #box physical dimensions
    self.meta['yl'] = float(infos.readline().split('=')[-1])
    self.meta['zl'] = float(infos.readline().split('=')[-1])
    self.meta['nx'] = int(infos.readline().split('=')[-1]) #box grid dimensions
    self.meta['ny'] = int(infos.readline().split('=')[-1])
    self.meta['nz'] = int(infos.readline().split('=')[-1])
    if not old_vers :
      infos.readline() #---------------------------
      self.meta['xc'] = float(infos.readline().split('=')[-1]) #planet position
      self.meta['yc'] = float(infos.readline().split('=')[-1])
      self.meta['zc'] = float(infos.readline().split('=')[-1])
      self.meta['R'] = float(infos.readline().split('=')[-1]) #planet radius
      self.meta['Doff'] = float(infos.readline().split('=')[-1])
      infos.readline() #---------------------------
      self.meta['SAL'] = int(infos.readline().split('=')[-1])
      self.meta['Nsal'] = int(infos.readline().split('=')[-1])
    infos.readline() #---------------------------
    self.meta['dt'] = float(infos.readline().split('=')[-1]) #timestep
    self.meta['nsteps'] = int(infos.readline().split('=')[-1]) #number of steps
    infos.readline() #---------------------------
    for i in range(self.meta['nss']):
      infos.readline() #rho init species
      infos.readline() #rho inject species
    infos.readline() #current sheet thickness
    self.meta['Bx0'] = float(infos.readline().split('=')[-1]) #Bx0
    self.meta['By0'] = float(infos.readline().split('=')[-1]) #By0
    self.meta['Bz0'] = float(infos.readline().split('=')[-1]) #Bz0
    if not old_vers :
      infos.readline() #---------------------------
      self.meta['Vx0'] = float(infos.readline().split('=')[-1]) #Vx0
      self.meta['Vy0'] = float(infos.readline().split('=')[-1]) #Vy0
      self.meta['Vz0'] = float(infos.readline().split('=')[-1]) #Vz0
      self.meta['vths'] = []
      for i in range(self.meta['nss']):
        self.meta['vths'].append(float(infos.readline().split('=')[-1])) #vth species
    infos.readline() #---------------------------
    infos.readline() #Smooth
    infos.readline() #2D smoothing
    infos.readline() #nvolte ?
    infos.readline() #GMRES error tolerance
    infos.readline() #CG error toletance
    infos.readline() # Mover error tolerance

    infos.readline() #---------------------------
    infos.readline() #Results saved in:
    infos.readline() #Restart saved in:

    infos.readline() #---------------------------
    infos.close()

  #------------------------------------------------------------
  def get_scal(self,
      tar_file,
      seg,
      fibo_obj = None,  
      tar_var = None,  
      double_y = False,
      silent = True):
    """ 
    Reads scalar from .vtk file
    
    Parameters :
      - tar_file           [str] target file to read (don't include '.vtk')
      - seg                [str] cycle of the simulation
      - fibo_obj = None    [None or fibo] fibo object you want to fill, else returns values 
      - tar_var = None     [None or str] name the.variable will be given
      - double_y = False   [bool] was your file printed twice in y?
      - silent = True      [bool] print status at the end?
    Returns :
      - scal               [fibo_var] 
    
    """  

    #create data vector, fill it!
    data_file = open(os.path.join(self.address,tar_file+'.vtk'),'r')

    if tar_var == None : 
      tar_var = data_file.readline().split()[0]+'%.8i'%int(seg)
    else : 
      tar_var = tar_var+'%.8i'%int(seg)
      data_file.readline()

    data_file.readline()
    data_format = data_file.readline()
    data_structure = data_file.readline().split()[1]
    self.meta['nx'], self.meta['ny'], self.meta['nz'] = map(int, data_file.readline().split()[1:4])
    data_file.readline()
    self.meta['dx'], self.meta['dy'], self.meta['dz'] = map(float, data_file.readline().split()[1:4])
    data_file.readline()
    data_file.readline()  #NB here you have the nx*ny*nz preduct
    data_file.readline()
    data_file.readline()

    data_file.close()

    if double_y : self.meta['ny'] = self.meta['ny']/2 #NB here you divide by two the box in y!

    if data_structure == 'STRUCTURED_POINTS': reader = vtk.vtkStructuredPointsReader() #here you can add other readers in case
    
    reader.SetFileName(os.path.join(self.address,tar_file+'.vtk'))
    reader.ReadAllScalarsOn()
    reader.Update()
    vtk_output = reader.GetOutput()

    if vtk_output.GetDimensions()[0] != self.meta['nx'] : print('ERROR: wrong number of cells along x (Nx)')
    if vtk_output.GetDimensions()[2] != self.meta['nz'] : print('ERROR: wrong number of cells along z (Nz)')
    if not double_y and vtk_output.GetDimensions()[1] != self.meta['ny'] :  print('ERROR: wrong number of cells along y (Ny) ; double_y=False')
    if double_y and vtk_output.GetDimensions()[1] != self.meta['ny']*2 :    print('ERROR: wrong number of cells along y (Ny) ; double_y=True')
        
    scal = VN.vtk_to_numpy(vtk_output.GetPointData().GetScalars())

    if double_y :     scal = scal.reshape(self.meta['nz'],2*self.meta['ny'],self.meta['nx']).transpose(2,1,0) #recast flatten array to 3D array
    else :            scal = scal.reshape(self.meta['nz'],self.meta['ny'],self.meta['nx']) .transpose(2,1,0)

    if double_y : scal = scal[:,:self.meta['ny'],:]

    if (fibo_obj != None) :
      fibo_obj.data[tar_var] = scal
    else: 
      return scal

    if not silent:
      print('get_scal_from_VTK> data format             :  ', data_format)
      print('get_scal_from_VTK> data structure          :  ', data_structure)
      print('get_scal_from_VTK> grid dimensions         :  ', self.meta['nnn'])
      print('get_scal_from_VTK> grid size               :  ', self.meta['lll'])
      print('get_scal_from_VTK> grid spacing            :  ', self.meta['ddd'])
      if (fibo_obj != None) :
        print('get_scal_from_VTK> created fibo_obj.data['+tar_var+']')
  #------------------------------------------------------------
  def get_vect(self,
      tar_file,
      seg,
      fibo_obj = None,
      tar_var = None,
      double_y = False,
      silent=True):
    """ 
    Reads vector from .vtk file
    
    Parameters :
      - tar_file           [str] target file to read (don't include '.vtk')
      - fibo_obj = None    [None or fibo] fibo object you want to fill, else returns values 
      - tar_var = None     [None or str] name the.variable will be given
      - double_y = False   [bool] was your file printed twice in y?
    
    Returns :
      - scal               [fibo_var] 
    
    """

    #create data vector, fill it!
    data_file = open(os.path.join(self.address,tar_file+'.vtk'),'r')

    if tar_var == None : 
      tar_var_x,tar_var_y,tar_var_z = data_file.readline().split()[0][1:-1].split(',')
      tar_var_x = tar_var_x+'%.8i'%int(seg)
      tar_var_y = tar_var_y+'%.8i'%int(seg)
      tar_var_z = tar_var_z+'%.8i'%int(seg)
    else : 
      tar_var_x = tar_var+'_x'+'%.8i'%int(seg)
      tar_var_y = tar_var+'_y'+'%.8i'%int(seg)
      tar_var_z = tar_var+'_z'+'%.8i'%int(seg)
      data_file.readline()

    data_file.readline()
    data_format = data_file.readline()
    data_structure = data_file.readline().split()[1]
    self.meta['nx'], self.meta['ny'], self.meta['nz'] = map(int, data_file.readline().split()[1:4])
    data_file.readline()
    self.meta['dx'], self.meta['dy'], self.meta['dz'] = map(float, data_file.readline().split()[1:4])
    data_file.readline()
    data_file.readline()  #NB here you have the nx*ny*nz preduct
    data_file.readline()
    
    data_file.close()

    if double_y : self.meta['ny'] = self.meta['ny']/2  #NB here you divide by two the box in y!

    if data_structure == 'STRUCTURED_POINTS': reader = vtk.vtkStructuredPointsReader() #here you can add other readers in case
    
    reader.SetFileName(os.path.join(self.address,tar_file+'.vtk'))
    reader.ReadAllVectorsOn()
    reader.Update()
    vtk_output = reader.GetOutput()

    if vtk_output.GetDimensions()[0] != self.meta['nx'] : print('ERROR: wrong number of cells along x (Nx)')
    if vtk_output.GetDimensions()[2] != self.meta['nz'] : print('ERROR: wrong number of cells along z (Nz)')
    if not double_y and vtk_output.GetDimensions()[1] != self.meta['ny'] :  print('ERROR: wrong number of cells along y (Ny) ; double_y=False')
    if double_y and vtk_output.GetDimensions()[1] != self.meta['ny']*2 :    print('ERROR: wrong number of cells along y (Ny) ; double_y=True')
        
    vect = VN.vtk_to_numpy(vtk_output.GetPointData().GetArray(tar_var))

    if double_y :
      vect_x = vect[:,0].reshape(self.meta['nz'],self.meta['ny']*2,self.meta['nx']).transpose(2,1,0)
      vect_y = vect[:,1].reshape(self.meta['nz'],self.meta['ny']*2,self.meta['nx']).transpose(2,1,0)
      vect_z = vect[:,2].reshape(self.meta['nz'],self.meta['ny']*2,self.meta['nx']).transpose(2,1,0)
    else :
      vect_x = vect[:,0].reshape(self.meta['nz'],self.meta['ny'],self.meta['nx']).transpose(2,1,0)
      vect_y = vect[:,1].reshape(self.meta['nz'],self.meta['ny'],self.meta['nx']).transpose(2,1,0)
      vect_z = vect[:,2].reshape(self.meta['nz'],self.meta['ny'],self.meta['nx']).transpose(2,1,0)
 
    if double_y : 
      vect_x = vect_x[:,:self.meta['ny'],:]
      vect_y = vect_y[:,:self.meta['ny'],:]
      vect_z = vect_z[:,:self.meta['ny'],:]

    if (fibo_obj != None) :
      fibo_obj.data[tar_var_x] = vect_x
      fibo_obj.data[tar_var_y] = vect_y
      fibo_obj.data[tar_var_z] = vect_z
    else: return np.array([vect_x, vect_y, vect_z])

    if not silent:
      print('get_vect_from_VTK> data format             :  ', data_format)
      print('get_vect_from_VTK> data structure          :  ', data_structure)
      print('get_vect_from_VTK> grid dimensions         :  ', self.meta['nnn'])
      print('get_vect_from_VTK> grid size               :  ', self.meta['lll'])
      print('get_vect_from_VTK> grid spacing            :  ', self.meta['ddd'])
      if (fibo_obj != None) :
        print('get_vect_from_VTK> created fibo_obj.data['+tar_var_x+']') 
        print('get_vect_from_VTK> created fibo_obj.data['+tar_var_y+']') 
        print('get_vect_from_VTK> created fibo_obj.data['+tar_var_z+']') 

