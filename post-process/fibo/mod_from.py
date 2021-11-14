import collections 
import codecs
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

#---------------------------------------------------------------------------------------
#------fill-fibo-objects-from-various-sources-------------------------------------------
#-------or-simply-get-your-data---------------------------------------------------------
#---------------------------------------------------------------------------------------

#this class contains all functions to read outputs from 2d simulations with EB, Ion, Press, Q and calculate remaining quantities
#this class could be expanded so to read also outputs of codes with chew-golberg-low and landau closures (explicit electron pressures)
class from_HVM (object):  

  def __init__(self, 
      address):
    """
    Creates the object to retrieve data from Hybrid-Vlasov-Maxwell codes
    
    Parameters :
      - address      [address] where your data are (folder with segs inside)
    
    """

    self.address = os.path.normpath(address)
    self.segs = {}
    self.meta = {}

  #------------------------------------------------------------
  def help(self):
    print('For further help, please shout:')
    print('!!!SIIIIIIIIIIIID!!!')

  #------------------------------------------------------------
  def get_meta(self,  #counts lines in file and calls the appropriate function to get the meta data
      extra_address = '',
      silent = True):
    """
    Fills the metadata list 
    
    Parameters :
      - extra_address = ''     [address] to reach any subfolder where your meta-data are
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    with open(os.path.join(self.address,extra_address,'input_parameters'),'r') as foo:
      line_number = len(foo.readlines())

    if line_number == 75 : self.get_meta_A(extra_address)
    elif line_number == 81 : self.get_meta_B(extra_address)
    elif line_number == 91 : self.get_meta_C(extra_address)
    elif line_number == 92 : self.get_meta_C(extra_address) #cm on nobody cares about satellite number ...
    else : print('FDP : unknown input_parameter file! write a new from_HVM.get_meta() for it ...')



    #---------get-dimensions-of-your-simulation----------

    self.meta['dx'] = self.meta['xl']/self.meta['nx']
    self.meta['dy'] = self.meta['yl']/self.meta['ny']
    self.meta['dz'] = self.meta['zl']/self.meta['nz']

    #nx_original = nx    ??? ASK THE BOSS ABOUT ALL THESE !!!
    #nx = int(nx / nwx)    ??? ASK THE BOSS ABOUT ALL THESE !!!
    #ny = int(ny / nwy)    ??? ASK THE BOSS ABOUT ALL THESE !!!
    #nz = int(nz / nwz)    ??? ASK THE BOSS ABOUT ALL THESE !!!
    
    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz']) 
    self.meta['ppp'] = (True, True, True)                                      #THIS IS HARDCODED! 3-PERIODICITY IS HARDCODED

    self.meta['ts'] = self.meta['dt']      #this is just for jeremy :)
    self.meta['fields'] = {}               #this is for all the world 



    #----------get-time-segment-infos-from all-subdirectories--------- 

    segments = [seg for seg in os.listdir(self.address) if os.path.isdir(os.path.join(self.address,seg))]
    segments = [seg for seg in segments if os.path.isfile(os.path.join(self.address,seg,'tempo.dat'))]
    segments = sorted(segments)

    for seg in segments:
      infot = open(os.path.join(self.address,seg,'tempo.dat'),'r')
      nexits = len(infot.readlines())
      self.segs[seg] = []
      infot.seek(0)
      for time_exit_num in range(nexits):
        time_exit = '%.3f' %np.around(float(infot.readline().split()[0]),decimals=1)  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        self.segs[seg].append(time_exit)
        self.segs = collections.OrderedDict(sorted(self.segs.items())) #Fra's ordering :)
      infot.close()
    
    #Fra's addition :)
    self.meta['time'] = np.concatenate( list(self.segs.values()) ).astype(float)
    self.meta['times'] = np.concatenate( list(self.segs.values()) )
    
    self.meta['time2seg'] = []
    for i in range(len(self.segs.keys())):
        self.meta['time2seg'].append(np.full((len( list(self.segs.values())[i] )), list(self.segs.keys())[i] ))
    self.meta['time2seg'] = np.concatenate(self.meta['time2seg'])
    
    self.meta['time2exit'] = []
    for i in range(len(self.segs.keys())):
        self.meta['time2exit'].append(np.arange(0,len( list(self.segs.values())[i] ),1))
    self.meta['time2exit'] = np.concatenate(self.meta['time2exit'])



    #----------add-informations-on-species-----------------
    #----ACHTUNG: these are hard-coded - change them in future!----------

    self.meta['nss'] = 2     #number of species

    species  = []
    species.append('ion     ')
    species.append('electron')

    charges = np.zeros([self.meta['nss']])
    charges[0] = 1.
    charges[1] = -1.

    masses = np.zeros([self.meta['nss']])
    masses[0] = 1.
    masses[1] = 1./self.meta['mime']

    self.meta['species']  = species
    self.meta['charges']  = { kk:charges[k] for k,kk in enumerate(species)}
    self.meta['masses']   = { kk:masses[k] for k,kk in enumerate(species)}


    #----------print-summary-----------------

    if not silent : 
      print('HVM_'+self.meta['space_dim']+'> cell number               :  ', self.meta['nnn'])
      print('HVM_'+self.meta['space_dim']+'> domain size               :  ', self.meta['lll'])
      print('HVM_'+self.meta['space_dim']+'> mesh spacing              :  ', self.meta['ddd'])
      print('HVM_'+self.meta['space_dim']+'> periodicity               :  ', self.meta['ppp'])
      print('HVM_'+self.meta['space_dim']+'> time step                 :  ', self.meta['ts'])
      print('HVM_'+self.meta['space_dim']+'> species                   :  ', self.meta['species'])
      for i in range(self.meta['nss']):
        print('          '+species[i]+' charge                :  ', self.meta['charges'][species[i]])
        print('          '+species[i]+' mass                  :  ', self.meta['masses'][species[i]])
      print('HVM_'+self.meta['space_dim']+'> teti                      :  ', self.meta['teti'])


  #------------------------------------------------------------
  def get_meta_A(self,
      extra_address = ''):
    """
    Extra routine, version A 
    
    """

    #get mesh infos from input_parameters (I take the input_parameters from subfolder 01)
    infos = open(os.path.join(self.address,extra_address,'input_parameters'),'r')
    infos.readline()
    self.meta['np_row'] = int(infos.readline().split()[2]) #number of mpi processeses in x direction
    self.meta['np_col'] = int(infos.readline().split()[2]) #number of mpi processeses in y direction
    self.meta['np_pla'] = int(infos.readline().split()[2]) #number of mpi processeses in z direction
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['nx'] = int(infos.readline().split()[2])
    self.meta['ny'] = int(infos.readline().split()[2])
    self.meta['nz'] = int(infos.readline().split()[2])
    self.meta['lvx'] = int(infos.readline().split()[2])
    self.meta['lvy'] = int(infos.readline().split()[2])
    self.meta['lvz'] = int(infos.readline().split()[2])
    self.meta['space_dim'] = infos.readline().split("'")[1]  #
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['model'] = int(infos.readline().split()[2]) #not totally sure: in case substitute this line with an innocent infos.readline()
    self.meta['dt'] = float(infos.readline().split()[2])
    nstep = infos.readline() #
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['xl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['yl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['zl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['vxpmax'] = float(infos.readline().split()[2])
    self.meta['vypmax'] = float(infos.readline().split()[2])
    self.meta['vzpmax'] = float(infos.readline().split()[2])
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    
    #infos.readline()
    noutt = infos.readline() #
    noutx = infos.readline() #
    noutf = infos.readline() #
    noutcheck = infos.readline() #
    outformat = infos.readline() #
    ifd = infos.readline() #
    nwx = infos.readline() #
    nwy = infos.readline() #
    nwz = infos.readline() #
    nwvx = infos.readline() #
    nwvy = infos.readline() #
    nwvz = infos.readline() #
    
    infos.readline()
    infos.readline()
    infos.readline()

    initcond = str(infos.readline().split("'")[1])      
    with_forcing = infos.readline() #
    incompressible_forcing = infos.readline() #
    self.meta['Bx0'] = float(infos.readline().split()[2]) #
    self.meta['By0'] = float(infos.readline().split()[2]) #
    self.meta['Bz0'] = float(infos.readline().split()[2]) #
    rhop = infos.readline() #
    self.meta['beta'] = float(infos.readline().split()[2])
    self.meta['beta_para'] = float(infos.readline().split()[2])
    self.meta['beta_perp'] = float(infos.readline().split()[2])
    self.meta['mime'] = float(infos.readline().split()[2])
    self.meta['teti'] = float(infos.readline().split()[2])
    gamma = float(infos.readline().split()[2])
    amp = infos.readline() #
    hk0 = infos.readline() #
    ux0 = infos.readline() #
    uy0 = infos.readline() #
    uz0 = infos.readline() #

    infos.close()

  #------------------------------------------------------------
  def get_meta_B(self,
      extra_address = ''):
    """
    Extra routine, version B 
    
    """

    #get mesh infos from input_parameters (I take the input_parameters from subfolder 01)
    infos = open(os.path.join(self.address,extra_address,'input_parameters'),'r')
    infos.readline()
    self.meta['np_row'] = int(infos.readline().split()[2]) #number of mpi processeses in x direction
    self.meta['np_col'] = int(infos.readline().split()[2]) #number of mpi processeses in y direction
    self.meta['np_pla'] = int(infos.readline().split()[2]) #number of mpi processeses in z direction
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['nx'] = int(infos.readline().split()[2])
    self.meta['ny'] = int(infos.readline().split()[2])
    self.meta['nz'] = int(infos.readline().split()[2])
    self.meta['lvx'] = int(infos.readline().split()[2])
    self.meta['lvy'] = int(infos.readline().split()[2])
    self.meta['lvz'] = int(infos.readline().split()[2])
    self.meta['space_dim'] = infos.readline().split("'")[1] #

    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['model'] = infos.readline() #int(infos.readline().split()[2]) #not totally sure: in case substitute this line with an innocent infos.readline()
    self.meta['dt'] = float(infos.readline().split()[2])
    nstep = infos.readline() #
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['xl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['yl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['zl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['vxpmax'] = float(infos.readline().split()[2])
    self.meta['vypmax'] = float(infos.readline().split()[2])
    self.meta['vzpmax'] = float(infos.readline().split()[2])
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    
    noutt = infos.readline() #
    noutx = infos.readline() #
    noutf = infos.readline() #
    noutcheck = infos.readline() #
    outformat = infos.readline() #
    ifd = infos.readline() #
    nwx = infos.readline() #
    nwy = infos.readline() #
    nwz = infos.readline() #
    nwvx = infos.readline() #
    nwvy = infos.readline() #
    nwvz = infos.readline() #
    
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['Bx0'] = float(infos.readline().split()[2]) 
    self.meta['By0'] = float(infos.readline().split()[2]) #
    self.meta['Bz0'] = float(infos.readline().split()[2]) #
    rhop = infos.readline() #
    self.meta['beta'] = float(infos.readline().split()[2])
    self.meta['mime'] = float(infos.readline().split()[2])
    self.meta['teti'] = float(infos.readline().split()[2])
    gamma = float(infos.readline().split()[2])
    amp = infos.readline() #
    hk0 = infos.readline() #
    ux0 = infos.readline() #
    uy0 = infos.readline() #
    uz0 = infos.readline() #
    initcond = infos.readline() #
    self.meta['beta_para'] = float(infos.readline().split()[2])
    self.meta['beta_perp'] = float(infos.readline().split()[2])
    with_forcing = infos.readline() #
    incompressible_forcing = infos.readline() #
    mx_max = infos.readline()
    my_max = infos.readline()
    mz_max = infos.readline()
    mf_min = infos.readline()
    mf_max = infos.readline()
    amp_forcing = infos.readline()

    infos.close()

  #------------------------------------------------------------
  def get_meta_C(self,  
      extra_address = ''):    #extra address to get to input_parameters
    """
    Extra routine, version C 
    
    """

    #get mesh infos from input_parameters (I take the input_parameters from subfolder 01)
    infos = open(os.path.join(self.address,extra_address,'input_parameters'),'r')
    infos.readline()
    self.meta['np_row'] = int(infos.readline().split()[2]) #number of mpi processeses in x direction
    self.meta['np_col'] = int(infos.readline().split()[2]) #number of mpi processeses in y direction
    self.meta['np_pla'] = int(infos.readline().split()[2]) #number of mpi processeses in z direction
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['nx'] = int(infos.readline().split()[2])
    self.meta['ny'] = int(infos.readline().split()[2])
    self.meta['nz'] = int(infos.readline().split()[2])
    self.meta['lvx'] = int(infos.readline().split()[2])
    self.meta['lvy'] = int(infos.readline().split()[2])
    self.meta['lvz'] = int(infos.readline().split()[2])
    self.meta['space_dim'] = infos.readline().split("'")[1] #
    #
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['model'] = int(infos.readline().split()[2])
    self.meta['dt'] = float(infos.readline().split()[2])
    nstep = infos.readline() #
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    self.meta['xl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['yl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['zl'] = float(infos.readline().split()[2]) * 2 * np.pi
    self.meta['vxpmax'] = float(infos.readline().split()[2])
    self.meta['vypmax'] = float(infos.readline().split()[2])
    self.meta['vzpmax'] = float(infos.readline().split()[2])
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    infos.readline()
    
    noutt = infos.readline() #
    noutx = infos.readline() #
    noutf = infos.readline() #
    infos.readline()
    infos.readline()
    infos.readline()
    noutcheck = infos.readline() #
    outformatx = infos.readline() #
    outformatf = infos.readline()
    iclock = infos.readline() #
    nwx = infos.readline() #
    nwy = infos.readline() #
    nwz = infos.readline() #
    nwvx = infos.readline() #
    nwvy = infos.readline() #
    nwvz = infos.readline() #
    
    infos.readline()
    infos.readline()
    infos.readline()
    initcond = infos.readline() #
    with_forcing = infos.readline() #
    incompressible_forcing = infos.readline() #
    self.meta['Bx0'] = float(infos.readline().split()[2])  #
    self.meta['By0'] = float(infos.readline().split()[2])  #
    self.meta['Bz0'] = float(infos.readline().split()[2])  #
    rhop = infos.readline() #
    self.meta['beta'] = float(infos.readline().split()[2])
    self.meta['mime'] = float(infos.readline().split()[2])
    self.meta['alpha_i'] = float(infos.readline().split()[2])
    self.meta['teti'] = float(infos.readline().split()[2])
    self.meta['alpha_e'] = float(infos.readline().split()[2])
    gamma = float(infos.readline().split()[2])
    amp = infos.readline() #
    hk0 = infos.readline() #
    ux0 = infos.readline() #
    uy0 = infos.readline() #
    uz0 = infos.readline() #
    mx_max = infos.readline()
    my_max = infos.readline()
    mz_max = infos.readline()
    mf_min = infos.readline()
    mf_max = infos.readline()

    infos.close()

  #-----routines-for-fields------------------------------------
  #------------------------------------------------------------
  def get_EB(self,
      seg,
      exit_num,
      fibo_obj = None,
      silent = True):
    """
    Gets the E,B fields at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    E_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    E_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    E_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    # open binary file - its structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 6*8*nx*ny*nz data, 4 empty) for each time exit
    if os.path.isfile(os.path.join(self.address,seg,'EB.bin')) :    
      rf = open(os.path.join(self.address,seg,'EB.bin'), 'rb') # version like jeremys - maybe ok with python3
      #rf = codecs.open(os.path.join(self.address,seg,'EB.bin'), encoding='quopri') # i tried this line to make it work with python3 - yet it does not work
      #rf = open(os.path.join(self.address,seg,'EB.bin'), 'r') # old version - worked with python2 
      #rfr = rf.read() # old version
      
      # jump to the correct line in the file
      offset = 0
      for l in range(exit_num):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        #time_exit = '%.3f' %float(struct.unpack('d',rfr[offset+4:offset+12])[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        offset += 32 + 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

      # check that the data exit corresponds with the correct one
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      #time_exit = '%.3f' %float(struct.unpack('d',rfr[offset+4:offset+12])[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      #if time_exit != self.segs[seg][exit_num] : print('=====WRONG=EXIT=====')
      if not silent: print('reading time:' , time_exit)
      offset += 32

      # fill data vectors
      rf.seek(offset,os.SEEK_SET) # like jeremy
      flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*6)
      #flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*6)+'d',rfr[offset:offset+6*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
      arr = np.reshape(flat_arr,(6,self.meta['nz'],self.meta['ny'],self.meta['nx']))
      E_x = np.transpose(arr[0,:,:,:],(2,1,0))
      E_y = np.transpose(arr[1,:,:,:],(2,1,0))
      E_z = np.transpose(arr[2,:,:,:],(2,1,0))
      B_x = np.transpose(arr[3,:,:,:],(2,1,0)) + self.meta['Bx0']
      B_y = np.transpose(arr[4,:,:,:],(2,1,0)) + self.meta['By0']
      B_z = np.transpose(arr[5,:,:,:],(2,1,0)) + self.meta['Bz0']      #I am putting this here since in the codes I know we have this asymmetry ...

    else :
      rf = open(os.path.join(self.address,seg,'EB.dat'), 'r')
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        for il in range(self.meta['nx']*self.meta['ny']*self.meta['nz']) : 
          rf.readline()

      #check that the data exit corresponds with the correct one
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)

      #fill data vectors
      for iz in range(self.meta['nz']):
        for iy in range(self.meta['ny']):
          for ix in range(self.meta['nx']):
            EB_line = rf.readline().split()
            E_x[ix,iy,iz] = float(EB_line[0])
            E_y[ix,iy,iz] = float(EB_line[1])
            E_z[ix,iy,iz] = float(EB_line[2])
            B_x[ix,iy,iz] = float(EB_line[3])
            B_y[ix,iy,iz] = float(EB_line[4])
            B_z[ix,iy,iz] = float(EB_line[5])

    rf.close()

    #copy all these in fibo, or return them! 
    if (fibo_obj != None) :
      time_exit = self.segs[seg][exit_num]
      fibo_obj.data['E_x_'+time_exit] = E_x
      fibo_obj.data['E_y_'+time_exit] = E_y
      fibo_obj.data['E_z_'+time_exit] = E_z
      fibo_obj.data['B_x_'+time_exit] = B_x
      fibo_obj.data['B_y_'+time_exit] = B_y
      fibo_obj.data['B_z_'+time_exit] = B_z
      fibo_obj.meta['fields']['E_'+time_exit] = ('E_x_'+time_exit,'E_y_'+time_exit,'E_z_'+time_exit)
      fibo_obj.meta['fields']['B_'+time_exit] = ('B_x_'+time_exit,'B_y_'+time_exit,'B_z_'+time_exit)
    else: return np.array([E_x, E_y, E_z]), np.array([B_x, B_y, B_z])

    if not silent: print('done with reading E, B!')

  #------------------------------------------------------------
  def get_Ion(self,
      seg,
      exit_num,
      fibo_obj = None,
      silent = True):
    """
    Gets the ni and ui fields at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    ui_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    ui_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    ui_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    ni = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    # binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 4*8*nx*ny*nz data, 4 empty) for each time exit
    if os.path.isfile(os.path.join(self.address,seg,'Ion.bin')) :   
      rf = open(os.path.join(self.address,seg,'Ion.bin'), 'rb') # version like jeremys - maybe ok with python3

      offset = 0
      #jump to the correct line in the file
      for l in range(exit_num):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        offset += 32 + 4*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] + 4

      # check that the data exit corresponds with the correct one
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)
      offset += 32

      # fill data vectors
      rf.seek(offset,os.SEEK_SET)
      flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*4)
      arr = np.reshape(flat_arr,(4,self.meta['nz'],self.meta['ny'],self.meta['nx']))
      ui_x = np.transpose(arr[0,:,:,:],(2,1,0))
      ui_y = np.transpose(arr[1,:,:,:],(2,1,0))
      ui_z = np.transpose(arr[2,:,:,:],(2,1,0))
      ni = np.transpose(arr[3,:,:,:],(2,1,0))

    else :    
      rf = open(os.path.join(self.address,seg,'Ion.dat'), 'r')
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        for il in range(self.meta['nx']*self.meta['ny']*self.meta['nz']) : 
          rf.readline()

      #check that the data exit corresponds with the correct one
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)

      #fill data vectors
      for iz in range(self.meta['nz']):
        for iy in range(self.meta['ny']):
          for ix in range(self.meta['nx']):
            Ion_line = rf.readline().split()
            ui_x[ix,iy,iz] = float(Ion_line[0])
            ui_y[ix,iy,iz] = float(Ion_line[1])
            ui_z[ix,iy,iz] = float(Ion_line[2])
            ni[ix,iy,iz] = float(Ion_line[3])

    rf.close()
    #ni[np.isnan(ni)] = 0.0

    ni = ni + 1.
    ui_x = np.divide(ui_x , ni)
    ui_y = np.divide(ui_y , ni)
    ui_z = np.divide(ui_z , ni)
    
    if (fibo_obj != None) :
      time_exit = self.segs[seg][exit_num]
      fibo_obj.data['ui_x_'+time_exit] = ui_x 
      fibo_obj.data['ui_y_'+time_exit] = ui_y
      fibo_obj.data['ui_z_'+time_exit] = ui_z
      fibo_obj.data['n_'+time_exit] = ni
      fibo_obj.meta['fields']['ui_'+time_exit] = ('ui_x_'+time_exit,'ui_y_'+time_exit,'ui_z_'+time_exit)
      fibo_obj.meta['fields']['n_'+time_exit]  = ('n_'+time_exit,)
    else: return ni, np.array([ui_x, ui_y, ui_z])

    if not silent: print('done with reading ni, ui!')

  #------------------------------------------------------------
  def get_Press(self,
      seg,            #segment number
      exit_num,      #number of the exit you choose
      fibo_obj = None,  #fibo object you are considering - if None, Ion values will just be returned 
      silent = True): 
    """
    Gets the Pi tensor field at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    Pi_xx = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_yy = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_zz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_xy = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_xz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_yz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    # binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 6*8*nx*ny*nz data, 4 empty) for each time exit
    if os.path.isfile(os.path.join(self.address,seg,'Press.bin')) :   
      rf = open(os.path.join(self.address,seg,'Press.bin'), 'rb') # version like jeremys - maybe ok with python3

      offset = 0
      #jump to the correct line in the file
      for l in range(exit_num):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        offset += 32 + 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

      # check that the data exit corresponds with the correct one
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)
      offset += 32

      # fill data vectors
      rf.seek(offset,os.SEEK_SET)
      flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*6)
      arr = np.reshape(flat_arr,(6,self.meta['nz'],self.meta['ny'],self.meta['nx']))
      Pi_xx = np.transpose(arr[0,:,:,:],(2,1,0))
      Pi_yy = np.transpose(arr[1,:,:,:],(2,1,0))
      Pi_zz = np.transpose(arr[2,:,:,:],(2,1,0))
      Pi_xy = np.transpose(arr[3,:,:,:],(2,1,0))
      Pi_xz = np.transpose(arr[4,:,:,:],(2,1,0))
      Pi_yz = np.transpose(arr[5,:,:,:],(2,1,0))

    else :
      rf = open(os.path.join(self.address,seg,'Press.dat'), 'r')
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        for il in range(self.meta['nx']*self.meta['ny']*self.meta['nz']) : 
          rf.readline()

      #check that the data exit corresponds with the correct one
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)

      #fill data vectors
      for iz in range(self.meta['nz']):
        for iy in range(self.meta['ny']):
          for ix in range(self.meta['nx']):
            Press_line = rf.readline().split()
            Pi_xx[ix,iy,iz] = float(Press_line[0])
            Pi_yy[ix,iy,iz] = float(Press_line[1])
            Pi_zz[ix,iy,iz] = float(Press_line[2])
            Pi_xy[ix,iy,iz] = float(Press_line[3])
            Pi_xz[ix,iy,iz] = float(Press_line[4])
            Pi_yz[ix,iy,iz] = float(Press_line[5])

    rf.close()

    if (fibo_obj != None) :
      time_exit = self.segs[seg][exit_num]
      fibo_obj.data['Pi_xx_'+time_exit] = Pi_xx[:,:,:]
      fibo_obj.data['Pi_yy_'+time_exit] = Pi_yy[:,:,:]
      fibo_obj.data['Pi_zz_'+time_exit] = Pi_zz[:,:,:]
      fibo_obj.data['Pi_xy_'+time_exit] = Pi_xy[:,:,:]
      fibo_obj.data['Pi_xz_'+time_exit] = Pi_xz[:,:,:]
      fibo_obj.data['Pi_yz_'+time_exit] = Pi_yz[:,:,:]
      fibo_obj.meta['fields']['Pi_'+time_exit]  = ('Pi_xx_'+time_exit,'Pi_yy_'+time_exit,'Pi_zz_'+time_exit)
      fibo_obj.meta['fields']['Pi_'+time_exit] += ('Pi_xy_'+time_exit,'Pi_xz_'+time_exit,'Pi_yz_'+time_exit)
    else: return np.array([[Pi_xx, Pi_xy, Pi_xz],[Pi_xy, Pi_yy, Pi_yz],[Pi_xz, Pi_yz, Pi_zz]])

    if not silent: print('done with reading Pi!')

  #------------------------------------------------------------
  def get_Q(self,    #achtung that two different outputs are given, depending on ux - so I will need to post-process it
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the Qi_per,Qi_par fields OR sQi field at the nth exit of the data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    if os.path.isfile(os.path.join(self.address,seg,'Q.bin')) :

      data_size = 8*self.meta['nx']*self.meta['ny']*self.meta['nz']*len(self.segs[seg])

      #case with Q_par and Q_per saved
      if os.path.getsize(os.path.join(self.address,seg,'Q.bin')) / data_size == 6 :
        #create data vectors
        Qi_par_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_par_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_par_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_per_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_per_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_per_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

        #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 6*8*nx*ny*nz data, 4 empty) for each time exit
        rf = open(os.path.join(self.address,seg,'Q.bin'), 'rb')
        
        offset = 0
        #jump to the correct line in the file
        for l in range(exit_num):
          rf.seek(offset+4,os.SEEK_SET)
          time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
          time_exit = time_exit.zfill(8)        # ..and three decimal digits
          if not silent: print('jumping time:' , time_exit)
          offset += 32 + 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4
        
        # check that the data exit corresponds with the correct one
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        offset += 32
        
        # fill data vectors
        rf.seek(offset,os.SEEK_SET)
        flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*6)
        arr = np.reshape(flat_arr,(6,self.meta['nz'],self.meta['ny'],self.meta['nx']))
        Qi_par_x = np.transpose(arr[0,:,:,:],(2,1,0))
        Qi_par_y = np.transpose(arr[1,:,:,:],(2,1,0))
        Qi_par_z = np.transpose(arr[2,:,:,:],(2,1,0))
        Qi_per_x = np.transpose(arr[3,:,:,:],(2,1,0))
        Qi_per_y = np.transpose(arr[4,:,:,:],(2,1,0))
        Qi_per_z = np.transpose(arr[5,:,:,:],(2,1,0))
        
        rf.close()
        
        if (fibo_obj != None) :
          time_exit = self.segs[seg][exit_num]
          fibo_obj.data['Qi_par_x_'+time_exit] = Qi_par_x[:,:,:]
          fibo_obj.data['Qi_par_y_'+time_exit] = Qi_par_y[:,:,:]
          fibo_obj.data['Qi_par_z_'+time_exit] = Qi_par_z[:,:,:]
          fibo_obj.data['Qi_per_x_'+time_exit] = Qi_per_x[:,:,:]
          fibo_obj.data['Qi_per_y_'+time_exit] = Qi_per_y[:,:,:]
          fibo_obj.data['Qi_per_z_'+time_exit] = Qi_per_z[:,:,:]
          fibo_obj.meta['fields']['Qi_par_'+time_exit] = ('Qi_par_x_'+time_exit,'Qi_par_y_'+time_exit,'Qi_par_z_'+time_exit)
          fibo_obj.meta['fields']['Qi_per_'+time_exit] = ('Qi_per_x_'+time_exit,'Qi_per_y_'+time_exit,'Qi_per_z_'+time_exit)
        else: return np.array([Qi_par_x, Qi_par_y, Qi_par_z]), np.array([Qi_per_x, Qi_per_y, Qi_per_z])
        if not silent: print('done with reading Qi_par and Qi_per!')
      
      #case with all 10 components of Q
      else :
        Qi_xxx = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_yyy = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_zzz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_xxy = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_xxz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_yyx = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])      
        Qi_yyz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_zzx = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_zzy = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        Qi_xyz = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']]) 

        #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 10*8*nx*ny*nz data, 4 empty) for each time exit
        rf = open(os.path.join(self.address,seg,'Q.bin'), 'rb')
        
        offset = 0
        #jump to the correct line in the file
        for l in range(exit_num):
          rf.seek(offset+4,os.SEEK_SET)
          time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
          time_exit = time_exit.zfill(8)        # ..and three decimal digits
          if not silent: print('jumping time:' , time_exit)
          offset += 32 + 10*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4
        
        # check that the data exit corresponds with the correct one
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        offset += 32
        
        # fill data vectors
        rf.seek(offset,os.SEEK_SET)
        flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*10)
        arr = np.reshape(flat_arr,(10,self.meta['nz'],self.meta['ny'],self.meta['nx']))
        Qi_xxx = np.transpose(arr[0,:,:,:],(2,1,0))
        Qi_yyy = np.transpose(arr[1,:,:,:],(2,1,0))
        Qi_zzz = np.transpose(arr[2,:,:,:],(2,1,0))
        Qi_xxy = np.transpose(arr[3,:,:,:],(2,1,0))
        Qi_xxz = np.transpose(arr[4,:,:,:],(2,1,0))
        Qi_yyx = np.transpose(arr[5,:,:,:],(2,1,0))
        Qi_yyz = np.transpose(arr[6,:,:,:],(2,1,0))
        Qi_zzx = np.transpose(arr[7,:,:,:],(2,1,0)) 
        Qi_zzy = np.transpose(arr[8,:,:,:],(2,1,0))
        Qi_xyz = np.transpose(arr[9,:,:,:],(2,1,0))
        rf.close()
        
        if (fibo_obj != None) :
          time_exit = self.segs[seg][exit_num]
          fibo_obj.data['Qi_xxx_'+time_exit] =  Qi_xxx[:,:,:]
          fibo_obj.data['Qi_yyy_'+time_exit] =  Qi_yyy[:,:,:]
          fibo_obj.data['Qi_zzz_'+time_exit] =  Qi_zzz[:,:,:]
          fibo_obj.data['Qi_xxy_'+time_exit] =  Qi_xxy[:,:,:]
          fibo_obj.data['Qi_xxz_'+time_exit] =  Qi_xxz[:,:,:]
          fibo_obj.data['Qi_yyx_'+time_exit] =  Qi_yyx[:,:,:]
          fibo_obj.data['Qi_yyz_'+time_exit] =  Qi_yyz[:,:,:]
          fibo_obj.data['Qi_zzx_'+time_exit] =  Qi_zzx[:,:,:]
          fibo_obj.data['Qi_zzy_'+time_exit] =  Qi_zzy[:,:,:]
          fibo_obj.data['Qi_xyz_'+time_exit] =  Qi_xyz[:,:,:]
          fibo_obj.meta['fields']['Qi_'+time_exit]  = ('Qi_xxx_'+time_exit,'Qi_yyy_'+time_exit,'Qi_zzz_'+time_exit)
          fibo_obj.meta['fields']['Qi_'+time_exit] += ('Qi_xxy_'+time_exit,'Qi_xxz_'+time_exit,'Qi_yyx_'+time_exit)
          fibo_obj.meta['fields']['Qi_'+time_exit] += ('Qi_yyz_'+time_exit,'Qi_zzx_'+time_exit,'Qi_zzy_'+time_exit,'Qi_xyz_'+time_exit)
        else: return np.array([Qi_xxx, Qi_yyy, Qi_zzz, Qi_xxy, Qi_xxz, Qi_yyx, Qi_yyz, Q_zzx, Q_zzy, Q_xyz])
        if not silent: print('done with reading Qi!')
        #if not silent: print('Qi_xxx, Qi_yyy, Qi_zzz, Qi_xxy, Qi_xxz, Qi_yyx, Qi_yyz, Q_zzx, Q_zzy, Q_xyz')

    else :
      #create data vectors
      sQi_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
      sQi_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
      sQi_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

      rf = open(os.path.join(self.address,seg,'Q.dat'), 'r')
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('jumping time:' , time_exit)
        for il in range(self.meta['nx']*self.meta['ny']*self.meta['nz']) : 
          rf.readline()

      #check that the data exit corresponds with the correct one
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)

      #fill data vectors
      for iz in range(self.meta['nz']):
        for iy in range(self.meta['ny']):
          for ix in range(self.meta['nx']):
            Q_line = rf.readline().split()
            sQi_x[ix,iy,iz] = float(Q_line[0])
            sQi_y[ix,iy,iz] = float(Q_line[1])
            sQi_z[ix,iy,iz] = float(Q_line[2])


      rf.close()

      if (fibo_obj != None) :
        time_exit = self.segs[seg][exit_num]
        fibo_obj.data['sQi_x_'+time_exit] = sQi_x[:,:,:]
        fibo_obj.data['sQi_y_'+time_exit] = sQi_y[:,:,:]
        fibo_obj.data['sQi_z_'+time_exit] = sQi_z[:,:,:]
        fibo_obj.meta['fields']['sQi_'+time_exit] = ('sQi_x_'+time_exit,'sQi_y_'+time_exit,'sQi_z_'+time_exit)
      else: return np.array([sQi_x, sQi_y, sQi_z])
      if not silent: print('done with reading sQi!')

  #------------------------------------------------------------
  def get_Te(self,   
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the Te_per,Te_par fields at the nth exit of the data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    Te_par = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Te_per = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 2*8*nx*ny*nz data, 4 empty) for each time exit
    rf = open(os.path.join(self.address,seg,'Te.bin'), 'rb')
    
    offset = 0
    #jump to the correct line in the file
    for l in range(exit_num):
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('jumping time:' , time_exit)
      offset += 32 + 2*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

    # check that the data exit corresponds with the correct one
    rf.seek(offset+4,os.SEEK_SET)
    time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
    time_exit = time_exit.zfill(8)        # ..and three decimal digits
    if not silent: print('reading time:' , time_exit)
    offset += 32

    # fill data vectors
    rf.seek(offset,os.SEEK_SET)
    flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*2)
    arr = np.reshape(flat_arr,(2,self.meta['nz'],self.meta['ny'],self.meta['nx']))
    Te_par = np.transpose(arr[0,:,:,:],(2,1,0)) + 0.5 * self.meta['beta'] * self.meta['teti']
    Te_per = np.transpose(arr[1,:,:,:],(2,1,0)) + 0.5 * self.meta['beta'] * self.meta['teti'] * self.meta['alpha_e']

    rf.close()

    if (fibo_obj != None) :
      time_exit = self.segs[seg][exit_num]
      fibo_obj.data['Te_par_'+time_exit] = Te_par[:,:,:]
      fibo_obj.data['Te_per_'+time_exit] = Te_per[:,:,:]
      fibo_obj.meta['fields']['Te_par_'+time_exit] = ('Te_par_'+time_exit,)
      fibo_obj.meta['fields']['Te_per_'+time_exit] = ('Te_per_'+time_exit,)
    else: return np.array([Te_par, Te_per])
    if not silent: print('done with reading Te_par and Te_per!')

  #------------------------------------------------------------
  def get_Qe(self,   
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the qe_per,qe_par fields at the nth exit of the data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    qe_par = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    qe_per = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 2*8*nx*ny*nz data, 4 empty) for each time exit
    rf = open(os.path.join(self.address,seg,'Qe.bin'), 'rb')
    
    offset = 0
    #jump to the correct line in the file
    for l in range(exit_num):
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('jumping time:' , time_exit)
      offset += 32 + 2*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

    # check that the data exit corresponds with the correct one
    rf.seek(offset+4,os.SEEK_SET)
    time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
    time_exit = time_exit.zfill(8)        # ..and three decimal digits
    if not silent: print('reading time:' , time_exit)
    offset += 32

    # fill data vectors
    rf.seek(offset,os.SEEK_SET)
    flat_arr = np.fromfile(rf,'float64',self.meta['nx']*self.meta['ny']*self.meta['nz']*2)
    arr = np.reshape(flat_arr,(2,self.meta['nz'],self.meta['ny'],self.meta['nx']))
    qe_par = np.transpose(arr[0,:,:,:],(2,1,0))
    qe_per = np.transpose(arr[1,:,:,:],(2,1,0))

    rf.close()

    if (fibo_obj != None) :
      time_exit = self.segs[seg][exit_num]
      fibo_obj.data['qe_par_'+time_exit] = qe_par[:,:,:]
      fibo_obj.data['qe_per_'+time_exit] = qe_per[:,:,:]
      fibo_obj.meta['fields']['qe_per_'+time_exit] = ('qe_per_'+time_exit,)
      fibo_obj.meta['fields']['qe_par_'+time_exit] = ('qe_par_'+time_exit,)
    else: return np.array([qe_par, qe_per])
    if not silent: print('done with reading qe_par and qe_per!')

  #------------------------------------------------------------
  #------------------------------------------------------------
  def get_seg_EB(self,
      seg,
      fibo_obj = None,
      silent = True): 
    """
    Gets the E,B fields in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nexits = len(self.segs[seg])
    E_x = np.zeros([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    E_y = np.zeros([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    E_z = np.zeros([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_x = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_y = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_z = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])

    if os.path.isfile(os.path.join(self.address,seg,'EB.bin')) :    
      #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 6*8*nx*ny*nz data, 4 empty) for each time exit
      rf = open(os.path.join(self.address,seg,'EB.bin'), 'r')
      rfr = rf.read()
      if len(rfr) != (36 + 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'])*nexits : print('===!!!==wrong=segs=:=correct==!!!===')
      offset = 0

      #jump to the correct line in the file and read data vectors
      for l in range(nexits):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        offset += 32 

        flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*6)+'d',rfr[offset:offset+6*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
        arr = np.reshape(flat_arr,(6,self.meta['nz'],self.meta['ny'],self.meta['nx']))
        E_x[l,:,:,:] = np.transpose(arr[0,:,:,:],(2,1,0))
        E_y[l,:,:,:] = np.transpose(arr[1,:,:,:],(2,1,0))
        E_z[l,:,:,:] = np.transpose(arr[2,:,:,:],(2,1,0))
        B_x[l,:,:,:] = np.transpose(arr[3,:,:,:],(2,1,0)) + self.meta['Bx0']
        B_y[l,:,:,:] = np.transpose(arr[4,:,:,:],(2,1,0)) + self.meta['By0']
        B_z[l,:,:,:] = np.transpose(arr[5,:,:,:],(2,1,0)) + self.meta['Bz0']      #I am putting this here since in the codes I know we have this asymmetry ...
        offset += 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

    else :
      rf = open(os.path.join(self.address,seg,'EB.dat'), 'r')

      for l in range(nexits):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        for iz in range(self.meta['nz']):
          for iy in range(self.meta['ny']):
            for ix in range(self.meta['nx']):
              EB_line = rf.readline().split()
              E_x[l,ix,iy,iz] = float(EB_line[0])
              E_y[l,ix,iy,iz] = float(EB_line[1])
              E_z[l,ix,iy,iz] = float(EB_line[2])
              B_x[l,ix,iy,iz] = float(EB_line[3])
              B_y[l,ix,iy,iz] = float(EB_line[4])
              B_z[l,ix,iy,iz] = float(EB_line[5])

    rf.close()

    if (fibo_obj != None) :
      for l in range(nexits):
        time_exit = self.segs[seg][l]
        fibo_obj.data['E_x_'+time_exit] = E_x[l,:,:,:]
        fibo_obj.data['E_y_'+time_exit] = E_y[l,:,:,:]
        fibo_obj.data['E_z_'+time_exit] = E_z[l,:,:,:]
        fibo_obj.data['B_x_'+time_exit] = B_x[l,:,:,:]
        fibo_obj.data['B_y_'+time_exit] = B_y[l,:,:,:]
        fibo_obj.data['B_z_'+time_exit] = B_z[l,:,:,:]
        fibo_obj.meta['fields']['E_'+time_exit] = ('E_x_'+time_exit,'E_y_'+time_exit,'E_z_'+time_exit)
        fibo_obj.meta['fields']['B_'+time_exit] = ('B_x_'+time_exit,'B_y_'+time_exit,'B_z_'+time_exit)
    else: return np.array([E_x, E_y, E_z]), np.array([B_x, B_y, B_z])
    if not silent: print('done with reading E, B!')

  #------------------------------------------------------------
  def get_seg_Ion(self,
      seg,        #segment considered (str)
      fibo_obj = None,  #fibo object you are considering - if None, Ion values will just be returned 
      silent = True): 
    """
    Gets the ni,ui fields in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nexits = len(self.segs[seg])
    ui_x = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    ui_y = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    ui_z = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    ni = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])

    if os.path.isfile(os.path.join(self.address,seg,'Ion.bin')) :    
      #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 4*8*nx*ny*nz data, 4 empty) for each time exit
      rf = open(os.path.join(self.address,seg,'Ion.bin'), 'r')
      rfr = rf.read()
      if len(rfr) != (36 + 4*8*self.meta['nx']*self.meta['ny']*self.meta['nz'])*nexits : print('===!!!==wrong=segs=:=correct==!!!===')
      offset = 0

      #jump to the correct line in the file and read data vectors
      for l in range(nexits):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        offset += 32 

        flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*4)+'d',rfr[offset:offset+4*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
        arr = np.reshape(flat_arr,(4,self.meta['nz'],self.meta['ny'],self.meta['nx']))
        ui_x[l,:,:,:] = np.transpose(arr[0,:,:,:],(2,1,0))
        ui_y[l,:,:,:] = np.transpose(arr[1,:,:,:],(2,1,0))
        ui_z[l,:,:,:] = np.transpose(arr[2,:,:,:],(2,1,0))
        ni[l,:,:,:]   = np.transpose(arr[3,:,:,:],(2,1,0))
        offset += 4*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

    else :
      rf = open(os.path.join(self.address,seg,'Ion.dat'), 'r')
  
      for l in range(nexits):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        for iz in range(self.meta['nz']):
          for iy in range(self.meta['ny']):
            for ix in range(self.meta['nx']):
              Ion_line = rf.readline().split()
              ui_x[l,ix,iy,iz] = float(Ion_line[0])
              ui_y[l,ix,iy,iz] = float(Ion_line[1])
              ui_z[l,ix,iy,iz] = float(Ion_line[2])
              ni[l,ix,iy,iz] = float(Ion_line[3])

    rf.close()

    ni = ni + 1.
    ui_x = np.divide(ui_x , ni)
    ui_y = np.divide(ui_y , ni)
    ui_z = np.divide(ui_z , ni)

    if (fibo_obj != None) :
      for l in range(nexits):
        time_exit = self.segs[seg][l]
        fibo_obj.data['ui_x_'+time_exit] = ui_x[l,:,:,:]
        fibo_obj.data['ui_y_'+time_exit] = ui_y[l,:,:,:]
        fibo_obj.data['ui_z_'+time_exit] = ui_z[l,:,:,:]
        fibo_obj.data['n_'+time_exit] = ni[l,:,:,:] 
        fibo_obj.meta['fields']['ui_'+time_exit] = ('ui_x_'+time_exit,'ui_y_'+time_exit,'ui_z_'+time_exit)
        fibo_obj.meta['fields']['n_'+time_exit]  = ('n_'+time_exit,)
    else: return ni, np.array([ui_x, ui_y, ui_z])
    if not silent: print('done with reading ni, ui!')

  #------------------------------------------------------------
  def get_seg_Press(self,
      seg,        #segment considered (str)
      fibo_obj = None,  #fibo object you are considering - if None, Ion values will just be returned 
      silent = True): 
    """
    Gets the Pi tensor field in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nexits = len(self.segs[seg])
    Pi_xx = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_yy = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_zz = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_xy = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_xz = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pi_yz = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])

    if os.path.isfile(os.path.join(self.address,seg,'Press.bin')) :  
    #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 6*8*nx*ny*nz data, 4 empty) for each time exit
      rf = open(os.path.join(self.address,seg,'Press.bin'), 'r')
      rfr = rf.read()
      if len(rfr) != (36 + 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'])*nexits : print('===!!!==wrong=segs=:=correct==!!!===')
      offset = 0
      #jump to the correct line in the file and read data vectors
      for l in range(nexits):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        offset += 32 

        flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*6)+'d',rfr[offset:offset+6*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
        arr = np.reshape(flat_arr,(6,self.meta['nz'],self.meta['ny'],self.meta['nx']))
        Pi_xx[l,:,:,:] = np.transpose(arr[0,:,:,:],(2,1,0))
        Pi_yy[l,:,:,:] = np.transpose(arr[1,:,:,:],(2,1,0))
        Pi_zz[l,:,:,:] = np.transpose(arr[2,:,:,:],(2,1,0))
        Pi_xy[l,:,:,:] = np.transpose(arr[3,:,:,:],(2,1,0))
        Pi_xz[l,:,:,:] = np.transpose(arr[4,:,:,:],(2,1,0))
        Pi_yz[l,:,:,:] = np.transpose(arr[5,:,:,:],(2,1,0)) 
        offset += 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

    else :
      rf = open(os.path.join(self.address,seg,'Press.dat'), 'r')
    
      for l in range(nexits):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        for iz in range(self.meta['nz']):
          for iy in range(self.meta['ny']):
            for ix in range(self.meta['nx']):
              Press_line = rf.readline().split()
              Pi_xx[l,ix,iy,iz] = float(Press_line[0])
              Pi_yy[l,ix,iy,iz] = float(Press_line[1])
              Pi_zz[l,ix,iy,iz] = float(Press_line[2])
              Pi_xy[l,ix,iy,iz] = float(Press_line[3])
              Pi_xz[l,ix,iy,iz] = float(Press_line[4])
              Pi_yz[l,ix,iy,iz] = float(Press_line[5])

    rf.close()

    if (fibo_obj != None) :
      for l in range(nexits):
        time_exit = self.segs[seg][l]
        fibo_obj.data['Pi_xx_'+time_exit] = Pi_xx[l,:,:,:]
        fibo_obj.data['Pi_yy_'+time_exit] = Pi_yy[l,:,:,:]
        fibo_obj.data['Pi_zz_'+time_exit] = Pi_zz[l,:,:,:]
        fibo_obj.data['Pi_xy_'+time_exit] = Pi_xy[l,:,:,:]
        fibo_obj.data['Pi_xz_'+time_exit] = Pi_xz[l,:,:,:]
        fibo_obj.data['Pi_yz_'+time_exit] = Pi_yz[l,:,:,:]
        fibo_obj.meta['fields']['Pi_'+time_exit]  = ('Pi_xx_'+time_exit,'Pi_yy_'+time_exit,'Pi_zz_'+time_exit)
        fibo_obj.meta['fields']['Pi_'+time_exit] += ('Pi_xy_'+time_exit,'Pi_xz_'+time_exit,'Pi_yz_'+time_exit)
    else: return np.array([[Pi_xx, Pi_xy, Pi_xz],[Pi_xy, Pi_yy, Pi_yz],[Pi_xz, Pi_yz, Pi_zz]])
    if not silent: print('done with reading Pi!')

  #------------------------------------------------------------
  def get_seg_Q(self,  #achtung that silvio's Q is actually more than the heat flux - so I will need to post-process it
      seg,      
      fibo_obj = None,
      silent = True): 
    """
    Gets the Qi_par,Qi_per fields OR the sQi vector in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nexits = len(self.segs[seg])

    if os.path.isfile(os.path.join(self.address,seg,'Q.bin')) :  
      Qi_par_x = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      Qi_par_y = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      Qi_par_z = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      Qi_per_x = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      Qi_per_y = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      Qi_per_z = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])

      #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 6*8*nx*ny*nz data, 4 empty) for each time exit
      rf = open(os.path.join(self.address,seg,'Q.bin'), 'r')
      rfr = rf.read()
      if len(rfr) != (36 + 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'])*nexits : print('===!!!==wrong=segs=:=correct==!!!===')
      offset = 0

      #jump to the correct line in the file and read data vectors
      for l in range(nexits):
        rf.seek(offset+4,os.SEEK_SET)
        time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        offset += 32 

        flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*6)+'d',rfr[offset:offset+6*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
        arr = np.reshape(flat_arr,(6,self.meta['nz'],self.meta['ny'],self.meta['nx']))
        Qi_par_x[l,:,:,:] = np.transpose(arr[0,:,:,:],(2,1,0))
        Qi_par_y[l,:,:,:] = np.transpose(arr[1,:,:,:],(2,1,0))
        Qi_par_z[l,:,:,:] = np.transpose(arr[2,:,:,:],(2,1,0))
        Qi_per_x[l,:,:,:] = np.transpose(arr[3,:,:,:],(2,1,0))
        Qi_per_y[l,:,:,:] = np.transpose(arr[4,:,:,:],(2,1,0))
        Qi_per_z[l,:,:,:] = np.transpose(arr[5,:,:,:],(2,1,0)) 
        offset += 6*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

      rf.close()

      if (fibo_obj != None) :  
        fibo_obj.xl = self.meta['xl']
        fibo_obj.yl = self.meta['yl']
        fibo_obj.zl = self.meta['zl']
        fibo_obj.nx = self.meta['nx']
        fibo_obj.ny = self.meta['ny']
        fibo_obj.nz = self.meta['nz']
        for l in range(nexits):
          time_exit = self.segs[seg][l]
          fibo_obj.data['Qi_par_x_'+time_exit] = Qi_par_x[l,:,:,:]
          fibo_obj.data['Qi_par_y_'+time_exit] = Qi_par_y[l,:,:,:]
          fibo_obj.data['Qi_par_z_'+time_exit] = Qi_par_z[l,:,:,:]
          fibo_obj.data['Qi_per_x_'+time_exit] = Qi_per_x[l,:,:,:]
          fibo_obj.data['Qi_per_y_'+time_exit] = Qi_per_y[l,:,:,:]
          fibo_obj.data['Qi_per_z_'+time_exit] = Qi_per_z[l,:,:,:]
          fibo_obj.meta['fields']['Qi_par_'+time_exit] = ('Qi_par_x_'+time_exit,'Qi_par_y_'+time_exit,'Qi_par_z_'+time_exit)
          fibo_obj.meta['fields']['Qi_per_'+time_exit] = ('Qi_per_x_'+time_exit,'Qi_per_y_'+time_exit,'Qi_per_z_'+time_exit)
      else: return np.array([Qi_par_x, Qi_par_y, Qi_par_z]), np.array([Qi_per_x, Qi_per_y, Qi_per_z])
      if not silent: print('done with reading Qi_par and Qi_per!')

    else :
      #create data vectors
      sQi_x = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      sQi_y = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
      sQi_z = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])

      rf = open(os.path.join(self.address,seg,'Q.dat'), 'r')
  
      for l in range(nexits):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        if not silent: print('reading time:' , time_exit)
        for iz in range(self.meta['nz']):
          for iy in range(self.meta['ny']):
            for ix in range(self.meta['nx']):
              Q_line = rf.readline().split()
              sQi_x[l,ix,iy,iz] = float(Q_line[0])
              sQi_y[l,ix,iy,iz] = float(Q_line[1])
              sQi_z[l,ix,iy,iz] = float(Q_line[2])

      rf.close()

      if (fibo_obj != None) :
        for l in range(nexits):
          time_exit = self.segs[seg][l]
          fibo_obj.data['sQi_x_'+time_exit] = sQi_x[l,:,:,:]
          fibo_obj.data['sQi_y_'+time_exit] = sQi_y[l,:,:,:]
          fibo_obj.data['sQi_z_'+time_exit] = sQi_z[l,:,:,:]
          fibo_obj.meta['fields']['sQi_'+time_exit]  = ('sQi_x_'+time_exit,'sQi_y_'+time_exit,'sQi_z_'+time_exit)
      else: return np.array([sQi_x, sQi_y, sQi_z])
      if not silent: print('done with reading sQi!')

  #------------------------------------------------------------
  def get_seg_Te(self,   
      seg,
      fibo_obj = None,
      silent = True): 
    """
    Gets the Te_per,Te_par fields in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nexits = len(self.segs[seg])

    Te_par = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Te_per = np.empty([nexits,self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 2*8*nx*ny*nz data, 4 empty) for each time exit
    rf = open(os.path.join(self.address,seg,'Te.bin'), 'r')
    rfr = rf.read()
    if len(rfr) != (36 + 2*8*self.meta['nx']*self.meta['ny']*self.meta['nz'])*nexits : print('===!!!==wrong=segs=:=correct==!!!===')
    offset = 0

    #jump to the correct line in the file and read data vectors
    for l in range(nexits):
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)
      offset += 32 

      flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*2)+'d',rfr[offset:offset+2*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
      arr = np.reshape(flat_arr,(2,self.meta['nz'],self.meta['ny'],self.meta['nx']))
      Te_par[l,:,:,:] = np.transpose(arr[0,:,:,:],(2,1,0)) + 0.5 * self.meta['beta'] * self.meta['teti']
      Te_per[l,:,:,:] = np.transpose(arr[1,:,:,:],(2,1,0)) + 0.5 * self.meta['beta'] * self.meta['teti'] * self.meta['alpha_e']
      offset += 2*8*self.meta['nx']*self.meta['ny']*self.meta['nz'] +4

    rf.close()

    if (fibo_obj != None) :
      for l in range(nexits): 
        time_exit = self.segs[seg][l]
        fibo_obj.data['Te_par_'+time_exit] = Te_par[l,:,:,:]
        fibo_obj.data['Te_per_'+time_exit] = Te_per[l,:,:,:]
        fibo_obj.meta['fields']['Te_par_'+time_exit]  = ('Te_par_'+time_exit,)
        fibo_obj.meta['fields']['Te_per_'+time_exit]  = ('Te_per_'+time_exit,)
    else: return np.array([Te_par, Te_per])
    if not silent: print('done with reading Te_par and Te_per!')

  #------------------------------------------------------------
  def get_seg_Qe(self,   
      seg,
      fibo_obj = None,
      silent = True): 
    """
    Gets the qe_per,qe_par fields at the nth exit of the data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nexits = len(self.segs[seg])

    qe_par = np.empty([nexits, self.meta['nx'],self.meta['ny'],self.meta['nz']])
    qe_per = np.empty([nexits, self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB binary file structure: (4 empty, 8 time_exit, 4+4+4 nx ny nz, 8 empty, 2*8*nx*ny*nz data, 4 empty) for each time exit
    rf = open(os.path.join(self.address,seg,'Qe.bin'), 'r')
    rfr = rf.read()
    if len(rfr) != (36 + 2*8*self.meta['nx']*self.meta['ny']*self.meta['nz'])*nexits : print('===!!!==wrong=segs=:=correct==!!!===')
    offset = 0

    #jump to the correct line in the file
    for l in range(nexits):
      rf.seek(offset+4,os.SEEK_SET)
      time_exit = '%.3f' %float(np.fromfile(rf,'float64',8)[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      if not silent: print('reading time:' , time_exit)
      offset += 32 

      flat_arr = struct.unpack(str(self.meta['nx']*self.meta['ny']*self.meta['nz']*2)+'d',rfr[offset:offset+2*8*self.meta['nx']*self.meta['ny']*self.meta['nz']])
      arr = np.reshape(flat_arr,(2,self.meta['nz'],self.meta['ny'],self.meta['nx']))
      qe_par[l,:,:,:] = np.transpose(arr[0,:,:,:],(2,1,0))
      qe_per[l,:,:,:] = np.transpose(arr[1,:,:,:],(2,1,0))

    rf.close()

    if (fibo_obj != None) :
      for l in range(nexits): 
        time_exit = self.segs[seg][exit_num]
        fibo_obj.data['qe_par_'+time_exit] = qe_par[l,:,:,:]
        fibo_obj.data['qe_per_'+time_exit] = qe_per[l,:,:,:]
        fibo_obj.meta['fields']['qe_par_'+time_exit]  = ('qe_par_'+time_exit,)
        fibo_obj.meta['fields']['qe_per_'+time_exit]  = ('qe_per_'+time_exit,)
    else: return np.array([qe_par, qe_per])
    if not silent: print('done with reading qe_par and qe_per!')


  #------------------------------------------------------------
  #------------------------------------------------------------
  def calc_all(self,
      time_exit,
      fibo_obj = None,
      list_terms = ['mom_1','mom_2','mom_3'],
      also_par = False,
      also_one = True,
      also_ene = True,
      silent = True): 
    """
    Calcs all from E,B,ni,ui,Pi,sQi OR Qi_par,Qi_per ((Te_par,Te_per,qe_par,qe_per)) 
    
    Parameters :
      - time_exit              [str] time exit (usually format '0000.000')
      - fibo_obj = None        [fibo_obj] fibo to fill - if None, you will be insulted
      - list_terms = ['mom_1','mom_2','mom_3'] choose how much you want ... 
      - also_par = False       [bool] for the parallel component also ...
      - also_one = True        [bool] if you also want the one-fluid stuff
      - also_ene = True        [bool] if you also want the energies
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """


    t = time_exit
    if fibo_obj == None : print('FDP: creer un object fibo SVP')

    #rename ni as n
    #if 'ni_'+time_exit in fibo_obj.data.keys() :
    #  fibo_obj.data['n_'+time_exit] = fibo_obj.data.pop('ni_'+time_exit) 

    #create current density and magnetic field curvature
    iB = np.sqrt(np.reciprocal(fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,'B_z_'+t,'B_z_'+t)))
    bx, by, bz = iB*fibo_obj.data['B_x_'+t], iB*fibo_obj.data['B_y_'+t], iB*fibo_obj.data['B_z_'+t]

    fibo_obj.data['J_x_'+t], fibo_obj.data['J_y_'+t], fibo_obj.data['J_z_'+t] =  fibo_obj.calc_curl('B_x_'+t,'B_y_'+t,'B_z_'+t)
    fibo_obj.meta['fields']['J_'+t]  = ('J_x_'+t,'J_y_'+t,'J_z_'+t)

    fibo_obj.data['curvB_x_'+t] = bx * fibo_obj.calc_gradx(bx) + by * fibo_obj.calc_grady(bx) + bz * fibo_obj.calc_gradz(bx)
    fibo_obj.data['curvB_y_'+t] = bx * fibo_obj.calc_gradx(by) + by * fibo_obj.calc_grady(by) + bz * fibo_obj.calc_gradz(by)
    fibo_obj.data['curvB_z_'+t] = bx * fibo_obj.calc_gradx(bz) + by * fibo_obj.calc_grady(bz) + bz * fibo_obj.calc_gradz(bz)
    fibo_obj.meta['fields']['curvB_'+t]  = ('curvB_x_'+t,'curvB_y_'+t,'curvB_z_'+t)

    #local electromagnetic energy density
    if also_ene : 
      fibo_obj.data['enB_'+t]  = np.square(fibo_obj.data['B_x_'+t])
      fibo_obj.data['enB_'+t] += np.square(fibo_obj.data['B_y_'+t])
      fibo_obj.data['enB_'+t] += np.square(fibo_obj.data['B_z_'+t])
      fibo_obj.data['enB_'+t] *= 0.5
      fibo_obj.meta['fields']['enB_'+t]  = ('enB_'+t,)
  
      fibo_obj.data['enE_'+t]  = np.square(fibo_obj.data['E_x_'+t])
      fibo_obj.data['enE_'+t] += np.square(fibo_obj.data['E_y_'+t])
      fibo_obj.data['enE_'+t] += np.square(fibo_obj.data['E_z_'+t])
      fibo_obj.data['enE_'+t] *= 0.5
      fibo_obj.meta['fields']['enE_'+t]  = ('enE_'+t,)

    if also_par :
      e1, e2, e3, e4, e5, e6 = fibo_obj.calc_par_per('E_x_'+t,'B_x_'+t,'E_y_'+t,'B_y_'+t,'E_z_'+t,'B_z_'+t)
      fibo_obj.data['E_par_x_'+t] = e1
      fibo_obj.data['E_par_y_'+t] = e2
      fibo_obj.data['E_par_z_'+t] = e3
      fibo_obj.meta['fields']['E_par_'+t]  = ('E_par_x_'+t,'E_par_y_'+t,'E_par_z_'+t)

    if also_par and also_ene:
      fibo_obj.data['enE_par_'+t]  = fibo_obj.calc_scalr(e1,e1,e2,e2,e3,e3) 
      fibo_obj.data['enE_par_'+t] *= 0.5
      fibo_obj.meta['fields']['enE_par_'+t] = ('enE_par_'+t,)

    # now call the other routines to finish this dirty job
    if np.any(['mom_1' in list_terms,'mom_2' in list_terms,'mom_3' in list_terms]) : 
      self.calc_mom_1(t,fibo_obj,also_par,also_one,also_one)
      if not silent: print('done with calculating mom_1!')
    if np.any(['mom_2' in list_terms,'mom_3' in list_terms]) :
      self.calc_mom_2(t,fibo_obj,also_par,also_one)
      if not silent: print('done with calculating mom_2!')
    if ('mom_3' in list_terms) :
      self.calc_mom_3(t,fibo_obj,also_par)
      if not silent: print('done with calculating mom_3!')
    if not silent: print('done with calculating all!')


  #------------------------------------------------------------
  def calc_mom_1(self,
      time_exit,
      fibo_obj,
      also_par = False,
      also_one = True,
      also_ene = True): 
    """
      calcs first moments
    
    Parameters :
      - time_exit              [str] time exit (usually format '0000.000')
      - fibo_obj               [fibo_obj] fibo to fill
      - also_par = False       [bool] if you also want the parallel component
      - also_one = True        [bool] if you also want the one-fluid stuff
      - also_ene = True        [bool] if you also want the kinetic energies
    """

    t = time_exit

    #--bulk-velocities--------------------
    fibo_obj.data['ue_x_'+t] = fibo_obj.data['ui_x_'+t] - np.divide(fibo_obj.data['J_x_'+t],fibo_obj.data['n_'+t])
    fibo_obj.data['ue_y_'+t] = fibo_obj.data['ui_y_'+t] - np.divide(fibo_obj.data['J_y_'+t],fibo_obj.data['n_'+t])
    fibo_obj.data['ue_z_'+t] = fibo_obj.data['ui_z_'+t] - np.divide(fibo_obj.data['J_z_'+t],fibo_obj.data['n_'+t])
    fibo_obj.meta['fields']['ue_'+t] = ('ue_x_'+t,'ue_y_'+t,'ue_z_'+t)

    if also_one : 
      fibo_obj.data['u_x_'+t] = (self.meta['mime'] * fibo_obj.data['ui_x_'+t] + fibo_obj.data['ue_x_'+t]) / (1. + self.meta['mime'])
      fibo_obj.data['u_y_'+t] = (self.meta['mime'] * fibo_obj.data['ui_y_'+t] + fibo_obj.data['ue_y_'+t]) / (1. + self.meta['mime'])
      fibo_obj.data['u_z_'+t] = (self.meta['mime'] * fibo_obj.data['ui_z_'+t] + fibo_obj.data['ue_z_'+t]) / (1. + self.meta['mime'])
      fibo_obj.meta['fields']['u_'+t] = ('u_x_'+t,'u_y_'+t,'u_z_'+t)

    #--kinetic-energy-densities-----------
    if also_ene :
      fibo_obj.data['Ki_'+t]  = np.square(fibo_obj.data['ui_x_'+t])
      fibo_obj.data['Ki_'+t] += np.square(fibo_obj.data['ui_y_'+t])
      fibo_obj.data['Ki_'+t] += np.square(fibo_obj.data['ui_z_'+t])
      fibo_obj.data['Ki_'+t]  = np.multiply(fibo_obj.data['n_'+t],fibo_obj.data['Ki_'+t])/2.
      fibo_obj.meta['fields']['Ki_'+t] = ('Ki_'+t,)
  
      fibo_obj.data['Ke_'+t]  = np.square(fibo_obj.data['ue_x_'+t])
      fibo_obj.data['Ke_'+t] += np.square(fibo_obj.data['ue_y_'+t])
      fibo_obj.data['Ke_'+t] += np.square(fibo_obj.data['ue_z_'+t])
      fibo_obj.data['Ke_'+t]  = np.multiply(fibo_obj.data['n_'+t],fibo_obj.data['Ke_'+t])/(2.*self.meta['mime'])
      fibo_obj.meta['fields']['Ke_'+t] = ('Ke_'+t,)

    if also_ene and also_one :
      fibo_obj.data['K_'+t]  = np.square(fibo_obj.data['u_x_'+t])
      fibo_obj.data['K_'+t] += np.square(fibo_obj.data['u_y_'+t])
      fibo_obj.data['K_'+t] += np.square(fibo_obj.data['u_z_'+t])
      fibo_obj.data['K_'+t]  = np.multiply(fibo_obj.data['n_'+t],fibo_obj.data['K_'+t])*(self.meta['mime']+1)/(2.*self.meta['mime'])
      fibo_obj.meta['fields']['K_'+t] = ('K_'+t,)

    if also_par : 
      e1, e2, e3, e4, e5, e6 = fibo_obj.calc_par_per('ui_x_'+t,'B_x_'+t,'ui_y_'+t,'B_y_'+t,'ui_z_'+t,'B_z_'+t)
      fibo_obj.data['ui_par_x_'+t] = e1
      fibo_obj.data['ui_par_y_'+t] = e2
      fibo_obj.data['ui_par_z_'+t] = e3

    if also_par and also_ene : 
      fibo_obj.data['Ki_par_'+t]  = fibo_obj.calc_scalr(e1,e1,e2,e2,e3,e3)
      fibo_obj.data['Ki_par_'+t] *= fibo_obj.data['n_'+t] / 2.
      fibo_obj.meta['fields']['Ki_par_'+t] = ('Ki_par_'+t,)

    if also_par:
      e1, e2, e3, e4, e5, e6 = fibo_obj.calc_par_per('ue_x_'+t,'B_x_'+t,'ue_y_'+t,'B_y_'+t,'ue_z_'+t,'B_z_'+t)
      fibo_obj.data['ue_par_x_'+t] = e1
      fibo_obj.data['ue_par_y_'+t] = e2
      fibo_obj.data['ue_par_z_'+t] = e3

    if also_par and also_ene :
      fibo_obj.data['Ke_par_'+t]  = fibo_obj.calc_scalr(e1,e1,e2,e2,e3,e3)
      fibo_obj.data['Ke_par_'+t] *= fibo_obj.data['n_'+t] / (2.*fibo_obj.meta['mime'])
      fibo_obj.meta['fields']['Ke_par_'+t] = ('Ke_par_'+t,)

    if also_par and also_one :
      e1, e2, e3, e4, e5, e6 = fibo_obj.calc_par_per('u_x_'+t,'B_x_'+t,'u_y_'+t,'B_y_'+t,'u_z_'+t,'B_z_'+t)
      fibo_obj.data['u_par_x_'+t] = e1
      fibo_obj.data['u_par_y_'+t] = e2
      fibo_obj.data['u_par_z_'+t] = e3

    if also_par and also_ene and also_one :
      fibo_obj.data['K_par_'+t]  = fibo_obj.calc_scalr(e1,e1,e2,e2,e3,e3)
      fibo_obj.data['K_par_'+t] *= fibo_obj.data['n_'+t] * (fibo_obj.meta['mime']+1) / (2.*fibo_obj.meta['mime'])
      fibo_obj.meta['fields']['K_par_'+t] = ('K_par_'+t,)

  #------------------------------------------------------------
  def calc_mom_2(self,
      time_exit,
      fibo_obj = None,
      also_par = False,
      also_one = True,
      also_ene = True): 
    """
      calcs second moments
    
    Parameters :
      - time_exit              [str] time exit (usually format '0000.000')
      - fibo_obj               [fibo_obj] fibo to fill
      - also_par = False       [bool] if you also want the parallel component
      - also_one = True        [bool] if you also want the one-fluid stuff
      - also_ene = True        [bool] if you also want the internal energies
   
    """

    t = time_exit
    iB2 = np.reciprocal(fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,'B_z_'+t,'B_z_'+t))

    #--electrons---------------------------
    fibo_obj.data['Pe_xx_'+t]  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Pe_yy_'+t]  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Pe_zz_'+t]  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Pe_xy_'+t]  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Pe_xz_'+t]  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Pe_yz_'+t]  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    if 'Te_par_'+t in fibo_obj.data.keys() :
      myfac = iB2*(fibo_obj.data['Te_par_'+t] - fibo_obj.data['Te_per_'+t]) 

      fibo_obj.data['Pe_xx_'+t] += fibo_obj.data['B_x_'+t]*fibo_obj.data['B_x_'+t]*myfac + fibo_obj.data['Te_per_'+t]
      fibo_obj.data['Pe_xx_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['Pe_yy_'+t] += fibo_obj.data['B_y_'+t]*fibo_obj.data['B_y_'+t]*myfac + fibo_obj.data['Te_per_'+t]
      fibo_obj.data['Pe_yy_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['Pe_zz_'+t] += fibo_obj.data['B_z_'+t]*fibo_obj.data['B_z_'+t]*myfac + fibo_obj.data['Te_per_'+t]
      fibo_obj.data['Pe_zz_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['Pe_xy_'+t] += fibo_obj.data['B_x_'+t]*fibo_obj.data['B_y_'+t]*myfac
      fibo_obj.data['Pe_xy_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['Pe_xz_'+t] += fibo_obj.data['B_x_'+t]*fibo_obj.data['B_z_'+t]*myfac
      fibo_obj.data['Pe_xz_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['Pe_yz_'+t] += fibo_obj.data['B_y_'+t]*fibo_obj.data['B_z_'+t]*myfac
      fibo_obj.data['Pe_yz_'+t] *= fibo_obj.data['n_'+t]

    else : 
      Pe = 0.5 * self.meta['beta'] * self.meta['teti'] * fibo_obj.data['n_'+t] #- 0.5d0*beta*TeTi 

      fibo_obj.data['Pe_xx_'+t] += Pe
      fibo_obj.data['Pe_yy_'+t] += Pe
      fibo_obj.data['Pe_zz_'+t] += Pe

    #--single-fluid------------------------
    if also_one :
      fibo_obj.data['P_xx_'+t]  = + np.square(fibo_obj.data['ui_x_'+t]) 
      fibo_obj.data['P_xx_'+t] += + np.square(fibo_obj.data['ue_x_'+t]) / self.meta['mime'] 
      fibo_obj.data['P_xx_'+t] -= + np.square(fibo_obj.data['u_x_'+t]) * (1. + self.meta['mime']) / self.meta['mime'] 
      fibo_obj.data['P_yy_'+t]  = + np.square(fibo_obj.data['ui_y_'+t]) 
      fibo_obj.data['P_yy_'+t] += + np.square(fibo_obj.data['ue_y_'+t]) / self.meta['mime']
      fibo_obj.data['P_yy_'+t] -= + np.square(fibo_obj.data['u_y_'+t]) * (1. + self.meta['mime']) / self.meta['mime']
      fibo_obj.data['P_zz_'+t]  = + np.square(fibo_obj.data['ui_z_'+t]) 
      fibo_obj.data['P_zz_'+t] += + np.square(fibo_obj.data['ue_z_'+t]) / self.meta['mime']
      fibo_obj.data['P_zz_'+t] -= + np.square(fibo_obj.data['u_z_'+t]) * (1. + self.meta['mime']) / self.meta['mime']
      fibo_obj.data['P_xy_'+t]  = + np.multiply(fibo_obj.data['ui_x_'+t] , fibo_obj.data['ui_y_'+t])
      fibo_obj.data['P_xy_'+t] += + np.multiply(fibo_obj.data['ue_x_'+t] , fibo_obj.data['ue_y_'+t]) / self.meta['mime']
      fibo_obj.data['P_xy_'+t] -= + np.multiply(fibo_obj.data['u_x_'+t] ,  fibo_obj.data['u_y_'+t]) * (1. + self.meta['mime']) / self.meta['mime']
      fibo_obj.data['P_xz_'+t]  = + np.multiply(fibo_obj.data['ui_x_'+t] , fibo_obj.data['ui_z_'+t])
      fibo_obj.data['P_xz_'+t] += + np.multiply(fibo_obj.data['ue_x_'+t] , fibo_obj.data['ue_z_'+t]) / self.meta['mime']
      fibo_obj.data['P_xz_'+t] -= + np.multiply(fibo_obj.data['u_x_'+t] ,  fibo_obj.data['u_z_'+t]) * (1. + self.meta['mime']) / self.meta['mime']
      fibo_obj.data['P_yz_'+t]  = + np.multiply(fibo_obj.data['ui_y_'+t] , fibo_obj.data['ui_z_'+t])
      fibo_obj.data['P_yz_'+t] += + np.multiply(fibo_obj.data['ue_y_'+t] , fibo_obj.data['ue_z_'+t]) / self.meta['mime']
      fibo_obj.data['P_yz_'+t] -= + np.multiply(fibo_obj.data['u_y_'+t] ,  fibo_obj.data['u_z_'+t]) * (1. + self.meta['mime']) / self.meta['mime']
  
      fibo_obj.data['P_xx_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['P_yy_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['P_zz_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['P_xy_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['P_xz_'+t] *= fibo_obj.data['n_'+t]
      fibo_obj.data['P_yz_'+t] *= fibo_obj.data['n_'+t]
  
      fibo_obj.data['P_xx_'+t] += fibo_obj.data['Pi_xx_'+t] + fibo_obj.data['Pe_xx_'+t]  
      fibo_obj.data['P_yy_'+t] += fibo_obj.data['Pi_yy_'+t] + fibo_obj.data['Pe_yy_'+t]
      fibo_obj.data['P_zz_'+t] += fibo_obj.data['Pi_zz_'+t] + fibo_obj.data['Pe_zz_'+t]
      fibo_obj.data['P_xy_'+t] += fibo_obj.data['Pi_xy_'+t] + fibo_obj.data['Pe_xy_'+t] 
      fibo_obj.data['P_xz_'+t] += fibo_obj.data['Pi_xz_'+t] + fibo_obj.data['Pe_xz_'+t] 
      fibo_obj.data['P_yz_'+t] += fibo_obj.data['Pi_yz_'+t] + fibo_obj.data['Pe_yz_'+t] 

    #local internal energy density
    if also_ene :
      fibo_obj.data['Ui_'+t] = (fibo_obj.data['Pi_xx_'+t] + fibo_obj.data['Pi_yy_'+t] + fibo_obj.data['Pi_zz_'+t]) /2.
      fibo_obj.data['Ue_'+t] = (fibo_obj.data['Pe_xx_'+t] + fibo_obj.data['Pe_yy_'+t] + fibo_obj.data['Pe_zz_'+t]) /2.

    if also_ene and also_one :  
      fibo_obj.data['U_'+t]  = (fibo_obj.data['P_xx_'+t] + fibo_obj.data['P_yy_'+t] + fibo_obj.data['P_zz_'+t]) /2.

    if also_par and also_ene :
      fibo_obj.data['Ui_par_'+t]  = fibo_obj.data['B_x_'+t] * fibo_obj.calc_scalr('Pi_xx_'+t,'B_x_'+t,'Pi_xy_'+t,'B_y_'+t,'Pi_xz_'+t,'B_z_'+t)
      fibo_obj.data['Ui_par_'+t] += fibo_obj.data['B_y_'+t] * fibo_obj.calc_scalr('Pi_xy_'+t,'B_x_'+t,'Pi_yy_'+t,'B_y_'+t,'Pi_yz_'+t,'B_z_'+t)
      fibo_obj.data['Ui_par_'+t] += fibo_obj.data['B_z_'+t] * fibo_obj.calc_scalr('Pi_xz_'+t,'B_x_'+t,'Pi_yz_'+t,'B_y_'+t,'Pi_zz_'+t,'B_z_'+t)
      fibo_obj.data['Ui_par_'+t] *= iB2 / 2.

      fibo_obj.data['Ue_par_'+t]  = fibo_obj.data['B_x_'+t] * fibo_obj.calc_scalr('Pe_xx_'+t,'B_x_'+t,'Pe_xy_'+t,'B_y_'+t,'Pe_xz_'+t,'B_z_'+t)
      fibo_obj.data['Ue_par_'+t] += fibo_obj.data['B_y_'+t] * fibo_obj.calc_scalr('Pe_xy_'+t,'B_x_'+t,'Pe_yy_'+t,'B_y_'+t,'Pe_yz_'+t,'B_z_'+t)
      fibo_obj.data['Ue_par_'+t] += fibo_obj.data['B_z_'+t] * fibo_obj.calc_scalr('Pe_xz_'+t,'B_x_'+t,'Pe_yz_'+t,'B_y_'+t,'Pe_zz_'+t,'B_z_'+t)
      fibo_obj.data['Ue_par_'+t] *= iB2 / 2.

    if also_par and also_ene and also_one : 
      fibo_obj.data['U_par_'+t]   = fibo_obj.data['B_x_'+t] * fibo_obj.calc_scalr('P_xx_'+t,'B_x_'+t,'P_xy_'+t,'B_y_'+t,'P_xz_'+t,'B_z_'+t)
      fibo_obj.data['U_par_'+t]  += fibo_obj.data['B_y_'+t] * fibo_obj.calc_scalr('P_xy_'+t,'B_x_'+t,'P_yy_'+t,'B_y_'+t,'P_yz_'+t,'B_z_'+t)
      fibo_obj.data['U_par_'+t]  += fibo_obj.data['B_z_'+t] * fibo_obj.calc_scalr('P_xz_'+t,'B_x_'+t,'P_yz_'+t,'B_y_'+t,'P_zz_'+t,'B_z_'+t)
      fibo_obj.data['U_par_'+t]  *= iB2 / 2.

  #------------------------------------------------------------
  def calc_mom_3(self,
      time_exit,
      fibo_obj,
      also_par = False,
      also_one = True): 
    """
      calcs third moments
    
    Parameters :
      - time_exit              [str] time exit (usually format '0000.000')
      - fibo_obj               [fibo_obj] fibo to fill
      - also_par = False       [bool] if you also want the parallel component
      - also_one = True        [bool] if you also want the one-fluid stuff
   
    """

    t = time_exit
    iB = np.sqrt(np.reciprocal(fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,'B_z_'+t,'B_z_'+t)))
    bx, by, bz = iB*fibo_obj.data['B_x_'+t], iB*fibo_obj.data['B_y_'+t], iB*fibo_obj.data['B_z_'+t]


    #--electrons---------------------------
    fibo_obj.data['Qe_x_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Qe_y_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    fibo_obj.data['Qe_z_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    if also_par :
      fibo_obj.data['Qe_par_x_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
      fibo_obj.data['Qe_par_y_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
      fibo_obj.data['Qe_par_z_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    if 'qe_par_'+t in fibo_obj.data.keys() :

      myfac = iB * fibo_obj.data['n_'+t] * fibo_obj.data['Te_per_'+t]
      myfactor = 2 * fibo_obj.data['Te_par_'+t] * (1 - fibo_obj.data['Te_par_'+t] / fibo_obj.data['Te_per_'+t])

      #calculate the parallel heat flux
      addx = fibo_obj.calc_gradx('Te_par_'+t) - fibo_obj.data['curvB_x_'+t] * myfactor
      addy = fibo_obj.calc_grady('Te_par_'+t) - fibo_obj.data['curvB_y_'+t] * myfactor
      addz = fibo_obj.calc_gradz('Te_par_'+t) - fibo_obj.data['curvB_z_'+t] * myfactor
      
      fibo_obj.data['Qe_par_x_'+t], fibo_obj.data['Qe_par_y_'+t], fibo_obj.data['Qe_par_z_'+t] = fibo_obj.calc_cross(addx,bx,addy,by,addz,bz)

      fibo_obj.data['Qe_par_x_'+t] *= myfac 
      fibo_obj.data['Qe_par_y_'+t] *= myfac 
      fibo_obj.data['Qe_par_z_'+t] *= myfac 
      
      fibo_obj.data['Qe_par_x_'+t] += bx * fibo_obj.data['qe_par_'+t]
      fibo_obj.data['Qe_par_y_'+t] += by * fibo_obj.data['qe_par_'+t]
      fibo_obj.data['Qe_par_z_'+t] += bz * fibo_obj.data['qe_par_'+t]
      
      #calculate the perpendicular heat flux
      addx = 2 * fibo_obj.calc_gradx('Te_per_'+t) 
      addy = 2 * fibo_obj.calc_grady('Te_per_'+t) 
      addz = 2 * fibo_obj.calc_gradz('Te_per_'+t) 
      
      fibo_obj.data['Qe_x_'+t], fibo_obj.data['Qe_y_'+t], fibo_obj.data['Qe_z_'+t] = fibo_obj.calc_cross(addx,bx,addy,by,addz,bz)

      fibo_obj.data['Qe_x_'+t] *= myfac 
      fibo_obj.data['Qe_y_'+t] *= myfac 
      fibo_obj.data['Qe_z_'+t] *= myfac 
      
      fibo_obj.data['Qe_x_'+t] += bx * fibo_obj.data['qe_per_'+t]
      fibo_obj.data['Qe_y_'+t] += by * fibo_obj.data['qe_per_'+t]
      fibo_obj.data['Qe_z_'+t] += bz * fibo_obj.data['qe_per_'+t]

      #and add twice it to the parallel to get the total
      fibo_obj.data['Qe_x_'+t] = fibo_obj.data['Qe_par_x_'+t] +  fibo_obj.data['Qe_x_'+t] 
      fibo_obj.data['Qe_y_'+t] = fibo_obj.data['Qe_par_y_'+t] +  fibo_obj.data['Qe_y_'+t] 
      fibo_obj.data['Qe_z_'+t] = fibo_obj.data['Qe_par_z_'+t] +  fibo_obj.data['Qe_z_'+t] 

      if not also_par : 
        del(fibo_obj.data['Qe_par_x_'+t])
        del(fibo_obj.data['Qe_par_y_'+t])
        del(fibo_obj.data['Qe_par_z_'+t])

    fibo_obj.data['sQe_x_'+t]  = + fibo_obj.data['Qe_x_'+t] 
    fibo_obj.data['sQe_x_'+t] +=   2 * np.multiply(fibo_obj.data['Ue_'+t] + fibo_obj.data['Ke_'+t], fibo_obj.data['ue_x_'+t])
    fibo_obj.data['sQe_x_'+t] +=   2 * fibo_obj.calc_scalr('Pe_xx_'+t,'ue_x_'+t,'Pe_xy_'+t,'ue_y_'+t,'Pe_xz_'+t,'ue_z_'+t)
    fibo_obj.data['sQe_y_'+t]  = + fibo_obj.data['Qe_y_'+t] 
    fibo_obj.data['sQe_y_'+t] +=   2 * np.multiply(fibo_obj.data['Ue_'+t] + fibo_obj.data['Ke_'+t], fibo_obj.data['ue_y_'+t])
    fibo_obj.data['sQe_y_'+t] +=   2 * fibo_obj.calc_scalr('Pe_xy_'+t,'ue_x_'+t,'Pe_yy_'+t,'ue_y_'+t,'Pe_yz_'+t,'ue_z_'+t)
    fibo_obj.data['sQe_z_'+t]  = + fibo_obj.data['Qe_z_'+t] 
    fibo_obj.data['sQe_z_'+t] +=   2 * np.multiply(fibo_obj.data['Ue_'+t] + fibo_obj.data['Ke_'+t], fibo_obj.data['ue_z_'+t])
    fibo_obj.data['sQe_z_'+t] +=   2 * fibo_obj.calc_scalr('Pe_xz_'+t,'ue_x_'+t,'Pe_yz_'+t,'ue_y_'+t,'Pe_zz_'+t,'ue_z_'+t)

    if also_par : 
      upar = np.sqrt(fibo_obj.data['Ke_par_'+t] * 2.)
      fibo_obj.data['sQe_par_x_'+t]  = + fibo_obj.data['Qe_par_x_'+t] 
      fibo_obj.data['sQe_par_x_'+t] +=   2 * np.multiply(fibo_obj.data['Ue_par_'+t] + fibo_obj.data['Ke_par_'+t], fibo_obj.data['ue_x_'+t])
      fibo_obj.data['sQe_par_x_'+t] +=   2 * fibo_obj.calc_scalr('Pe_xx_'+t,bx,'Pe_xy_'+t,by,'Pe_xz_'+t,bz) * upar
      fibo_obj.data['sQe_par_y_'+t]  = + fibo_obj.data['Qe_par_y_'+t] 
      fibo_obj.data['sQe_par_y_'+t] +=   2 * np.multiply(fibo_obj.data['Ue_par_'+t] + fibo_obj.data['Ke_par_'+t], fibo_obj.data['ue_y_'+t])
      fibo_obj.data['sQe_par_y_'+t] +=   2 * fibo_obj.calc_scalr('Pe_xy_'+t,bx,'Pe_yy_'+t,by,'Pe_yz_'+t,bz) * upar
      fibo_obj.data['sQe_par_z_'+t]  = + fibo_obj.data['Qe_par_z_'+t] 
      fibo_obj.data['sQe_par_z_'+t] +=   2 * np.multiply(fibo_obj.data['Ue_par_'+t] + fibo_obj.data['Ke_par_'+t], fibo_obj.data['ue_z_'+t])
      fibo_obj.data['sQe_par_z_'+t] +=   2 * fibo_obj.calc_scalr('Pe_xz_'+t,bx,'Pe_yz_'+t,by,'Pe_zz_'+t,bz) * upar

    #--ions--------------------------------
    if 'sQi_x_'+t in fibo_obj.data.keys() : #from sQi to Qi - careful that in this case there is no way to get the Q_par
      fibo_obj.data['Qi_x_'+t]  = + fibo_obj.data['sQi_x_'+t] 
      fibo_obj.data['Qi_x_'+t] -= 2 * np.multiply(fibo_obj.data['Ui_'+t] + fibo_obj.data['Ki_'+t], fibo_obj.data['ui_x_'+t])
      fibo_obj.data['Qi_x_'+t] -= 2 * fibo_obj.calc_scalr('Pi_xx_'+t,'ui_x_'+t,'Pi_xy_'+t,'ui_y_'+t,'Pi_xz_'+t,'ui_z_'+t)
      fibo_obj.data['Qi_y_'+t]  = + fibo_obj.data['sQi_y_'+t] 
      fibo_obj.data['Qi_y_'+t] -= 2 * np.multiply(fibo_obj.data['Ui_'+t] + fibo_obj.data['Ki_'+t], fibo_obj.data['ui_y_'+t])
      fibo_obj.data['Qi_y_'+t] -= 2 * fibo_obj.calc_scalr('Pi_xy_'+t,'ui_x_'+t,'Pi_yy_'+t,'ui_y_'+t,'Pi_yz_'+t,'ui_z_'+t)
      fibo_obj.data['Qi_z_'+t]  = + fibo_obj.data['sQi_z_'+t]
      fibo_obj.data['Qi_z_'+t] -= 2 * np.multiply(fibo_obj.data['Ui_'+t] + fibo_obj.data['Ki_'+t], fibo_obj.data['ui_y_'+t])
      fibo_obj.data['Qi_z_'+t] -= 2 * fibo_obj.calc_scalr('Pi_xz_'+t,'ui_x_'+t,'Pi_yz_'+t,'ui_y_'+t,'Pi_zz_'+t,'ui_z_'+t)

      if also_par :
        fibo_obj.data['Qi_par_x_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        fibo_obj.data['Qi_par_y_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        fibo_obj.data['Qi_par_z_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])

        fibo_obj.data['sQi_par_x_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        fibo_obj.data['sQi_par_y_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])
        fibo_obj.data['sQi_par_z_'+t] = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    else : #reconstruct Qi, then from Qi to sQi
      fibo_obj.data['Qi_x_'+t] = fibo_obj.data['Qi_par_x_'+t] + 2.*fibo_obj.data['Qi_per_x_'+t]
      fibo_obj.data['Qi_y_'+t] = fibo_obj.data['Qi_par_y_'+t] + 2.*fibo_obj.data['Qi_per_y_'+t]
      fibo_obj.data['Qi_z_'+t] = fibo_obj.data['Qi_par_z_'+t] + 2.*fibo_obj.data['Qi_per_z_'+t]

      fibo_obj.data['sQi_x_'+t]  = + fibo_obj.data['Qi_x_'+t] 
      fibo_obj.data['sQi_x_'+t] += 2 * np.multiply(fibo_obj.data['Ui_'+t] + fibo_obj.data['Ki_'+t], fibo_obj.data['ui_x_'+t])
      fibo_obj.data['sQi_x_'+t] += 2 * fibo_obj.calc_scalr('Pi_xx_'+t,'ui_x_'+t,'Pi_xy_'+t,'ui_y_'+t,'Pi_xz_'+t,'ui_z_'+t)
      fibo_obj.data['sQi_y_'+t]  = + fibo_obj.data['Qi_y_'+t] 
      fibo_obj.data['sQi_y_'+t] += 2 * np.multiply(fibo_obj.data['Ui_'+t] + fibo_obj.data['Ki_'+t], fibo_obj.data['ui_y_'+t])
      fibo_obj.data['sQi_y_'+t] += 2 * fibo_obj.calc_scalr('Pi_xy_'+t,'ui_x_'+t,'Pi_yy_'+t,'ui_y_'+t,'Pi_yz_'+t,'ui_z_'+t)
      fibo_obj.data['sQi_z_'+t]  = + fibo_obj.data['Qi_z_'+t]
      fibo_obj.data['sQi_z_'+t] += 2 * np.multiply(fibo_obj.data['Ui_'+t] + fibo_obj.data['Ki_'+t], fibo_obj.data['ui_z_'+t])
      fibo_obj.data['sQi_z_'+t] += 2 * fibo_obj.calc_scalr('Pi_xz_'+t,'ui_x_'+t,'Pi_yz_'+t,'ui_y_'+t,'Pi_zz_'+t,'ui_z_'+t)

      if also_par : 
        upar = np.sqrt(fibo_obj.data['Ki_par_'+t] * 2.)
        fibo_obj.data['sQi_par_x_'+t]  = + fibo_obj.data['Qi_par_x_'+t] 
        fibo_obj.data['sQi_par_x_'+t] +=   2 * np.multiply(fibo_obj.data['Ui_par_'+t] + fibo_obj.data['Ki_par_'+t], fibo_obj.data['ui_x_'+t])
        fibo_obj.data['sQi_par_x_'+t] +=   2 * fibo_obj.calc_scalr('Pi_xx_'+t,bx,'Pi_xy_'+t,by,'Pi_xz_'+t,bz) * upar
        fibo_obj.data['sQi_par_y_'+t]  = + fibo_obj.data['Qi_par_y_'+t] 
        fibo_obj.data['sQi_par_y_'+t] +=   2 * np.multiply(fibo_obj.data['Ui_par_'+t] + fibo_obj.data['Ki_par_'+t], fibo_obj.data['ui_y_'+t])
        fibo_obj.data['sQi_par_y_'+t] +=   2 * fibo_obj.calc_scalr('Pi_xy_'+t,bx,'Pi_yy_'+t,by,'Pi_yz_'+t,bz) * upar
        fibo_obj.data['sQi_par_z_'+t]  = + fibo_obj.data['Qi_par_z_'+t] 
        fibo_obj.data['sQi_par_z_'+t] +=   2 * np.multiply(fibo_obj.data['Ui_par_'+t] + fibo_obj.data['Ki_par_'+t], fibo_obj.data['ui_z_'+t])
        fibo_obj.data['sQi_par_z_'+t] +=   2 * fibo_obj.calc_scalr('Pi_xz_'+t,bx,'Pi_yz_'+t,by,'Pi_zz_'+t,bz) * upar

    #--single-fluid------------------------
    if also_one : 
      fibo_obj.data['sQ_x_'+t] = fibo_obj.data['sQi_x_'+t] + fibo_obj.data['sQe_x_'+t]
      fibo_obj.data['sQ_y_'+t] = fibo_obj.data['sQi_y_'+t] + fibo_obj.data['sQe_y_'+t]
      fibo_obj.data['sQ_z_'+t] = fibo_obj.data['sQi_z_'+t] + fibo_obj.data['sQe_z_'+t]
  
      fibo_obj.data['Q_x_'+t]  = + fibo_obj.data['sQ_x_'+t] 
      fibo_obj.data['Q_x_'+t] -= 2 * np.multiply(fibo_obj.data['U_'+t] + fibo_obj.data['K_'+t], fibo_obj.data['u_x_'+t])
      fibo_obj.data['Q_x_'+t] -= 2 * fibo_obj.calc_scalr('P_xx_'+t,'u_x_'+t,'P_xy_'+t,'u_y_'+t,'P_xz_'+t,'u_z_'+t)
      fibo_obj.data['Q_y_'+t]  = + fibo_obj.data['sQ_y_'+t] 
      fibo_obj.data['Q_y_'+t] -= 2 * np.multiply(fibo_obj.data['U_'+t] + fibo_obj.data['K_'+t], fibo_obj.data['u_y_'+t])
      fibo_obj.data['Q_y_'+t] -= 2 * fibo_obj.calc_scalr('P_xy_'+t,'u_x_'+t,'P_yy_'+t,'u_y_'+t,'P_yz_'+t,'u_z_'+t)
      fibo_obj.data['Q_z_'+t]  = + fibo_obj.data['sQ_z_'+t] 
      fibo_obj.data['Q_z_'+t] -= 2 * np.multiply(fibo_obj.data['U_'+t] + fibo_obj.data['K_'+t], fibo_obj.data['u_z_'+t])
      fibo_obj.data['Q_z_'+t] -= 2 * fibo_obj.calc_scalr('P_xz_'+t,'u_x_'+t,'P_yz_'+t,'u_y_'+t,'P_zz_'+t,'u_z_'+t)

    if also_one and also_par : 
      fibo_obj.data['sQ_par_x_'+t] = fibo_obj.data['sQi_par_x_'+t] + fibo_obj.data['sQe_par_x_'+t]
      fibo_obj.data['sQ_par_y_'+t] = fibo_obj.data['sQi_par_y_'+t] + fibo_obj.data['sQe_par_y_'+t]
      fibo_obj.data['sQ_par_z_'+t] = fibo_obj.data['sQi_par_z_'+t] + fibo_obj.data['sQe_par_z_'+t]

      upar = np.sqrt(fibo_obj.data['K_par_'+t] * 2.)
      fibo_obj.data['Q_par_x_'+t]  = + fibo_obj.data['sQ_par_x_'+t] 
      fibo_obj.data['Q_par_x_'+t] -=   2 * np.multiply(fibo_obj.data['U_par_'+t] + fibo_obj.data['K_par_'+t], fibo_obj.data['u_x_'+t])
      fibo_obj.data['Q_par_x_'+t] -=   2 * fibo_obj.calc_scalr('P_xx_'+t,bx,'P_xy_'+t,by,'P_xz_'+t,bz) * upar
      fibo_obj.data['Q_par_y_'+t]  = + fibo_obj.data['sQ_par_y_'+t] 
      fibo_obj.data['Q_par_y_'+t] -=   2 * np.multiply(fibo_obj.data['U_par_'+t] + fibo_obj.data['K_par_'+t], fibo_obj.data['u_y_'+t])
      fibo_obj.data['Q_par_y_'+t] -=   2 * fibo_obj.calc_scalr('P_xy_'+t,bx,'P_yy_'+t,by,'P_yz_'+t,bz) * upar
      fibo_obj.data['Q_par_z_'+t]  = + fibo_obj.data['sQ_par_z_'+t] 
      fibo_obj.data['Q_par_z_'+t] -=   2 * np.multiply(fibo_obj.data['U_par_'+t] + fibo_obj.data['K_par_'+t], fibo_obj.data['u_z_'+t])
      fibo_obj.data['Q_par_z_'+t] -=   2 * fibo_obj.calc_scalr('P_xz_'+t,bx,'P_yz_'+t,by,'P_zz_'+t,bz) * upar



  #------------------------------------------------------------
  def clean_old_vars(self,
      time_exit,
      fibo_obj = None,
      silent = True): 
    """
    Deletes sQi,sQe,sQ ((Te_par, Te_per, qe_par, qe_per)) 
    
    Parameters :
      - time_exit              [str] time exit (usually format '0000.000')
      - fibo_obj = None        [fibo_obj] fibo to fill - if None, you will be insulted
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    if fibo_obj == None : print('FDP: creer un object fibo SVP')

    if 'sQi_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['sQi_x_'+time_exit]
      del fibo_obj.data['sQi_y_'+time_exit]
      del fibo_obj.data['sQi_z_'+time_exit]
    if 'sQe_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['sQe_x_'+time_exit]
      del fibo_obj.data['sQe_y_'+time_exit]
      del fibo_obj.data['sQe_z_'+time_exit]
    if 'sQ_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['sQ_x_'+time_exit]
      del fibo_obj.data['sQ_y_'+time_exit]
      del fibo_obj.data['sQ_z_'+time_exit]

    if 'sQi_par_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['sQi_par_x_'+time_exit]
      del fibo_obj.data['sQi_par_y_'+time_exit]
      del fibo_obj.data['sQi_par_z_'+time_exit]
    if 'sQe_par_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['sQe_par_x_'+time_exit]
      del fibo_obj.data['sQe_par_y_'+time_exit]
      del fibo_obj.data['sQe_par_z_'+time_exit]
    if 'sQ_par_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['sQ_par_x_'+time_exit]
      del fibo_obj.data['sQ_par_y_'+time_exit]
      del fibo_obj.data['sQ_par_z_'+time_exit]

    if 'Te_par_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['Te_par_'+time_exit]
      del fibo_obj.data['Te_per_'+time_exit]
    if 'qe_par_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['qe_par_'+time_exit]
      del fibo_obj.data['qe_per_'+time_exit]

    if not silent: print('done with cleaning sQi,sQe,sQ ((Te_par, Te_per, qe_par, qe_per)) ')

  #------------------------------------------------------------
  def clean_new_vars(self,
      time_exit,
      fibo_obj = None,
      silent = True): 
    """
    Deletes ue,Pe,Qe,sQe, u,P,Q,sQ, enE,enB, Ki,Ke,K,Ui,Ue,U 
    
    Parameters :
      - time_exit              [str] time exit (usually format '0000.000')
      - fibo_obj = None        [fibo_obj] fibo to fill - if None, you will be insulted
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    if 'ui_par_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['ui_par_x_'+time_exit]
      del fibo_obj.data['ui_par_y_'+time_exit]
      del fibo_obj.data['ui_par_z_'+time_exit]
    if 'sQi_par_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['sQi_par_x_'+time_exit]
      del fibo_obj.data['sQi_par_y_'+time_exit]
      del fibo_obj.data['sQi_par_z_'+time_exit]

    if 'ue_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['ue_x_'+time_exit]
      del fibo_obj.data['ue_y_'+time_exit]
      del fibo_obj.data['ue_z_'+time_exit]
    if 'ue_par_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['ue_par_x_'+time_exit]
      del fibo_obj.data['ue_par_y_'+time_exit]
      del fibo_obj.data['ue_par_z_'+time_exit]
    if 'Pe_xx_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Pe_xx_'+time_exit]
      del fibo_obj.data['Pe_yy_'+time_exit]
      del fibo_obj.data['Pe_zz_'+time_exit]
      del fibo_obj.data['Pe_xy_'+time_exit]
      del fibo_obj.data['Pe_xz_'+time_exit]
      del fibo_obj.data['Pe_yz_'+time_exit]
    if 'Qe_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Qe_x_'+time_exit]
      del fibo_obj.data['Qe_y_'+time_exit]
      del fibo_obj.data['Qe_z_'+time_exit]
    if 'sQe_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['sQe_x_'+time_exit]
      del fibo_obj.data['sQe_y_'+time_exit]
      del fibo_obj.data['sQe_z_'+time_exit]
    if 'Qe_par_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Qe_par_x_'+time_exit]
      del fibo_obj.data['Qe_par_y_'+time_exit]
      del fibo_obj.data['Qe_par_z_'+time_exit]
    if 'sQe_par_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['sQe_par_x_'+time_exit]
      del fibo_obj.data['sQe_par_y_'+time_exit]
      del fibo_obj.data['sQe_par_z_'+time_exit]

    if 'u_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['u_x_'+time_exit]
      del fibo_obj.data['u_y_'+time_exit]
      del fibo_obj.data['u_z_'+time_exit]
    if 'u_par_x_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['u_par_x_'+time_exit]
      del fibo_obj.data['u_par_y_'+time_exit]
      del fibo_obj.data['u_par_z_'+time_exit]
    if 'P_xx_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['P_xx_'+time_exit]
      del fibo_obj.data['P_yy_'+time_exit]
      del fibo_obj.data['P_zz_'+time_exit]
      del fibo_obj.data['P_xy_'+time_exit]
      del fibo_obj.data['P_xz_'+time_exit]
      del fibo_obj.data['P_yz_'+time_exit]
    if 'Q_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Q_x_'+time_exit]
      del fibo_obj.data['Q_y_'+time_exit]
      del fibo_obj.data['Q_z_'+time_exit]
    if 'sQ_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['sQ_x_'+time_exit]
      del fibo_obj.data['sQ_y_'+time_exit]
      del fibo_obj.data['sQ_z_'+time_exit]
    if 'Q_par_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Q_par_x_'+time_exit]
      del fibo_obj.data['Q_par_y_'+time_exit]
      del fibo_obj.data['Q_par_z_'+time_exit]
    if 'sQ_par_x_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['sQ_par_x_'+time_exit]
      del fibo_obj.data['sQ_par_y_'+time_exit]
      del fibo_obj.data['sQ_par_z_'+time_exit]
  
  
    del fibo_obj.data['enE_'+time_exit]
    del fibo_obj.data['enB_'+time_exit]
    if 'K_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Ki_'+time_exit]
      del fibo_obj.data['Ke_'+time_exit]
      del fibo_obj.data['K_'+time_exit]
    if 'K_par_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Ki_par_'+time_exit]
      del fibo_obj.data['Ke_par_'+time_exit]
      del fibo_obj.data['K_par_'+time_exit]
    if 'U_'+time_exit in fibo_obj.data.keys() :
      del fibo_obj.data['Ui_'+time_exit]
      del fibo_obj.data['Ue_'+time_exit]
      del fibo_obj.data['U_'+time_exit]
    if 'U_par_'+time_exit in fibo_obj.data.keys() :  
      del fibo_obj.data['Ui_par_'+time_exit]
      del fibo_obj.data['Ue_par_'+time_exit]
      del fibo_obj.data['U_par_'+time_exit]

    if not silent: print('done with cleaning ue,Pe,Qe,sQe, u,P,Q,sQ, enE,enB, Ki,Ke,K,Ui,Ue,U ')

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#this should read Multi-FLuid codes outputs 
#for now it reads the 2fluid_3d outputs only, but might be expanded to get also those from 3_fluid etc.
class from_MFL (object):

  def __init__(self, 
      address,
      prefix,
      nproc):
    """
    Creates the object to retrieve data from Multi-FLuid codes
    
    Parameters :
      - address      [address] where your data are (folder with segs inside)
      - prefix       [str] name of the simulation run you are using
      - nproc        [int] number of processors in the run you are using
    
    """

    self.address = address
    self.prefix = prefix
    self.nproc = nproc
    self.segs = {}
    self.meta = {}

  #------------------------------------------------------------
  def help(self):
    print('For further help, please shout:')
    print('!!!SIIIIIIIIIIIID!!!')

  #------------------------------------------------------------
  def get_meta(self,
      extra_address = '',
      silent = True):
    """
    Fills the metadata list 
    
    Parameters :
      - extra_address = ''     [address] to reach any subfolder where your meta-data are
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """
    
    #get mesh infos from input_parameters (I take the input_parameters from subfolder 01)
    infos = open(os.path.join(self.address,extra_address,self.prefix+'_cpu.dat'),'r')
    
    nx, ny, nz, nwx, nwy, nwz = map(float, infos.readline().split())
    self.meta['nx'] = int(nx / nwx)
    self.meta['ny'] = int(ny / nwy)
    self.meta['nz'] = int(nz / nwz)

    self.meta['xl'], self.meta['yl'], self.meta['zl'] = map(float, infos.readline().split())
    self.meta['dx'] = self.meta['xl'] / self.meta['nx'] #(nx_original - 1)
    self.meta['dy'] = self.meta['yl'] / self.meta['ny']
    self.meta['dz'] = self.meta['zl'] / self.meta['nz']

    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz']) 

    self.meta['dt'] = 'BOH WHO KNOWS'
    self.meta['ts'] = self.meta['dt']      #this is just for jeremy :)
    
    infos.close()
    self.meta['space_dim'] = str(np.count_nonzero(np.array(self.meta['nnn'])))+'D'

    #get time segment infos from all subdirectories: 
    segments = sorted([d for d in os.listdir(self.address) if os.path.isdir(os.path.join(self.address,d))])

    for seg in segments:
      infot = open(os.path.join(self.address,seg,self.prefix+'_cpu.dat'),'r')
      nexits = len(infot.readlines())
      self.segs[seg] = []
      infot.seek(0)
      for time_exit_num in range(nexits):
        time_exit = '%.3f' %float(infot.readline().split()[2])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        self.segs[seg].append(time_exit)
      infot.close()
    
    if not silent: print(sorted(self.segs))

    #add informations on species (sorry, some of these are hard-coded - change them in future!)
    self.meta['teti'] = 1. #HORRIBLE! PLEASE CHANGE THIS !! 

    self.meta['nss'] = 2     #number of species

    species  = []
    species.append('ion     ')
    species.append('electron')

    charges = np.zeros([self.meta['nss']])
    charges[0] = 1.
    charges[1] = -1.

    masses = np.zeros([self.meta['nss']])
    masses[0] = 1.
    masses[1] = 0. #1./self.meta['mime']

    self.meta['species']  = species
    self.meta['charges']  = { kk:charges[k] for k,kk in enumerate(species)}
    self.meta['masses']   = { kk:masses[k] for k,kk in enumerate(species)}

    if not silent : 
      print('MFL_'+self.meta['space_dim']+'> cell number               :  ', self.meta['nnn'])
      print('MFL_'+self.meta['space_dim']+'> domain size               :  ', self.meta['lll'])
      print('MFL_'+self.meta['space_dim']+'> mesh spacing              :  ', self.meta['ddd'])
      print('MFL_'+self.meta['space_dim']+'> time step                 :  ', self.meta['dt'])
      print('MFL_'+self.meta['space_dim']+'> species                   :  ', self.meta['species'])
      for i in range(self.meta['nss']):
        print('          '+species[i]+' charge                :  ', self.meta['charges'][species[i]])
        print('          '+species[i]+' mass                  :  ', self.meta['masses'][species[i]])
      print('MFL_'+self.meta['space_dim']+'> teti                      :  ', self.meta['teti'])

  #------------------------------------------------------------
  def get_B(self, 
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the B field at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nx,ny,nz = self.meta['nnn']
    Bx = np.empty([nx,ny,nz])
    By = np.empty([nx,ny,nz])
    Bz = np.empty([nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_EB.da'), 'r')
    
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        #if time_exit != self.segs[seg][exit_num] : print('=====WRONG=EXIT=====')
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('jumping time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)

        for il in range(nx*ny*nz) : 
          rf.readline()

      #fill data vectors
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      #if time_exit != self.segs[seg][exit_num] : print('=====WRONG=EXIT=====')
      iyi, iyf, izi, izf = map(int, rf.readline().split())
      if not silent: print('reading time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)

      for iz in range(-1,izf-izi):
        for iy in range(-1,iyf-iyi):
          for ix in range(nx):
            EB_line = rf.readline().split()
            Bx[ix,iy,iz] = float(EB_line[3])
            By[ix,iy,iz] = float(EB_line[4])
            Bz[ix,iy,iz] = float(EB_line[5])

      rf.close()

    if (fibo_obj != None) :
      fibo_obj.data['B_x_'+time_exit] = Bx
      fibo_obj.data['B_y_'+time_exit] = By
      fibo_obj.data['B_z_'+time_exit] = Bz
    else: return np.array([Bx, By, Bz])

    if not silent: print('done with the reading!')

  #------------------------------------------------------------
  def get_E(self,
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the B field at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns valuesd
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """
    
    #create data vectors
    nx,ny,nz = self.meta['nnn']
    Ex = np.empty([nx,ny,nz])
    Ey = np.empty([nx,ny,nz])
    Ez = np.empty([nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_EB.da'), 'r')
    
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('jumping time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)
        for il in range(nx*ny*nz) : 
          rf.readline()

      #fill data vectors
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      iyi, iyf, izi, izf = map(int, rf.readline().split())
      if not silent: print('reading time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)
      for iz in range(-1,izf-izi):
        for iy in range(-1,iyf-iyi):
          for ix in range(nx):
            EB_line = rf.readline().split()
            Ex[ix,iy,iz] = float(EB_line[0])
            Ey[ix,iy,iz] = float(EB_line[1])
            Ez[ix,iy,iz] = float(EB_line[2])

      rf.close()

    if (fibo_obj != None) :
      fibo_obj.data['E_x_'+time_exit] = Ex
      fibo_obj.data['E_y_'+time_exit] = Ey
      fibo_obj.data['E_z_'+time_exit] = Ez
    else: return np.array([Ex, Ey, Ez])

    if not silent: print('done with the reading!')

  #------------------------------------------------------------
  def get_Ui(self,
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the B field at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """

    #create data vectors
    nx,ny,nz = self.meta['nnn']
    Uix = np.empty([nx,ny,nz])
    Uiy = np.empty([nx,ny,nz])
    Uiz = np.empty([nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_U.da'), 'r')
    
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('jumping time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)
        for il in range(nx*ny*nz) : 
          rf.readline()

      #fill data vectors
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      iyi, iyf, izi, izf = map(int, rf.readline().split())
      if not silent: print('reading time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)
      for iz in range(-1,izf-izi):
        for iy in range(-1,iyf-iyi):
          for ix in range(nx):
            U_line = rf.readline().split()
            Uix[ix,iy,iz] = float(U_line[3])
            Uiy[ix,iy,iz] = float(U_line[4])
            Uiz[ix,iy,iz] = float(U_line[5])

      rf.close()

    if (fibo_obj != None) :
      fibo_obj.data['ui_x_'+time_exit] = Uix
      fibo_obj.data['ui_y_'+time_exit] = Uiy
      fibo_obj.data['ui_z_'+time_exit] = Uiz
    else: return np.array([Uix, Uiy, Uiz])

    if not silent: print('done with the reading!')

  #------------------------------------------------------------
  def get_Ue(self,
      seg,
      exit_num,
      fibo_obj = None,
      silent = True): 
    """
    Gets the B field at the nth exit of specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - exit_num               [int] number of time exit (0,1,...)
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """
    
    #create data vectors
    nx,ny,nz = self.meta['nnn']
    Uex = np.empty([nx,ny,nz])
    Uey = np.empty([nx,ny,nz])
    Uez = np.empty([nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_U.da'), 'r')
    
      #jump to the correct line in the file
      for l in range(exit_num):
        time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
        time_exit = time_exit.zfill(8)        # ..and three decimal digits
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('jumping time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)
        for il in range(nx*ny*nz) : 
          rf.readline()

      #fill data vectors
      time_exit = '%.3f' %float(rf.readline().split()[0])  #writes time exit with four integer .. 
      time_exit = time_exit.zfill(8)        # ..and three decimal digits
      iyi, iyf, izi, izf = map(int, rf.readline().split())
      if not silent: print('reading time:' , time_exit, 'procs', ip, iyi, iyf, izi, izf)
      for iz in range(-1,izf-izi):
        for iy in range(-1,iyf-iyi):
          for ix in range(nx):
            U_line = rf.readline().split()
            Uex[ix,iy,iz] = float(U_line[0])
            Uey[ix,iy,iz] = float(U_line[1])
            Uez[ix,iy,iz] = float(U_line[2])

      rf.close()

    if (fibo_obj != None) :
      fibo_obj.data['ue_x_'+time_exit] = Uex
      fibo_obj.data['ue_y_'+time_exit] = Uey
      fibo_obj.data['ue_z_'+time_exit] = Uez
    else: return np.array([Uex, Uey, Uez])

    if not silent: print('done with the reading!')

  #------------------------------------------------------------
  def get_seg_U(self,
      seg,
      fibo_obj = None,
      silent = True): 
    """
    Gets the ue,ui vector fields in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """
    
    nexits = len(self.segs[seg])
    nx,ny,nz = self.meta['nnn']

    #create data vectors
    Uex = np.empty([nexits,nx,ny,nz])
    Uey = np.empty([nexits,nx,ny,nz])
    Uez = np.empty([nexits,nx,ny,nz])
    Uix = np.empty([nexits,nx,ny,nz])
    Uiy = np.empty([nexits,nx,ny,nz])
    Uiz = np.empty([nexits,nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_U.da'), 'r')
    
      #read all data in the file
      for l in range(nexits):
        t = float(rf.readline().split()[0])
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('reading time:' , t, 'procs', ip, iyi, iyf, izi, izf)
        for iz in range(-1,izf-izi):
          for iy in range(-1,iyf-iyi):
            for ix in range(nx):
              U_line = rf.readline().split()
              Uex[l,ix,iy,iz] = float(U_line[0])
              Uey[l,ix,iy,iz] = float(U_line[1])
              Uez[l,ix,iy,iz] = float(U_line[2])
              Uix[l,ix,iy,iz] = float(U_line[3])
              Uiy[l,ix,iy,iz] = float(U_line[4])
              Uiz[l,ix,iy,iz] = float(U_line[5])

      rf.close()

    if (fibo_obj != None) :
      for l in range(nexits):
        time_exit = self.segs[seg][l]
        fibo_obj.data['ue_x_'+time_exit] = Uex[l,:,:,:]
        fibo_obj.data['ue_y_'+time_exit] = Uey[l,:,:,:]
        fibo_obj.data['ue_z_'+time_exit] = Uez[l,:,:,:]
        fibo_obj.data['ui_x_'+time_exit] = Uix[l,:,:,:]
        fibo_obj.data['ui_y_'+time_exit] = Uiy[l,:,:,:]
        fibo_obj.data['ui_z_'+time_exit] = Uiz[l,:,:,:]
    else: return np.array([Uex, Uey, Uez]), np.array([Uix, Uiy, Uiz])

    if not silent: print('done with the reading!')



  #------------------------------------------------------------
  def get_seg_EB(self,
      seg,        #segment considered (str)
      fibo_obj = None,  #fibo object you are considering - if None, EB values will just be returned 
      silent = True): 
    """
    Gets the E,B vector fields in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """
    
    nexits = len(self.segs[seg])
    nx,ny,nz = self.meta['nnn']

    #create data vectors
    Ex = np.empty([nexits,nx,ny,nz])
    Ey = np.empty([nexits,nx,ny,nz])
    Ez = np.empty([nexits,nx,ny,nz])
    Bx = np.empty([nexits,nx,ny,nz])
    By = np.empty([nexits,nx,ny,nz])
    Bz = np.empty([nexits,nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_EB.da'), 'r')
    
      #read all data in the file
      for l in range(nexits):
        t = float(rf.readline().split()[0])
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('reading time:' , t, 'procs', ip, iyi, iyf, izi, izf)
        for iz in range(-1,izf-izi):
          for iy in range(-1,iyf-iyi):
            for ix in range(nx):
              EB_line = rf.readline().split()
              Ex[l,ix,iy,iz] = float(EB_line[0])
              Ey[l,ix,iy,iz] = float(EB_line[1])
              Ez[l,ix,iy,iz] = float(EB_line[2])
              Bx[l,ix,iy,iz] = float(EB_line[3])
              By[l,ix,iy,iz] = float(EB_line[4])
              Bz[l,ix,iy,iz] = float(EB_line[5])

      rf.close()

    if (fibo_obj != None) :
      for l in range(nexits):
        time_exit = self.segs[seg][l]
        fibo_obj.data['E_x_'+time_exit] = Ex[l,:,:,:]
        fibo_obj.data['E_y_'+time_exit] = Ey[l,:,:,:]
        fibo_obj.data['E_z_'+time_exit] = Ez[l,:,:,:]
        fibo_obj.data['B_x_'+time_exit] = Bx[l,:,:,:]
        fibo_obj.data['B_y_'+time_exit] = By[l,:,:,:]
        fibo_obj.data['B_z_'+time_exit] = Bz[l,:,:,:]
    else: return np.array([Ex, Ey, Ez]), np.array([Bx, By, Bz])

    if not silent: print('done with the reading!')

  #------------------------------------------------------------
  def get_seg_DPJ(self,
      seg,
      fibo_obj = None,
      silent = True): 
    """
    Gets the n,Pe,Pi,Trac scalar fields in the specified data segment
    
    Parameters :
      - seg                    [str] segment name
      - fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
      - silent = True          [bool] don't you want to see all infos printed on shell?
    
    """    

    nexits = len(self.segs[seg])
    nx,ny,nz = self.meta['nnn']

    #create data vectors
    Den  = np.empty([nexits,nx,ny,nz])
    Pe   = np.empty([nexits,nx,ny,nz])
    Pi   = np.empty([nexits,nx,ny,nz])
    Trac = np.empty([nexits,nx,ny,nz])

    for ip in range (self.nproc):

      rf = open(os.path.join(self.address,seg,self.prefix+'_'+str(ip).zfill(4)+'_DPJ.da'), 'r')
    
      #read all data in the file
      for l in range(nexits):
        t = float(rf.readline().split()[0])
        iyi, iyf, izi, izf = map(int, rf.readline().split())
        if not silent: print('reading time:' , t, 'procs', ip, iyi, iyf, izi, izf)
        for iz in range(-1,izf-izi):
          for iy in range(-1,iyf-iyi):
            for ix in range(nx):
              DPJ_line = rf.readline().split()
              Den[l,ix,iy,iz] = float(DPJ_line[0])
              Pe[l,ix,iy,iz] = float(DPJ_line[1])
              Pi[l,ix,iy,iz] = float(DPJ_line[2])
              Trac[l,ix,iy,iz] = float(DPJ_line[3])

      rf.close()

    if (fibo_obj != None) :
      for l in range(nexits):
        time_exit = self.segs[seg][l]
        fibo_obj.data['Den_'+time_exit] = Den[l,:,:,:]
        fibo_obj.data['Pe_'+time_exit] = Pe[l,:,:,:]
        fibo_obj.data['Pi_'+time_exit] = Pi[l,:,:,:]
        fibo_obj.data['Trac_'+time_exit] = Trac[l,:,:,:]
    else: return Den, Pe, Pi, Trac

    if not silent: print('done with the reading!')


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
    self.segs = collections.OrderedDict()
    self.meta = {}

  #------------------------------------------------------------
  def get_meta(self,  #counts lines in file and calls the appropriate function to get the meta data
      extra_address = '',
      silent = True):
    """
    ------------------------------------------------------------------------------------
      fills the metadata list 
    ------------------------------------------------------------------------------------
    extra_address = ''     [address] to reach any subfolder where your meta-data are
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """

    with open(os.path.join(self.address,extra_address,'SimulationData.txt'),'r') as foo:
      line_number = len(foo.readlines())

    self.get_meta_A(extra_address)

    #---------get-dimensions-of-your-simulation----------

    self.meta['dx'] = self.meta['xl']/self.meta['nx']
    self.meta['dy'] = self.meta['yl']/self.meta['ny']
    self.meta['dz'] = self.meta['zl']/self.meta['nz']

    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz']) 
    self.meta['ppp'] = (False, False, False)             #THIS IS HARDCODED! 3-non-PERIODICITY IS HARDCODED

    self.meta['ts'] = self.meta['dt']      #this is just for jeremy :)

    self.meta['x'] = np.arange(0.,self.meta['xl'],self.meta['dx'])
    try:
        self.meta['y'] = np.arange(0.,self.meta['yl'],self.meta['dy'])
    except:
        self.meta['y'] = np.array([0.])
    try:
        self.meta['z'] = np.arange(0.,self.meta['zl'],self.meta['dz'])
    except:
        self.meta['z'] = np.array([0.])


    #----------get-time-infos-from all-h5-files----------------- 

    segments = [f for f in os.listdir(self.address) if f.split('.')[-1]=='vtk']
    for i in range(len(segments)):
        if i == 0:
          self.meta['name'] = segments[0].split('_')[0]
        segments[i] = segments[i].split('_')[-1].split('.')[0]
    segments = list(set(segments))
    segments = map(str, sorted(map(int, segments)))

    for seg in segments:
      self.segs[seg] = []
      self.segs[seg].append(float(seg)*self.meta['dt']) #PLASMA FREQ., NOT CYCLOTRON TIME!!!

    self.meta['time'] = np.concatenate(self.segs.values()).astype(float)
    self.meta['times'] = np.concatenate(self.segs.values())
    
    self.meta['time2seg'] = []
    for i in range(len(self.segs.keys())):
        self.meta['time2seg'].append(np.full((len(self.segs.values()[i])),self.segs.keys()[i]))
    self.meta['time2seg'] = np.concatenate(self.meta['time2seg'])
    
    self.meta['time2exit'] = []
    for i in range(len(self.segs.keys())):
        self.meta['time2exit'].append(np.arange(0,len(self.segs.values()[i]),1))
    self.meta['time2exit'] = np.concatenate(self.meta['time2exit'])


    #----------add-informations-on-species-----------------
    #----ACHTUNG: these are hard-coded - change them in future!----------

    self.meta['model'] = -1 #let's use negative values to indicate full kinetics models

    species  = [] #THIS BIT IS HARD-CODED, MAYBE USE stag (SPECIES TAG) INSTEAD
    species.append('ele_sw')
    species.append('ion_sw')
    species.append('ele_exo')
    species.append('ion_exo')

    msQOM = list(np.unique(self.meta['sQOM']))
    macro_species = [] #THIS BIT IS HARD-CODED
    macro_species.append('ele_sw')
    macro_species.append('ion_sw')
    macro_species.append('ele_exo')
    macro_species.append('ion_exo')

    self.meta['species']  = species
    self.meta['charges']  = {}
    self.meta['masses']   = {}
    self.meta['msQOM'] = msQOM
    self.meta['macro_species'] = macro_species

    if self.meta['ny'] == 1:
      self.meta['space_dim'] = '1D'
    elif self.meta['nz'] == 1:
      self.meta['space_dim'] = '2D'
    else:
      self.meta['space_dim'] = '3D'

    #----------print-summary-----------------

    if not silent : 
      print('iPIC_'+self.meta['space_dim']+'> cell number               :  ', self.meta['nnn'])
      print('iPIC_'+self.meta['space_dim']+'> domain size               :  ', self.meta['lll'])
      print('iPIC_'+self.meta['space_dim']+'> mesh spacing              :  ', self.meta['ddd'])
      print('iPIC_'+self.meta['space_dim']+'> periodicity               :  ', self.meta['ppp'])
      print('iPIC_'+self.meta['space_dim']+'> time step                 :  ', self.meta['ts'])
      print('iPIC_'+self.meta['space_dim']+'> species                   :  ', self.meta['species'])
      for i in range(self.meta['nss']):
        print('          '+species[i]+' charge-over-mass           :  ', self.meta['sQOM'][i])

  #------------------------------------------------------------
  def get_meta_A(self,
      extra_address = ''):
    """
    ------------------------------------------------------------------------------------
      extra routine, version A 
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
    infos.readline() #---------------------------
    self.meta['xc'] = float(infos.readline().split('=')[-1]) #box physical dimensions
    self.meta['yc'] = float(infos.readline().split('=')[-1])
    self.meta['zc'] = float(infos.readline().split('=')[-1])
    self.meta['R'] = float(infos.readline().split('=')[-1]) #box grid dimensions
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
    infos.readline() #---------------------------
    self.meta['Vx0'] = float(infos.readline().split('=')[-1]) #Vx0
    self.meta['Vy0'] = float(infos.readline().split('=')[-1]) #Vy0
    self.meta['Vz0'] = float(infos.readline().split('=')[-1]) #Vz0
    self.meta['vths'] = []
    for i in range(self.meta['nss']):
      self.meta['vths'].append(float(infos.readline().split('=')[-1])) #vth species
    print(self.meta['vths'])
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
      - tar_name           [str] target file to read (don't include '.vtk')
      - fibo_obj = None    [None or fibo] fibo object you want to fill, else returns values 
      - tar_var = None         [None or str] name the.variable will be given
      - double_y = False   [bool] was your file printed twice in y?
    
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

    self.meta['xl'] = self.meta['nx']*self.meta['dx'] #(self.n#x-1)*dx
    self.meta['yl'] = self.meta['ny']*self.meta['dy']
    self.meta['zl'] = self.meta['nz']*self.meta['dz']

    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz'])

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

    self.meta['xl'] = self.meta['nx']*self.meta['dx'] #(self.n#x-1)*dx
    self.meta['yl'] = self.meta['ny']*self.meta['dy']
    self.meta['zl'] = self.meta['nz']*self.meta['dz']

    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz'])

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

