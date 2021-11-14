# coding: utf-8

#---------------------------------------------------------------------------------------
#-(1)-----16.11.19----F-----------------------------------------------------------------
#---------------------------------------------------------------------------------------

import numpy as np
import os
import collections
import h5py

#This class read data from iPIC simulations
class from_iPIC (object):  

  def __init__(self, 
      address):
    """
    ------------------------------------------------------------------------------------
      creates the object to retrieve data from iPIC codes
    ------------------------------------------------------------------------------------
    address      [address] where your data are (folder with segs inside)
    ------------------------------------------------------------------------------------
    """

    self.address = os.path.normpath(address)
    self.segs = {}
    self.meta = {}

  #------------------------------------------------------------
  def help(self):
    print('For further help, please shout:')
    print('!!!FRAAAAAAAAAAAAAAAAAAAAAAA!!!')

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

    if line_number == 48 : self.get_meta_A(extra_address)
    else : print('FDP : unknown input_parameter file! write a new from_iPIC.get_meta() for it ...')


    #---------get-dimensions-of-your-simulation----------

    self.meta['dx'] = self.meta['xl']/self.meta['nx']
    self.meta['dy'] = self.meta['yl']/self.meta['ny']
    self.meta['dz'] = self.meta['zl']/self.meta['nz']

    self.meta['nnn'] = (self.meta['nx'], self.meta['ny'], self.meta['nz'])
    self.meta['lll'] = (self.meta['xl'], self.meta['yl'], self.meta['zl'])
    self.meta['ddd'] = (self.meta['dx'], self.meta['dy'], self.meta['dz']) 
    self.meta['ppp'] = (True, True, True)             #THIS IS HARDCODED! 3-PERIODICITY IS HARDCODED

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

    segments = [f for f in os.listdir(self.address) if f.split('.')[-1]=='h5']
    for i in range(len(segments)):
        if i == 0:
          self.meta['name'] = ''.join(segments[0].split('_')[:-1])
        segments[i] = segments[i].split('_')[-1].split('.')[0]
    segments = sorted(segments)

    for seg in segments:
      self.segs[seg] = []
      self.segs[seg].append(float(seg)*self.meta['dt']) #PLASMA FREQ., NOT CYCLOTRON TIME!!!
    self.segs = collections.OrderedDict(sorted(self.segs.items())) 
 
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

    self.meta['model'] = -1 #let's use negatuve values to indicate full kinetics models

    species  = [] #THIS BIT IS HARD-CODED, MYBE USE stag (SPECIES TAG) INSTEAD
    species.append('bg_ele')
    species.append('hs1_ele')
    species.append('hs2_ele')
    species.append('bg_ion')
    species.append('hs1_ion')
    species.append('hs2_ion')

    msQOM = list(np.unique(self.meta['sQOM']))
    macro_species = [] #THIS BIT IS HARD-CODED
    macro_species.append('ele')
    macro_species.append('ion')

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
    stag = [] #species tag
    sppp = [] #Number of particles per proc of species
    sQOM = [] #charge-over-mass ratio for species
    for i in range(self.meta['nss']):
      tmp = infos.readline().split('species')[1].split('=')
      stag.append(tmp[0].split()[0])
      sppp.append(int(tmp[1].split()[0]))
      sQOM.append(float(tmp[3]))
    self.meta['stag'] = stag
    self.meta['sppp'] = sppp
    self.meta['sQOM'] = sQOM

    infos.readline() #---------------------------
    self.meta['xl'] = float(infos.readline().split('=')[-1]) #box physical dimensions
    self.meta['yl'] = float(infos.readline().split('=')[-1])
    self.meta['zl'] = float(infos.readline().split('=')[-1])
    self.meta['nx'] = int(infos.readline().split('=')[-1]) #box grid dimensions
    self.meta['ny'] = int(infos.readline().split('=')[-1])
    self.meta['nz'] = int(infos.readline().split('=')[-1])

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

  #-----routines-for-fields------------------------------------
  #------------------------------------------------------------
  def get_EB(self,
      t_,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the E,B fields at the specified data segment
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #create data vectors
    E_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    E_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    E_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_x = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_y = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    B_z = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
      
    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])

    #fill data vectors
    E_x = np.transpose(rf['/Step#0/Block/Ex/0'][...])[:-1,:-1,:-1]
    E_y = np.transpose(rf['/Step#0/Block/Ey/0'][...])[:-1,:-1,:-1]
    E_z = np.transpose(rf['/Step#0/Block/Ez/0'][...])[:-1,:-1,:-1]
    B_x = np.transpose(rf['/Step#0/Block/Bx/0'][...])[:-1,:-1,:-1]
    B_y = np.transpose(rf['/Step#0/Block/By/0'][...])[:-1,:-1,:-1]
    B_z = np.transpose(rf['/Step#0/Block/Bz/0'][...])[:-1,:-1,:-1]

    rf.close()

    #copy all these in fibo, or return them!
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['E_x_'+time_exit] = E_x
      fibo_obj.data['E_y_'+time_exit] = E_y
      fibo_obj.data['E_z_'+time_exit] = E_z
      fibo_obj.data['B_x_'+time_exit] = B_x
      fibo_obj.data['B_y_'+time_exit] = B_y
      fibo_obj.data['B_z_'+time_exit] = B_z
    else: return np.array([E_x, E_y, E_z]), np.array([B_x, B_y, B_z])
    if not silent: print('done with reading E, B!')

  #------------------------------------------------------------
  def get_Spec(self,
      t_,
      st,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the ns and us fields at the specified data segment (for s a specified species)
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    st                     [str] species tag
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #create data vectors
    Jx_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jy_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jz_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    rho_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])
    #for 6 species with tags 0,1,2,3,4,5,
    #rf['/Step#0/Block'].keys() is equal to:
    #[u'Bx', u'By', u'Bz', 
    # u'EFx_0', u'EFx_1', u'EFx_2', u'EFx_3', u'EFx_4', u'EFx_5', 
    # u'EFy_0', u'EFy_1', u'EFy_2', u'EFy_3', u'EFy_4', u'EFy_5', 
    # u'EFz_0', u'EFz_1', u'EFz_2', u'EFz_3', u'EFz_4', u'EFz_5', 
    # u'Ex', u'Ey', u'Ez', 
    # u'Jx_0', u'Jx_1', u'Jx_2', u'Jx_3', u'Jx_4', u'Jx_5', 
    # u'Jy_0', u'Jy_1', u'Jy_2', u'Jy_3', u'Jy_4', u'Jy_5', 
    # u'Jz_0', u'Jz_1', u'Jz_2', u'Jz_3', u'Jz_4', u'Jz_5', 
    # u'Pxx_0', u'Pxx_1', u'Pxx_2', u'Pxx_3', u'Pxx_4', u'Pxx_5', 
    # u'Pxy_0', u'Pxy_1', u'Pxy_2', u'Pxy_3', u'Pxy_4', u'Pxy_5', 
    # u'Pxz_0', u'Pxz_1', u'Pxz_2', u'Pxz_3', u'Pxz_4', u'Pxz_5', 
    # u'Pyy_0', u'Pyy_1', u'Pyy_2', u'Pyy_3', u'Pyy_4', u'Pyy_5', 
    # u'Pyz_0', u'Pyz_1', u'Pyz_2', u'Pyz_3', u'Pyz_4', u'Pyz_5', 
    # u'Pzz_0', u'Pzz_1', u'Pzz_2', u'Pzz_3', u'Pzz_4', u'Pzz_5', 
    # u'rho_0', u'rho_1', u'rho_2', u'rho_3', u'rho_4', u'rho_5']

    #fill data vectors
    Jx_s  = np.transpose(rf['/Step#0/Block/Jx_' +st+'/0'][...])[:-1,:-1,:-1]
    Jy_s  = np.transpose(rf['/Step#0/Block/Jy_' +st+'/0'][...])[:-1,:-1,:-1]
    Jz_s  = np.transpose(rf['/Step#0/Block/Jz_' +st+'/0'][...])[:-1,:-1,:-1]
    rho_s = np.transpose(rf['/Step#0/Block/rho_'+st+'/0'][...])[:-1,:-1,:-1]

    #get velocities from currents
    ind = self.meta['stag'].index(st)
    Jx_s  = np.divide(Jx_s,rho_s)
    Jy_s  = np.divide(Jy_s,rho_s)
    Jz_s  = np.divide(Jz_s,rho_s)
    rho_s = np.sign(self.meta['sQOM'][ind])*rho_s

    rf.close()

    #copy all these in fibo, or return them!
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['ux_' +st+'_'+time_exit] = Jx_s
      fibo_obj.data['uy_' +st+'_'+time_exit] = Jy_s
      fibo_obj.data['uz_' +st+'_'+time_exit] = Jz_s
      fibo_obj.data['rho_'+st+'_'+time_exit] = rho_s
    else: return rho_s, np.array([Jx_s, Jy_s, Jz_s])
    if not silent: print('done with reading n_'+st+', u_'+st+'!')

  #------------------------------------------------------------
  def get_Ion(self,
      t_,
      qom = 1.0,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the n and u fields at the specified data segment (for ions OR a species
      specified via the qom)
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    qom                    [float] species charge.over-mass
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #find all species which share the given qom
    ss = [i for i,x in enumerate(self.meta['sQOM']) if x == qom]

    #create data vectors
    Jx_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jy_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jz_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    rho_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)

    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])

    #fill data vectors
    for st in np.array(self.meta['stag'])[ss]:
      Jx_qom  += np.transpose(rf['/Step#0/Block/Jx_' +st+'/0'][...])[:-1,:-1,:-1]
      Jy_qom  += np.transpose(rf['/Step#0/Block/Jy_' +st+'/0'][...])[:-1,:-1,:-1]
      Jz_qom  += np.transpose(rf['/Step#0/Block/Jz_' +st+'/0'][...])[:-1,:-1,:-1]
      rho_qom += np.transpose(rf['/Step#0/Block/rho_'+st+'/0'][...])[:-1,:-1,:-1]

    #get velocities from currents
    Jx_qom  = np.divide(Jx_qom,rho_qom)
    Jy_qom  = np.divide(Jy_qom,rho_qom)
    Jz_qom  = np.divide(Jz_qom,rho_qom)
    rho_qom = np.sign(qom)*rho_qom

    rf.close()

    #copy all these in fibo, or return them!
    ind = self.meta['msQOM'].index(qom)
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['ux_' +self.meta['macro_species'][ind]+'_'+time_exit] = Jx_qom
      fibo_obj.data['uy_' +self.meta['macro_species'][ind]+'_'+time_exit] = Jy_qom
      fibo_obj.data['uz_' +self.meta['macro_species'][ind]+'_'+time_exit] = Jz_qom
      fibo_obj.data['rho_'+self.meta['macro_species'][ind]+'_'+time_exit] = rho_qom
    else: return rho_qom, np.array([Jx_qom, Jy_qom, Jz_qom])
    if not silent: print('done with reading n_'+self.meta['macro_species'][ind]+', u_'+self.meta['macro_species'][ind]+'!')

  #------------------------------------------------------------
  def get_Pspec(self,
      t_,
      st,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the Ps fields at the specified data segment (for s a specified species)
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    st                     [str] species tag
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #create data vectors
    Pxx_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pyy_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pzz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pxy_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pxz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pyz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jx_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jy_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jz_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    rho_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])

    #fill data vectors
    Pxx_s = np.transpose(rf['/Step#0/Block/Pxx_'+st+'/0'][...])[:-1,:-1,:-1]
    Pyy_s = np.transpose(rf['/Step#0/Block/Pyy_'+st+'/0'][...])[:-1,:-1,:-1]
    Pzz_s = np.transpose(rf['/Step#0/Block/Pzz_'+st+'/0'][...])[:-1,:-1,:-1]
    Pxy_s = np.transpose(rf['/Step#0/Block/Pxy_'+st+'/0'][...])[:-1,:-1,:-1]
    Pxz_s = np.transpose(rf['/Step#0/Block/Pxz_'+st+'/0'][...])[:-1,:-1,:-1]
    Pyz_s = np.transpose(rf['/Step#0/Block/Pyz_'+st+'/0'][...])[:-1,:-1,:-1]
    Jx_s  = np.transpose(rf['/Step#0/Block/Jx_' +st+'/0'][...])[:-1,:-1,:-1]
    Jy_s  = np.transpose(rf['/Step#0/Block/Jy_' +st+'/0'][...])[:-1,:-1,:-1]
    Jz_s  = np.transpose(rf['/Step#0/Block/Jz_' +st+'/0'][...])[:-1,:-1,:-1]
    rho_s = np.transpose(rf['/Step#0/Block/rho_'+st+'/0'][...])[:-1,:-1,:-1]

    #get actual pressure: remove bulk speed: Pij = 1/qom * (Pij - JiJj/rho)
    ind = self.meta['stag'].index(st)
    qom = self.meta['sQOM'][ind]
    Pxx_s -= np.divide(Jx_s*Jx_s,rho_s)
    Pxx_s /= qom
    Pyy_s -= np.divide(Jy_s*Jy_s,rho_s)
    Pyy_s /= qom
    Pzz_s -= np.divide(Jz_s*Jz_s,rho_s)
    Pzz_s /= qom
    Pxy_s -= np.divide(Jx_s*Jy_s,rho_s)
    Pxy_s /= qom
    Pxz_s -= np.divide(Jx_s*Jz_s,rho_s)
    Pxz_s /= qom
    Pyz_s -= np.divide(Jy_s*Jz_s,rho_s)
    Pyz_s /= qom

    del Jx_s
    del Jy_s
    del Jz_s
    del rho_s

    rf.close()

    #copy all these in fibo, or return them!
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['P'+st+'_xx_'+time_exit] = Pxx_s[:,:,:]
      fibo_obj.data['P'+st+'_yy_'+time_exit] = Pyy_s[:,:,:]
      fibo_obj.data['P'+st+'_zz_'+time_exit] = Pzz_s[:,:,:]
      fibo_obj.data['P'+st+'_xy_'+time_exit] = Pxy_s[:,:,:]
      fibo_obj.data['P'+st+'_xz_'+time_exit] = Pxz_s[:,:,:]
      fibo_obj.data['P'+st+'_yz_'+time_exit] = Pyz_s[:,:,:]
    else: return np.array([[Pxx_s,Pxy_s,Pxz_s],[Pxy_s,Pyy_s,Pyz_s],[Pxz_s,Pyz_s,Pzz_s]])
    if not silent: print('done with reading P_'+st+'!')

  #------------------------------------------------------------
  def get_Press(self,
      t_,
      qom = 1.0,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the P fields at the specified data segment (for ions OR a species specified
      by qom)
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    qom                    [float] species charge-over-mass
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #find all species which share the given qom
    ss = [i for i,x in enumerate(self.meta['sQOM']) if x == qom]

    #create data vectors
    Pxx_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pyy_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pzz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pxy_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pxz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pyz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jx_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jy_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jz_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    rho_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)

    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])

    #fill data vectors
    for st in np.array(self.meta['stag'])[ss]:
      Pxx_qom += np.transpose(rf['/Step#0/Block/Pxx_'+st+'/0'][...])[:-1,:-1,:-1]
      Pyy_qom += np.transpose(rf['/Step#0/Block/Pyy_'+st+'/0'][...])[:-1,:-1,:-1]
      Pzz_qom += np.transpose(rf['/Step#0/Block/Pzz_'+st+'/0'][...])[:-1,:-1,:-1]
      Pxy_qom += np.transpose(rf['/Step#0/Block/Pxy_'+st+'/0'][...])[:-1,:-1,:-1]
      Pxz_qom += np.transpose(rf['/Step#0/Block/Pxz_'+st+'/0'][...])[:-1,:-1,:-1]
      Pyz_qom += np.transpose(rf['/Step#0/Block/Pyz_'+st+'/0'][...])[:-1,:-1,:-1]
      Jx_qom  += np.transpose(rf['/Step#0/Block/Jx_' +st+'/0'][...])[:-1,:-1,:-1]
      Jy_qom  += np.transpose(rf['/Step#0/Block/Jy_' +st+'/0'][...])[:-1,:-1,:-1]
      Jz_qom  += np.transpose(rf['/Step#0/Block/Jz_' +st+'/0'][...])[:-1,:-1,:-1]
      rho_qom += np.transpose(rf['/Step#0/Block/rho_'+st+'/0'][...])[:-1,:-1,:-1]

    #get actual pressure: remove bulk speed: Pij = 1/qom * (Pij - JiJj/rho)
    Pxx_qom -= np.divide(Jx_qom*Jx_qom,rho_qom)
    Pxx_qom /= qom                            
    Pyy_qom -= np.divide(Jy_qom*Jy_qom,rho_qom)
    Pyy_qom /= qom                            
    Pzz_qom -= np.divide(Jz_qom*Jz_qom,rho_qom)
    Pzz_qom /= qom                            
    Pxy_qom -= np.divide(Jx_qom*Jy_qom,rho_qom)
    Pxy_qom /= qom                            
    Pxz_qom -= np.divide(Jx_qom*Jz_qom,rho_qom)
    Pxz_qom /= qom                            
    Pyz_qom -= np.divide(Jy_qom*Jz_qom,rho_qom)
    Pyz_qom /= qom

    del Jx_qom
    del Jy_qom
    del Jz_qom
    del rho_qom

    rf.close()

    #copy all these in fibo, or return them!
    ind = self.meta['msQOM'].index(qom)
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['P'+self.meta['macro_species'][ind]+'_xx_'+time_exit] = Pxx_qom[:,:,:]
      fibo_obj.data['P'+self.meta['macro_species'][ind]+'_yy_'+time_exit] = Pyy_qom[:,:,:]
      fibo_obj.data['P'+self.meta['macro_species'][ind]+'_zz_'+time_exit] = Pzz_qom[:,:,:]
      fibo_obj.data['P'+self.meta['macro_species'][ind]+'_xy_'+time_exit] = Pxy_qom[:,:,:]
      fibo_obj.data['P'+self.meta['macro_species'][ind]+'_xz_'+time_exit] = Pxz_qom[:,:,:]
      fibo_obj.data['P'+self.meta['macro_species'][ind]+'_yz_'+time_exit] = Pyz_qom[:,:,:]
    else: return np.array([[Pxx_qom,Pxy_qom,Pxz_qom],[Pxy_qom,Pyy_qom,Pyz_qom],[Pxz_qom,Pyz_qom,Pzz_qom]])
    if not silent: print('done with reading P'+self.meta['macro_species'][ind]+'!')

  #------------------------------------------------------------
  def get_Qspec(self,
      t_,
      st,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the Qs fields at the specified data segment (for s a specified species)
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    st                     [str] species tag
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #create data vectors
    EFx_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    EFy_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    EFz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pxx_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pyy_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pzz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pxy_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pxz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Pyz_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jx_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jy_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    Jz_s  = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])
    rho_s = np.empty([self.meta['nx'],self.meta['ny'],self.meta['nz']])

    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])

    #fill data vectors
    EFx_s = np.transpose(rf['/Step#0/Block/EFx_'+st+'/0'][...])[:-1,:-1,:-1]
    EFy_s = np.transpose(rf['/Step#0/Block/EFy_'+st+'/0'][...])[:-1,:-1,:-1]
    EFz_s = np.transpose(rf['/Step#0/Block/EFz_'+st+'/0'][...])[:-1,:-1,:-1]
    Pxx_s = np.transpose(rf['/Step#0/Block/Pxx_'+st+'/0'][...])[:-1,:-1,:-1]
    Pyy_s = np.transpose(rf['/Step#0/Block/Pyy_'+st+'/0'][...])[:-1,:-1,:-1]
    Pzz_s = np.transpose(rf['/Step#0/Block/Pzz_'+st+'/0'][...])[:-1,:-1,:-1]
    Pxy_s = np.transpose(rf['/Step#0/Block/Pxy_'+st+'/0'][...])[:-1,:-1,:-1]
    Pxz_s = np.transpose(rf['/Step#0/Block/Pxz_'+st+'/0'][...])[:-1,:-1,:-1]
    Pyz_s = np.transpose(rf['/Step#0/Block/Pyz_'+st+'/0'][...])[:-1,:-1,:-1]
    Jx_s  = np.transpose(rf['/Step#0/Block/Jx_' +st+'/0'][...])[:-1,:-1,:-1]
    Jy_s  = np.transpose(rf['/Step#0/Block/Jy_' +st+'/0'][...])[:-1,:-1,:-1]
    Jz_s  = np.transpose(rf['/Step#0/Block/Jz_' +st+'/0'][...])[:-1,:-1,:-1]
    rho_s = np.transpose(rf['/Step#0/Block/rho_'+st+'/0'][...])[:-1,:-1,:-1]

    #get actual pressure: remove bulk speed: Pij = 1/qom * (Pij - JiJj/rho)
    ind = self.meta['stag'].index(st)
    qom = self.meta['sQOM'][ind]
    Pxx_s -= np.divide(Jx_s*Jx_s,rho_s)
    Pxx_s /= qom
    Pyy_s -= np.divide(Jy_s*Jy_s,rho_s)
    Pyy_s /= qom
    Pzz_s -= np.divide(Jz_s*Jz_s,rho_s)
    Pzz_s /= qom
    Pxy_s -= np.divide(Jx_s*Jy_s,rho_s)
    Pxy_s /= qom
    Pxz_s -= np.divide(Jx_s*Jz_s,rho_s)
    Pxz_s /= qom
    Pyz_s -= np.divide(Jy_s*Jz_s,rho_s)
    Pyz_s /= qom
    J2_s = Jx_s*Jx_s+Jy_s*Jy_s+Jz_s*Jz_s
    TrP_s = Pxx_s+Pyy_s+Pzz_s
    EFx_s = (EFx_s - np.divide(TrP_s*Jx_s,2.*rho_s) - np.divide(J2_s*Jx_s,2.*qom*rho_s*rho_s) 
            - np.divide(Jx_s*Pxx_s+Jy_s*Pxy_s+Jz_s*Pxz_s,rho_s))
    EFy_s = (EFy_s - np.divide(TrP_s*Jy_s,2.*rho_s) - np.divide(J2_s*Jy_s,2.*qom*rho_s*rho_s) 
            - np.divide(Jx_s*Pxy_s+Jy_s*Pyy_s+Jz_s*Pyz_s,rho_s))
    EFz_s = (EFz_s - np.divide(TrP_s*Jz_s,2.*rho_s) - np.divide(J2_s*Jz_s,2.*qom*rho_s*rho_s) 
            - np.divide(Jx_s*Pxz_s+Jy_s*Pyz_s+Jz_s*Pzz_s,rho_s))

    del Pxx_s
    del Pyy_s
    del Pzz_s
    del Pxy_s
    del Pxz_s
    del Pyz_s
    del TrP_s
    del Jx_s
    del Jy_s
    del Jz_s
    del J2_s
    del rho_s

    rf.close()

    #copy all these in fibo, or return them!
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['Q'+st+'_x_'+time_exit] = EFx_s[:,:,:]
      fibo_obj.data['Q'+st+'_y_'+time_exit] = EFy_s[:,:,:]
      fibo_obj.data['Q'+st+'_z_'+time_exit] = EFz_s[:,:,:]
    else: return np.array([EFx_s,EFy_s,EFz_s])
    if not silent: print('done with reading Q'+st+'!')

  #------------------------------------------------------------
  def get_Q(self,
      t_,
      qom = 1.0,
      fibo_obj = None,
      silent = True):
    """
    ------------------------------------------------------------------------------------
      gets the Q fields at the specified data segment (for ions OR a species specified
      by qom)
    ------------------------------------------------------------------------------------
    t_                     [int/long/float/double] time exit identifier. If int or 
                           long, will be used as time index, or else as time value
    qom                    [float] species charge-over-mass
    fibo_obj = None        [None OR fibo_obj] fibo to fill - if None, returns values
    silent = True          [bool] don't you want to see all infos printed on shell?
    ------------------------------------------------------------------------------------
    """
    #processing time input
    if type(t_) != float:
        ind = t_ #if t_ is intended as the (global) time index
    else:
        ind = np.argmin(np.abs(self.meta['time'] - t_)) #if t_ is intended as a time value. 
                                                        #The closest time will be used
    if not silent:
        print('Closest time is: '+str(self.meta['time'][ind]))
    seg = self.meta['time2seg'][ind]

    #find all species which share the given qom
    ss = [i for i,x in enumerate(self.meta['sQOM']) if x == qom]

    #create data vectors
    EFx_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    EFy_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    EFz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pxx_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pyy_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pzz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pxy_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pxz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Pyz_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jx_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jy_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    Jz_qom  = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)
    rho_qom = np.zeros([self.meta['nx'],self.meta['ny'],self.meta['nz']],dtype=np.float)

    #NB hdf5 file
    rf = h5py.File(os.path.join(self.address,self.meta['name']+'_'+seg+'.h5'), 'r')
    if not silent: print('reading time:' , self.segs[seg][0])

    #fill data vectors
    for st in np.array(self.meta['stag'])[ss]:
      EFx_qom += np.transpose(rf['/Step#0/Block/EFx_'+st+'/0'][...])[:-1,:-1,:-1]
      EFy_qom += np.transpose(rf['/Step#0/Block/EFy_'+st+'/0'][...])[:-1,:-1,:-1]
      EFz_qom += np.transpose(rf['/Step#0/Block/EFz_'+st+'/0'][...])[:-1,:-1,:-1]
      Pxx_qom += np.transpose(rf['/Step#0/Block/Pxx_'+st+'/0'][...])[:-1,:-1,:-1]
      Pyy_qom += np.transpose(rf['/Step#0/Block/Pyy_'+st+'/0'][...])[:-1,:-1,:-1]
      Pzz_qom += np.transpose(rf['/Step#0/Block/Pzz_'+st+'/0'][...])[:-1,:-1,:-1]
      Pxy_qom += np.transpose(rf['/Step#0/Block/Pxy_'+st+'/0'][...])[:-1,:-1,:-1]
      Pxz_qom += np.transpose(rf['/Step#0/Block/Pxz_'+st+'/0'][...])[:-1,:-1,:-1]
      Pyz_qom += np.transpose(rf['/Step#0/Block/Pyz_'+st+'/0'][...])[:-1,:-1,:-1]
      Jx_qom  += np.transpose(rf['/Step#0/Block/Jx_' +st+'/0'][...])[:-1,:-1,:-1]
      Jy_qom  += np.transpose(rf['/Step#0/Block/Jy_' +st+'/0'][...])[:-1,:-1,:-1]
      Jz_qom  += np.transpose(rf['/Step#0/Block/Jz_' +st+'/0'][...])[:-1,:-1,:-1]
      rho_qom += np.transpose(rf['/Step#0/Block/rho_'+st+'/0'][...])[:-1,:-1,:-1]

    #get actual pressure: remove bulk speed: Pij = 1/qom * (Pij - JiJj/rho)
    Pxx_qom -= np.divide(Jx_qom*Jx_qom,rho_qom)
    Pxx_qom /= qom                            
    Pyy_qom -= np.divide(Jy_qom*Jy_qom,rho_qom)
    Pyy_qom /= qom                            
    Pzz_qom -= np.divide(Jz_qom*Jz_qom,rho_qom)
    Pzz_qom /= qom                            
    Pxy_qom -= np.divide(Jx_qom*Jy_qom,rho_qom)
    Pxy_qom /= qom                            
    Pxz_qom -= np.divide(Jx_qom*Jz_qom,rho_qom)
    Pxz_qom /= qom                            
    Pyz_qom -= np.divide(Jy_qom*Jz_qom,rho_qom)
    Pyz_qom /= qom
    J2_qom = Jx_qom*Jx_qom+Jy_qom*Jy_qom+Jz_qom*Jz_qom
    TrP_qom = Pxx_qom+Pyy_qom+Pzz_qom
    EFx_qom = (EFx_qom - np.divide(TrP_qom*Jx_qom,2.*rho_qom) 
            - np.divide(J2_qom*Jx_qom,2.*qom*rho_qom*rho_qom) 
            - np.divide(Jx_qom*Pxx_qom+Jy_qom*Pxy_qom+Jz_qom*Pxz_qom,rho_qom))
    EFy_qom = (EFy_qom - np.divide(TrP_qom*Jy_qom,2.*rho_qom) 
            - np.divide(J2_qom*Jy_qom,2.*qom*rho_qom*rho_qom) 
            - np.divide(Jx_qom*Pxy_qom+Jy_qom*Pyy_qom+Jz_qom*Pyz_qom,rho_qom))
    EFz_qom = (EFz_qom - np.divide(TrP_qom*Jz_qom,2.*rho_qom) 
            - np.divide(J2_qom*Jz_qom,2.*qom*rho_qom*rho_qom) 
            - np.divide(Jx_qom*Pxz_qom+Jy_qom*Pyz_qom+Jz_qom*Pzz_qom,rho_qom))

    del Pxx_qom
    del Pyy_qom
    del Pzz_qom
    del Pxy_qom
    del Pxz_qom
    del Pyz_qom
    del TrP_qom
    del Jx_qom
    del Jy_qom
    del Jz_qom
    del rho_qom

    rf.close()

    #copy all these in fibo, or return them!
    ind = self.meta['msQOM'].index(qom)
    if (fibo_obj != None) :
      time_exit = self.segs[seg][0]
      fibo_obj.data['Q'+self.meta['macro_species'][ind]+'_x_'+time_exit] = EFx_qom[:,:,:]
      fibo_obj.data['Q'+self.meta['macro_species'][ind]+'_y_'+time_exit] = EFy_qom[:,:,:]
      fibo_obj.data['Q'+self.meta['macro_species'][ind]+'_z_'+time_exit] = EFz_qom[:,:,:]
    else: return np.array([EFx_qom,EFy_qom,EFz_qom])
    if not silent: print('done with reading Q'+self.meta['macro_species'][ind]+'!')

