#---------------------------------------------------------------------------------------
#------do-some-physics-with-fibo-objects------------------------------------------------
#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

import numpy as np


class phybo (object):

  def __init__(self, 
      fibo_obj): 
    """
    Basically a repository for all your scientific routines 
    
    Parameters :
      - fibo_obj    [fibo] whose data and meta-data will be used in computations 
    
    """

    self.fibo_obj = fibo_obj

  #------------------------------------------------------------
  def calc_ref_change(self,
      time_str,
      ref_point,
      ref_vel_x,
      ref_vel_y,
      ref_vel_z):
    """
    Transforms reference frame so that in ref_point single-fluid is still
    
    Parameters :
      - time_str     [str] to complete the var names 
      - ref_point    [int,int,int] choose the reference point 
      - ref_vel_x    [fibo.data] choode the reference velocity, x comp.
      - ref_vel_y    [fibo.data] choode the reference velocity, y comp.
      - ref_vel_z    [fibo.data] choode the reference velocity, z comp.
    
    """
    
    nx,ny,nz = self.fibo_obj.meta['nnn']
    t = time_str
    
    
    #calculate barycenter plasma velocity at the ref_point
    du_x = np.ones([nx,ny,nz]) * self.fibo_obj.data[ref_vel_x][ref_point[0],ref_point[1],ref_point[2]]
    du_y = np.ones([nx,ny,nz]) * self.fibo_obj.data[ref_vel_y][ref_point[0],ref_point[1],ref_point[2]]
    du_z = np.ones([nx,ny,nz]) * self.fibo_obj.data[ref_vel_z][ref_point[0],ref_point[1],ref_point[2]]
    
    #transform ui
    self.fibo_obj.data['ui_x_'+t] = self.fibo_obj.data['ui_x_'+t] - du_x
    self.fibo_obj.data['ui_y_'+t] = self.fibo_obj.data['ui_y_'+t] - du_y
    self.fibo_obj.data['ui_z_'+t] = self.fibo_obj.data['ui_z_'+t] - du_z
    
    #transform E
    myd_x,myd_y,myd_z = self.fibo_obj.calc_cross(du_x,'B_x_'+t,du_y,'B_y_'+t,du_z,'B_z_'+t)
    
    self.fibo_obj.data['E_x_'+t] = self.fibo_obj.data['E_x_'+t] + myd_x 
    self.fibo_obj.data['E_y_'+t] = self.fibo_obj.data['E_y_'+t] + myd_y 
    self.fibo_obj.data['E_z_'+t] = self.fibo_obj.data['E_z_'+t] + myd_z 
    


  #------------------------------------------------------------
  def calc_Psi(self,
      time_str,
      cut_z = 0):
    """
    Calculates magnetic flux function 
    
    Parameters :
      - time_str     [str] to complete the var names 
      - cut_z = 0    [int] choose the level for this computation
    
    """

    xl,yl,zl = self.fibo_obj.meta['lll']
    nx,ny,nz = self.fibo_obj.meta['nnn']
    nym = ny//2
    t = time_str

    mtrc1 = yl / (2. * np.pi) 
    mtrc2 = xl / (nx) #substitute nx wit add -1 ????

    #a1 = np.zeros([ny],dtype=complex)
    a2 = np.zeros([ny],dtype=complex)
    #a3 = np.zeros([ny],dtype=complex)
    self.fibo_obj.data['Psi_'+t] = np.zeros([nx,ny,1],dtype=float)

    a2[0] = 0.0
    for ix in range(nx):
      #transform the array
      a1 = np.fft.fft(self.fibo_obj.data['B_x_'+t][ix,:,cut_z]) #* 2 / (nx*ny*nz)

      #process all elements in the array least for a2[0]
      for j in range(1,nym+1):
        a2[j] = -1j * a1[j] / float(j) * mtrc1
      for j in range(nym+1,ny):
        a2[j] =  1j * a1[j] / float(ny-j) * mtrc1

      #process a2[0]
      a3 = np.fft.fft(self.fibo_obj.data['B_y_'+t][ix,:,cut_z]) #* 2 / (nx*ny*nz)
      a2[0] -= a3[0]*mtrc2
      
      #anti-transform
      self.fibo_obj.data['Psi_'+t][ix,:,cut_z] = np.real(np.fft.ifft(a2)) #* 2 / (nx*ny*nz)


  #------------------------------------------------------------
  def calc_lengths(self,
      time_str):
    """
    Calculates inertial and gyration lengths for ions and electrons
    
    Parameters :
      - time_str     [str] to complete the var names 
    
    """

    t = time_str
    iB = np.reciprocal(np.sqrt( self.fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t, 'B_y_'+t,'B_y_'+t, 'B_z_'+t,'B_z_'+t) ))

    if 'Ue_par_'+t not in self.fibo_obj.data.keys() : 
      self.fibo_obj.data['Ue_par_'+t]  = self.fibo_obj.data['B_x_'+t] * self.fibo_obj.calc_scalr('Pe_xx_'+t,'B_x_'+t,'Pe_xy_'+t,'B_y_'+t,'Pe_xz_'+t,'B_z_'+t)
      self.fibo_obj.data['Ue_par_'+t] += self.fibo_obj.data['B_y_'+t] * self.fibo_obj.calc_scalr('Pe_xy_'+t,'B_x_'+t,'Pe_yy_'+t,'B_y_'+t,'Pe_yz_'+t,'B_z_'+t)
      self.fibo_obj.data['Ue_par_'+t] += self.fibo_obj.data['B_z_'+t] * self.fibo_obj.calc_scalr('Pe_xz_'+t,'B_x_'+t,'Pe_yz_'+t,'B_y_'+t,'Pe_zz_'+t,'B_z_'+t)
      self.fibo_obj.data['Ue_par_'+t] *= iB*iB / 2.
    
    if 'Ui_par_'+t not in self.fibo_obj.data.keys() : 
      self.fibo_obj.data['Ui_par_'+t]  = self.fibo_obj.data['B_x_'+t] * self.fibo_obj.calc_scalr('Pi_xx_'+t,'B_x_'+t,'Pi_xy_'+t,'B_y_'+t,'Pi_xz_'+t,'B_z_'+t)
      self.fibo_obj.data['Ui_par_'+t] += self.fibo_obj.data['B_y_'+t] * self.fibo_obj.calc_scalr('Pi_xy_'+t,'B_x_'+t,'Pi_yy_'+t,'B_y_'+t,'Pi_yz_'+t,'B_z_'+t)
      self.fibo_obj.data['Ui_par_'+t] += self.fibo_obj.data['B_z_'+t] * self.fibo_obj.calc_scalr('Pi_xz_'+t,'B_x_'+t,'Pi_yz_'+t,'B_y_'+t,'Pi_zz_'+t,'B_z_'+t)
      self.fibo_obj.data['Ui_par_'+t] *= iB*iB / 2.

  
    #thermal velocity of ions and electrons
    self.fibo_obj.data['cTi_'+t] = 2. * np.sqrt(2./3. * self.fibo_obj.data['Ui_'+t] / self.fibo_obj.data['n_'+t]) 
    self.fibo_obj.data['cTe_'+t] = 2. * np.sqrt(2./3. * self.fibo_obj.data['Ue_'+t] / self.fibo_obj.data['n_'+t]) * self.fibo_obj.meta['mime']
    #perpendicular thermal velocity of ions and electrons
    self.fibo_obj.data['cTi_per_'+t] = self.fibo_obj.data['Ui_'+t] - self.fibo_obj.data['Ui_par_'+t]
    self.fibo_obj.data['cTe_per_'+t] = self.fibo_obj.data['Ue_'+t] - self.fibo_obj.data['Ue_par_'+t]
    self.fibo_obj.data['cTi_per_'+t] = 2. * np.sqrt(2./3. * self.fibo_obj.data['cTi_per_'+t] / self.fibo_obj.data['n_'+t]) 
    self.fibo_obj.data['cTe_per_'+t] = 2. * np.sqrt(2./3. * self.fibo_obj.data['cTe_per_'+t] / self.fibo_obj.data['n_'+t] * self.fibo_obj.meta['mime']) 
    #perpendicular gyration lengths of ions and electrons (Larmor radia)
    self.fibo_obj.data['li_g_per_'+t] = self.fibo_obj.data['cTi_per_'+t] * iB
    self.fibo_obj.data['le_g_per_'+t] = self.fibo_obj.data['cTe_per_'+t] * iB / self.fibo_obj.meta['mime']
    #parallel gyration lengths would give an idea of how much interaction along the line is present (maybe to add)

    #inertial lengths of ions and electrons (I fear that when we have quasineutrality the two are not that interesting together)
    self.fibo_obj.data['li_d_'+t] = np.sqrt(np.reciprocal(self.fibo_obj.data['n_'+t]))
    self.fibo_obj.data['le_d_'+t] = self.fibo_obj.data['li_d_'+t] * np.sqrt(1./self.fibo_obj.meta['mime'])
    



  #------------------------------------------------------------  
  def calc_energy_conv(self,
      time_str,
      par_per = False, 
      iso_Te = True):
    """
    Calculates all energy conversion rates
    
    Parameters :
      - time_str     [str] to complete the var names 
      - par_per      [bool] do you want also par ((and per))? 
      - iso_Te       [bool] isothermal code?
    
    """

    nx,ny,nz = self.fibo_obj.meta['nnn']
    t = time_str

    #divergence of the velocity field
    self.fibo_obj.data['div_ui_'+t] = self.fibo_obj.calc_divr('ui_x_'+t,'ui_y_'+t,'ui_z_'+t)
    self.fibo_obj.data['div_ue_'+t] = self.fibo_obj.calc_divr('ue_x_'+t,'ue_y_'+t,'ue_z_'+t)
    self.fibo_obj.data['div_u_'+t]  = self.fibo_obj.calc_divr('u_x_'+t,'u_y_'+t,'u_z_'+t)

    self.fibo_obj.data['uXB_x_'+t], self.fibo_obj.data['uXB_y_'+t], self.fibo_obj.data['uXB_z_'+t] = self.fibo_obj.calc_cross('u_x_'+t,'B_x_'+t,'u_y_'+t,'B_y_'+t,'u_z_'+t,'B_z_'+t)

    if par_per:
      iB2 = np.reciprocal(self.fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,'B_z_'+t,'B_z_'+t))
      iB = np.sqrt(iB2)
      ui_par = np.sqrt(self.fibo_obj.calc_scalr('ui_par_x_'+t,'ui_par_x_'+t,'ui_par_y_'+t,'ui_par_y_'+t,'ui_par_z_'+t,'ui_par_z_'+t))
      ue_par = np.sqrt(self.fibo_obj.calc_scalr('ue_par_x_'+t,'ue_par_x_'+t,'ue_par_y_'+t,'ue_par_y_'+t,'ue_par_z_'+t,'ue_par_z_'+t))
      u_par  = np.sqrt(self.fibo_obj.calc_scalr('u_par_x_'+t,'u_par_x_'+t,'u_par_y_'+t,'u_par_y_'+t,'u_par_z_'+t,'u_par_z_'+t))

      self.fibo_obj.data['-dB_t_x_'+t],self.fibo_obj.data['-dB_t_y_'+t],self.fibo_obj.data['-dB_t_z_'+t] = self.fibo_obj.calc_curl('E_x_'+t,'E_y_'+t,'E_z_'+t)

      self.fibo_obj.data['dB_ic_x_'+t]  = self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.calc_gradx('B_x_'+t)
      self.fibo_obj.data['dB_ic_x_'+t] += self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.calc_grady('B_x_'+t)
      self.fibo_obj.data['dB_ic_x_'+t] += self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.calc_gradz('B_x_'+t)
      self.fibo_obj.data['dB_ic_x_'+t] -= self.fibo_obj.data['-dB_t_x_'+t]
      self.fibo_obj.data['dB_ic_y_'+t]  = self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.calc_gradx('B_y_'+t)
      self.fibo_obj.data['dB_ic_y_'+t] += self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.calc_grady('B_y_'+t)
      self.fibo_obj.data['dB_ic_y_'+t] += self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.calc_gradz('B_y_'+t)
      self.fibo_obj.data['dB_ic_x_'+t] -= self.fibo_obj.data['-dB_t_y_'+t]
      self.fibo_obj.data['dB_ic_z_'+t]  = self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.calc_gradx('B_z_'+t)
      self.fibo_obj.data['dB_ic_z_'+t] += self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.calc_grady('B_z_'+t)
      self.fibo_obj.data['dB_ic_z_'+t] += self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.calc_gradz('B_z_'+t)
      self.fibo_obj.data['dB_ic_x_'+t] -= self.fibo_obj.data['-dB_t_z_'+t]

      self.fibo_obj.data['dB_ec_x_'+t]  = self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.calc_gradx('B_x_'+t)
      self.fibo_obj.data['dB_ec_x_'+t] += self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.calc_grady('B_x_'+t)
      self.fibo_obj.data['dB_ec_x_'+t] += self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.calc_gradz('B_x_'+t)
      self.fibo_obj.data['dB_ec_x_'+t] -= self.fibo_obj.data['-dB_t_x_'+t]
      self.fibo_obj.data['dB_ec_y_'+t]  = self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.calc_gradx('B_y_'+t)
      self.fibo_obj.data['dB_ec_y_'+t] += self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.calc_grady('B_y_'+t)
      self.fibo_obj.data['dB_ec_y_'+t] += self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.calc_gradz('B_y_'+t)
      self.fibo_obj.data['dB_ec_x_'+t] -= self.fibo_obj.data['-dB_t_y_'+t]
      self.fibo_obj.data['dB_ec_z_'+t]  = self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.calc_gradx('B_z_'+t)
      self.fibo_obj.data['dB_ec_z_'+t] += self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.calc_grady('B_z_'+t)
      self.fibo_obj.data['dB_ec_z_'+t] += self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.calc_gradz('B_z_'+t)
      self.fibo_obj.data['dB_ec_x_'+t] -= self.fibo_obj.data['-dB_t_z_'+t]

      self.fibo_obj.data['dB_c_x_'+t]  = self.fibo_obj.data['u_x_'+t] * self.fibo_obj.calc_gradx('B_x_'+t)
      self.fibo_obj.data['dB_c_x_'+t] += self.fibo_obj.data['u_y_'+t] * self.fibo_obj.calc_grady('B_x_'+t)
      self.fibo_obj.data['dB_c_x_'+t] += self.fibo_obj.data['u_z_'+t] * self.fibo_obj.calc_gradz('B_x_'+t)
      self.fibo_obj.data['dB_c_x_'+t] -= self.fibo_obj.data['-dB_t_x_'+t]
      self.fibo_obj.data['dB_c_y_'+t]  = self.fibo_obj.data['u_x_'+t] * self.fibo_obj.calc_gradx('B_y_'+t)
      self.fibo_obj.data['dB_c_y_'+t] += self.fibo_obj.data['u_y_'+t] * self.fibo_obj.calc_grady('B_y_'+t)
      self.fibo_obj.data['dB_c_y_'+t] += self.fibo_obj.data['u_z_'+t] * self.fibo_obj.calc_gradz('B_y_'+t)
      self.fibo_obj.data['dB_c_x_'+t] -= self.fibo_obj.data['-dB_t_y_'+t]
      self.fibo_obj.data['dB_c_z_'+t]  = self.fibo_obj.data['u_x_'+t] * self.fibo_obj.calc_gradx('B_z_'+t)
      self.fibo_obj.data['dB_c_z_'+t] += self.fibo_obj.data['u_y_'+t] * self.fibo_obj.calc_grady('B_z_'+t)
      self.fibo_obj.data['dB_c_z_'+t] += self.fibo_obj.data['u_z_'+t] * self.fibo_obj.calc_gradz('B_z_'+t)   
      self.fibo_obj.data['dB_c_x_'+t] -= self.fibo_obj.data['-dB_t_z_'+t]

      #calculation of the Yi, Ye and Y tensors:
      self.fibo_obj.data['Yi_xx_'+t] = self.fibo_obj.data['dB_ic_x_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Yi_xx_'+t] *= iB2
      self.fibo_obj.data['Yi_yy_'+t] = self.fibo_obj.data['dB_ic_y_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Yi_yy_'+t] *= iB2
      self.fibo_obj.data['Yi_zz_'+t] = self.fibo_obj.data['dB_ic_z_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Yi_zz_'+t] *= iB2
      self.fibo_obj.data['Yi_xy_'+t]  = self.fibo_obj.data['dB_ic_x_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Yi_xy_'+t] += self.fibo_obj.data['dB_ic_y_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Yi_xy_'+t] *= 0.5 * iB2
      self.fibo_obj.data['Yi_xz_'+t]  = self.fibo_obj.data['dB_ic_x_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Yi_xz_'+t] += self.fibo_obj.data['dB_ic_z_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Yi_xz_'+t] *= 0.5 * iB2
      self.fibo_obj.data['Yi_yz_'+t]  = self.fibo_obj.data['dB_ic_y_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Yi_yz_'+t] += self.fibo_obj.data['dB_ic_z_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Yi_yz_'+t] *= 0.5 * iB2

      myfac = iB2 * iB2 * self.fibo_obj.calc_scalr('B_x_'+t,'dB_ic_x_'+t,'B_y_'+t,'dB_ic_y_'+t,'B_z_'+t,'dB_ic_z_'+t)
      self.fibo_obj.data['Yi_xx_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Yi_yy_'+t] -= myfac * self.fibo_obj.data['B_y_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Yi_zz_'+t] -= myfac * self.fibo_obj.data['B_z_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Yi_xy_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Yi_xz_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Yi_yz_'+t] -= myfac * self.fibo_obj.data['B_y_'+t] * self.fibo_obj.data['B_z_'+t]

      self.fibo_obj.data['Ye_xx_'+t] = self.fibo_obj.data['dB_ec_x_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Ye_xx_'+t] *= iB2 
      self.fibo_obj.data['Ye_yy_'+t] = self.fibo_obj.data['dB_ec_y_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Ye_yy_'+t] *= iB2
      self.fibo_obj.data['Ye_zz_'+t] = self.fibo_obj.data['dB_ec_z_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Ye_zz_'+t] *= iB2
      self.fibo_obj.data['Ye_xy_'+t]  = self.fibo_obj.data['dB_ec_x_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Ye_xy_'+t] += self.fibo_obj.data['dB_ec_y_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Ye_xy_'+t] *= 0.5 * iB2
      self.fibo_obj.data['Ye_xz_'+t]  = self.fibo_obj.data['dB_ec_x_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Ye_xz_'+t] += self.fibo_obj.data['dB_ec_z_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Ye_xz_'+t] *= 0.5 * iB2
      self.fibo_obj.data['Ye_yz_'+t]  = self.fibo_obj.data['dB_ec_y_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Ye_yz_'+t] += self.fibo_obj.data['dB_ec_z_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Ye_yz_'+t] *= 0.5 * iB2

      myfac = iB2 * iB2 * self.fibo_obj.calc_scalr('B_x_'+t,'dB_ec_x_'+t,'B_y_'+t,'dB_ec_y_'+t,'B_z_'+t,'dB_ec_z_'+t)
      self.fibo_obj.data['Ye_xx_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Ye_yy_'+t] -= myfac * self.fibo_obj.data['B_y_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Ye_zz_'+t] -= myfac * self.fibo_obj.data['B_z_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Ye_xy_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Ye_xz_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Ye_yz_'+t] -= myfac * self.fibo_obj.data['B_y_'+t] * self.fibo_obj.data['B_z_'+t]

      self.fibo_obj.data['Y_xx_'+t] = self.fibo_obj.data['dB_c_x_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Y_xx_'+t] *= iB2
      self.fibo_obj.data['Y_yy_'+t] = self.fibo_obj.data['dB_c_y_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Y_yy_'+t] *= iB2
      self.fibo_obj.data['Y_zz_'+t] = self.fibo_obj.data['dB_c_z_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Y_zz_'+t] *= iB2
      self.fibo_obj.data['Y_xy_'+t]  = self.fibo_obj.data['dB_c_x_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Y_xy_'+t] += self.fibo_obj.data['dB_c_y_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Y_xy_'+t] *= 0.5 * iB2
      self.fibo_obj.data['Y_xz_'+t]  = self.fibo_obj.data['dB_c_x_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Y_xz_'+t] += self.fibo_obj.data['dB_c_z_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Y_xz_'+t] *= 0.5 * iB2
      self.fibo_obj.data['Y_yz_'+t]  = self.fibo_obj.data['dB_c_y_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Y_yz_'+t] += self.fibo_obj.data['dB_c_z_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Y_yz_'+t] *= 0.5 * iB2

      myfac = iB2 * iB2 * self.fibo_obj.calc_scalr('B_x_'+t,'dB_c_x_'+t,'B_y_'+t,'dB_c_y_'+t,'B_z_'+t,'dB_c_z_'+t)
      self.fibo_obj.data['Y_xx_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_x_'+t]
      self.fibo_obj.data['Y_yy_'+t] -= myfac * self.fibo_obj.data['B_y_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Y_zz_'+t] -= myfac * self.fibo_obj.data['B_z_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Y_xy_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_y_'+t]
      self.fibo_obj.data['Y_xz_'+t] -= myfac * self.fibo_obj.data['B_x_'+t] * self.fibo_obj.data['B_z_'+t]
      self.fibo_obj.data['Y_yz_'+t] -= myfac * self.fibo_obj.data['B_y_'+t] * self.fibo_obj.data['B_z_'+t]


    #-------------------------------------------
    #energy conversion between local electromagnetic and local kinetic energy (qnu dot E) or (u dot (J cross B))
    self.fibo_obj.data['dKi_F_'+t]  = self.fibo_obj.calc_scalr('ui_x_'+t,'E_x_'+t,'ui_y_'+t,'E_y_'+t,'ui_z_'+t,'E_z_'+t)
    self.fibo_obj.data['dKi_F_'+t]  = np.multiply(self.fibo_obj.data['dKi_F_'+t],self.fibo_obj.data['n_'+t])

    self.fibo_obj.data['dKe_F_'+t]  = -self.fibo_obj.calc_scalr('ue_x_'+t,'E_x_'+t,'ue_y_'+t,'E_y_'+t,'ue_z_'+t,'E_z_'+t)
    self.fibo_obj.data['dKe_F_'+t]  = np.multiply(self.fibo_obj.data['dKe_F_'+t],self.fibo_obj.data['n_'+t])

    self.fibo_obj.data['dK_F_'+t] = -self.fibo_obj.calc_scalr('J_x_'+t,'uXB_x_'+t,'J_y_'+t,'uXB_y_'+t,'J_z_'+t,'uXB_z_'+t)

    if par_per:
      self.fibo_obj.data['dK_F_par_'+t]  = iB * self.fibo_obj.calc_scalr('E_x_'+t,'B_x_'+t,'E_y_'+t,'B_y_'+t,'E_z_'+t,'B_z_'+t)
      self.fibo_obj.data['dKi_F_par_'+t] =  ui_par * self.fibo_obj.data['n_'+t] * self.fibo_obj.data['dK_F_par_'+t]
      self.fibo_obj.data['dKe_F_par_'+t] = -ue_par * self.fibo_obj.data['n_'+t] * self.fibo_obj.data['dK_F_par_'+t]
      self.fibo_obj.data['dK_F_par_'+t] *= 0.

    #energy conversion between surrounding internal and local kinetic -(u dot div P)
    divP_x = self.fibo_obj.calc_divr('Pi_xx_'+t,'Pi_xy_'+t,'Pi_xz_'+t)
    divP_y = self.fibo_obj.calc_divr('Pi_xy_'+t,'Pi_yy_'+t,'Pi_yz_'+t)
    divP_z = self.fibo_obj.calc_divr('Pi_xz_'+t,'Pi_yz_'+t,'Pi_zz_'+t)
    self.fibo_obj.data['dKi_U_'+t]  = -np.multiply(self.fibo_obj.data['ui_x_'+t],divP_x)
    self.fibo_obj.data['dKi_U_'+t] += -np.multiply(self.fibo_obj.data['ui_y_'+t],divP_y)
    self.fibo_obj.data['dKi_U_'+t] += -np.multiply(self.fibo_obj.data['ui_z_'+t],divP_z)
    if par_per:
      self.fibo_obj.data['dKi_U_par_'+t]  = self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],divP_x,self.fibo_obj.data['B_y_'+t],divP_y,self.fibo_obj.data['B_z_'+t],divP_z)
      self.fibo_obj.data['dKi_U_par_'+t] *= - iB * ui_par

    divP_x = self.fibo_obj.calc_divr('Pe_xx_'+t,'Pe_xy_'+t,'Pe_xz_'+t)
    divP_y = self.fibo_obj.calc_divr('Pe_xy_'+t,'Pe_yy_'+t,'Pe_yz_'+t)
    divP_z = self.fibo_obj.calc_divr('Pe_xz_'+t,'Pe_yz_'+t,'Pe_zz_'+t)
    self.fibo_obj.data['dKe_U_'+t]  = -np.multiply(self.fibo_obj.data['ue_x_'+t],divP_x)
    self.fibo_obj.data['dKe_U_'+t] += -np.multiply(self.fibo_obj.data['ue_y_'+t],divP_y)
    self.fibo_obj.data['dKe_U_'+t] += -np.multiply(self.fibo_obj.data['ue_z_'+t],divP_z)
    if par_per:
      self.fibo_obj.data['dKe_U_par_'+t]  = self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],divP_x,self.fibo_obj.data['B_y_'+t],divP_y,self.fibo_obj.data['B_z_'+t],divP_z)
      self.fibo_obj.data['dKe_U_par_'+t] *= - iB * ue_par

    divP_x = self.fibo_obj.calc_divr('P_xx_'+t,'P_xy_'+t,'P_xz_'+t)
    divP_y = self.fibo_obj.calc_divr('P_xy_'+t,'P_yy_'+t,'P_yz_'+t)
    divP_z = self.fibo_obj.calc_divr('P_xz_'+t,'P_yz_'+t,'P_zz_'+t)
    self.fibo_obj.data['dK_U_'+t]  = -np.multiply(self.fibo_obj.data['u_x_'+t],divP_x)
    self.fibo_obj.data['dK_U_'+t] += -np.multiply(self.fibo_obj.data['u_y_'+t],divP_y)
    self.fibo_obj.data['dK_U_'+t] += -np.multiply(self.fibo_obj.data['u_z_'+t],divP_z)
    if par_per:
      self.fibo_obj.data['dK_U_par_'+t]  = self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],divP_x,self.fibo_obj.data['B_y_'+t],divP_y,self.fibo_obj.data['B_z_'+t],divP_z)
      self.fibo_obj.data['dK_U_par_'+t] *= - iB * u_par

    #local kinetic energy accumulation/dispersion due to convection -K(div u)
    self.fibo_obj.data['dKi_c_'+t] = -np.multiply(self.fibo_obj.data['Ki_'+t],self.fibo_obj.data['div_ui_'+t])
    self.fibo_obj.data['dKe_c_'+t] = -np.multiply(self.fibo_obj.data['Ke_'+t],self.fibo_obj.data['div_ue_'+t])
    self.fibo_obj.data['dK_c_'+t] = -np.multiply(self.fibo_obj.data['K_'+t],self.fibo_obj.data['div_u_'+t])
    if par_per:
      self.fibo_obj.data['dKi_c_par_'+t] = -np.multiply(self.fibo_obj.data['Ki_par_'+t],self.fibo_obj.data['div_ui_'+t])
      self.fibo_obj.data['dKe_c_par_'+t] = -np.multiply(self.fibo_obj.data['Ke_par_'+t],self.fibo_obj.data['div_ue_'+t])
      self.fibo_obj.data['dK_c_par_'+t] = -np.multiply(self.fibo_obj.data['K_par_'+t],self.fibo_obj.data['div_u_'+t])

    #local kinetic energy accumulation/dispersion due to convection -K(div u)
    self.fibo_obj.data['dKi_a_'+t]  = np.multiply(self.fibo_obj.data['ui_x_'+t],self.fibo_obj.calc_gradx(self.fibo_obj.data['Ki_'+t]))
    self.fibo_obj.data['dKi_a_'+t] += np.multiply(self.fibo_obj.data['ui_y_'+t],self.fibo_obj.calc_grady(self.fibo_obj.data['Ki_'+t]))
    self.fibo_obj.data['dKi_a_'+t] += np.multiply(self.fibo_obj.data['ui_z_'+t],self.fibo_obj.calc_gradz(self.fibo_obj.data['Ki_'+t]))
    
    self.fibo_obj.data['dKe_a_'+t]  = np.multiply(self.fibo_obj.data['ue_x_'+t],self.fibo_obj.calc_gradx(self.fibo_obj.data['Ke_'+t]))
    self.fibo_obj.data['dKe_a_'+t] += np.multiply(self.fibo_obj.data['ue_y_'+t],self.fibo_obj.calc_grady(self.fibo_obj.data['Ke_'+t]))
    self.fibo_obj.data['dKe_a_'+t] += np.multiply(self.fibo_obj.data['ue_z_'+t],self.fibo_obj.calc_gradz(self.fibo_obj.data['Ke_'+t]))
    
    self.fibo_obj.data['dK_a_'+t]  = np.multiply(self.fibo_obj.data['u_x_'+t],self.fibo_obj.calc_gradx(self.fibo_obj.data['K_'+t]))
    self.fibo_obj.data['dK_a_'+t] += np.multiply(self.fibo_obj.data['u_y_'+t],self.fibo_obj.calc_grady(self.fibo_obj.data['K_'+t]))
    self.fibo_obj.data['dK_a_'+t] += np.multiply(self.fibo_obj.data['u_z_'+t],self.fibo_obj.calc_gradz(self.fibo_obj.data['K_'+t]))

    if par_per:
      self.fibo_obj.data['dKi_a_par_'+t]  = np.multiply(self.fibo_obj.data['ui_x_'+t],self.fibo_obj.calc_gradx(self.fibo_obj.data['Ki_par_'+t]))
      self.fibo_obj.data['dKi_a_par_'+t] += np.multiply(self.fibo_obj.data['ui_y_'+t],self.fibo_obj.calc_grady(self.fibo_obj.data['Ki_par_'+t]))
      self.fibo_obj.data['dKi_a_par_'+t] += np.multiply(self.fibo_obj.data['ui_z_'+t],self.fibo_obj.calc_gradz(self.fibo_obj.data['Ki_par_'+t]))
      
      self.fibo_obj.data['dKe_a_par_'+t]  = np.multiply(self.fibo_obj.data['ue_x_'+t],self.fibo_obj.calc_gradx(self.fibo_obj.data['Ke_par_'+t]))
      self.fibo_obj.data['dKe_a_par_'+t] += np.multiply(self.fibo_obj.data['ue_y_'+t],self.fibo_obj.calc_grady(self.fibo_obj.data['Ke_par_'+t]))
      self.fibo_obj.data['dKe_a_par_'+t] += np.multiply(self.fibo_obj.data['ue_z_'+t],self.fibo_obj.calc_gradz(self.fibo_obj.data['Ke_par_'+t]))
      
      self.fibo_obj.data['dK_a_par_'+t]  = np.multiply(self.fibo_obj.data['u_x_'+t],self.fibo_obj.calc_gradx(self.fibo_obj.data['K_par_'+t]))
      self.fibo_obj.data['dK_a_par_'+t] += np.multiply(self.fibo_obj.data['u_y_'+t],self.fibo_obj.calc_grady(self.fibo_obj.data['K_par_'+t]))
      self.fibo_obj.data['dK_a_par_'+t] += np.multiply(self.fibo_obj.data['u_z_'+t],self.fibo_obj.calc_gradz(self.fibo_obj.data['K_par_'+t]))


    #parallel energy gain due to change in b direction along fluid motions
    if par_per: 
      self.fibo_obj.data['dKi_Y_par_'+t]  = self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.calc_scalr('ui_x_'+t,'Yi_xx_'+t,'ui_y_'+t,'Yi_xy_'+t,'ui_z_'+t,'Yi_xz_'+t)
      self.fibo_obj.data['dKi_Y_par_'+t] += self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.calc_scalr('ui_x_'+t,'Yi_xy_'+t,'ui_y_'+t,'Yi_yy_'+t,'ui_z_'+t,'Yi_yz_'+t)
      self.fibo_obj.data['dKi_Y_par_'+t] += self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.calc_scalr('ui_x_'+t,'Yi_xz_'+t,'ui_y_'+t,'Yi_yz_'+t,'ui_z_'+t,'Yi_zz_'+t)
      self.fibo_obj.data['dKi_Y_par_'+t] *= self.fibo_obj.data['n_'+t]
      
      self.fibo_obj.data['dKe_Y_par_'+t]  = self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.calc_scalr('ue_x_'+t,'Ye_xx_'+t,'ue_y_'+t,'Ye_xy_'+t,'ue_z_'+t,'Ye_xz_'+t)
      self.fibo_obj.data['dKe_Y_par_'+t] += self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.calc_scalr('ue_x_'+t,'Ye_xy_'+t,'ue_y_'+t,'Ye_yy_'+t,'ue_z_'+t,'Ye_yz_'+t)
      self.fibo_obj.data['dKe_Y_par_'+t] += self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.calc_scalr('ue_x_'+t,'Ye_xz_'+t,'ue_y_'+t,'Ye_yz_'+t,'ue_z_'+t,'Ye_zz_'+t)
      self.fibo_obj.data['dKe_Y_par_'+t] *= self.fibo_obj.data['n_'+t] / (self.fibo_obj.meta['mime'])
      
      self.fibo_obj.data['dK_Y_par_'+t]  = self.fibo_obj.data['u_x_'+t] * self.fibo_obj.calc_scalr('u_x_'+t,'Y_xx_'+t,'u_y_'+t,'Y_xy_'+t,'u_z_'+t,'Y_xz_'+t)
      self.fibo_obj.data['dK_Y_par_'+t] += self.fibo_obj.data['u_y_'+t] * self.fibo_obj.calc_scalr('u_x_'+t,'Y_xy_'+t,'u_y_'+t,'Y_yy_'+t,'u_z_'+t,'Y_yz_'+t)
      self.fibo_obj.data['dK_Y_par_'+t] += self.fibo_obj.data['u_z_'+t] * self.fibo_obj.calc_scalr('u_x_'+t,'Y_xz_'+t,'u_y_'+t,'Y_yz_'+t,'u_z_'+t,'Y_zz_'+t)
      self.fibo_obj.data['dK_Y_par_'+t] *= self.fibo_obj.data['n_'+t] * (self.fibo_obj.meta['mime']+1) / (self.fibo_obj.meta['mime'])

    #total energy conversion to local kinetic energy (electromagnetic + surrounding internal + convection term)
    #self.fibo_obj.data['dKi_'+t] = self.fibo_obj.data['dKi_F_'+t] + self.fibo_obj.data['dKi_U_'+t] #+ self.fibo_obj.data['dKi_c_'+t]
    #self.fibo_obj.data['dKe_'+t] = self.fibo_obj.data['dKe_F_'+t] #+ self.fibo_obj.data['dKe_c_'+t]
    #self.fibo_obj.data['dK_'+t] = self.fibo_obj.data['dK_F_'+t] + self.fibo_obj.data['dK_U_'+t] #+ self.fibo_obj.data['dK_c_'+t]

    #if par_per:
    #  self.fibo_obj.data['dKi_F_per_'+t] = self.fibo_obj.data['dKi_F_'+t] - self.fibo_obj.data['dKi_F_par_'+t]
    #  self.fibo_obj.data['dKi_U_per_'+t] = self.fibo_obj.data['dKi_U_'+t] - self.fibo_obj.data['dKi_U_par_'+t]
    #  self.fibo_obj.data['dKi_c_per_'+t] = self.fibo_obj.data['dKi_c_'+t] - self.fibo_obj.data['dKi_c_par_'+t]
    #  self.fibo_obj.data['dKi_Y_per_'+t] = - self.fibo_obj.data['dKi_Y_par_'+t]
    #  self.fibo_obj.data['dKe_F_per_'+t] = self.fibo_obj.data['dKe_F_'+t] - self.fibo_obj.data['dKe_F_par_'+t]
    #  self.fibo_obj.data['dKe_U_per_'+t] = self.fibo_obj.data['dKe_U_'+t] - self.fibo_obj.data['dKe_U_par_'+t]
    #  self.fibo_obj.data['dKe_c_per_'+t] = self.fibo_obj.data['dKe_c_'+t] - self.fibo_obj.data['dKe_c_par_'+t]
    #  self.fibo_obj.data['dKe_Y_per_'+t] = - self.fibo_obj.data['dKe_Y_par_'+t]
    #  self.fibo_obj.data['dK_F_per_'+t] = self.fibo_obj.data['dK_F_'+t] - self.fibo_obj.data['dK_F_par_'+t]
    #  self.fibo_obj.data['dK_U_per_'+t] = self.fibo_obj.data['dK_U_'+t] - self.fibo_obj.data['dK_U_par_'+t]
    #  self.fibo_obj.data['dK_c_per_'+t] = self.fibo_obj.data['dK_c_'+t] - self.fibo_obj.data['dK_c_par_'+t]
    #  self.fibo_obj.data['dK_Y_per_'+t] = - self.fibo_obj.data['dK_Y_par_'+t]

    #-------------------------------------------
    #energy conversion between fields and internal (J dot (E + u cross B)) ... only for barycentre! 
    self.fibo_obj.data['dU_F_'+t]  = self.fibo_obj.calc_scalr('J_x_'+t,'E_x_'+t,'J_y_'+t,'E_y_'+t,'J_z_'+t,'E_z_'+t)
    self.fibo_obj.data['dU_F_'+t] += self.fibo_obj.calc_scalr('J_x_'+t,'uXB_x_'+t,'J_y_'+t,'uXB_y_'+t,'J_z_'+t,'uXB_z_'+t)
    if par_per:
      self.fibo_obj.data['dU_F_par_'+t] = self.fibo_obj.calc_scalr('J_x_'+t,'E_par_x_'+t,'J_y_'+t,'E_par_y_'+t,'J_z_'+t,'E_par_z_'+t)

    #energy conversion between surrounding kinetic and local internal -(P ddot grad u)
    self.fibo_obj.data['dUi_K_'+t]  = -np.multiply(self.fibo_obj.data['Pi_xx_'+t],self.fibo_obj.calc_gradx('ui_x_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_xy_'+t],self.fibo_obj.calc_gradx('ui_y_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_xz_'+t],self.fibo_obj.calc_gradx('ui_z_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_xy_'+t],self.fibo_obj.calc_grady('ui_x_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_yy_'+t],self.fibo_obj.calc_grady('ui_y_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_yz_'+t],self.fibo_obj.calc_grady('ui_z_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_xz_'+t],self.fibo_obj.calc_gradz('ui_x_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_yz_'+t],self.fibo_obj.calc_gradz('ui_y_'+t))
    self.fibo_obj.data['dUi_K_'+t] += -np.multiply(self.fibo_obj.data['Pi_zz_'+t],self.fibo_obj.calc_gradz('ui_z_'+t))

    if par_per:
      tens_xx  = self.fibo_obj.data['Pi_xx_'+t] * self.fibo_obj.calc_gradx('ui_x_'+t)
      tens_xx += self.fibo_obj.data['Pi_xy_'+t] * self.fibo_obj.calc_grady('ui_x_'+t)
      tens_xx += self.fibo_obj.data['Pi_xz_'+t] * self.fibo_obj.calc_gradz('ui_x_'+t)
      tens_yy  = self.fibo_obj.data['Pi_xy_'+t] * self.fibo_obj.calc_gradx('ui_y_'+t)
      tens_yy += self.fibo_obj.data['Pi_yy_'+t] * self.fibo_obj.calc_grady('ui_y_'+t)
      tens_yy += self.fibo_obj.data['Pi_yz_'+t] * self.fibo_obj.calc_gradz('ui_y_'+t)
      tens_zz  = self.fibo_obj.data['Pi_xz_'+t] * self.fibo_obj.calc_gradx('ui_z_'+t)
      tens_zz += self.fibo_obj.data['Pi_yz_'+t] * self.fibo_obj.calc_grady('ui_z_'+t)
      tens_zz += self.fibo_obj.data['Pi_zz_'+t] * self.fibo_obj.calc_gradz('ui_z_'+t)
      tens_xy  = self.fibo_obj.data['Pi_xx_'+t] * self.fibo_obj.calc_gradx('ui_y_'+t)
      tens_xy += self.fibo_obj.data['Pi_xy_'+t] * self.fibo_obj.calc_grady('ui_y_'+t)
      tens_xy += self.fibo_obj.data['Pi_xz_'+t] * self.fibo_obj.calc_gradz('ui_y_'+t)
      tens_xy += self.fibo_obj.data['Pi_xy_'+t] * self.fibo_obj.calc_gradx('ui_x_'+t)
      tens_xy += self.fibo_obj.data['Pi_yy_'+t] * self.fibo_obj.calc_grady('ui_x_'+t)
      tens_xy += self.fibo_obj.data['Pi_yz_'+t] * self.fibo_obj.calc_gradz('ui_x_'+t)
      tens_xy *= 0.5
      tens_xz  = self.fibo_obj.data['Pi_xx_'+t] * self.fibo_obj.calc_gradx('ui_z_'+t)
      tens_xz += self.fibo_obj.data['Pi_xy_'+t] * self.fibo_obj.calc_grady('ui_z_'+t)
      tens_xz += self.fibo_obj.data['Pi_xz_'+t] * self.fibo_obj.calc_gradz('ui_z_'+t)
      tens_xz += self.fibo_obj.data['Pi_xz_'+t] * self.fibo_obj.calc_gradx('ui_x_'+t)
      tens_xz += self.fibo_obj.data['Pi_yz_'+t] * self.fibo_obj.calc_grady('ui_x_'+t)
      tens_xz += self.fibo_obj.data['Pi_zz_'+t] * self.fibo_obj.calc_gradz('ui_x_'+t)
      tens_xz *= 0.5
      tens_yz  = self.fibo_obj.data['Pi_xy_'+t] * self.fibo_obj.calc_gradx('ui_z_'+t)
      tens_yz += self.fibo_obj.data['Pi_yy_'+t] * self.fibo_obj.calc_grady('ui_z_'+t)
      tens_yz += self.fibo_obj.data['Pi_yz_'+t] * self.fibo_obj.calc_gradz('ui_z_'+t)
      tens_yz += self.fibo_obj.data['Pi_xz_'+t] * self.fibo_obj.calc_gradx('ui_y_'+t)
      tens_yz += self.fibo_obj.data['Pi_yz_'+t] * self.fibo_obj.calc_grady('ui_y_'+t)
      tens_yz += self.fibo_obj.data['Pi_zz_'+t] * self.fibo_obj.calc_gradz('ui_y_'+t)
      tens_yz *= 0.5
      self.fibo_obj.data['dUi_K_par_'+t]  = self.fibo_obj.data['B_x_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xx,self.fibo_obj.data['B_y_'+t],tens_xy,self.fibo_obj.data['B_z_'+t],tens_xz)
      self.fibo_obj.data['dUi_K_par_'+t] += self.fibo_obj.data['B_y_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xy,self.fibo_obj.data['B_y_'+t],tens_yy,self.fibo_obj.data['B_z_'+t],tens_yz)
      self.fibo_obj.data['dUi_K_par_'+t] += self.fibo_obj.data['B_z_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xz,self.fibo_obj.data['B_y_'+t],tens_yz,self.fibo_obj.data['B_z_'+t],tens_zz)
      self.fibo_obj.data['dUi_K_par_'+t] *= -iB2

    self.fibo_obj.data['dUe_K_'+t]  = -np.multiply(self.fibo_obj.data['Pe_xx_'+t],self.fibo_obj.calc_gradx('ue_x_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_xy_'+t],self.fibo_obj.calc_gradx('ue_y_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_xz_'+t],self.fibo_obj.calc_gradx('ue_z_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_xy_'+t],self.fibo_obj.calc_grady('ue_x_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_yy_'+t],self.fibo_obj.calc_grady('ue_y_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_yz_'+t],self.fibo_obj.calc_grady('ue_z_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_xz_'+t],self.fibo_obj.calc_gradz('ue_x_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_yz_'+t],self.fibo_obj.calc_gradz('ue_y_'+t))
    self.fibo_obj.data['dUe_K_'+t] += -np.multiply(self.fibo_obj.data['Pe_zz_'+t],self.fibo_obj.calc_gradz('ue_z_'+t))

    if par_per:
      tens_xx  = self.fibo_obj.data['Pe_xx_'+t] * self.fibo_obj.calc_gradx('ue_x_'+t)
      tens_xx += self.fibo_obj.data['Pe_xy_'+t] * self.fibo_obj.calc_grady('ue_x_'+t)
      tens_xx += self.fibo_obj.data['Pe_xz_'+t] * self.fibo_obj.calc_gradz('ue_x_'+t)
      tens_yy  = self.fibo_obj.data['Pe_xy_'+t] * self.fibo_obj.calc_gradx('ue_y_'+t)
      tens_yy += self.fibo_obj.data['Pe_yy_'+t] * self.fibo_obj.calc_grady('ue_y_'+t)
      tens_yy += self.fibo_obj.data['Pe_yz_'+t] * self.fibo_obj.calc_gradz('ue_y_'+t)
      tens_zz  = self.fibo_obj.data['Pe_xz_'+t] * self.fibo_obj.calc_gradx('ue_z_'+t)
      tens_zz += self.fibo_obj.data['Pe_yz_'+t] * self.fibo_obj.calc_grady('ue_z_'+t)
      tens_zz += self.fibo_obj.data['Pe_zz_'+t] * self.fibo_obj.calc_gradz('ue_z_'+t)
      tens_xy  = self.fibo_obj.data['Pe_xx_'+t] * self.fibo_obj.calc_gradx('ue_y_'+t)
      tens_xy += self.fibo_obj.data['Pe_xy_'+t] * self.fibo_obj.calc_grady('ue_y_'+t)
      tens_xy += self.fibo_obj.data['Pe_xz_'+t] * self.fibo_obj.calc_gradz('ue_y_'+t)
      tens_xy += self.fibo_obj.data['Pe_xy_'+t] * self.fibo_obj.calc_gradx('ue_x_'+t)
      tens_xy += self.fibo_obj.data['Pe_yy_'+t] * self.fibo_obj.calc_grady('ue_x_'+t)
      tens_xy += self.fibo_obj.data['Pe_yz_'+t] * self.fibo_obj.calc_gradz('ue_x_'+t)
      tens_xy *= 0.5
      tens_xz  = self.fibo_obj.data['Pe_xx_'+t] * self.fibo_obj.calc_gradx('ue_z_'+t)
      tens_xz += self.fibo_obj.data['Pe_xy_'+t] * self.fibo_obj.calc_grady('ue_z_'+t)
      tens_xz += self.fibo_obj.data['Pe_xz_'+t] * self.fibo_obj.calc_gradz('ue_z_'+t)
      tens_xz += self.fibo_obj.data['Pe_xz_'+t] * self.fibo_obj.calc_gradx('ue_x_'+t)
      tens_xz += self.fibo_obj.data['Pe_yz_'+t] * self.fibo_obj.calc_grady('ue_x_'+t)
      tens_xz += self.fibo_obj.data['Pe_zz_'+t] * self.fibo_obj.calc_gradz('ue_x_'+t)
      tens_xz *= 0.5
      tens_yz  = self.fibo_obj.data['Pe_xy_'+t] * self.fibo_obj.calc_gradx('ue_z_'+t)
      tens_yz += self.fibo_obj.data['Pe_yy_'+t] * self.fibo_obj.calc_grady('ue_z_'+t)
      tens_yz += self.fibo_obj.data['Pe_yz_'+t] * self.fibo_obj.calc_gradz('ue_z_'+t)
      tens_yz += self.fibo_obj.data['Pe_xz_'+t] * self.fibo_obj.calc_gradx('ue_y_'+t)
      tens_yz += self.fibo_obj.data['Pe_yz_'+t] * self.fibo_obj.calc_grady('ue_y_'+t)
      tens_yz += self.fibo_obj.data['Pe_zz_'+t] * self.fibo_obj.calc_gradz('ue_y_'+t)
      tens_yz *= 0.5
      self.fibo_obj.data['dUe_K_par_'+t]  = self.fibo_obj.data['B_x_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xx,self.fibo_obj.data['B_y_'+t],tens_xy,self.fibo_obj.data['B_z_'+t],tens_xz)
      self.fibo_obj.data['dUe_K_par_'+t] += self.fibo_obj.data['B_y_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xy,self.fibo_obj.data['B_y_'+t],tens_yy,self.fibo_obj.data['B_z_'+t],tens_yz)
      self.fibo_obj.data['dUe_K_par_'+t] += self.fibo_obj.data['B_z_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xz,self.fibo_obj.data['B_y_'+t],tens_yz,self.fibo_obj.data['B_z_'+t],tens_zz)
      self.fibo_obj.data['dUe_K_par_'+t] *= -iB2

    self.fibo_obj.data['dU_K_'+t]  = -np.multiply(self.fibo_obj.data['P_xx_'+t],self.fibo_obj.calc_gradx('u_x_'+t))    
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_xy_'+t],self.fibo_obj.calc_gradx('u_y_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_xz_'+t],self.fibo_obj.calc_gradx('u_z_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_xy_'+t],self.fibo_obj.calc_grady('u_x_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_yy_'+t],self.fibo_obj.calc_grady('u_y_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_yz_'+t],self.fibo_obj.calc_grady('u_z_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_xz_'+t],self.fibo_obj.calc_gradz('u_x_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_yz_'+t],self.fibo_obj.calc_gradz('u_y_'+t))
    self.fibo_obj.data['dU_K_'+t] += -np.multiply(self.fibo_obj.data['P_zz_'+t],self.fibo_obj.calc_gradz('u_z_'+t))

    if par_per:
      tens_xx  = self.fibo_obj.data['P_xx_'+t] * self.fibo_obj.calc_gradx('u_x_'+t)
      tens_xx += self.fibo_obj.data['P_xy_'+t] * self.fibo_obj.calc_grady('u_x_'+t)
      tens_xx += self.fibo_obj.data['P_xz_'+t] * self.fibo_obj.calc_gradz('u_x_'+t)
      tens_yy  = self.fibo_obj.data['P_xy_'+t] * self.fibo_obj.calc_gradx('u_y_'+t)
      tens_yy += self.fibo_obj.data['P_yy_'+t] * self.fibo_obj.calc_grady('u_y_'+t)
      tens_yy += self.fibo_obj.data['P_yz_'+t] * self.fibo_obj.calc_gradz('u_y_'+t)
      tens_zz  = self.fibo_obj.data['P_xz_'+t] * self.fibo_obj.calc_gradx('u_z_'+t)
      tens_zz += self.fibo_obj.data['P_yz_'+t] * self.fibo_obj.calc_grady('u_z_'+t)
      tens_zz += self.fibo_obj.data['P_zz_'+t] * self.fibo_obj.calc_gradz('u_z_'+t)
      tens_xy  = self.fibo_obj.data['P_xx_'+t] * self.fibo_obj.calc_gradx('u_y_'+t)
      tens_xy += self.fibo_obj.data['P_xy_'+t] * self.fibo_obj.calc_grady('u_y_'+t)
      tens_xy += self.fibo_obj.data['P_xz_'+t] * self.fibo_obj.calc_gradz('u_y_'+t)
      tens_xy += self.fibo_obj.data['P_xy_'+t] * self.fibo_obj.calc_gradx('u_x_'+t)
      tens_xy += self.fibo_obj.data['P_yy_'+t] * self.fibo_obj.calc_grady('u_x_'+t)
      tens_xy += self.fibo_obj.data['P_yz_'+t] * self.fibo_obj.calc_gradz('u_x_'+t)
      tens_xy *= 0.5
      tens_xz  = self.fibo_obj.data['P_xx_'+t] * self.fibo_obj.calc_gradx('u_z_'+t)
      tens_xz += self.fibo_obj.data['P_xy_'+t] * self.fibo_obj.calc_grady('u_z_'+t)
      tens_xz += self.fibo_obj.data['P_xz_'+t] * self.fibo_obj.calc_gradz('u_z_'+t)
      tens_xz += self.fibo_obj.data['P_xz_'+t] * self.fibo_obj.calc_gradx('u_x_'+t)
      tens_xz += self.fibo_obj.data['P_yz_'+t] * self.fibo_obj.calc_grady('u_x_'+t)
      tens_xz += self.fibo_obj.data['P_zz_'+t] * self.fibo_obj.calc_gradz('u_x_'+t)
      tens_xz *= 0.5
      tens_yz  = self.fibo_obj.data['P_xy_'+t] * self.fibo_obj.calc_gradx('u_z_'+t)
      tens_yz += self.fibo_obj.data['P_yy_'+t] * self.fibo_obj.calc_grady('u_z_'+t)
      tens_yz += self.fibo_obj.data['P_yz_'+t] * self.fibo_obj.calc_gradz('u_z_'+t)
      tens_yz += self.fibo_obj.data['P_xz_'+t] * self.fibo_obj.calc_gradx('u_y_'+t)
      tens_yz += self.fibo_obj.data['P_yz_'+t] * self.fibo_obj.calc_grady('u_y_'+t)
      tens_yz += self.fibo_obj.data['P_zz_'+t] * self.fibo_obj.calc_gradz('u_y_'+t)
      tens_yz *= 0.5
      self.fibo_obj.data['dU_K_par_'+t]  = self.fibo_obj.data['B_x_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xx,self.fibo_obj.data['B_y_'+t],tens_xy,self.fibo_obj.data['B_z_'+t],tens_xz)
      self.fibo_obj.data['dU_K_par_'+t] += self.fibo_obj.data['B_y_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xy,self.fibo_obj.data['B_y_'+t],tens_yy,self.fibo_obj.data['B_z_'+t],tens_yz)
      self.fibo_obj.data['dU_K_par_'+t] += self.fibo_obj.data['B_z_'+t] * self.fibo_obj.calc_scalr(self.fibo_obj.data['B_x_'+t],tens_xz,self.fibo_obj.data['B_y_'+t],tens_yz,self.fibo_obj.data['B_z_'+t],tens_zz)
      self.fibo_obj.data['dU_K_par_'+t] *= -iB2

    #energy conversion between surrounding internal and local internal -(div Q)/2
    self.fibo_obj.data['dUi_U_'+t]  =  -0.5 *self.fibo_obj.calc_divr('Qi_x_'+t,'Qi_y_'+t,'Qi_z_'+t) 
    self.fibo_obj.data['dUe_U_'+t]  =  -0.5 *self.fibo_obj.calc_divr('Qe_x_'+t,'Qe_y_'+t,'Qe_z_'+t) 
    self.fibo_obj.data['dU_U_'+t]  =   -0.5 *self.fibo_obj.calc_divr('Q_x_'+t,'Q_y_'+t,'Q_z_'+t)

    if iso_Te: 
      self.fibo_obj.data['dUe_U_'+t] -= self.fibo_obj.data['dUe_K_'+t]
      self.fibo_obj.data['dU_U_'+t]  += self.fibo_obj.data['dUe_U_'+t]

    if par_per: #achtung: this is approximate - yet the best we can have when the 10 elements of rank-3 Q are not provided ... 
      self.fibo_obj.data['dUi_U_par_'+t]  = self.fibo_obj.calc_gradx('Qi_par_x_'+t)
      self.fibo_obj.data['dUi_U_par_'+t] += self.fibo_obj.calc_grady('Qi_par_y_'+t)
      self.fibo_obj.data['dUi_U_par_'+t] += self.fibo_obj.calc_gradz('Qi_par_z_'+t)
      self.fibo_obj.data['dUi_U_par_'+t] *= -0.5 * iB2

      self.fibo_obj.data['dUe_U_par_'+t]  = self.fibo_obj.calc_gradx('Qe_par_x_'+t)
      self.fibo_obj.data['dUe_U_par_'+t] += self.fibo_obj.calc_grady('Qe_par_y_'+t)
      self.fibo_obj.data['dUe_U_par_'+t] += self.fibo_obj.calc_gradz('Qe_par_z_'+t)
      self.fibo_obj.data['dUe_U_par_'+t] *= -0.5 * iB2

      self.fibo_obj.data['dU_U_par_'+t]  = self.fibo_obj.calc_gradx('Q_par_x_'+t)
      self.fibo_obj.data['dU_U_par_'+t] += self.fibo_obj.calc_grady('Q_par_y_'+t)
      self.fibo_obj.data['dU_U_par_'+t] += self.fibo_obj.calc_gradz('Q_par_z_'+t)
      self.fibo_obj.data['dU_U_par_'+t] *= -0.5 * iB2

      if iso_Te: 
        self.fibo_obj.data['dUe_U_par_'+t] -= self.fibo_obj.data['dUe_K_par_'+t]
        self.fibo_obj.data['dU_U_par_'+t]  += self.fibo_obj.data['dUe_U_par_'+t]

    #local internal energy accumulation/dispersion due to convection -U(div u)
    self.fibo_obj.data['dUi_c_'+t] = -np.multiply(self.fibo_obj.data['Ui_'+t],self.fibo_obj.data['div_ui_'+t])
    self.fibo_obj.data['dUe_c_'+t] = -np.multiply(self.fibo_obj.data['Ue_'+t],self.fibo_obj.data['div_ue_'+t])
    self.fibo_obj.data['dU_c_'+t]  = -np.multiply(self.fibo_obj.data['U_'+t],self.fibo_obj.data['div_u_'+t])

    if par_per:
      self.fibo_obj.data['dUi_c_par_'+t] = -np.multiply(self.fibo_obj.data['Ui_par_'+t],self.fibo_obj.data['div_ui_'+t])
      self.fibo_obj.data['dUe_c_par_'+t] = -np.multiply(self.fibo_obj.data['Ue_par_'+t],self.fibo_obj.data['div_ue_'+t])
      self.fibo_obj.data['dU_c_par_'+t]  = -np.multiply(self.fibo_obj.data['U_par_'+t],self.fibo_obj.data['div_u_'+t])

    #parallel energy gain due to change in b direction along fluid motions
    if par_per:
      self.fibo_obj.data['dUi_Y_par_'+t]  = self.fibo_obj.calc_scalr('Pi_xx_'+t,'Yi_xx_'+t,'Pi_xy_'+t,'Yi_xy_'+t,'Pi_xz_'+t,'Yi_xz_'+t)
      self.fibo_obj.data['dUi_Y_par_'+t] += self.fibo_obj.calc_scalr('Pi_xy_'+t,'Yi_xy_'+t,'Pi_yy_'+t,'Yi_yy_'+t,'Pi_yz_'+t,'Yi_yz_'+t)
      self.fibo_obj.data['dUi_Y_par_'+t] += self.fibo_obj.calc_scalr('Pi_xz_'+t,'Yi_xz_'+t,'Pi_yz_'+t,'Yi_yz_'+t,'Pi_zz_'+t,'Yi_zz_'+t)
      
      self.fibo_obj.data['dUe_Y_par_'+t]  = self.fibo_obj.calc_scalr('Pe_xx_'+t,'Ye_xx_'+t,'Pe_xy_'+t,'Ye_xy_'+t,'Pe_xz_'+t,'Ye_xz_'+t)
      self.fibo_obj.data['dUe_Y_par_'+t] += self.fibo_obj.calc_scalr('Pe_xy_'+t,'Ye_xy_'+t,'Pe_yy_'+t,'Ye_yy_'+t,'Pe_yz_'+t,'Ye_yz_'+t)
      self.fibo_obj.data['dUe_Y_par_'+t] += self.fibo_obj.calc_scalr('Pe_xz_'+t,'Ye_xz_'+t,'Pe_yz_'+t,'Ye_yz_'+t,'Pe_zz_'+t,'Ye_zz_'+t)
      
      self.fibo_obj.data['dU_Y_par_'+t]  = self.fibo_obj.calc_scalr('P_xx_'+t,'Y_xx_'+t,'P_xy_'+t,'Y_xy_'+t,'P_xz_'+t,'Y_xz_'+t)
      self.fibo_obj.data['dU_Y_par_'+t] += self.fibo_obj.calc_scalr('P_xy_'+t,'Y_xy_'+t,'P_yy_'+t,'Y_yy_'+t,'P_yz_'+t,'Y_yz_'+t)
      self.fibo_obj.data['dU_Y_par_'+t] += self.fibo_obj.calc_scalr('P_xz_'+t,'Y_xz_'+t,'P_yz_'+t,'Y_yz_'+t,'P_zz_'+t,'Y_zz_'+t)


    #total energy conversion to local internal energy (electromagnetic + surrounding internal + convection term)
    #self.fibo_obj.data['dUi_'+t] = self.fibo_obj.data['dUi_K_'+t] + self.fibo_obj.data['dUi_U_'+t] #+ self.fibo_obj.data['dUi_c_'+t]
    #self.fibo_obj.data['dUe_'+t] = self.fibo_obj.data['dUe_K_'+t] #+ self.fibo_obj.data['dUe_c_'+t]
    #self.fibo_obj.data['dU_'+t]  = self.fibo_obj.data['dU_F_'+t] + self.fibo_obj.data['dU_K_'+t] + self.fibo_obj.data['dU_U_'+t] #+ self.fibo_obj.data['dU_c_'+t]  
    
    #if par_per:
    #  self.fibo_obj.data['dUi_K_per_'+t] = self.fibo_obj.data['dUi_K_'+t] - self.fibo_obj.data['dUi_K_par_'+t]
    #  self.fibo_obj.data['dUi_U_per_'+t] = self.fibo_obj.data['dUi_U_'+t] - self.fibo_obj.data['dUi_U_par_'+t]
    #  self.fibo_obj.data['dUi_c_per_'+t] = self.fibo_obj.data['dUi_c_'+t] - self.fibo_obj.data['dUi_c_par_'+t]
    #  self.fibo_obj.data['dUe_K_per_'+t] = self.fibo_obj.data['dUe_K_'+t] - self.fibo_obj.data['dUe_K_par_'+t]
    #  self.fibo_obj.data['dUe_U_per_'+t] = self.fibo_obj.data['dUe_U_'+t] - self.fibo_obj.data['dUe_U_par_'+t]
    #  self.fibo_obj.data['dUe_c_per_'+t] = self.fibo_obj.data['dUe_c_'+t] - self.fibo_obj.data['dUe_c_par_'+t]
    #  self.fibo_obj.data['dU_F_per_'+t] = self.fibo_obj.data['dU_F_'+t] - self.fibo_obj.data['dU_F_par_'+t]
    #  self.fibo_obj.data['dU_K_per_'+t] = self.fibo_obj.data['dU_K_'+t] - self.fibo_obj.data['dU_K_par_'+t]
    #  self.fibo_obj.data['dU_U_per_'+t] = self.fibo_obj.data['dU_U_'+t] - self.fibo_obj.data['dU_U_par_'+t]
    #  self.fibo_obj.data['dU_c_per_'+t] = self.fibo_obj.data['dU_c_'+t] - self.fibo_obj.data['dU_c_par_'+t]

  #------------------------------------------------------------
  def calc_E_dec(self,
      time_str):
    """
    Decomposes E according to the generalized Ohm's law 
    (most general form: two Euler eqs summed, each multiplied by charge/mass)
    
    E = E_m + E_P + E_t + E_u
    
    Parameters :
      - time_str      [str] to indicate temporal exit ..
    
    """

    t = time_str
    twoden = np.reciprocal(self.fibo_obj.data['n_'+t] * (1. + 1./self.fibo_obj.meta['mime']))
    
    #motional field component
    self.fibo_obj.data['E_mi_x_'+t],self.fibo_obj.data['E_mi_y_'+t],self.fibo_obj.data['E_mi_z_'+t] = self.fibo_obj.calc_cross('ui_x_'+t,'B_x_'+t,'ui_y_'+t,'B_y_'+t,'ui_z_'+t,'B_z_'+t) 
    self.fibo_obj.data['E_mi_x_'+t] *= -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_mi_y_'+t] *= -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_mi_z_'+t] *= -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_me_x_'+t],self.fibo_obj.data['E_me_y_'+t],self.fibo_obj.data['E_me_z_'+t] = self.fibo_obj.calc_cross('ue_x_'+t,'B_x_'+t,'ue_y_'+t,'B_y_'+t,'ue_z_'+t,'B_z_'+t) 
    self.fibo_obj.data['E_me_x_'+t] *= -twoden
    self.fibo_obj.data['E_me_y_'+t] *= -twoden
    self.fibo_obj.data['E_me_z_'+t] *= -twoden 
    self.fibo_obj.data['E_m_x_'+t] = self.fibo_obj.data['E_mi_x_'+t] + self.fibo_obj.data['E_me_x_'+t]
    self.fibo_obj.data['E_m_y_'+t] = self.fibo_obj.data['E_mi_y_'+t] + self.fibo_obj.data['E_me_y_'+t]
    self.fibo_obj.data['E_m_z_'+t] = self.fibo_obj.data['E_mi_z_'+t] + self.fibo_obj.data['E_me_z_'+t]
    
    #pressure divergence term
    self.fibo_obj.data['E_Pi_x_'+t] = self.fibo_obj.calc_divr('Pi_xx_'+t,'Pi_xy_'+t,'Pi_xz_'+t) * twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_Pi_y_'+t] = self.fibo_obj.calc_divr('Pi_xy_'+t,'Pi_yy_'+t,'Pi_yz_'+t) * twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_Pi_z_'+t] = self.fibo_obj.calc_divr('Pi_xz_'+t,'Pi_yz_'+t,'Pi_zz_'+t) * twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_Pe_x_'+t] = - self.fibo_obj.calc_divr('Pe_xx_'+t,'Pe_xy_'+t,'Pe_xz_'+t) * twoden
    self.fibo_obj.data['E_Pe_y_'+t] = - self.fibo_obj.calc_divr('Pe_xy_'+t,'Pe_yy_'+t,'Pe_yz_'+t) * twoden
    self.fibo_obj.data['E_Pe_z_'+t] = - self.fibo_obj.calc_divr('Pe_xz_'+t,'Pe_yz_'+t,'Pe_zz_'+t) * twoden
    self.fibo_obj.data['E_P_x_'+t] = self.fibo_obj.data['E_Pi_x_'+t] + self.fibo_obj.data['E_Pe_x_'+t]
    self.fibo_obj.data['E_P_y_'+t] = self.fibo_obj.data['E_Pi_y_'+t] + self.fibo_obj.data['E_Pe_y_'+t]
    self.fibo_obj.data['E_P_z_'+t] = self.fibo_obj.data['E_Pi_z_'+t] + self.fibo_obj.data['E_Pe_z_'+t]
    
    #inertial term: temporal derivative
    self.fibo_obj.data['E_t_x_'+t],self.fibo_obj.data['E_t_y_'+t],self.fibo_obj.data['E_t_z_'+t] = self.fibo_obj.calc_curl('E_x_'+t,'E_y_'+t,'E_z_'+t)
    self.fibo_obj.data['E_t_x_'+t],self.fibo_obj.data['E_t_y_'+t],self.fibo_obj.data['E_t_z_'+t] = self.fibo_obj.calc_curl('E_t_x_'+t,'E_t_y_'+t,'E_t_z_'+t)
    self.fibo_obj.data['E_t_x_'+t] *=  -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_t_y_'+t] *=  -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_t_z_'+t] *=  -twoden / self.fibo_obj.meta['mime']
    
    #this other term is there if you want to confront this with Valentini's Ohm's law
    self.fibo_obj.data['E_tt_x_'+t],self.fibo_obj.data['E_tt_y_'+t],self.fibo_obj.data['E_tt_z_'+t] = self.fibo_obj.calc_grad(self.fibo_obj.calc_divr('E_x_'+t,'E_y_'+t,'E_z_'+t))
    self.fibo_obj.data['E_tt_x_'+t] *=  -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_tt_y_'+t] *=  -twoden / self.fibo_obj.meta['mime']
    self.fibo_obj.data['E_tt_z_'+t] *=  -twoden / self.fibo_obj.meta['mime']

    #inertial term: spatial derivative
    nuu_x = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.data['ui_x_'+t]
    nuu_y = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.data['ui_x_'+t]
    nuu_z = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.data['ui_x_'+t]
    self.fibo_obj.data['E_ui_x_'+t] = (self.fibo_obj.calc_divr(nuu_x,nuu_y,nuu_z)) * twoden / self.fibo_obj.meta['mime']
    nuu_x = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.data['ui_y_'+t]
    nuu_y = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.data['ui_y_'+t]
    nuu_z = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.data['ui_y_'+t]
    self.fibo_obj.data['E_ui_y_'+t] = (self.fibo_obj.calc_divr(nuu_x,nuu_y,nuu_z)) * twoden / self.fibo_obj.meta['mime']
    nuu_x = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_x_'+t] * self.fibo_obj.data['ui_z_'+t]
    nuu_y = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_y_'+t] * self.fibo_obj.data['ui_z_'+t]
    nuu_z = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ui_z_'+t] * self.fibo_obj.data['ui_z_'+t]
    self.fibo_obj.data['E_ui_z_'+t] = (self.fibo_obj.calc_divr(nuu_x,nuu_y,nuu_z)) * twoden / self.fibo_obj.meta['mime']
    
    nuu_x = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.data['ue_x_'+t]
    nuu_y = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.data['ue_x_'+t]
    nuu_z = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.data['ue_x_'+t]
    self.fibo_obj.data['E_ue_x_'+t] = (self.fibo_obj.calc_divr(nuu_x,nuu_y,nuu_z)) * twoden / self.fibo_obj.meta['mime']
    nuu_x = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.data['ue_y_'+t]
    nuu_y = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.data['ue_y_'+t]
    nuu_z = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.data['ue_y_'+t]
    self.fibo_obj.data['E_ue_y_'+t] = (self.fibo_obj.calc_divr(nuu_x,nuu_y,nuu_z)) * twoden / self.fibo_obj.meta['mime']
    nuu_x = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_x_'+t] * self.fibo_obj.data['ue_z_'+t]
    nuu_y = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_y_'+t] * self.fibo_obj.data['ue_z_'+t]
    nuu_z = self.fibo_obj.data['n_'+t] * self.fibo_obj.data['ue_z_'+t] * self.fibo_obj.data['ue_z_'+t]
    self.fibo_obj.data['E_ue_z_'+t] = (self.fibo_obj.calc_divr(nuu_x,nuu_y,nuu_z)) * twoden / self.fibo_obj.meta['mime']
    
    self.fibo_obj.data['E_u_x_'+t] = self.fibo_obj.data['E_ui_x_'+t] - self.fibo_obj.data['E_ue_x_'+t]
    self.fibo_obj.data['E_u_y_'+t] = self.fibo_obj.data['E_ui_y_'+t] - self.fibo_obj.data['E_ue_y_'+t]
    self.fibo_obj.data['E_u_z_'+t] = self.fibo_obj.data['E_ui_z_'+t] - self.fibo_obj.data['E_ue_z_'+t]
    
  #------------------------------------------------------------
  def calc_u_dec(self,
      time_str,
      par_per = False, 
      cE = True,
      cP = False,
      cn = False,
      cK = False,
      cX = False,
      cB = False):
    """
    Decomposes u parallel and perpendicular, and calculates fluid drifts
    
    Parameters :
      - time_str      [str] to indicate temporal exit ..
      - par_per       [bool] do you want also par and per? 
      - cE            [bool] do you want also cE?
      - cP            [bool] do you want also cP?
      - cn            [bool] do you want also cn?
      - cK            [bool] do you want also cK?
      - cX            [bool] do you want also cX?
      - cB            [bool] do you want also cB?
    
    """
    
    t = time_str
    
    inv_n = np.reciprocal(self.fibo_obj.data['n_'+t])
    inv_B = np.sqrt(np.reciprocal(  self.fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,'B_z_'+t,'B_z_'+t) ))
    
    self.fibo_obj.data['eB_x_'+t] = np.multiply(self.fibo_obj.data['B_x_'+t], inv_B)
    self.fibo_obj.data['eB_y_'+t] = np.multiply(self.fibo_obj.data['B_y_'+t], inv_B)
    self.fibo_obj.data['eB_z_'+t] = np.multiply(self.fibo_obj.data['B_z_'+t], inv_B)
    self.fibo_obj.data['invB_x_'+t] = np.multiply(self.fibo_obj.data['eB_x_'+t], inv_B)
    self.fibo_obj.data['invB_y_'+t] = np.multiply(self.fibo_obj.data['eB_y_'+t], inv_B)
    self.fibo_obj.data['invB_z_'+t] = np.multiply(self.fibo_obj.data['eB_z_'+t], inv_B)
    
    #----parallel-perpendicular-velocity--------------------------------------
    if par_per : 
      par_x,par_y,par_z, per_x,per_y,per_z = self.fibo_obj.calc_par_per('ui_x_'+t,'B_x_'+t,'ui_y_'+t,'B_y_'+t,'ui_z_'+t,'B_z_'+t)
      self.fibo_obj.data['ui_par_x_'+t], self.fibo_obj.data['ui_par_y_'+t], self.fibo_obj.data['ui_par_z_'+t] = par_x, par_y, par_z
      self.fibo_obj.data['ui_per_x_'+t], self.fibo_obj.data['ui_per_y_'+t], self.fibo_obj.data['ui_per_z_'+t] = per_x, per_y, per_z
      par_x,par_y,par_z, per_x,per_y,per_z = self.fibo_obj.calc_par_per('ue_x_'+t,'B_x_'+t,'ue_y_'+t,'B_y_'+t,'ue_z_'+t,'B_z_'+t)
      self.fibo_obj.data['ue_par_x_'+t],self.fibo_obj.data['ue_par_y_'+t],self.fibo_obj.data['ue_par_z_'+t] = par_x, par_y, par_z
      self.fibo_obj.data['ue_per_x_'+t],self.fibo_obj.data['ue_per_y_'+t],self.fibo_obj.data['ue_per_z_'+t] = per_x, per_y, per_z
      #par_x,par_y,par_z, per_x,per_y,per_z = self.fibo_obj.calc_par_per('u_x_'+t,'B_x_'+t,'u_y_'+t,'B_y_'+t,'u_z_'+t,'B_z_'+t)
      #self.fibo_obj.data['u_par_x_'+t],self.fibo_obj.data['u_par_y_'+t],self.fibo_obj.data['u_par_z_'+t] = par_x, par_y, par_z
      #self.fibo_obj.data['u_per_x_'+t],self.fibo_obj.data['u_per_y_'+t],self.fibo_obj.data['u_per_z_'+t] = per_x, per_y, per_z 
    
    
    
    #----E-cross-B-drift---------------------------------------
    if cE :     
      self.fibo_obj.data['cE_x_'+t], self.fibo_obj.data['cE_y_'+t], self.fibo_obj.data['cE_z_'+t] = self.fibo_obj.calc_cross('E_x_'+t,'invB_x_'+t,'E_y_'+t,'invB_y_'+t,'E_z_'+t,'invB_z_'+t)

    #----diamagnetic-(pressure)-drift--------------------------------------
    if cn : 

      gn_x = np.multiply( np.multiply(self.fibo_obj.calc_gradx('n_'+t), inv_n), inv_n)
      gn_y = np.multiply( np.multiply(self.fibo_obj.calc_grady('n_'+t), inv_n), inv_n)
      gn_z = np.multiply( np.multiply(self.fibo_obj.calc_gradz('n_'+t), inv_n), inv_n)

      self.fibo_obj.data['cni_x_'+t] = -self.fibo_obj.calc_scalr('Pi_xx_'+t,gn_x, 'Pi_xy_'+t,gn_y, 'Pi_xz_'+t,gn_z)
      self.fibo_obj.data['cni_y_'+t] = -self.fibo_obj.calc_scalr('Pi_xy_'+t,gn_x, 'Pi_yy_'+t,gn_y, 'Pi_yz_'+t,gn_z)
      self.fibo_obj.data['cni_z_'+t] = -self.fibo_obj.calc_scalr('Pi_xz_'+t,gn_x, 'Pi_yz_'+t,gn_y, 'Pi_zz_'+t,gn_z)
      self.fibo_obj.data['cni_x_'+t], self.fibo_obj.data['cni_y_'+t], self.fibo_obj.data['cni_z_'+t] = self.fibo_obj.calc_cross('cni_x_'+t,'invB_x_'+t,'cni_y_'+t,'invB_y_'+t,'cni_z_'+t,'invB_z_'+t)

      self.fibo_obj.data['cne_x_'+t] = self.fibo_obj.calc_scalr('Pe_xx_'+t,gn_x, 'Pe_xy_'+t,gn_y, 'Pe_xz_'+t,gn_z)
      self.fibo_obj.data['cne_y_'+t] = self.fibo_obj.calc_scalr('Pe_xy_'+t,gn_x, 'Pe_yy_'+t,gn_y, 'Pe_yz_'+t,gn_z)
      self.fibo_obj.data['cne_z_'+t] = self.fibo_obj.calc_scalr('Pe_xz_'+t,gn_x, 'Pe_yz_'+t,gn_y, 'Pe_zz_'+t,gn_z)
      self.fibo_obj.data['cne_x_'+t], self.fibo_obj.data['cne_y_'+t], self.fibo_obj.data['cne_z_'+t] = self.fibo_obj.calc_cross('cne_x_'+t,'invB_x_'+t,'cne_y_'+t,'invB_y_'+t,'cne_z_'+t,'invB_z_'+t)
    
    #----diamagnetic-(pressure)-drift--------------------------------------
    if cP : 

      self.fibo_obj.data['cPi_x_'+t] = -np.multiply(self.fibo_obj.calc_divr('Pi_xx_'+t,'Pi_xy_'+t,'Pi_xz_'+t), inv_n)
      self.fibo_obj.data['cPi_y_'+t] = -np.multiply(self.fibo_obj.calc_divr('Pi_xy_'+t,'Pi_yy_'+t,'Pi_yz_'+t), inv_n)
      self.fibo_obj.data['cPi_z_'+t] = -np.multiply(self.fibo_obj.calc_divr('Pi_xz_'+t,'Pi_yz_'+t,'Pi_zz_'+t), inv_n)
      self.fibo_obj.data['cPi_x_'+t], self.fibo_obj.data['cPi_y_'+t], self.fibo_obj.data['cPi_z_'+t] = self.fibo_obj.calc_cross('cPi_x_'+t,'invB_x_'+t,'cPi_y_'+t,'invB_y_'+t,'cPi_z_'+t,'invB_z_'+t)
      
      self.fibo_obj.data['cPe_x_'+t] = np.multiply(self.fibo_obj.calc_divr('Pe_xx_'+t,'Pe_xy_'+t,'Pe_xz_'+t), inv_n)
      self.fibo_obj.data['cPe_y_'+t] = np.multiply(self.fibo_obj.calc_divr('Pe_xy_'+t,'Pe_yy_'+t,'Pe_yz_'+t), inv_n)
      self.fibo_obj.data['cPe_z_'+t] = np.multiply(self.fibo_obj.calc_divr('Pe_xz_'+t,'Pe_yz_'+t,'Pe_zz_'+t), inv_n)
      self.fibo_obj.data['cPe_x_'+t], self.fibo_obj.data['cPe_y_'+t], self.fibo_obj.data['cPe_z_'+t] = self.fibo_obj.calc_cross('cPe_x_'+t,'invB_x_'+t,'cPe_y_'+t,'invB_y_'+t,'cPe_z_'+t,'invB_z_'+t)
      
      #self.fibo_obj.data['JP_x_'+t] = -self.fibo_obj.calc_divr('P_xx_'+t,'P_xy_'+t,'P_xz_'+t)
      #self.fibo_obj.data['JP_y_'+t] = -self.fibo_obj.calc_divr('P_xy_'+t,'P_yy_'+t,'P_yz_'+t)
      #self.fibo_obj.data['JP_z_'+t] = -self.fibo_obj.calc_divr('P_xz_'+t,'P_yz_'+t,'P_zz_'+t)
      #self.fibo_obj.data['JP_x_'+t], self.fibo_obj.data['JP_y_'+t], self.fibo_obj.data['JP_z_'+t] = self.fibo_obj.calc_cross('JP_x_'+t,'invB_x_'+t,'JP_y_'+t,'invB_y_'+t,'JP_z_'+t,'invB_z_'+t)
    
    
    
    #----curvature-inertial-drift-------------------------------------
    if cK :

      if not 'curvB_x_'+t in self.fibo_obj.data.keys() :
        self.fibo_obj.data['curvB_x_'+t]  = self.fibo_obj.data['eB_x_'+t] * self.fibo_obj.calc_gradx('eB_x_'+t) 
        self.fibo_obj.data['curvB_x_'+t] += self.fibo_obj.data['eB_y_'+t] * self.fibo_obj.calc_grady('eB_x_'+t) 
        self.fibo_obj.data['curvB_x_'+t] += self.fibo_obj.data['eB_z_'+t] * self.fibo_obj.calc_gradz('eB_x_'+t)
      if not 'curvB_y_'+t in self.fibo_obj.data.keys() :
        self.fibo_obj.data['curvB_y_'+t]  = self.fibo_obj.data['eB_x_'+t] * self.fibo_obj.calc_gradx('eB_y_'+t) 
        self.fibo_obj.data['curvB_y_'+t] += self.fibo_obj.data['eB_y_'+t] * self.fibo_obj.calc_grady('eB_y_'+t) 
        self.fibo_obj.data['curvB_y_'+t] += self.fibo_obj.data['eB_z_'+t] * self.fibo_obj.calc_gradz('eB_y_'+t)
      if not 'curvB_z_'+t in self.fibo_obj.data.keys() :
        self.fibo_obj.data['curvB_z_'+t]  = self.fibo_obj.data['eB_x_'+t] * self.fibo_obj.calc_gradx('eB_z_'+t) 
        self.fibo_obj.data['curvB_z_'+t] += self.fibo_obj.data['eB_y_'+t] * self.fibo_obj.calc_grady('eB_z_'+t) 
        self.fibo_obj.data['curvB_z_'+t] += self.fibo_obj.data['eB_z_'+t] * self.fibo_obj.calc_gradz('eB_z_'+t)
      
      self.fibo_obj.data['cCK_x_'+t],self.fibo_obj.data['cCK_y_'+t],self.fibo_obj.data['cCK_z_'+t] = self.fibo_obj.calc_cross('curvB_x_'+t,'invB_x_'+t,'curvB_y_'+t,'invB_y_'+t,'curvB_z_'+t,'invB_z_'+t)
      
      ui2_par = self.fibo_obj.calc_scalr('ui_par_x_'+t,'ui_par_x_'+t,'ui_par_y_'+t,'ui_par_y_'+t,'ui_par_z_'+t,'ui_par_z_'+t) 
      self.fibo_obj.data['cCKi_x_'+t] = self.fibo_obj.data['cCK_x_'+t] * ui2_par
      self.fibo_obj.data['cCKi_y_'+t] = self.fibo_obj.data['cCK_y_'+t] * ui2_par
      self.fibo_obj.data['cCKi_z_'+t] = self.fibo_obj.data['cCK_z_'+t] * ui2_par
      
      ue2_par = self.fibo_obj.calc_scalr('ue_par_x_'+t,'ue_par_x_'+t,'ue_par_y_'+t,'ue_par_y_'+t,'ue_par_z_'+t,'ue_par_z_'+t) / self.fibo_obj.meta['mime']
      self.fibo_obj.data['cCKe_x_'+t] = self.fibo_obj.data['cCK_x_'+t] * ue2_par
      self.fibo_obj.data['cCKe_y_'+t] = self.fibo_obj.data['cCK_y_'+t] * ue2_par
      self.fibo_obj.data['cCKe_z_'+t] = self.fibo_obj.data['cCK_z_'+t] * ue2_par
      
      #nu2_par = self.fibo_obj.data['n_'+t] * self.fibo_obj.calc_scalr('u_par_x_'+t,'u_par_x_'+t,'u_par_y_'+t,'u_par_y_'+t,'u_par_z_'+t,'u_par_z_'+t) 
      #self.fibo_obj.data['cCK_x_'+t] = self.fibo_obj.data['cCK_x_'+t] * nu2_par
      #self.fibo_obj.data['cCK_y_'+t] = self.fibo_obj.data['cCK_y_'+t] * nu2_par
      #self.fibo_obj.data['cCK_z_'+t] = self.fibo_obj.data['cCK_z_'+t] * nu2_par


    #----magnetic-structure-velocity-------------------------------------
    if cX :
    
      gradB = np.zeros((3,3)+self.fibo_obj.meta['nnn'])
      curlE = np.zeros((3,)+self.fibo_obj.meta['nnn'])
      
      gradB[0,0,:,:,:] = self.fibo_obj.calc_gradx('B_x_'+t)
      gradB[0,1,:,:,:] = self.fibo_obj.calc_gradx('B_y_'+t)
      gradB[0,2,:,:,:] = self.fibo_obj.calc_gradx('B_z_'+t)  
    
      gradB[1,0,:,:,:] = self.fibo_obj.calc_grady('B_x_'+t)
      gradB[1,1,:,:,:] = self.fibo_obj.calc_grady('B_y_'+t)
      gradB[1,2,:,:,:] = self.fibo_obj.calc_grady('B_z_'+t)  
    
      gradB[2,0,:,:,:] = self.fibo_obj.calc_gradz('B_x_'+t)
      gradB[2,1,:,:,:] = self.fibo_obj.calc_gradz('B_y_'+t)
      gradB[2,2,:,:,:] = self.fibo_obj.calc_gradz('B_z_'+t)    
      
      curlE[0,:,:,:], curlE[1,:,:,:], curlE[2,:,:,:] = self.fibo_obj.calc_curl('E_x_'+t,'E_y_'+t,'E_z_'+t)
      
      self.fibo_obj.data['cX_x_'+t] = np.zeros(self.fibo_obj.meta['nnn'])
      self.fibo_obj.data['cX_y_'+t] = np.zeros(self.fibo_obj.meta['nnn'])
      self.fibo_obj.data['cX_z_'+t] = np.zeros(self.fibo_obj.meta['nnn'])  

      if self.fibo_obj.meta['space_dim'] == '2D' : 
        inv_gradB = np.zeros((2,2))      
        
        for ix in range( 0, self.fibo_obj.meta['nx'] ) : 
          for iy in range( 0, self.fibo_obj.meta['ny'] ) : 
            for iz in range( 0, self.fibo_obj.meta['nz'] ) :   
              
              inv_tra_gradB = np.linalg.inv(np.transpose(gradB[0:2,0:2,ix,iy,iz]))
              self.fibo_obj.data['cX_x_'+t][ix,iy,iz], self.fibo_obj.data['cX_y_'+t][ix,iy,iz] = np.dot(inv_tra_gradB,curlE[0:2,ix,iy,iz])
    
      if self.fibo_obj.meta['space_dim'] == '3D' : 
        inv_gradB = np.zeros((3,3))      
        
        for ix in range( 0, self.fibo_obj.meta['nx'] ) : 
          for iy in range( 0, self.fibo_obj.meta['ny'] ) : 
            for iz in range( 0, self.fibo_obj.meta['nz'] ) :   
              
              inv_tra_gradB = np.linalg.inv(np.transpose(gradB[:,:,ix,iy,iz]))
              self.fibo_obj.data['cX_x_'+t][ix,iy,iz], self.fibo_obj.data['cX_y_'+t][ix,iy,iz], self.fibo_obj.data['cX_z_'+t][ix,iy,iz] = np.dot(inv_tra_gradB,curlE[:,ix,iy,iz])

    #----magnetic-field-line-velocity-------------------------------------
    if cB :
      
      if self.fibo_obj.meta['space_dim'] == '2D' : 
        
        invB_pla = np.reciprocal( self.fibo_obj.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,None,None) )
        
        self.fibo_obj.data['cB_x_'+t] = - invB_pla * np.multiply(self.fibo_obj.data['E_z_'+t], self.fibo_obj.data['B_y_'+t])
        self.fibo_obj.data['cB_y_'+t] =   invB_pla * np.multiply(self.fibo_obj.data['E_z_'+t], self.fibo_obj.data['B_x_'+t])
        self.fibo_obj.data['cB_z_'+t] = np.zeros(self.fibo_obj.meta['nnn'])
        
      if self.fibo_obj.meta['space_dim'] == '3D' : 
        
        print('ha ha still missing')