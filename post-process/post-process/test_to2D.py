#--------------------------------------------------------------------
#  Code based on fibo library that (1) READS, (2) CUT and (3) PRINT
#  the three-dimensional output of an iPIC3D simu. The cuts are done
#  in Lz/2 and Ly/2. The cuts are saved in the same dir as the data,
#  named as e.g. Beq, Bdp for equatorial and dipolar cut.
#  Add global path to data dir in command line.
#  Only parameters to change are:
#
#  tars         = [str] target variables you want to use (in order!)
#  cycle_min    = [int] minimum cycle to read
#  cycle_max    = [int] maximum cycle to read
#  diffRho      = [bool] do you want to compute total charge rhoi+rhoe?
#
# Job 11/2021
#-------------------------------------------------------------------

######################################################################################
#======regulate=parameters=============================================================
tars = ['B','E','Je','Ji','PXX','PXY','PXZ','PYY','PYZ','PZZ','rho']
cycle_min = 0
cycle_max  = 200000
offset = 0
diffRho = True
#=======================================================================================
########################################################################################

########################################################################################
####################### DO NOT TOUCH BELOW! ############################################
########################################################################################
#----import-modules---------------
import fibo as fb 
import time
import sys
import numpy as np
from scipy.ndimage.filters import gaussian_filter as gf

#----path-to-run-------------------
data_address = sys.argv[1]

#----create-your-objects-----------
fibo_name = 'ipic3d' 
data3D = fb.fibo(fibo_name)
from_VTK = fb.from_VTK(data_address)

#----load-metadata------------------
from_VTK.get_meta(silent=False)
data3D.meta = from_VTK.meta
nx,ny,nz = from_VTK.meta['nnn']

#---define-cuts-you-want----------
cuts  = []
fibos = []
cuts.append(['eq',[0,nx],[0,ny],[int(nz/2),int(nz/2)+1]])
cuts.append(['dp',[0,nx],[int(ny/2),int(ny/2)+1],[0,nz]])
for cut in cuts:
  data2D = fb.fibo(fibo_name)
  data2D.meta = from_VTK.meta
  fibos.append(data2D)

#----load-everything-in-code-units-------------------
for seg in data3D.meta['segcycles']:
  seg=seg+offset
  t0 = time.time()
  if cycle_min<=seg<=cycle_max :
    print('LOAD->',seg)
    for tar in tars:
      if (('P' in tar) or ('rho' in tar)) :
        for sp in data3D.meta['species']:
          from_VTK.get_scal(data3D.meta['name']+'_'+tar+sp+'_'+str(seg),seg,fibo_obj=data3D,tar_var=tar+sp,silent=True)
      else:
        from_VTK.get_vect(data3D.meta['name']+'_'+tar+'_'+str(seg),seg,fibo_obj=data3D,tar_var=tar,silent=True)

#----cut-to-2D------------------
    print('CUT->',seg)
    for ic,cut in enumerate(cuts):
      for tar in tars:
        if (('P' in tar) or ('rho' in tar)) :
          for sp in data3D.meta['species']:
            fibos[ic].data[tar+sp+'%.8i'%int(seg)], fibos[ic].meta = data3D.extract_range(tar+sp+'%.8i'%int(seg),cut[1],cut[2],cut[3])
        else:
          fibos[ic].data[tar+'_x%.8i'%int(seg)], fibos[ic].meta = data3D.extract_range(tar+'_x%.8i'%int(seg),cut[1],cut[2],cut[3])
          fibos[ic].data[tar+'_y%.8i'%int(seg)], fibos[ic].meta = data3D.extract_range(tar+'_y%.8i'%int(seg),cut[1],cut[2],cut[3])
          fibos[ic].data[tar+'_z%.8i'%int(seg)], fibos[ic].meta = data3D.extract_range(tar+'_z%.8i'%int(seg),cut[1],cut[2],cut[3])
      if diffRho :
        fibos[ic].data['Drho%.8i'%int(seg)] = np.zeros(np.shape(fibos[ic].data['rhoi1%.8i'%int(seg)][:,:,:]))
        for isp,sp in enumerate(data3D.meta['species']):
          fibos[ic].data['Drho%.8i'%int(seg)] += fibos[ic].data['rho'+sp+'%.8i'%int(seg)]
        fibos[ic].data['Drho%.8i'%int(seg)] = gf(fibos[ic].data['Drho%.8i'%int(seg)],3)

#----print-derived-fields-----------------------
    print('PRINT->',seg)
    for ic,cut in enumerate(cuts):
      for tar in tars:
        if (('P' in tar) or ('rho' in tar)) :
          for sp in data3D.meta['species']:
            fibos[ic].print_vtk_scal(tar+sp,seg,data_address,data3D.meta['name']+'_'+tar+sp+cut[0]+'_'+str(seg),silent=True)
        else:
          fibos[ic].print_vtk_vect(tar+'_x',tar+'_y',tar+'_z',seg,data_address,data3D.meta['name']+'_'+tar+cut[0]+'_'+str(seg),silent=True)
      if diffRho :
        fibos[ic].print_vtk_scal('Drho',seg,data_address,data3D.meta['name']+'_Drho'+cut[0]+'_'+str(seg),silent=True)

#----del-3D-fields----------------------------
    print('DELETE->',seg)
    for ic,cut in enumerate(cuts):
      for key in fibos[ic].data.keys():
        del fibos[ic].data[key]

#----print-comp-time--------------------
    print('COMP.TIME [seg%i] = %.3f s\n'%(seg,time.time()-t0))
