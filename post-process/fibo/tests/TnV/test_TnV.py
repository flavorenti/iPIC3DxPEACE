#--------------------------------------------------------------------
#  Code based on fibo library that (1) READS, (2) COMPUTE n,T,V and (3) PRINT
#  the derived two-dimensional and three-dimensional output of an iPIC3D simu.
#  The computation starts form the pristine output of the simu (either 2D or 3D)
#  The new fields are saved in the same dir as the data, with the same dimensions,
#  named as follows PXXi1 -> TXXi1, rhoi1 -> ni1, Ji -> Vi etc.
#  Add global path to data dir in command line.
#  Only parameters to change are:
#
#  inp_tars     = [str] target variables you want to use as input (in order!)
#  out_tars     = [str] name variables you want to create as output (same order)
#  cuts         = [str] cuts to use if the data are 2D (eg data1), if [''] data are 3D (eg data2)
#  cycle_min    = [int] minimum cycle to read
#  cycle_max    = [int] maximum cycle to read
#
# Job 12/2021
#-------------------------------------------------------------------

######################################################################################
#======regulate=parameters=============================================================
inp_tars = ['Je','Ji','PXX','PXY','PXZ','PYY','PYZ','PZZ','rho']
out_tars = ['Ve','Vi','TXX','TXY','TXZ','TYY','TYZ','TZZ','den']
cuts = ['eq','dp']
cycle_min = 0
cycle_max  = 2000000
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

#----path-to-run-------------------
data_address = sys.argv[1]

#----create-your-objects-----------
fibo_name = 'ipic3d' 
alldata = fb.fibo(fibo_name)
from_VTK = fb.from_VTK(data_address)

#----loop-on-types-of-data--------
for cut in cuts:
  
  #----load-metadata------------------
  from_VTK.get_meta(silent=False)
  if cut == 'eq':
    from_VTK.meta['nz']=1
    from_VTK.meta['nnn']=(from_VTK.meta['nx'],from_VTK.meta['ny'],1)
  if cut == 'dp':
    from_VTK.meta['ny']=1
    from_VTK.meta['nnn']=(from_VTK.meta['nx'],1,from_VTK.meta['nz'])
  alldata.meta = from_VTK.meta

  #----load-everything-in-code-units-------------------
  for seg in alldata.meta['segcycles']:
    t0 = time.time()
    if cycle_min<=seg<=cycle_max :
      print('LOAD->',seg)
      for tar in inp_tars:
        if len(tar)<3 :
          from_VTK.get_vect(alldata.meta['name']+'_'+tar+cut+'_'+str(seg),seg,fibo_obj=alldata,tar_var=tar,silent=True)
        else :
          for sp in alldata.meta['species']:
            from_VTK.get_scal(alldata.meta['name']+'_'+tar+sp+cut+'_'+str(seg),seg,fibo_obj=alldata,tar_var=tar+sp,silent=True)

      print(alldata.data.keys())

#----change-units------------------------------
      alldata.calc_units('iPIC', fibo_obj=alldata, silent=False)

#----compute-derived-fields------------------
      print('COMPUTE->',seg)
      for tar in out_tars:
        if tar=='Vi' :
          alldata.calc_divid('Ji','rhoi1',new_tar=tar,seg=seg)
        if tar=='Ve' :
          alldata.calc_divid('Je','rhoe0',new_tar=tar,seg=seg)
        if 'T' in tar :
          for sp in alldata.meta['species']:
            alldata.calc_temp(tar.replace('T','P')+sp,new_tar=tar+sp,seg=seg)
        if tar=='den' :
          for isp,sp in enumerate(alldata.meta['species']):
            alldata.data[tar+sp+'%.8i'%int(seg)] = alldata.data['rho'+sp+'%.8i'%int(seg)] / np.sign(alldata.meta['sQOM'][isp])
    
      print(alldata.data.keys())

#----print-derived-fields-----------------------
      print('PRINT->',seg)
      for tar in out_tars:
        if len(tar)<3 :
          alldata.print_vtk_vect(tar+'_x',tar+'_y',tar+'_z',seg,data_address,alldata.meta['name']+'_'+tar+cut+'_'+str(seg),silent=True)
        else :
          for sp in alldata.meta['species']:
            alldata.print_vtk_scal(tar+sp,seg,data_address,alldata.meta['name']+'_'+tar+sp+cut+'_'+str(seg),silent=True)

#----del-3D-fields----------------------------
      print('DELETE->',seg)
      for tar in (inp_tars+out_tars):
        if len(tar)<3 :
          del alldata.data[tar+'_x'+'%.8i'%int(seg)], alldata.data[tar+'_y'+'%.8i'%int(seg)], alldata.data[tar+'_z'+'%.8i'%int(seg)]
        else :
          for sp in alldata.meta['species']:
            del alldata.data[tar+sp+'%.8i'%int(seg)]

#----print-comp-time--------------------
      print('COMP.TIME [seg%i] = %.3f s\n'%(seg,time.time()-t0))
