#--------------------------------------------------------------------
#  Code based on fibo library that (1) READS, (2) COMPUTE Anis,Gyro
#  and (3) PRINT the new fields. Uses the derived two-dimensional 
#  and three-dimensional output of an iPIC3D simu.
#  The new fields are saved in the same dir as the data, with the same dimensions.
#  Add global path to data dir in command line.
#  Only parameters to change are:
#
#  inp_tars     = [str] target variables you want to use (in order!)
#  out_tars     = [str] name of the variables creates to print (only scalars!) must follow the order!
#  cuts         = [str] cuts to use if the data are 2D (eg data1), if [''] data are 3D (eg data2)
#  cycle_min    = [int] minimum cycle to read
#  cycle_max    = [int] maximum cycle to read
#
#  Extra rules for out_tars: Epar require Eper, Eprime require Eper, Tper require Tiso, 
#                            Tpar require Tper
#
# Job 12/2021
#-------------------------------------------------------------------

######################################################################################
#======regulate=parameters=============================================================
inp_tars = ['B','TXX','TXY','TXZ','TYY','TYZ','TZZ','Tpar','Tper']
out_tars = ['Anis','Gyro']
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

#----compute-derived-fields------------------
      print('COMPUTE->',seg)
      for sp in alldata.meta['species']:
        I1  = alldata.data['TXX'+sp+'%.8i'%int(seg)]+alldata.data['TYY'+sp+'%.8i'%int(seg)]+alldata.data['TZZ'+sp+'%.8i'%int(seg)]
        I2  = alldata.calc_scalr('TXX'+sp,'TYY'+sp,'TXX'+sp,'TZZ'+sp,'TYY'+sp,'TZZ'+sp,seg=seg,fill=False)
        I2 -= alldata.calc_scalr('TXY'+sp,'TXY'+sp,'TXZ'+sp,'TXZ'+sp,'TYZ'+sp,'TYZ'+sp,seg=seg,fill=False)
        for tar in out_tars:
          if tar=='Anis' :
            alldata.data[tar+sp+'%.8i'%int(seg)] = 2.*alldata.data['Tpar'+sp+'%.8i'%int(seg)]/alldata.data['Tper'+sp+'%.8i'%int(seg)]
          if tar=='Gyro' :
            alldata.data[tar+sp+'%.8i'%int(seg)] = 1. - (4.*I2/(I1-alldata.data['Tpar'+sp+'%.8i'%int(seg)])/(I1+3.*alldata.data['Tpar'+sp+'%.8i'%int(seg)]))

      print(alldata.data.keys())
      

#----print-derived-fields-----------------------
      print('PRINT->',seg)
      for tar in out_tars:
        for sp in alldata.meta['species']:
          alldata.print_vtk_scal(tar+sp,seg,data_address,alldata.meta['name']+'_'+tar+sp+cut+'_'+str(seg),silent=True)

#----del-3D-fields----------------------------
      print('DELETE->',seg)
      del I1,I2
      for key in alldata.data.keys():
          del alldata.data[key]

#----print-comp-time--------------------
      print('COMP.TIME [seg%i] = %.3f s\n'%(seg,time.time()-t0))
