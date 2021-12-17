#--------------------------------------------------------------------
#  Code based on fibo library that (1) READS inp_tars, (2) SELECT the
#  values of these fields in a given region, (3) PRINT these values on .txt
#  The .txt are saved in the dir Closure_files in the same dir as python script.
#  Add global path to data dir in command line.
#  Only parameters to change are:
#
#  inp_tars     = [str] target variables you want to use (in order!)
#  limits       = [list] limits to compute subset [[xmin,xmax],[ymin,ymax],[zmin,zmax]] 
#  cuts         = [str] cuts to use if the data are 2D (eg data1), if [''] data are 3D (eg data2)
#  cycle_min    = [int] minimum cycle to read
#  cycle_max    = [int] maximum cycle to read
#
#  Extra rules: only scalars are admitted in inp_tars!
#
# Job 12/2021
#-------------------------------------------------------------------

######################################################################################
#======regulate=parameters=============================================================
inp_tars = ['Tiso','Tpar','Tper','den']
limits = [[22.5,24.],[19.25,20.75],[0.,0.]]
cuts = ['eq']
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

#----compute-limits--------------- 
  ix1 = int(limits[0][0]/from_VTK.meta['dx'])
  ix2 = int(limits[0][1]/from_VTK.meta['dx'])
  iy1 = int(limits[1][0]/from_VTK.meta['dy'])
  iy2 = int(limits[1][1]/from_VTK.meta['dy'])
  iz1 = int(limits[2][0]/from_VTK.meta['dz'])
  iz2 = int(limits[2][1]/from_VTK.meta['dz']) 

#----load-everything-in-code-units-------------------
  for seg in alldata.meta['segcycles']:
    t0 = time.time()
    if cycle_min<=seg<=cycle_max :
      print('LOAD->',seg)
      for tar in inp_tars:
        for sp in alldata.meta['species']:
          from_VTK.get_scal(alldata.meta['name']+'_'+tar+sp+cut+'_'+str(seg),seg,fibo_obj=alldata,tar_var=tar+sp,silent=True)

      print(alldata.data.keys())

#----compute-derived-fields------------------
      print('SELECT+PRINT->',seg)
      for tar in inp_tars:
        for sp in alldata.meta['species']:
          file_out = open('closure_files/closure_'+tar+sp+cut+'.txt','a')
          if int(seg)==0 : file_out.write('#seg\t value\n')
          if ix1>=ix2 :
              out = np.ravel(alldata.data[tar+sp+'%.8i'%int(seg)][:,iy1:iy2,iz1:iz2])
          if iy1>=iy2 :
            out = np.ravel(alldata.data[tar+sp+'%.8i'%int(seg)][ix1:ix2,:,iz1:iz2])
          if iz1>=iz2 :
            out = np.ravel(alldata.data[tar+sp+'%.8i'%int(seg)][ix1:ix2,iy1:iy2,:])
          out = out[~np.isnan(out)]
          out = out[~np.isinf(out)]
          file_out.write('%i\t%.9f\n'%(int(seg),np.mean(out)))
          file_out.close()

#----del-3D-fields----------------------------
      print('DELETE->',seg)
      for key in alldata.data.keys():
          del alldata.data[key]

#----print-comp-time--------------------
      print('COMP.TIME [seg%i] = %.3f s\n'%(seg,time.time()-t0))
