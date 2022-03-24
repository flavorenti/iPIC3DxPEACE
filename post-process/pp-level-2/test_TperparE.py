#--------------------------------------------------------------------
#  Code based on fibo library that (1) READS, (2) COMPUTE Tper,Tpar,Epar,Eper,
#  |Eper+VixB|,|Eper+VexB|, and (3) PRINT the new fields. Uses the derived two-dimensional 
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
inp_tars = ['B','E','Je','Ji','Ve','Vi','TXX','TXY','TXZ','TYY','TYZ','TZZ']
out_tars = ['Jper','Jpar','Eper','Epar','Eprimei','Eprimee','JdotE','JdotEe','JdotEi']#['Tiso','Tper','Tpar','Eper','Epar','Eprimei','Eprimee','Jper','Jpar','Tiso','Tper','Tpar']
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
      for tar in out_tars:
        if tar=='Eper' :
          alldata.calc_par_per('E_x','B_x','E_y','B_y','E_z','B_z',seg=seg)
          mod2Eper = alldata.calc_scalr('EperB_x','EperB_x','EperB_y','EperB_y','EperB_z','EperB_z',seg=seg,fill=False)
          alldata.data[tar+'%.8i'%int(seg)] = np.sqrt( mod2Eper )
        if tar=='Epar' :
          mod2Epar = alldata.calc_scalr('EparB_x','EparB_x','EparB_y','EparB_y','EparB_z','EparB_z',seg=seg,fill=False)
          alldata.data[tar+'%.8i'%int(seg)] = np.sqrt( mod2Epar )
        if tar=='Eprimei' :
          alldata.calc_cross('Vi_x','B_x','Vi_y','B_y','Vi_z','B_z',seg=seg)
          mod2VixB = alldata.calc_scalr('VixB_x','VixB_x','VixB_y','VixB_y','VixB_z','VixB_z',seg=seg,fill=False)
          alldata.data[tar+'%.8i'%int(seg)] = np.sqrt( mod2Eper+mod2VixB+2.*alldata.calc_scalr('EperB_x','VixB_x','EperB_y','VixB_y','EperB_z','VixB_z',seg=seg,fill=False) )
        if tar=='Eprimee' :
          alldata.calc_cross('Ve_x','B_x','Ve_y','B_y','Ve_z','B_z',seg=seg)
          mod2VexB = alldata.calc_scalr('VexB_x','VexB_x','VexB_y','VexB_y','VexB_z','VexB_z',seg=seg,fill=False)
          alldata.data[tar+'%.8i'%int(seg)] = np.sqrt( mod2Eper+mod2VexB+2.*alldata.calc_scalr('EperB_x','VexB_x','EperB_y','VexB_y','EperB_z','VexB_z',seg=seg,fill=False) )
        if tar=='Jper' or tar=='Jpar':
          alldata.data['J_x%.8i'%seg]=alldata.data['Je_x%.8i'%seg]+alldata.data['Ji_x%.8i'%seg]
          alldata.data['J_y%.8i'%seg]=alldata.data['Je_y%.8i'%seg]+alldata.data['Ji_y%.8i'%seg]
          alldata.data['J_z%.8i'%seg]=alldata.data['Je_z%.8i'%seg]+alldata.data['Ji_z%.8i'%seg]
          alldata.calc_par_per('J_x','B_x','J_y','B_y','J_z','B_z',seg=seg)
        if tar=='JdotE' :
          alldata.data[tar+'%.8i'%int(seg)] = alldata.calc_scalr('J_x','E_x','J_y','E_y','J_z','E_z',seg=seg,fill=False)
        if tar=='JdotEi' :
          alldata.data[tar+'%.8i'%int(seg)] = alldata.data['JdotE%.8i'%int(seg)] + alldata.calc_scalr('J_x','VixB_x','J_y','VixB_y','J_z','VixB_z',seg=seg,fill=False)
        if tar=='JdotEe' :
          alldata.data[tar+'%.8i'%int(seg)] = alldata.data['JdotE%.8i'%int(seg)] + alldata.calc_scalr('J_x','VexB_x','J_y','VexB_y','J_z','VexB_z',seg=seg,fill=False)
        for sp in alldata.meta['species']:
          if tar=='Tiso' :
            alldata.data[tar+sp+'%.8i'%int(seg)] = (alldata.data['TXX'+sp+'%.8i'%int(seg)]+alldata.data['TYY'+sp+'%.8i'%int(seg)]+alldata.data['TZZ'+sp+'%.8i'%int(seg)])/3.
          if tar=='Tper' :
            alldata.data['temp_x%.8i'%int(seg)] = alldata.calc_scalr('TXX'+sp,'B_x','TXY'+sp,'B_y','TXZ'+sp,'B_z',seg=seg,fill=False)
            alldata.data['temp_y%.8i'%int(seg)] = alldata.calc_scalr('TXY'+sp,'B_x','TYY'+sp,'B_y','TYZ'+sp,'B_z',seg=seg,fill=False)
            alldata.data['temp_z%.8i'%int(seg)] = alldata.calc_scalr('TXZ'+sp,'B_x','TYZ'+sp,'B_y','TZZ'+sp,'B_z',seg=seg,fill=False)
            mod2B = alldata.calc_scalr('B_x','B_x','B_y','B_y','B_z','B_z',seg=seg,fill=False)
            alldata.data[tar+sp+'%.8i'%int(seg)] = (3.*alldata.data['Tiso'+sp+'%.8i'%int(seg)]-(alldata.calc_scalr('temp_x','B_x','temp_y','B_y','temp_z','B_z',seg=seg,fill=False)/mod2B))
          if tar=='Tpar' :
            alldata.data[tar+sp+'%.8i'%int(seg)] = alldata.calc_scalr('temp_x','B_x','temp_y','B_y','temp_z','B_z',seg=seg,fill=False)/mod2B

      print(alldata.data.keys())
      

#----print-derived-fields-----------------------
      print('PRINT->',seg)
      for tar in out_tars:
        if 'Jp' in tar :
          alldata.print_vtk_vect(tar+'B_x',tar+'B_y',tar+'B_z',seg,data_address,alldata.meta['name']+'_'+tar+cut+'_'+str(seg),silent=True)
        elif 'E' in tar :
          alldata.print_vtk_scal(tar,seg,data_address,alldata.meta['name']+'_'+tar+cut+'_'+str(seg),silent=True)
        else :
          for sp in alldata.meta['species']:
            alldata.print_vtk_scal(tar+sp,seg,data_address,alldata.meta['name']+'_'+tar+sp+cut+'_'+str(seg),silent=True)

#----del-3D-fields----------------------------
      print('DELETE->',seg)
      for key in alldata.data.keys():
          del alldata.data[key]

#----print-comp-time--------------------
      print('COMP.TIME [seg%i] = %.3f s\n'%(seg,time.time()-t0))
