#------------------------------------------------------
# tests the basic functionality of from_VTK 
# 
#
# fada 19
#------------------------------------------------------
import fibo as fb

#======regulate=parameters=======================================================
data_address = '/home/flavorenti/Bureau/data_simu_iPIC3D/run5-0_512_3600/all_periodic_ipic-clean/data' 
fibo_name = 'ipic3d' 
#================================================================================

#----create-your-objects-----------
alldata = fb.fibo(fibo_name)
from_VTK = fb.from_VTK(data_address)

#----load-metadata------------------
from_VTK.get_meta(silent=False)
alldata.meta = from_VTK.meta
print('\n')

#----load-everything-in-code-units-------------------
for seg in alldata.meta['segcycles']:
  print('LOAD->',seg)
  from_VTK.get_vect(alldata.meta['name']+'_B_'    +str(seg),seg,fibo_obj=alldata,tar_var='B',silent=False)
  from_VTK.get_vect(alldata.meta['name']+'_E_'    +str(seg),seg,fibo_obj=alldata,tar_var='E')
  from_VTK.get_vect(alldata.meta['name']+'_Je_'   +str(seg),seg,fibo_obj=alldata,tar_var='Je')
  from_VTK.get_vect(alldata.meta['name']+'_Ji_'   +str(seg),seg,fibo_obj=alldata,tar_var='Ji')
  from_VTK.get_scal(alldata.meta['name']+'_PXXe0_'+str(seg),seg,fibo_obj=alldata,tar_var='PXXe0')
  from_VTK.get_scal(alldata.meta['name']+'_PXXi1_'+str(seg),seg,fibo_obj=alldata,tar_var='PXXi1')
  from_VTK.get_scal(alldata.meta['name']+'_PXYe0_'+str(seg),seg,fibo_obj=alldata,tar_var='PXYe0')
  from_VTK.get_scal(alldata.meta['name']+'_PXYi1_'+str(seg),seg,fibo_obj=alldata,tar_var='PXYi1')
  from_VTK.get_scal(alldata.meta['name']+'_PXZe0_'+str(seg),seg,fibo_obj=alldata,tar_var='PXZe0')
  from_VTK.get_scal(alldata.meta['name']+'_PXZi1_'+str(seg),seg,fibo_obj=alldata,tar_var='PXZi1')
  from_VTK.get_scal(alldata.meta['name']+'_PYYe0_'+str(seg),seg,fibo_obj=alldata,tar_var='PYYe0')
  from_VTK.get_scal(alldata.meta['name']+'_PYYi1_'+str(seg),seg,fibo_obj=alldata,tar_var='PYYi1')
  from_VTK.get_scal(alldata.meta['name']+'_PYZe0_'+str(seg),seg,fibo_obj=alldata,tar_var='PYZe0')
  from_VTK.get_scal(alldata.meta['name']+'_PYZi1_'+str(seg),seg,fibo_obj=alldata,tar_var='PYZi1')
  from_VTK.get_scal(alldata.meta['name']+'_PZZe0_'+str(seg),seg,fibo_obj=alldata,tar_var='PZZe0')
  from_VTK.get_scal(alldata.meta['name']+'_PZZi1_'+str(seg),seg,fibo_obj=alldata,tar_var='PZZi1')
  from_VTK.get_scal(alldata.meta['name']+'_rhoe0_'+str(seg),seg,fibo_obj=alldata,tar_var='rhoe0')
  from_VTK.get_scal(alldata.meta['name']+'_rhoi1_'+str(seg),seg,fibo_obj=alldata,tar_var='rhoi1')

#----change-units------------------------------
units = 'SW' # possible units are iPIC, SW, SI
alldata.calc_units(units, fibo_obj=alldata, silent=False)

#----calculate-derived-fields-with-aritm-operations------------------
for seg in alldata.meta['segcycles']:
  print('CALC->',seg)
  # Velocity = Je,rhoe0 --> Ve_x,Ve_y,Ve_z
  alldata.calc_divid('Je','rhoe0',new_tar='Ve',seg=seg)
  # Velocity = Ji,rhoi1 --> Vi_x,Vi_y,Vi_z
  alldata.calc_divid('Ji','rhoi1',new_tar='Vi',seg=seg)
  # module of B
  alldata.calc_scalr('B_x','B_x','B_y','B_y','B_z','B_z',seg=seg)
  # module of E
  alldata.calc_scalr('E_x','E_x','E_y','E_y','E_z','E_z',seg=seg)
  # VixB
  alldata.calc_cross('Vi_x','B_x','Vi_y','B_y','Vi_z','B_z',seg=seg)
  # E perpendicular to B --> EperB
  # E parallel to B --> EparB
  alldata.calc_par_per('E_x','B_x','E_y','B_y','E_z','B_z',seg=seg)


#----print-derived-fields-----------------------
for seg in alldata.meta['segcycles']:
  print('PRINT->',seg)
  alldata.print_vtk_vect('Ve_x','Ve_y','Ve_z',seg,data_address,alldata.meta['name']+'_Ve_'+str(seg))
  alldata.print_vtk_vect('Vi_x','Vi_y','Vi_z',seg,data_address,alldata.meta['name']+'_Vi_'+str(seg))

'''
#---extras-----------------------------------
cycle_ok = int( alldata.get_time_exit(from_VTK.segs,100.,silent=False)[0] )
print(cycle_ok)

#---plots------------------------------------
key='rhoi1'
xlab = {'text':'$x/d_i$','fontsize':12,'usetex':True}
ylab = {'text':'$y/d_i$','fontsize':12,'usetex':True}
import matplotlib.pyplot as plt
plot = alldata.draw_canvas([[key],[key],[key]],'line')
alldata.draw_spotted(plot[0],alldata.data[key+'00000400'], range_z=[30,31], label_xyb=[xlab,ylab,{'text':key,'usetex':True,'rotation':90}])
alldata.draw_spotted(plot[1],alldata.data[key+'00000400'], range_z=[30,31], label_xyb=[xlab,ylab,{'text':key,'usetex':True,'rotation':90}])
alldata.draw_spotted(plot[2],alldata.data[key+'00000400'], range_z=[30,31], label_xyb=[xlab,ylab,{'text':key,'usetex':True,'rotation':90}])
plt.show()
'''
