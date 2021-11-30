#------------------------------------------------------
# tests the basic functionality of from_VTK 
# 
#
# fada 19
#------------------------------------------------------
import fibo as fb

#======regulate=parameters=======================================================
data_address = '/home/flavorenti/Bureau/data_simu_iPIC3D/Mercury_SaeInit/PR1/run1/data' 
fibo_name = 'ipic3d' 
#================================================================================

#----create-your-objects-----------
data3D = fb.fibo(fibo_name)
data2D_eq = fb.fibo(fibo_name)
data2D_dp = fb.fibo(fibo_name)
from_VTK = fb.from_VTK(data_address)

#----load-metadata------------------
from_VTK.get_meta(silent=False)
data3D.meta = from_VTK.meta
data2D_eq.meta = from_VTK.meta
data2D_dp.meta = from_VTK.meta
nx,ny,nz = from_VTK.meta['nnn']
print('\n')

#----load-everything-in-code-units-------------------
for seg in data3D.meta['segcycles']:
  print('LOAD->',seg)
  from_VTK.get_vect(data3D.meta['name']+'_B_'    +str(seg),seg,fibo_obj=data3D,tar_var='B',silent=False)
  #from_VTK.get_vect(data3D.meta['name']+'_E_'    +str(seg),seg,fibo_obj=data3D,tar_var='E')
  #from_VTK.get_vect(data3D.meta['name']+'_Je_'   +str(seg),seg,fibo_obj=data3D,tar_var='Je')
  #from_VTK.get_vect(data3D.meta['name']+'_Ji_'   +str(seg),seg,fibo_obj=data3D,tar_var='Ji')
  #from_VTK.get_scal(data3D.meta['name']+'_PXXe0_'+str(seg),seg,fibo_obj=data3D,tar_var='PXXe0')
  #from_VTK.get_scal(data3D.meta['name']+'_PXXi1_'+str(seg),seg,fibo_obj=data3D,tar_var='PXXi1')
  #from_VTK.get_scal(data3D.meta['name']+'_PXYe0_'+str(seg),seg,fibo_obj=data3D,tar_var='PXYe0')
  #from_VTK.get_scal(data3D.meta['name']+'_PXYi1_'+str(seg),seg,fibo_obj=data3D,tar_var='PXYi1')
  #from_VTK.get_scal(data3D.meta['name']+'_PXZe0_'+str(seg),seg,fibo_obj=data3D,tar_var='PXZe0')
  #from_VTK.get_scal(data3D.meta['name']+'_PXZi1_'+str(seg),seg,fibo_obj=data3D,tar_var='PXZi1')
  #from_VTK.get_scal(data3D.meta['name']+'_PYYe0_'+str(seg),seg,fibo_obj=data3D,tar_var='PYYe0')
  #from_VTK.get_scal(data3D.meta['name']+'_PYYi1_'+str(seg),seg,fibo_obj=data3D,tar_var='PYYi1')
  #from_VTK.get_scal(data3D.meta['name']+'_PYZe0_'+str(seg),seg,fibo_obj=data3D,tar_var='PYZe0')
  #from_VTK.get_scal(data3D.meta['name']+'_PYZi1_'+str(seg),seg,fibo_obj=data3D,tar_var='PYZi1')
  #from_VTK.get_scal(data3D.meta['name']+'_PZZe0_'+str(seg),seg,fibo_obj=data3D,tar_var='PZZe0')
  #from_VTK.get_scal(data3D.meta['name']+'_PZZi1_'+str(seg),seg,fibo_obj=data3D,tar_var='PZZi1')
  #from_VTK.get_scal(data3D.meta['name']+'_rhoe0_'+str(seg),seg,fibo_obj=data3D,tar_var='rhoe0')
  #from_VTK.get_scal(data3D.meta['name']+'_rhoi1_'+str(seg),seg,fibo_obj=data3D,tar_var='rhoi1')

#----cut-to-2D------------------
for seg in data3D.meta['segcycles']:
  print('CUT->',seg)
  data2D_eq.data['B_x%.8i'%int(seg)], data2D_eq.meta = data3D.extract_range('B_x%.8i'%int(seg),[0,nx],[0,ny],[int(nz/2),int(nz/2)+1])
  data2D_eq.data['B_y%.8i'%int(seg)], data2D_eq.meta = data3D.extract_range('B_y%.8i'%int(seg),[0,nx],[0,ny],[int(nz/2),int(nz/2)+1])
  data2D_eq.data['B_z%.8i'%int(seg)], data2D_eq.meta = data3D.extract_range('B_z%.8i'%int(seg),[0,nx],[0,ny],[int(nz/2),int(nz/2)+1])

#----print-derived-fields-----------------------
for seg in data3D.meta['segcycles']:
  print('PRINT->',seg)
  data2D_eq.print_vtk_vect('B_x','B_y','B_z',seg,data_address,data3D.meta['name']+'_Beq_'+str(seg))


