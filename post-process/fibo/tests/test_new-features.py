#------------------------------------------------------
# tests the basic functionality of from_VTK 
# 
#
# fada 19
#------------------------------------------------------
import fibo as fb
from scipy.ndimage.filters import gaussian_filter as gf

#======regulate=parameters=======================================================
data_address = '/home/flavorenti/Bureau/data_simu_iPIC3D/PR0_tests-new-features/'
run_list = ['000','100','010','110','001','101','011','111']
fibo_name = 'ipic3d' 
#================================================================================


#----create-your-objects-----------
alldata = fb.fibo(fibo_name)
for rr in run_list:
  from_VTK = fb.from_VTK(data_address+rr+'/data/')

  #----load-metadata------------------
  from_VTK.get_meta(silent=False)
  alldata.meta = from_VTK.meta
  alldata.segs = from_VTK.segs
  print(alldata.segs)
  print('\n')

  #----load-everything-in-code-units-------------------
  for seg in alldata.meta['time2seg']:
    print(seg)
    from_VTK.get_scal(alldata.meta['name']+'_rhoe0_'+str(seg),seg,fibo_obj=alldata,tar_var='rhoe0'+rr)
    from_VTK.get_scal(alldata.meta['name']+'_rhoi1_'+str(seg),seg,fibo_obj=alldata,tar_var='rhoi1'+rr)
    from_VTK.get_scal(alldata.meta['name']+'_rhoe2_'+str(seg),seg,fibo_obj=alldata,tar_var='rhoe2'+rr)
    from_VTK.get_scal(alldata.meta['name']+'_rhoi3_'+str(seg),seg,fibo_obj=alldata,tar_var='rhoi3'+rr)
  print(alldata.data.keys())

  #----compute-derived-fields----------------------
  alldata.data['rho'+rr+seg.zfill(8)]=gf(alldata.data['rhoe0'+rr+seg.zfill(8)]+alldata.data['rhoi1'+rr+seg.zfill(8)],3)

#---plots------------------------------------
key  = 'rho'
time = '00004500'
xlab = {'text':'','fontsize':18,'usetex':True}
ylab = {'text':'','fontsize':18,'usetex':True}
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plot = alldata.draw_canvas([[key,key,key,key],[key,key,key,key]],'line')
i=0
for rr in run_list:
  alldata.draw_spotted(plot[i], key+rr+time, range_x=[0,256], range_y=[104,105], range_z=[0,208], args_fig={'extent':(0,25,0,20),'cmap':cm.seismic}, label_xyb=[xlab,ylab,None])
  plot[i].text(2.,15.,rr,fontsize=16)
  plot[i].set_xlim(0,20.)
  plot[i].set_ylim(0,25.)
  i+=1

#---add-planet-to-plots----------------------
for pp in plot:
  planet = plt.Circle((alldata.meta['zc']+alldata.meta['Doff'], alldata.meta['xc']), alldata.meta['R'], color='grey', fill=False, linewidth=2.)
  pp.add_patch(planet)

#---show-plots----------------------
plt.show()
#plt.savefig('/home/flavorenti/Bureau/Global-simulations-Mercury/global/Rhoi1SW_PR0_tests-new-features_cutXY_'+time+'.png',format='png')
