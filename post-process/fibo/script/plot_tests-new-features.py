'''
Script to plot the output of the 8 test runs with the same setup as PR0
but many new features:
1) Exosphere (only H+ and e- for now)
2) Reinject Particles entering planet to keep net charge zero
3) new particle initialization to avoid pcls inside the planet 

The script loads the data of 8 runs specified by which of these
three features is active (I for active, 0 for silent).
'''

import fibo as fb

data_path = ''

alldata = fb.fibo('')
from_VTK000 = fb.from_VTK(data_path+'/000/data')
from_VTK100 = fb.from_VTK(data_path+'/100/data')
from_VTK010 = fb.from_VTK(data_path+'/010/data')
from_VTK110 = fb.from_VTK(data_path+'/110/data')
from_VTK001 = fb.from_VTK(data_path+'/001/data')
from_VTK101 = fb.from_VTK(data_path+'/101/data')
from_VTK011 = fb.from_VTK(data_path+'/011/data')
from_VTK111 = fb.from_VTK(data_path+'/111/data')

alldata.meta = from_VTK000.meta
alldata.segs = from_VTK000.segs

for seg in alldata.meta['time2seg']:
  print(seg)
  from_VTK000.get_vect(alldata.meta['name']+'_E_'        +str(seg),seg,fibo_obj=alldata,tar_var='000E'    ,silent=False)
  from_VTK000.get_scal(alldata.meta['name']+'_rhoe0_'    +str(seg),seg,fibo_obj=alldata,tar_var='000rhoe0')
  from_VTK000.get_scal(alldata.meta['name']+'_rhoi1_'    +str(seg),seg,fibo_obj=alldata,tar_var='000rhoi1')

alldata.calc_units('SW', fibo_obj=alldata, silent=False)
