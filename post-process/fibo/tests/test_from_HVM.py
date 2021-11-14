#------------------------------------------------------
# tests the basic functionality of from_HVM 
# 
#
# fada 19
#------------------------------------------------------

import os 
import sys

sys.path.append('/home/s.fadanelli/Libraries') 
sys.path.append('/home/s.fadanelli/Libraries/fibo_gamma') 

import fibo_gamma as fb


#======regulate=parameters=======================================================
#data_address = '/work1/VLASOV_SIMS/PRACE_3D/Run_LF_Dimitri/'
#fibo_name = 'turbLF'
#seg = '06_24-28'
#exit_num = 1  #int, from zero onwards or None if you wanna get a full segment
#================================================================================

#======regulate=parameters=======================================================
data_address = '/work1/VLASOV_SIMS/HVM2D/HVM_2d_e-rec.a' 
fibo_name = 'eMRa' 
seg = '15'
exit_num = 4  #int, from zero onwards or None if you wanna get a full segment
#================================================================================

#======regulate=parameters=======================================================
#data_address = '/work1/VLASOV_SIMS/PRACE_3D/HVM_3d_e-rec_WT/'
#fibo_name = '3DWT'
#seg = '13'
#exit_num = 0  #int, from zero onwards or None if you wanna get a full segment
#================================================================================

#======regulate=parameters=======================================================
#data_address = '/work1/VLASOV_SIMS/HVM3D/AIDA_1/'
#fibo_name = 'simA'
#seg = '09'
#exit_num = 0  #int, from zero onwards or None if you wanna get a full segment
#================================================================================

#======regulate=parameters=======================================================
#data_address = '/work1/VLASOV_SIMS/HVM2D/HVM_3072' 
#fibo_name = 'H3072' 
#seg = '16'
#exit_num = 0  #int, from zero onwards or None if you wanna get a full segment
#================================================================================





#----create-your-objects-----------
alldata = fb.fibo(fibo_name)
from_HVM = fb.from_HVM(data_address)


#----load-everything-you-already-have--
from_HVM.get_meta(seg,silent=False)
alldata.meta = from_HVM.meta

if exit_num != None :#option 1: load one time only
  from_HVM.get_EB(seg,exit_num,fibo_obj=alldata,silent=False)
  from_HVM.get_Ion(seg,exit_num,fibo_obj=alldata,silent=False)
  from_HVM.get_Press(seg,exit_num,fibo_obj=alldata,silent=False)
  from_HVM.get_Q(seg,exit_num,fibo_obj=alldata,silent=False)
  if os.path.isfile(os.path.join(data_address,seg,'Te.bin')):
    from_HVM.get_Te(seg,exit_num,fibo_obj=alldata,silent=False)
    from_HVM.get_Qe(seg,exit_num,fibo_obj=alldata,silent=False)
  times = [from_HVM.segs[seg][exit_num]]

else : #option 2: load a full data  segment
  from_HVM.get_seg_EB(seg,fibo_obj=alldata,silent=False)
  from_HVM.get_seg_Ion(seg,fibo_obj=alldata,silent=False)
  from_HVM.get_seg_Press(seg,fibo_obj=alldata,silent=False)
  from_HVM.get_seg_Q(seg,fibo_obj=alldata,silent=False)
  if os.path.isfile(os.path.join(data_address,seg,'Te.bin')):
    from_HVM.get_seg_Te(seg,fibo_obj=alldata,silent=False)
    from_HVM.get_seg_Qe(seg,fibo_obj=alldata,silent=False)
  times = from_HVM.segs[seg]

#----calculate-everything-else---------
for time_exit in times :
  from_HVM.calc_all(time_exit,fibo_obj=alldata,silent=False)


