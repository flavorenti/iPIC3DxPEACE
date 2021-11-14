#------------------------------------------------------
# tests the basic functionality of from_MFL 
# 
#
# fada 19
#------------------------------------------------------


import sys
sys.path.append('/home/fadanelli/fibo_beta') 

import fibo_beta as fb



#======regulate=parameters=======================================================
data_address = '/work2/manuela/Run_mms_real_tris/' 
fibo_name = 'mmr3' 
nproc = 1024    #number of processors 
seg = '10'		#str  
exit_num = 0	#int, from zero onwards or None if you wanna get a full segment
#================================================================================

#======regulate=parameters=======================================================
data_address = '/work2/fadanelli/_mmsr/'
fibo_name = 'mmsr'
nproc = 256 #number of processors 
seg = '27' 		#str  
exit_num = 0	#int, from zero onwards or None if you wanna get a full segment
#================================================================================







#----create-your-objects-----------
alldata = fb.fibo(fibo_name)
from_MFL = fb.from_MFL(data_address,fibo_name,nproc)


#----load-everything-------------------
from_MFL.get_meta('',silent=False)
alldata.meta = from_MFL.meta

if exit_num != None :#option 1: load one time only
	from_MFL.get_B(seg,exit_num,fibo_obj=alldata,silent=False)
	from_MFL.get_E(seg,exit_num,fibo_obj=alldata,silent=False)
	from_MFL.get_Ui(seg,exit_num,fibo_obj=alldata,silent=False)
	from_MFL.get_Ue(seg,exit_num,fibo_obj=alldata,silent=False)
	times = [from_MFL.segs[seg][exit_num]]

else : #option 2: load a full data  segment
	from_MFL.get_seg_U(seg,fibo_obj=alldata,silent=False)
	from_MFL.get_seg_EB(seg,fibo_obj=alldata,silent=False)
	from_MFL.get_seg_DPJ(seg,fibo_obj=alldata,silent=False)
	times = from_MFL.segs[seg]



