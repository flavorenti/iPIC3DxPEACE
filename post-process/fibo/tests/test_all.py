#------------------------------------------------------
# tests fibo functions and plotting routines
#
# we should operate on fibo so that all these  
# work correctly all the time ...
# 
# 
# fada 19
#------------------------------------------------------

#=====================================================#
# plz run the command                                 #
exec(open('test_from_VTK.py').read())
# so to get all necessary variables ... then go       #
#   >>>exec(open('test_all.py').read())               #
#=====================================================#

#------import-modules--------------------------------
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndm
import os



#------choose-ranges--------------------------------------
nx,ny,nz = alldata.meta['nnn']

ranx = [150,550]#
rany = [120,700]#
if alldata.meta['space_dim'] == '3D' : cut_z = 80
if alldata.meta['space_dim'] == '2D' : cut_z = 0

t = str(times[0])

diag_len  = (alldata.meta['dx']*(ranx[1]-ranx[0]))**2
diag_len += (alldata.meta['dy']*(rany[1]-rany[0]))**2
diag_len = np.sqrt(diag_len)



#------create-some-other-variables-----------------------------------
allcalc = fb.phybo(alldata) 
allcalc.calc_Psi(t,cut_z)

alldata.data['|J|_'+t] = np.sqrt(alldata.calc_scalr('J_x_'+t,'J_x_'+t,'J_y_'+t,'J_y_'+t,'J_z_'+t,'J_z_'+t))
alldata.data['|B|_'+t] = np.sqrt(alldata.calc_scalr('B_x_'+t,'B_x_'+t,'B_y_'+t,'B_y_'+t,'B_z_'+t,'B_z_'+t))





#=============set=your=requests=from=the=code==================
plot_set = set(['03'])
show_me = True
print_me = True
print_address = '/home/flavorenti/Bureau/Global-simulations-Mercury'
#==============================================================

my_address = os.getcwd()


if '00' in plot_set: exec(open(my_address+'/test_00.py').read())
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      TEST ZERO 
# uses phybo to calculate magnetic flux function
# then finds coordinates of critical points with one
# method and levels of critical points with another
# plots all in three subplots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if '01' in plot_set: exec(open(my_address+'/test_01.py').read())
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      TEST ONE 
# plots energy densities: kinetic, electric, internal,
# magnetic. superplots electron, ion and barycenter 
# fluid velocities as streamplots, electric and magnetic
# fields as arrow plots and highlights zones of high 
# current density by dotting them
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if '02' in plot_set: exec(open(my_address+'/test_02.py').read())
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      TEST TWO 
# plots energy and current densities along a cut that 
# traverses diagonally the zone plotted in plot one
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if '03' in plot_set: exec(open(my_address+'/test_03.py').read())
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      TEST THREE 
# loads data & plots energy spectra for 
# B, E, ui_irr, ui_sol, ue_irr, ue_sol
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++



