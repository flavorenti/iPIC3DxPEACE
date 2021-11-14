#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      TEST ZERO 
# uses phybo to calculate magnetic flux function
# then finds coordinates of critical points with one
# method and levels of critical points with another.
# plots all in three subplots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#------find-coordinates-of-critical-points-in-Psi---------- 
alldata.pnts['Oa_A'], alldata.pnts['Ob_A'], alldata.pnts['X_A'] = alldata.find_critical_delaunay_2d('Psi_'+t)
print('calculated first time the O and X points')
alldata.pnts['Oa_B'], alldata.pnts['Ob_B'], alldata.pnts['X_B'] = alldata.find_critical_delaunay_2d('Psi_'+t,aa=1,bb=1)
print('calculated second time the O and X points')

coord_Oa = alldata.find_common('Oa_A','Oa_B')
coord_Ob = alldata.find_common('Ob_A','Ob_B')
coord_X = alldata.find_common('X_A','X_B')
print('calculated coomon coordinates')

#lvl_Oa = np.unique(alldata.data['Psi_'+t][coord_Oa[0],coord_Oa[1],0]) 
#lvl_Ob = np.unique(alldata.data['Psi_'+t][coord_Ob[0],coord_Ob[1],0])
#lvl_X  = np.unique(alldata.data['Psi_'+t][coord_X[0],coord_X[1],0])

##coord_X[0]  = coord_X[0]  * alldata.meta['dx']
##coord_Oa[0] = coord_Oa[0] * alldata.meta['dx'] 
##coord_Ob[0] = coord_Ob[0] * alldata.meta['dx'] 
##coord_X[1]  = coord_X[1]  * alldata.meta['dy'] 
##coord_Oa[1] = coord_Oa[1] * alldata.meta['dy'] 
##coord_Ob[1] = coord_Ob[1] * alldata.meta['dy']

#------find-interesting-values-for-levels--------------------
# (yes, I am doing this the difficult way, basically to show off)

mymax = np.max(alldata.data['Psi_'+t])
mymin = np.min(alldata.data['Psi_'+t])
mystep = (mymax-mymin)/200.
mylvl = np.arange(mymin+mystep/2.,mymax+mystep/2.,mystep) #mymax * np.exp(np.arange(-30,0))

reg_num_above, reg_num_below = alldata.comp_region_number('Psi_'+t,mylvl)
#print(reg_num_above) 
#print(reg_num_below)
#print(mylvl)
#cs = plt.contourf(np.transpose(Psi),levels=mylvl,alpha=0.8) #careful: the transpose is necessary for having correct contour-plots
#plt.show()

#find values of Psi at X and O points
change_above = np.ediff1d(reg_num_above)
change_below = np.ediff1d(reg_num_below)

lvl_X = np.argwhere(np.logical_or(change_above > 0.,change_below < 0.))
lvl_Oa = np.argwhere(np.logical_and(change_above < 0., change_below == 0.))
lvl_Ob = np.argwhere(np.logical_and(change_below > 0., change_above == 0.))

lvl_X = np.reshape(lvl_X,(np.shape(lvl_X)[0]))
lvl_Oa = np.reshape(lvl_Oa,(np.shape(lvl_Oa)[0]))
lvl_Ob = np.reshape(lvl_Ob,(np.shape(lvl_Ob)[0]))

#change_all = np.array([change_above,change_below])
#change_pos = np.nonzero(change_all)
#goodlvl, repeats = np.unique(change_pos[1],return_counts=True)

#lvl_X = goodlvl[np.nonzero(repeats-1)[0]]
#lvl_O = goodlvl[np.nonzero(repeats-2)[0]]

avlvl = ndm.uniform_filter1d(mylvl, size=2)

lvl_X = avlvl[lvl_X+1]
lvl_Oa = avlvl[lvl_Oa]
lvl_Ob = avlvl[lvl_Ob+1]
print('calculated levels of O and X points')

#----choose-variables-to-be-plotted--------------------------
tar_labs = [['X points','Oa points','Ob points']]
my_fig = alldata.draw_canvas(tar_labs,'cont_h')

alldata.draw_spotted(my_fig[0],'Psi_'+t,args_fig={'levels':lvl_X})
alldata.draw_spotted(my_fig[1],'Psi_'+t,args_fig={'levels':lvl_X})
alldata.draw_spotted(my_fig[2],'Psi_'+t,args_fig={'levels':lvl_X})

alldata.draw_scatter(my_fig[0],coord_X)
alldata.draw_scatter(my_fig[1],coord_Oa)
alldata.draw_scatter(my_fig[2],coord_Ob)

if print_me: plt.savefig(print_address+alldata.fibo_name+'_test_plot_zero.png',format='png')
if show_me: plt.show()




