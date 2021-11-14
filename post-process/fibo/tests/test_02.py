#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      PLOT TWO 
# plots energy and current densities along a cut that 
# traverses diagonally the zone plotted in plot one
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#----extract-the-diagonal------------------------------------
diadata = fb.fibo('diagonal')

line_pnts = 200
line_dims = np.sqrt(alldata.meta['dx']**2 + alldata.meta['yl']**2)

center = [alldata.meta['nx']//2,alldata.meta['ny']//2,0]
versor = [1.,1.,0.]

for var in alldata.data.keys(): 
  diadata.data[var], diadata.meta = alldata.extract_line(var,line_pnts,line_dims,center,versor)


#----choose-labels-------------------------------------------
tar_labs = [
  ['Kinetic en. dens.','Internal en. dens.'],
  [None,'Electromagnetic en. dens.']]
my_fig = diadata.draw_canvas(tar_labs,'line')

diadata.draw_lineplt(my_fig[0],'Ki_'+t,args_fig={'c':'b'})
diadata.draw_lineplt(my_fig[0],'Ke_'+t,args_fig={'c':'darkorange'})
diadata.draw_lineplt(my_fig[0],'K_'+t ,args_fig={'c':'m'})

diadata.draw_lineplt(my_fig[1],'Ui_'+t,args_fig={'c':'b'})
diadata.draw_lineplt(my_fig[1],'Ue_'+t,args_fig={'c':'darkorange'})
diadata.draw_lineplt(my_fig[1],'U_'+t ,args_fig={'c':'m'})

diadata.draw_lineplt(my_fig[2],'enE_'+t,args_fig={'c':'gold'})
diadata.draw_lineplt(my_fig[2],'enB_'+t,args_fig={'c':'darkblue'})
diadata.draw_lineplt(my_fig[2],'|J|_'+t,args_fig={'c':'k','ls':':'})

if print_me: plt.savefig(print_address+alldata.fibo_name+'_test_plot_two.png',format='png')
if show_me: plt.show()
