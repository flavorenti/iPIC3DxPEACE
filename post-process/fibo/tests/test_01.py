#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      PLOT ONE 
# plots energy densities: kinetic, electric, internal,
# magnetic. superplots electron, ion and barycenter 
# fluid velocities as streamplots, electric and magnetic
# fields as arrow plots and highlights zones of high 
# current density by dotting them
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#-----generate-lists-of-subplots--------------------------
tar_labs = [
  ['Ki_'+t,'Ke_'+t,'K_'+t,'enE_'+t],
  ['Ui_'+t,'Ue_'+t,'U_'+t,'enB_'+t]]
my_fig = alldata.draw_canvas(tar_labs,'cont_v')
print('got the canvas')

alldata.draw_spotted(my_fig[0],'Ki_'+t ,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[1],'Ke_'+t ,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[2],'K_'+t  ,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[3],'enE_'+t,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[4],'Ui_'+t ,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[5],'Ue_'+t ,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[6],'U_'+t  ,range_x=ranx,range_y=rany)
alldata.draw_spotted(my_fig[7],'enB_'+t,range_x=ranx,range_y=rany)
print('got the energies')

args_fig = {'color':'w','linewidth':'field-mag'} #,'density':[8,10]}
alldata.draw_streams(my_fig[0],'ui_x_'+t,'ui_y_'+t,range_x=ranx,range_y=rany,args_fig=args_fig)
alldata.draw_streams(my_fig[1],'ue_x_'+t,'ue_y_'+t,range_x=ranx,range_y=rany,args_fig=args_fig)
alldata.draw_streams(my_fig[2],'u_x_'+t ,'u_y_'+t ,range_x=ranx,range_y=rany,args_fig=args_fig)
alldata.draw_streams(my_fig[4],'ui_x_'+t,'ui_y_'+t,range_x=ranx,range_y=rany,args_fig=args_fig)
alldata.draw_streams(my_fig[5],'ue_x_'+t,'ue_y_'+t,range_x=ranx,range_y=rany,args_fig=args_fig)
alldata.draw_streams(my_fig[6],'u_x_'+t ,'u_y_'+t ,range_x=ranx,range_y=rany,args_fig=args_fig)
print('got the streams')

mymax = np.log10(np.percentile(alldata.data['|J|_'+t][ranx[0]:ranx[1],rany[0]:rany[1],0],99.95))
mymin = np.log10(np.percentile(alldata.data['|J|_'+t][ranx[0]:ranx[1],rany[0]:rany[1],0],0.005)) #max(np.min(alldata.data['|J|_'+t]),1e-6))
mylvc = np.logspace(mymin,mymax,10)[-3:]
args_fig={'colors':['r','r','r'],'levels':mylvc,'hatches':[None,'..','....'],'alpha':0.3}
alldata.draw_spotted(my_fig[3],'|J|_'+t,range_x=ranx,range_y=rany,args_fig=args_fig,args_bar=None)
alldata.draw_spotted(my_fig[7],'|J|_'+t,range_x=ranx,range_y=rany,args_fig=args_fig,args_bar=None)
print('got the current')

args_fig = {'units':'dots','density':[40,30],'color':'w','headwidth':2.,'headlength':5.,'headaxislength':1.,'width':1,'alpha':0.8}
alldata.draw_quivers(my_fig[3],'E_x_'+t,'E_y_'+t,range_x=ranx,range_y=rany,args_fig=args_fig)
alldata.draw_quivers(my_fig[7],'B_x_'+t,'B_y_'+t,range_x=ranx,range_y=rany,args_fig=args_fig)
#alldata.draw_multi_quivers([my_fig[3],my_fig[7]],['E_x_'+t,'B_x_'+t],['E_y_'+t,'B_y_'+t],ranx,rany,0,args_fig)
print('got the arrows')

if print_me: plt.savefig(print_address+alldata.fibo_name+'_test_plot_one.png',format='png')
if show_me: plt.show()
