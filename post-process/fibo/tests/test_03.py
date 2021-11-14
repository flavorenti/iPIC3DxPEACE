#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      TEST THREE 
# loads data & plots energy spectra for 
# B, E, ui_irr, ui_sol, ue_irr, ue_sol
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#----------------------------------------
#calculate velocity decompositions

ui_irr_x, ui_irr_y, ui_irr_z, ui_sol_x, ui_sol_y, ui_sol_z = alldata.calc_irr_sol('ui_x_'+t,'ui_y_'+t,'ui_z_'+t)
ue_irr_x, ue_irr_y, ue_irr_z, ue_sol_x, ue_sol_y, ue_sol_z = alldata.calc_irr_sol('ue_x_'+t,'ue_y_'+t,'ue_z_'+t)
print('done the decomposition!')

#----------------------------------------
#create the fibo object to keep transformed stuff ... 

foudata = fb.fibo('fout')
foudata.meta = alldata.meta

foudata.meta['xl'] = 2.*np.pi*nx/alldata.meta['xl']
foudata.meta['yl'] = 2.*np.pi*ny/alldata.meta['yl']
foudata.meta['zl'] = 2.*np.pi*nz/alldata.meta['zl']

foudata.meta['dx'] = alldata.meta['xl']/nx
foudata.meta['dy'] = alldata.meta['yl']/ny
foudata.meta['dz'] = alldata.meta['zl']/nz

#----------------------------------------
#... calculate transforms (no transform is taken for the z component of each field, since we want just to calculate energies in k_xy)

FT_x = alldata.make_spect('B_x_'+t)
FT_y = alldata.make_spect('B_y_'+t)
FT_z = None #alldata.make_spect('B_z_'+t)
foudata.data['FTe_B_'+t] = foudata.calc_scalr(FT_x,FT_x,FT_y,FT_y,FT_z,FT_z)

FT_x = alldata.make_spect('E_x_'+t)
FT_y = alldata.make_spect('E_y_'+t)
FT_z = None #alldata.make_spect('E_z_'+t)
foudata.data['FTe_E_'+t] = foudata.calc_scalr(FT_x,FT_x,FT_y,FT_y,FT_z,FT_z)

FT_x = alldata.make_spect(ui_irr_x)
FT_y = alldata.make_spect(ui_irr_y)
FT_z = None #alldata.make_spect(ui_irr_z)
foudata.data['FTe_ui_irr_'+t] = foudata.calc_scalr(FT_x,FT_x,FT_y,FT_y,FT_z,FT_z)

FT_x = alldata.make_spect(ue_irr_x)
FT_y = alldata.make_spect(ue_irr_y)
FT_z = None #alldata.make_spect(ue_irr_z)
foudata.data['FTe_ue_irr_'+t] = foudata.calc_scalr(FT_x,FT_x,FT_y,FT_y,FT_z,FT_z)

FT_x = alldata.make_spect(ui_sol_x)
FT_y = alldata.make_spect(ui_sol_y)
FT_z = None #alldata.make_spect(ui_sol_z)
foudata.data['FTe_ui_sol_'+t] = foudata.calc_scalr(FT_x,FT_x,FT_y,FT_y,FT_z,FT_z)

FT_x = alldata.make_spect(ue_sol_x)
FT_y = alldata.make_spect(ue_sol_y)
FT_z = None #alldata.make_spect(ue_sol_z)
foudata.data['FTe_ue_sol_'+t] = foudata.calc_scalr(FT_x,FT_x,FT_y,FT_y,FT_z,FT_z)

print('done the transforms!')

old_vars = set(foudata.data.keys())

#----------------------------------------
#calculate shell spectra
lvl_num = 550 #513 #1537 # #number of levels (one more, yeah)
lvl_max = 30. #90. #61.4 / 2. #      #highest level value (see foudata.meta['xl'] if you don't believe me) 


##version sid :/
##calculate average and summation on spherical shells 
#foudata.data['FTea_B_'+t], ww = foudata.comp_shell_pop('FTe_B_'+t,nx//2,ny//2,nz//2,lvl_num,lvl_max)
#foudata.data['FTes_B_'+t], ww = foudata.comp_shell_pop('FTe_B_'+t,nx//2,ny//2,nz//2,lvl_num,lvl_max,density=False)

##version cerri <3
##calculate summation on cilindrical shells (same that before in the case of two-dimensional arrays)

foudata.data['FTes_B_'+t], ww = foudata.comp_shell_pop('FTe_B_'+t,nx//2,ny//2,0,lvl_num,lvl_max,density=False)
foudata.data['FTes_E_'+t], ww = foudata.comp_shell_pop('FTe_E_'+t,nx//2,ny//2,0,lvl_num,lvl_max,density=False)
foudata.data['FTes_ui_irr_'+t], ww = foudata.comp_shell_pop('FTe_ui_irr_'+t,nx//2,ny//2,0,lvl_num,lvl_max,density=False)
foudata.data['FTes_ui_sol_'+t], ww = foudata.comp_shell_pop('FTe_ui_sol_'+t,nx//2,ny//2,0,lvl_num,lvl_max,density=False)
foudata.data['FTes_ue_irr_'+t], ww = foudata.comp_shell_pop('FTe_ue_irr_'+t,nx//2,ny//2,0,lvl_num,lvl_max,density=False)
foudata.data['FTes_ue_sol_'+t], ww = foudata.comp_shell_pop('FTe_ue_sol_'+t,nx//2,ny//2,0,lvl_num,lvl_max,density=False)

print('done the shells!')



#----------------------------------------
#find best fit in some sub-range
fit_rng = [10,80]
x_vals = np.log2(np.linspace(0,lvl_max,lvl_num)[fit_rng[0]:fit_rng[1]])
foudata.data['theo'] = np.minimum(np.power(np.linspace(0,lvl_max,lvl_num),-5./3.),np.power(np.linspace(0,lvl_max,lvl_num),-3.)) * 10**(-3.)

y_vals = np.log2(foudata.data['FTes_B_'+t][fit_rng[0]:fit_rng[1]])
aa, bb = np.poly1d(np.polyfit(x_vals,y_vals,1))  #find coefficients
foudata.data['fit_FTes_B_'+t] = np.power(np.linspace(0,lvl_max,lvl_num),aa)*(2**bb)
foudata.data['fitp_FTes_B_'+t] = [aa, bb]

new_vars = set(foudata.data.keys())


for var in list(new_vars-old_vars): 
  foudata.data[var] = np.reshape(foudata.data[var],(-1,1,1))


#----------------------------------------
#plot all stuff! 
tar_labs = [['t = '+t]]
my_fig = foudata.draw_canvas(tar_labs,[12.,8.,0.8,0.8,0.8,0.8])

foudata.draw_lineplt(my_fig[0],'FTes_B_'+t,range_x=[1,lvl_num-2],args_fig={'color':'k','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'FTes_E_'+t,range_x=[1,lvl_num-2],args_fig={'color':'b','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'FTes_ui_irr_'+t,range_x=[1,lvl_num-2],args_fig={'color':'g','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'FTes_ue_irr_'+t,range_x=[1,lvl_num-2],args_fig={'color':'c','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'FTes_ui_sol_'+t,range_x=[1,lvl_num-2],args_fig={'color':'y','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'FTes_ue_sol_'+t,range_x=[1,lvl_num-2],args_fig={'color':'r','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'fit_FTes_B_'+t,range_x=[1,lvl_num-2],args_fig={'color':'k','linestyle':':','set-xscale':'log','set-yscale':'log'})
foudata.draw_lineplt(my_fig[0],'theo',range_x=[1,lvl_num-2],args_fig={'color':'m','linestyle':':','set-xscale':'log','set-yscale':'log'})

if print_me: plt.savefig(print_address+'_test_03.png',format='png')
if show_me : plt.show()
