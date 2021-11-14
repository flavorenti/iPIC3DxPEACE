
import collections 
import cv2
import numpy as np
import os
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.ndimage as ndm
import itertools as itt
import struct
import re
import time

class fibo_draw : 

  #-------------------------------------------------------------
  #--------------------routines-for-plotting--------------------
  #-------------------------------------------------------------
  def draw_canvas(self,
        tar_labs,  
        tar_dims,
        tar_char=12):
    """ 
    Prepares canvas for drawing
    
    Parameters :
      - tar_labs    [list of list of str] matrix of labels
      - tar_dims    [6*[float] OR 'line' OR 'cont_h' OR 'cont_v'] dims of figs and spaces within
      - tar_char    [float] dimension of the ticks and labels
    Returns :
      - tar_plot   [blank figure] 
    
    """ 

    y_fig_num, x_fig_num = np.shape(tar_labs)

    if tar_dims is 'line'   : tar_dims = [7.2,3.0,0.8,0.6,0.6,0.6]
    if tar_dims is 'cont_h' : tar_dims = [5.2,5.2,0.8,0.6,0.6,0.6]
    if tar_dims is 'cont_v' : tar_dims = [6.0,6.0,0.6,0.6,0.6,0.6]

    dim_x = tar_dims[0]*x_fig_num + tar_dims[2]*(x_fig_num-1) + 2*tar_dims[4]
    dim_y = tar_dims[1]*y_fig_num + tar_dims[3]*(y_fig_num-1) + 2*tar_dims[5]

    multi_fig = plt.figure(figsize=(dim_x,dim_y))
    multi_fig.patch.set_facecolor('white')

    for ii in range(y_fig_num):
      for jj in range(x_fig_num):
        if tar_labs[ii][jj] is not None : 
        
          pos_x = (tar_dims[4] + (tar_dims[0] + tar_dims[2]) * jj)/dim_x
          pos_y = (tar_dims[5] + (tar_dims[1] + tar_dims[3]) * (y_fig_num-1-ii) )/dim_y
  
          multi_fig.add_axes([pos_x,pos_y,tar_dims[0]/dim_x,tar_dims[1]/dim_y])
          mpl.pyplot.ylabel(tar_labs[ii][jj])

    tar_plot = multi_fig.get_axes()

    for ax in tar_plot : 
      ax.xaxis.set_tick_params(labelsize=tar_char)
      ax.yaxis.set_tick_params(labelsize=tar_char)

    return tar_plot

  #------------------------------------------------------------  TO BE TESTED!!!!!!!!
  def draw_lineplt(self,
        tar_plot,
        tar_data,
        range_x = None,
        range_y = None,
        range_z = None,
        args_fig = {},
        ticks_xy = [None,None],
        label_xy = [None,None]):
    """ 
    Draws classic plots of one-dimensional variables
    
    Parameters :
      - tar_plot   [(sub)plot] to draw on
      - tar_data   [fibo_data] target variable
      + range_x    [None OR int,int] x range 
      + range_y    [None OR int,int] y range
      + range_z    [None OR int,int] z range
      + args_fig   [line-kwargs]
      + ticks_xy   [2-list with None OR set_tick_params-kwargs]
      + label_xy   [2-list with None OR set_label-kwargs]

    """

    nx,ny,nz = self.meta['nnn']
    
    if range_x == None : range_x = [0,nx]
    if range_y == None : range_y = [0,1]
    if range_z == None : range_z = [0,1]
    
    # check that the cut you extracted is actually 1D
    case_x = range_x[1]!=range_x[0]+1
    case_y = range_y[1]!=range_y[0]+1
    case_z = range_z[1]!=range_z[0]+1
    if case_x + case_y + case_z != 1 : print('two among range_x, range_y, range_z must be [k,k+1]')
    
    # determine the coordinates of selected data
    x_coord, y_coord, z_coord = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=False)
    if case_x : coord = x_coord
    if case_y : coord = y_coord
    if case_z : coord = z_coord
    
    
    if ticks_xy[0] is not None: tar_plot.xaxis.set_tick_parameters(**ticks_xy[0])
    if ticks_xy[1] is not None: tar_plot.yaxis.set_tick_parameters(**ticks_xy[1])
    if label_xy[0] is not None: tar_plot.set_xlabel('',**label_xy[0])
    if label_xy[1] is not None: tar_plot.set_ylabel('',**label_xy[1])
    
    args_fig_complete = {}
    for k in args_fig.keys():
      args_fig_complete[k] = args_fig[k]
    
    # a little extra: if there's the set-xscale or set-yscale in args_fig ...
    if 'set-xscale' in args_fig_complete.keys(): tar_plot.set_xscale(args_fig_complete.pop('set-xscale'))
    if 'set-yscale' in args_fig_complete.keys(): tar_plot.set_yscale(args_fig_complete.pop('set-yscale'))
    
    # fetch the data
    targ, ww = self.extract_range(tar_data,range_x,range_y,range_z)
    targ = targ.flatten()
    
    tar_plot.plot(coord,targ,**args_fig_complete)
    tar_plot.set_xlim(coord[0],coord[-1])


  #------------------------------------------------------------  
  def draw_contour(self,    
        tar_plot, 
        tar_data,
        range_x = None,
        range_y = None,
        range_z = None,
        args_fig = {'colors':'k','alpha':0.6},
        args_bar = {'aspect':30},
        ticks_xyb = [None,None,None], 
        label_xyb = [None,None,None]): #[{'labelsize':9},{'labelrotation':20,'labeltop':True,'top':True},{}]
    """ 
    Draws level contours of the chosen field for a 2D subsection of the box
    
    Parameters :
      - tar_plot   [(sub)plot] in which you are drawing
      - tar_data   [fibo_data] you want to plot
      + range_x    [None OR int,int] x range 
      + range_y    [None OR int,int] y range
      + range_z    [None OR int,int] z range
      + args_fig   [contour-kwargs]
      + args_bar   [None OR colorbar-kwargs]
      + ticks_xyb  [3-list with None OR set_tick_params-kwargs]
      + label_xyb  [3-list with None OR set_label-kwargs]
      
    """  

    nx,ny,nz = np.shape(self.get_data(tar_data))

    if range_x == None : range_x = [0,nx]
    if range_y == None : range_y = [0,ny]
    if range_z == None : range_z = [0,1]
    
    if label_xyb[0] is None : label_xyb[0] = {'text':''}
    if label_xyb[1] is None : label_xyb[1] = {'text':''}
    if label_xyb[2] is None : label_xyb[2] = {'text':''}

    #len_x = self.meta['xl']
    #len_y = self.meta['yl']
    #x_coord = np.linspace(len_x*range_x[0]/nx,len_x*range_x[1]/nx,range_x[1]-range_x[0])
    #y_coord = np.linspace(len_y*range_y[0]/ny,len_y*range_y[1]/ny,range_y[1]-range_y[0])
    #x_coord, y_coord = np.meshgrid(x_coord,y_coord)

    # check that the cut you extracted is actually 2D
    case_x = range_x[1]==range_x[0]+1
    case_y = range_y[1]==range_y[0]+1
    case_z = range_z[1]==range_z[0]+1
    if case_x + case_y + case_z != 1 : print('one among range_x, range_y, range_z must be [k,k+1]')

    # determine the coordinates of selected data
    x_coord, y_coord, z_coord = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=True)
    if case_x : a_coord, b_coord = y_coord[0,:,:], z_coord[0,:,:]
    if case_y : a_coord, b_coord = z_coord[:,0,:], x_coord[:,0,:]
    if case_z : a_coord, b_coord = x_coord[:,:,0], y_coord[:,:,0]
    
    # cut out selected data
    targ, meta = self.extract_range(tar_data,range_x,range_y,range_z)
    if case_x : targ = targ[0,:,:] #np.transpose( targ[:,:,0] )
    if case_y : targ = targ[:,0,:] #np.transpose( targ[:,:,0] )
    if case_z : targ = targ[:,:,0] #np.transpose( targ[:,:,0] )

    #if levlog: norm = mpl.colors.LogNorm(vmin=targ.min(), vmax=targ.max())
    #else : norm = mpl.colors.Normalize(vmin=targ.min(), vmax=targ.max())

    if ticks_xyb[0] is not None: tar_plot.xaxis.set_tick_params(**ticks_xyb[0])
    if ticks_xyb[1] is not None: tar_plot.yaxis.set_tick_params(**ticks_xyb[1])
    tar_plot.set_xlabel('',**label_xyb[0])
    tar_plot.set_ylabel('',**label_xyb[1])

    tar_plot.set_xlim( a_coord[0,0], a_coord[-1,0] )
    tar_plot.set_ylim( b_coord[0,0], b_coord[0,-1] )
     
    cpts = tar_plot.contour(a_coord,b_coord,targ,**args_fig)
    #levels=levels,colors=styles[0],linestyles=styles[1],alpha=alpha)#,norm=norm
    if args_bar is not None: 
      cbar = plt.colorbar(cpts,ax=tar_plot,**args_bar)
      if ticks_xyb[2] is not None: cbar.ax.xaxis.set_tick_params(**ticks_xyb[2])
      if ticks_xyb[2] is not None: cbar.ax.yaxis.set_tick_params(**ticks_xyb[2])
      cbar.set_label('',**label_xyb[2])

  #------------------------------------------------------------  
  def draw_spotted(self,    
        tar_plot,
        tar_data,
        range_x = None,
        range_y = None,
        range_z = None,
        args_fig = {'colors':None,'hatches':[None],'extend':'both','alpha':0.6},
        args_bar = {'aspect':30},
        ticks_xyb = [None,None,None],
        label_xyb = [None,None,None]):        
    """ 
    Draws filled contours of the chosen fields
    
    Parameters :
      - tar_plot    [list of (sub)plot] in you are drawing
      - tar_data    [list of fibo_var] you want to plot
      + range_x     [None OR int,int] x range 
      + range_y     [None OR int,int] y range 
      + range_z     [None OR int,int] z range 
      + args_fig    [contourf-kwargs]
      + args_bar    [None OR colorbar-kwargs]
      + ticks_xyb   [3-list with None OR set_tick_params-kwargs]
      + label_xyb   [3-list with None OR set_label-kwargs]
    
    """

    nx,ny,nz = self.meta['nnn']  
    if range_x == None : range_x = [0,nx]
    if range_y == None : range_y = [0,ny]
    if range_z == None : range_z = [0,1]
    
    if label_xyb[0] is None : label_xyb[0] = {'text':''}
    if label_xyb[1] is None : label_xyb[1] = {'text':''}
    if label_xyb[2] is None : label_xyb[2] = {'text':''}

    #len_x = self.meta['xl']
    #len_y = self.meta['yl']
    #x_coord = np.linspace(len_x*range_x[0]/nx,len_x*range_x[1]/nx,range_x[1]-range_x[0])
    #y_coord = np.linspace(len_y*range_y[0]/ny,len_y*range_y[1]/ny,range_y[1]-range_y[0])
    #x_coord, y_coord = np.meshgrid(x_coord,y_coord)

    # check that the cut you extracted is actually 2D
    case_x = range_x[1]==range_x[0]+1
    case_y = range_y[1]==range_y[0]+1
    case_z = range_z[1]==range_z[0]+1
    if case_x + case_y + case_z != 1 : print('one among range_x, range_y, range_z must be [k,k+1]')
    
    x_coord, y_coord, z_coord = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=True)
    if case_x : a_coord, b_coord = y_coord[0,:,:], z_coord[0,:,:]
    if case_y : a_coord, b_coord = z_coord[:,0,:], x_coord[:,0,:]
    if case_z : a_coord, b_coord = x_coord[:,:,0], y_coord[:,:,0]
    
    targ, meta = self.extract_range(tar_data,range_x,range_y,range_z)
    if case_x : targ = targ[0,:,:] #np.transpose( targ[:,:,0] )
    if case_y : targ = targ[:,0,:] #np.transpose( targ[:,:,0] )
    if case_z : targ = targ[:,:,0] #np.transpose( targ[:,:,0] )

    if ticks_xyb[0] is not None :  tar_plot.xaxis.set_tick_params(**ticks_xyb[0])
    if ticks_xyb[1] is not None :  tar_plot.yaxis.set_tick_params(**ticks_xyb[1])
    tar_plot.set_xlabel('',**label_xyb[0])
    tar_plot.set_ylabel('',**label_xyb[1])

    tar_plot.set_xlim( a_coord[0,0], a_coord[-1,0] )
    tar_plot.set_ylim( b_coord[0,0], b_coord[0,-1] )

    cpts = tar_plot.contourf(a_coord,b_coord,targ,**args_fig)
    if args_bar is not None: 
      cbar = mpl.pyplot.colorbar(cpts,ax=tar_plot,**args_bar)

      if ticks_xyb[2] is not None : cbar.ax.xaxis.set_tick_params(**ticks_xyb[2])
      if ticks_xyb[2] is not None : cbar.ax.yaxis.set_tick_params(**ticks_xyb[2])
      cbar.set_label('',**label_xyb[2])


  #------------------------------------------------------------  
  def draw_quivers(self,    
        tar_plot,  
        tar_var_x,
        tar_var_y, 
        range_x = None,
        range_y = None,
        range_z = None,
        args_fig = {'color':'k','width':1.,'headwidth':1.,'headlength':5.,'headaxislength':1.,'alpha':0.6},
        ticks_xy = [None,None],
        label_xy = [None,None]):
    """
    Draws arrows to represent the chosen fields
    
    Parameters :
      - tar_plot    [(sub)plot] in you are drawing
      - tar_var_x   [fibo_data] you want to plot x
      - tar_var_y   [fibo_data] you want to plot y
      + range_x     [None OR int,int] x range 
      + range_y     [None OR int,int] y range 
      + range_z     [None OR int,int] z range 
      + args_fig    [quiver-kwargs]
      + ticks_xy    [2-list with None OR set_tick_params-kwargs]
      + label_xy    [2-list with None OR set_label-kwargs]
    
    """
    
    nx,ny,nz = self.meta['nnn']
    if range_x == None : range_x = [0,nx]
    if range_y == None : range_y = [0,ny]
    if range_z == None : range_z = [0,1]
    
    if ticks_xy[0] is None : ticks_xy[0] = {'labelsize':9}
    if ticks_xy[1] is None : ticks_xy[1] = {'labelsize':9}
    
    if label_xy[0] is None : label_xy[0] = {'text':''}
    if label_xy[1] is None : label_xy[1] = {'text':''}
    
    # check that the cut you extracted is actually 2D
    case_x = range_x[1]==range_x[0]+1
    case_y = range_y[1]==range_y[0]+1
    case_z = range_z[1]==range_z[0]+1
    if case_x + case_y + case_z != 1 : print('one among range_x, range_y, range_z must be [k,k+1]')
    
    args_fig_complete = {}
    args_fig_complete.update(args_fig)
    
    #box on the reduced resolution data
    if 'density' in args_fig_complete.keys() : density = args_fig_complete.pop('density')
    else : density = [30,30]
    
    #axes for the reduced resolution box
    if case_x : ww, b_mesh, a_mesh = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=True)
    if case_y : a_mesh, ww, b_mesh = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=True)
    if case_z : b_mesh, a_mesh, ww = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=True)
    
    a_mesh = np.transpose(a_mesh[::density[0],::density[1],0])
    b_mesh = np.transpose(b_mesh[::density[0],::density[1],0])
    
    #now reduce data resolution and plot!
    tar_plot.xaxis.set_tick_params(**ticks_xy[0])
    tar_plot.yaxis.set_tick_params(**ticks_xy[1])
    
    tar_plot.set_xlabel('',**label_xy[0])
    tar_plot.set_ylabel('',**label_xy[1])
    
    if case_x : 
      tar_a, ww = self.extract_range(tar_var_y,range_x,range_y,range_z)
      tar_b, ww = self.extract_range(tar_var_z,range_x,range_y,range_z)
    if case_y : 
      tar_a, ww = self.extract_range(tar_var_z,range_x,range_y,range_z)
      tar_b, ww = self.extract_range(tar_var_x,range_x,range_y,range_z)
    if case_z : 
      tar_a, ww = self.extract_range(tar_var_x,range_x,range_y,range_z)
      tar_b, ww = self.extract_range(tar_var_y,range_x,range_y,range_z)
    
    tar_a = np.transpose(tar_a[::density[0],::density[1],0])
    tar_b = np.transpose(tar_b[::density[0],::density[1],0])
    
    cpts = tar_plot.quiver(b_mesh,a_mesh,tar_a,tar_b,**args_fig_complete) 
    

  #------------------------------------------------------------  
  def draw_streams(self,    
        tar_plot,
        tar_data_x,
        tar_data_y,
        range_x = None,
        range_y = None,
        range_z = None,
        args_fig = {'color':'w','linewidth':'field-mag','density':[5,5],'cmap':None},
        args_bar = None,
        ticks_xyb = [None,None,None],
        label_xyb = [None,None,None]):

    """ 
    Draws stream-map of some vector field
    
    Parameters :
      - tar_plot    [sub)plot] in you are drawing
      - tar_data_x  [fibo.data] the x component of field you want to plot
      - tar_data_y  [fibo.data] the y component of field you want to plot
      + range_x     [None OR int,int] x range 
      + range_y     [None OR int,int] y range 
      + range_z     [None OR int,int] z range 
      + args_fig    [contourf-kwargs]
      + args_bar    [None OR colorbar-kwargs]
      + ticks_xyb   [3-list with None OR set_tick_params-kwargs]
      + label_xyb   [3-list with None OR set_label-kwargs]   
    
    """
    
    # here we start
    nx,ny,nz = self.meta['nnn']
    if range_x == None : range_x = [0,nx]
    if range_y == None : range_y = [0,ny]
    if range_z == None : range_z = [0,1]
    
    if ticks_xyb[0] is None : ticks_xyb[0] = {'labelsize':9}
    if ticks_xyb[1] is None : ticks_xyb[1] = {'labelsize':9}
    
    if label_xyb[0] is None : label_xyb[0] = {'text':''}
    if label_xyb[1] is None : label_xyb[1] = {'text':''}
    
    # check that the cut you extracted is actually 2D
    case_x = range_x[1]==range_x[0]+1
    case_y = range_y[1]==range_y[0]+1
    case_z = range_z[1]==range_z[0]+1
    if case_x + case_y + case_z != 1 : print('one among range_x, range_y, range_z must be [k,k+1]')
    
    x_coord, y_coord, z_coord = self.axis_data(range_x,range_y,range_z,coor_like=True,mesh_like=True)
    if case_x : a_coord, b_coord = np.transpose(y_coord[0,:,:]), np.transpose(z_coord[0,:,:])
    if case_y : a_coord, b_coord = np.transpose(z_coord[:,0,:]), np.transpose(x_coord[:,0,:])
    if case_z : a_coord, b_coord = np.transpose(x_coord[:,:,0]), np.transpose(y_coord[:,:,0])
    
    # here it becomes serious 
    tar_plot.xaxis.set_tick_params(**ticks_xyb[0])
    tar_plot.yaxis.set_tick_params(**ticks_xyb[1])
    tar_plot.set_xlabel('',**label_xyb[0])
    tar_plot.set_ylabel('',**label_xyb[1])
    
    tar_plot.set_xlim( a_coord[0,0], a_coord[-1,-1] )
    tar_plot.set_ylim( b_coord[0,0], b_coord[-1,-1] )
    
    tar_a, meta = self.extract_range(tar_data_x,range_x,range_y,range_z)
    tar_b, meta = self.extract_range(tar_data_y,range_x,range_y,range_z)
    if case_x : tar_a, tar_b = np.transpose(tar_a[0,:,:]), np.transpose(tar_b[0,:,:])
    if case_y : tar_a, tar_b = np.transpose(tar_a[:,0,:]), np.transpose(tar_b[:,0,:])
    if case_z : tar_a, tar_b = np.transpose(tar_a[:,:,0]), np.transpose(tar_b[:,:,0])
    
    args_fig_complete = {}
    args_fig_complete.update(args_fig)
    
    # a little extra: if you set as 'field-mag' the 'color' argument 
    # then the color of streams is modulated as field magnitude
    if 'color' in args_fig.keys() and args_fig['color'] == 'field-mag' :
      args_fig_complete['color'] = np.sqrt( self.calc_scalr(tar_a,tar_a,tar_b,tar_b) )
      if not('cmap' in args_fig.keys()) : args_fig_complete['cmap'] = 'autumn'
    
    # a little extra: if you set as 'field-mag' the 'linewidth' argument 
    # then the width of streams is modulated as field magnitude
    if 'linewidth' in args_fig.keys() and args_fig['linewidth'] == 'field-mag' :
      args_fig_complete['linewidth'] = np.sqrt( self.calc_scalr(tar_a,tar_a,tar_b,tar_b) )
      args_fig_complete['linewidth'] = 5. * args_fig_complete['linewidth'] / np.max( args_fig_complete['linewidth'] )
    
    cpts = tar_plot.streamplot(a_coord,b_coord,tar_a,tar_b,**args_fig_complete) 
    
    if args_bar is not None: 
      cbar = mpl.pyplot.colorbar(cpts.lines,ax=tar_plot,**args_bar)
      cbar.ax.xaxis.set_tick_params(**ticks_xyb[2])
      cbar.ax.yaxis.set_tick_params(**ticks_xyb[2])
      cbar.set_label('',**label_xyb[2])

  #------------------------------------------------------------  
  def draw_scatter(self,    
        tar_plot,
        tar_pnts,
        range_x = None,
        range_y = None,
        range_z = None,
        args_fig = {'c':'r','s':8.,'alpha':0.6},
        ticks_xy = [None,None],  #[{'labelsize':9},{'labelsize':9}]
        label_xy = [None,None]):  
    """
    Draws scatterplot of the chosen points
    
    Parameters :
      - tar_plot   [(sub)plot] in which you are drawing
      - tar_pnts   [fibo.pnts] you want to plot
      + range_x    [None OR int,int] x range 
      + range_y    [None OR int,int] y range 
      + range_z    [None OR int,int] z range 
      + args_fig   [scatter-kwargs]
      + ticks_xy   [2-list with None OR set_tick_params-kwargs]
      + label_xy   [2-list with None OR set_label-kwargs]   
    """

    nx,ny,nz = self.meta['nnn']
    dx,dy,dz = self.meta['ddd']

    if range_x == None : range_x = [0,nx]# all_max]
    if range_y == None : range_y = [0,ny] #all_max]
    if range_z == None : range_z = [0,1] #all_max]
   
    # check that the cut you extracted is actually 2D
    case_x = range_x[1]==range_x[0]+1
    case_y = range_y[1]==range_y[0]+1
    case_z = range_z[1]==range_z[0]+1
    if case_x + case_y + case_z != 1 : print('one among range_x, range_y, range_z must be [k,k+1]')
    
    # now start setting the plot specifics
    if ticks_xy[0] is not None : tar_plot.xaxis.set_tick_params(**ticks_xy[0])
    if ticks_xy[1] is not None : tar_plot.yaxis.set_tick_params(**ticks_xy[1])

    if label_xy[0] is None : label_xy[0] = {'text':''}
    if label_xy[1] is None : label_xy[1] = {'text':''}

    tar_plot.set_xlabel('',**label_xy[0])
    tar_plot.set_ylabel('',**label_xy[1])
    
    # select the tar_points inside the plotting box 
    # there might be something really slow in the follofing lines - to check :/
    plt_pnts = self.find_range(tar_pnts,range_x,range_y,range_z)
    plt_pnts = plt_pnts.astype('float') 
    
    # give their coordinates the actual values
    plt_pnts[0,:] =  plt_pnts[0,:]*dx
    plt_pnts[1,:] =  plt_pnts[1,:]*dy
    plt_pnts[2,:] =  plt_pnts[2,:]*dz
    
    # and plot them! 
    if case_x : tar_plot.scatter(plt_pnts[1,:],plt_pnts[2,:],**args_fig) #,c=styles[0],s=styles[1],alpha=alpha) 
    if case_y : tar_plot.scatter(plt_pnts[2,:],plt_pnts[0,:],**args_fig) #,c=styles[0],s=styles[1],alpha=alpha) 
    if case_z : tar_plot.scatter(plt_pnts[0,:],plt_pnts[1,:],**args_fig) #,c=styles[0],s=styles[1],alpha=alpha) 
    
    # finish by specifying box settings
    x_coord, y_coord, z_coord = self.axis_data(range_x,range_y,range_z,coor_like=True)
    if case_x : a_coord, b_coord = y_coord, z_coord  
    if case_y : a_coord, b_coord = z_coord, x_coord  
    if case_z : a_coord, b_coord = x_coord, y_coord  

    tar_plot.set_xlim( a_coord[0], a_coord[-1] )
    tar_plot.set_ylim( b_coord[0], b_coord[-1] )

