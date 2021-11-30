#! /bin/env python
# coding: utf8

#---------------------------------------------------------------------------------------
#-(1)-----21.12.18----FDP:S,F-----------------------------------------------------------
#-(2)-----30.12.18----FDP:S,J,R---------------------------------------------------------
#-(3)-----02.04.19----FDP:S,P,L---------------------------------------------------------
#-(4)-----04.04.19----FDP:S,J,L---------------------------------------------------------
#---------------------------------------------------------------------------------------
#-(alpha)-19.07.19----FDP:S,J,L---------------------------------------------------------
#-(beta0)-12.11.19----FDP:S,F-----------------------------------------------------------
#-(beta1)-19.11.19----FDP:S-------------------------------------------------------------
#-(beta2)-25.11.19----FDP:S,L-----------------------------------------------------------
#-(beta3)-21.03.20----FDP:S-------------------------------------------------------------
#-(gamma)-03.06.21----FDP:S,J-----------------------------------------------------------
#---------------------------------------------------------------------------------------


from mod_from import  from_VTK
from mod_phybo import phybo


import mod_get
import mod_axis
import mod_extract
import mod_calc
import mod_find
import mod_comp
import mod_draw
import mod_print
import mod_extra


class fibo (mod_get.fibo_get,
            mod_axis.fibo_axis,
            mod_extract.fibo_extract,
            mod_calc.fibo_calc,
            mod_find.fibo_find,
            mod_comp.fibo_comp,
            mod_draw.fibo_draw,
            mod_print.fibo_print,
            mod_extra.fibo_extra):
  """
  fibo is a python object designed to contain your simulation data and perform authomatically simple operations on it 
  all functions are designed for data arrays of the form [nx,ny,nz] (no other indeces, please - this is done to keep routines light and free from for cycles)

  [fibo.data] means [str] in data or [np.ndarray(nx,ny,nz)]
  """
  
  def __init__(self,
      fibo_name):    #name of the set of data on which you work (part of simulation)

    self.fibo_name = str(fibo_name)

    self.data = {}  #dict of  np.ndarray(nx,ny,nz)
    self.pnts = {}  #dict of  np.ndarray(3,points)

    self.meta = {}  #dict of  meta_data - must be copied from the loaders
    self.stat = {}  #dict of  statistical values 

  #------------------------------------------------------------
  def help(self):

    print("Qu'est-ce que c'est que ce bordel ?")

    print('pippo = fb.fibo("pippo")')
    print('dir(fb)')
    print('dir(fb.fibo)')
    print('pippo.data.keys() --> list of data available')
    print('pippo.data["newname"] = pippo.data.pop("oldname")')
    print('np.unravel_index(np.argmax(...),np.shape(...))')
