"""
oxidation
======
source code to generate figures for oxidation paper (Deng et al., 2020)

"""
"""
Pre-requisite
======
pandas
lmfit
numpy
scipy : interp1d, opt
burnman : modified `burnman` with `lf` eos

"""
"""
Code List
======
/src
fig1   : dV, eos
         tool => cal_dV_this_study
         oxidation_5.xlsx
         
fig2   : redox vs. depth of MO
         tool => cal_dV_this_study
         tool2 => X_cal_vector, fo2_cal_vector

fig3   : redox vs. depth of MO base
         tool => cal_dV_this_study
         tool2 => X_cal_vector, fo2_cal_vector


figS1  : dV comparision
         tool : cal_dV_this_study, pt2v_prev
         tool2 : teos
         oxidation_4.xlsx

figS2  : cn
         oxidation_5.xlsx
         
figS3  : goetherm
         geotherm_MgSiO3_melt => geotherm

figS4  : dV along geotherm
         tool => cal_dV_this_study
         geotherm_MgSiO3_melt => geotherm
         
figS5  : compare G models
         tool => Gr_janaf, Gr5, Gr
         
figS6_fit_all_2_JANAF.py : calculate fit 1,2,3,4
                           tool => Gr_janaf
                           
figS6  : plot fit 1,2,3,4

figS7  : fit 1,2,3,4, fo2 vs. depth
         tool => cal_dV_this_study
         tool2 => X_cal_vector, fo2_cal_vector
         
figS8  : compare with exp at high P
         tool => Gr_janaf
         tool2 => unpack_par, W
         
figS9  : comparision of different geotherms
         tool => cal_dV_this_study
         tool2 => X_cal_vector, fo2_cal_vector
         geotherm_MgSiO3_melt => geotherm
         
figS10 : run after fig3

figS11 : effects of spin transition
         use tool.cal_PV and tool.cal_dV_this_study

---------
A recent paper by borisov et al., 2018 presented some more data at 1 bar, which 
we may use to construct better model for activity ratio. Also, we can gest our 
best-fitting model against these new data ( 195 new data).

ref: Ferric/ferrous ratio in silicate melts: a new model for 1 atm data

"""

#### import src folder 
from . import src
from .src import *


