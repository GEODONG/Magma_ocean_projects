#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:52:43 2019

@author: jiedeng
"""

import numpy as np
import matplotlib.pyplot as plt
from vatic import uplot as uplot
#from scipy.interpolate import interp1d 
import IW_buffer as iw
from geotherm_MgSiO3_melt import geotherm
from tools import read_par,BM4_TH_pt2v_vector_uct,BM4_TH_vt2p_vector_uct,cal_PV,cal_PV_uct
from tools import pt2v,pt2v_this_study,pt2v_prev
import tools as tl
from burnman import minerals
from uncertainties import umath
from melt_eos import MgFeO_pt2v_vector


R = 8.314
A3_to_cm3   = 1e-30*1e6*6.022e23
cm3_to_A3   = 1e30*1e-6/6.022e23

Tm0          = 2500
#Fe_number    = 2
#skiprow_read = 87


P    = np.linspace(1e5,80e9,100)/1e9
Tgeo = geotherm(P*1e9,Tm0=Tm0)
#vx  = np.linspace(587.07*0.90,1350,100)
#Tgeo = np.ones(Tgeo.shape)*4000
V_k18_FeO_l = MgFeO_pt2v_vector(100,P,T=Tgeo,spin=0)[-1]*A3_to_cm3
FeO_k14_l   = minerals.Komabayashi_2014.feo_l()
V_FeO_k14_l = FeO_k14_l.evaluate(['molar_volume'],P*1e9,Tgeo).T*1e6

V_FeO_1 = pt2v_prev(P,Tgeo,'PFeO_1')
V_FeO_2 = pt2v_prev(P,Tgeo,'PFeO_2')
V_FeO1p5 = pt2v_prev(P,Tgeo,'PFeO1p5')
V_FeO_low = pt2v_prev(P,Tgeo,'PFeO_low')
V_FeO_high = pt2v_prev(P,Tgeo,'PFeO_high')
V_FeO1p5_low = pt2v_prev(P,Tgeo,'PFeO1p5_low')
V_FeO1p5_high = pt2v_prev(P,Tgeo,'PFeO1p5_high')

v_feo_re = pt2v_this_study(P,Tgeo,'FeO_re',False)/32*A3_to_cm3
v_feo_ox = pt2v_this_study(P,Tgeo,'FeO_ox',False)/32*A3_to_cm3;
dv_feo   = v_feo_ox - v_feo_re

v_fe_25_re = pt2v_this_study(P,Tgeo,'Fe_25_re',False)/4*A3_to_cm3
v_fe_25_ox = pt2v_this_study(P,Tgeo,'Fe_25_ox',False)/4*A3_to_cm3
dv_fe_25 = v_fe_25_ox - v_fe_25_re

v_fe_12p5_re = pt2v_this_study(P,Tgeo,'Fe_12p5_re',False)/2*A3_to_cm3
v_fe_12p5_ox = pt2v_this_study(P,Tgeo,'Fe_12p5_ox',False)/2*A3_to_cm3
dv_fe_12p5   = v_fe_12p5_ox - v_fe_12p5_re

v_fe_6p25_re = pt2v_this_study(P,Tgeo,'Fe_6p25_re',False)/2*A3_to_cm3
v_fe_6p25_ox = pt2v_this_study(P,Tgeo,'Fe_6p25_ox',False)/2*A3_to_cm3
dv_fe_6p25   = v_fe_6p25_ox - v_fe_6p25_re

plt.figure(figsize=(8,6))
plt.plot(P,V_k18_FeO_l,'*',label='FeO l k18')
plt.plot(P,V_FeO_k14_l,'*',label='FeO l k14')
plt.plot(P,V_FeO_1,'--',label='re dVdT1 Zhang')
plt.plot(P,V_FeO_2,'--',label='re dVdT2 Zhang')
plt.plot(P,V_FeO1p5,label='ox Zhang')
plt.plot(P,V_FeO_low,'--',label='re low Schaefer')
plt.plot(P,V_FeO_high,'--',label='re high Schaefer')
plt.plot(P,V_FeO1p5_low,label='ox low Schaefer')
plt.plot(P,V_FeO1p5_high,label='ox high Schaefer')
plt.plot(P,v_feo_re,'k--',label='FeO this')
plt.plot(P,v_feo_ox,'k-',label='FeO1p5 this')
plt.xlabel('P (GPa)')
plt.ylabel('V (cm^3/mol)')
plt.xlim([0,80])
plt.legend()
plt.grid(True)

plt.figure(figsize=(8,6))
plt.plot(P,V_FeO1p5 - V_FeO_1,'-',label='dVdT1 Zhang')
plt.plot(P,V_FeO1p5 - V_FeO_2,'-',label='dVdT2 Zhang')
plt.plot(P,V_FeO1p5_low - V_FeO_low,'-',label='low Schaefer')
plt.plot(P,V_FeO1p5_high - V_FeO_high,label='high Schaefer')

plt.plot(P,dv_feo,'k-',label = 'FeO this')
plt.plot(P,dv_fe_6p25,'k-.',label = 'Fe_6p25 this')
plt.plot(P,dv_fe_12p5,'k--',label = 'Fe_12p5 this')
plt.plot(P,dv_fe_25,'k:',label = 'Fe_25 this')

plt.xlabel('P (GPa)')
plt.ylabel('dV (cm^3/mol)')
plt.xlim([0,80])
plt.legend()
plt.grid(True)


#
#plt.figure()
#uplot.uplot(P,dV,'r')
#plt.xlabel('P (GPa)')
#plt.ylabel('V (cm^3/mol)')
#plt.legend()
#plt.grid(True)
#plt.xlim([0,80])
#
dV = dv_feo
fg      = []
log10fg = []
fiw = np.zeros(dV.shape)
for i in range(len(dV)):
    PV = cal_PV(dV[range(i+1)],P[range(i+1)],Tgeo[i],P[i],minP=0)
    fg.append(umath.exp(PV*4 + 4*107257/R/Tgeo[i] ))
    log10fg.append(umath.log10(fg[i]))
    fiw[i] = iw.f3(P[i],Tgeo[i])

f_z17_1,logf_z17_1 = tl.fo2_PV(P,Tgeo,V_FeO1p5 - V_FeO_1)
f_z17_2,logf_z17_2 = tl.fo2_PV(P,Tgeo,V_FeO1p5 - V_FeO_2)
f_s19_1,logf_s19_1 = tl.fo2_PV(P,Tgeo,V_FeO1p5_low - V_FeO_low)
f_s19_2,logf_s19_2 = tl.fo2_PV(P,Tgeo,V_FeO1p5_high - V_FeO_high)
f_feo,  logf_feo   = tl.fo2_PV(P,Tgeo,dv_feo)
f_fe_25,logf_fe_25 = tl.fo2_PV(P,Tgeo,dv_fe_25)
f_fe_12p5,logf_fe_12p5 = tl.fo2_PV(P,Tgeo,dv_fe_12p5)
f_fe_6p25,logf_fe_6p25 = tl.fo2_PV(P,Tgeo,dv_fe_6p25)

for i in range(len(Tgeo)):
    fiw[i] = iw.f3(P[i],Tgeo[i])

plt.figure(figsize=(8,6))
plt.plot(P,logf_z17_1,'-',label='dVdT1 Zhang')
plt.plot(P,logf_z17_2,'-',label='dVdT2 Zhang')
plt.plot(P,logf_s19_1,'-',label='low Schaefer')
plt.plot(P,logf_s19_2,label='high Schaefer')

plt.plot(P,logf_feo,'k-',label = 'FeO this')
plt.plot(P,logf_fe_6p25,'k-.',label = 'Fe_6p25 this')
plt.plot(P,logf_fe_12p5,'k--',label = 'Fe_12p5 this')
plt.plot(P,dv_fe_25*A3_to_cm3,'r:',label = 'Fe_25 this')
plt.plot(P, np.log10(fiw)-np.log10(fiw[0]),label='IW')
plt.xlabel('P (GPa)')
plt.ylabel('log10(fo2) - log10(fo2(1 bar))')
plt.xlim([0,80])
plt.legend()
plt.grid(True)


########## PV/R/T + deltaG/R/T, detlaG is a constant taken from Zhang et al., 2017  ##########
f_z17_1,logf_z17_1 = tl.fo2_PV_G(P,Tgeo,V_FeO1p5 - V_FeO_1)
f_z17_2,logf_z17_2 = tl.fo2_PV_G(P,Tgeo,V_FeO1p5 - V_FeO_2)
f_s19_1,logf_s19_1 = tl.fo2_PV_G(P,Tgeo,V_FeO1p5_low - V_FeO_low)
f_s19_2,logf_s19_2 = tl.fo2_PV_G(P,Tgeo,V_FeO1p5_high - V_FeO_high)
f_feo,  logf_feo   = tl.fo2_PV_G(P,Tgeo,dv_feo)
f_fe_25,logf_fe_25 = tl.fo2_PV(P,Tgeo,dv_fe_25)
f_fe_12p5,logf_fe_12p5 = tl.fo2_PV_G(P,Tgeo,dv_fe_12p5)
f_fe_6p25,logf_fe_6p25 = tl.fo2_PV_G(P,Tgeo,dv_fe_6p25)


plt.figure(figsize=(8,6))
plt.plot(P,logf_z17_1,'-',label='dVdT1 Zhang')
plt.plot(P,logf_z17_2,'-',label='dVdT2 Zhang')
plt.plot(P,logf_s19_1,'-',label='low Schaefer')
plt.plot(P,logf_s19_2,label='high Schaefer')

plt.plot(P,logf_feo,'k-',label = 'FeO this')
plt.plot(P,logf_fe_6p25,'k-.',label = 'Fe_6p25 this')
plt.plot(P,logf_fe_12p5,'k--',label = 'Fe_12p5 this')
plt.plot(P,dv_fe_25*A3_to_cm3,'r:',label = 'Fe_25 this')
plt.plot(P, np.log10(fiw)-np.log10(fiw[0])+logf_z17_1[0],label='IW')
plt.xlabel('P (GPa)')
plt.ylabel('log10(fo2) - log10(fo2(1 bar))')
plt.xlim([0,25])
plt.ylim([-10,2.5])
plt.legend()
plt.grid(True)

########## PV/R/T + deltaG/R/T, detlaG is a constant taken from Zhang et al., 2017  ##########

f_z17_1,logf_z17_1 = tl.fo2_PV_Gr(P,Tgeo,V_FeO1p5 - V_FeO_1)
f_z17_2,logf_z17_2 = tl.fo2_PV_Gr(P,Tgeo,V_FeO1p5 - V_FeO_2)
f_s19_1,logf_s19_1 = tl.fo2_PV_Gr(P,Tgeo,V_FeO1p5_low - V_FeO_low)
f_s19_2,logf_s19_2 = tl.fo2_PV_Gr(P,Tgeo,V_FeO1p5_high - V_FeO_high)
f_feo,  logf_feo   = tl.fo2_PV_Gr(P,Tgeo,dv_feo)
f_fe_25,logf_fe_25 = tl.fo2_PV(P,Tgeo,dv_fe_25)
f_fe_12p5,logf_fe_12p5 = tl.fo2_PV_Gr(P,Tgeo,dv_fe_12p5)
f_fe_6p25,logf_fe_6p25 = tl.fo2_PV_Gr(P,Tgeo,dv_fe_6p25)


plt.figure(figsize=(8,6))
plt.plot(P,logf_z17_1,'-',label='dVdT1 Zhang')
plt.plot(P,logf_z17_2,'-',label='dVdT2 Zhang')
plt.plot(P,logf_s19_1,'-',label='low Schaefer')
plt.plot(P,logf_s19_2,label='high Schaefer')

plt.plot(P,logf_feo,'k-',label = 'FeO this')
plt.plot(P,logf_fe_6p25,'k-.',label = 'Fe_6p25 this')
plt.plot(P,logf_fe_12p5,'k--',label = 'Fe_12p5 this')
#plt.plot(P,dv_fe_25*A3_to_cm3,'r:',label = 'Fe_25 this')
plt.plot(P, np.log10(fiw)-np.log10(fiw[0])+logf_z17_1[0],label='IW')
plt.xlabel('P (GPa)')
plt.ylabel('log10(fo2) - log10(fo2(1 bar))')
plt.xlim([0,25])
plt.ylim([0,10])
plt.legend()
plt.grid(True)

########## PV/R/T + deltaG/R/T, detlaG is a constant taken from Zhang et al., 2017  ##########

f_z17_1,logf_z17_1 = tl.fo2_PV_Gr4(P,Tgeo,V_FeO1p5 - V_FeO_1)
f_z17_2,logf_z17_2 = tl.fo2_PV_Gr4(P,Tgeo,V_FeO1p5 - V_FeO_2)
f_s19_1,logf_s19_1 = tl.fo2_PV_Gr4(P,Tgeo,V_FeO1p5_low - V_FeO_low)
f_s19_2,logf_s19_2 = tl.fo2_PV_Gr4(P,Tgeo,V_FeO1p5_high - V_FeO_high)
f_feo,  logf_feo   = tl.fo2_PV_Gr4(P,Tgeo,dv_feo)
f_fe_25,logf_fe_25 = tl.fo2_PV_Gr4(P,Tgeo,dv_fe_25)
f_fe_12p5,logf_fe_12p5 = tl.fo2_PV_Gr4(P,Tgeo,dv_fe_12p5)
f_fe_6p25,logf_fe_6p25 = tl.fo2_PV_Gr4(P,Tgeo,dv_fe_6p25)


plt.figure(figsize=(8,6))
plt.plot(P,logf_z17_1,'-',label='dVdT1 Zhang')
plt.plot(P,logf_z17_2,'-',label='dVdT2 Zhang')
plt.plot(P,logf_s19_1,'-',label='low Schaefer')
plt.plot(P,logf_s19_2,label='high Schaefer')

plt.plot(P,logf_feo,'k-',label = 'FeO this')
plt.plot(P,logf_fe_6p25,'k-.',label = 'Fe_6p25 this')
plt.plot(P,logf_fe_12p5,'k--',label = 'Fe_12p5 this')
plt.plot(P,logf_fe_25,'r:',label = 'Fe_25 this')
plt.plot(P, np.log10(fiw)-np.log10(fiw[0])+logf_z17_1[0],label='IW')
plt.xlabel('P (GPa)')
plt.ylabel('log10(fo2) - log10(fo2(1 bar))')
plt.xlim([0,25])
#plt.ylim([0,10])
plt.legend()
plt.grid(True)

f_z17_1,logf_z17_1 = tl.fo2_PV_Gr4_w(P,Tgeo,V_FeO1p5 - V_FeO_1)
f_z17_2,logf_z17_2 = tl.fo2_PV_Gr4_w(P,Tgeo,V_FeO1p5 - V_FeO_2)
f_s19_1,logf_s19_1 = tl.fo2_PV_Gr4_w(P,Tgeo,V_FeO1p5_low - V_FeO_low)
f_s19_2,logf_s19_2 = tl.fo2_PV_Gr4_w(P,Tgeo,V_FeO1p5_high - V_FeO_high)
f_feo,  logf_feo   = tl.fo2_PV_Gr4_w(P,Tgeo,dv_feo)
f_fe_25,logf_fe_25 = tl.fo2_PV_Gr4_w(P,Tgeo,dv_fe_25)
f_fe_12p5,logf_fe_12p5 = tl.fo2_PV_Gr4_w(P,Tgeo,dv_fe_12p5)
f_fe_6p25,logf_fe_6p25 = tl.fo2_PV_Gr4_w(P,Tgeo,dv_fe_6p25)


plt.figure(figsize=(8,6))
plt.plot(P,logf_z17_1,'-',label='dVdT1 Zhang')
plt.plot(P,logf_z17_2,'-',label='dVdT2 Zhang')
plt.plot(P,logf_s19_1,'-',label='low Schaefer')
plt.plot(P,logf_s19_2,label='high Schaefer')

plt.plot(P,logf_feo,'k-',label = 'FeO this')
plt.plot(P,logf_fe_6p25,'k-.',label = 'Fe_6p25 this')
plt.plot(P,logf_fe_12p5,'k--',label = 'Fe_12p5 this')
plt.plot(P,logf_fe_25,'r:',label = 'Fe_25 this')
plt.plot(P, np.log10(fiw)-np.log10(fiw[0])+logf_z17_1[0],label='IW')
plt.xlabel('P (GPa)')
plt.ylabel('log10(fo2) - log10(fo2(1 bar))')
plt.xlim([0,25])
#plt.ylim([0,10])
plt.legend()
plt.grid(True)
#plt.figure()
#uplot.uplot(P,log10fg,label='magma')
##plt.plot(P,np.log10(fg2),label='magma adjusted')
#plt.plot(P,np.log10(fiw),label='IW')
##plt.plot(P,np.log10(fg2) - np.log10(fiw),label=r'$\Delta IW$')

#def adjust_base(P,fg,fiw,Pbase,dIWbase):
#    fiw_int = interp1d(P,fiw)
#    fg_int  = interp1d(P,fg)
#    fgbase_real  = fiw_int(Pbase)*10**(dIWbase)
#    fg_base_prev = fg_int(Pbase)# base fo2 before adjust
##    print('fg_base is', fg_base)
#    adjust  = fgbase_real/fg_base_prev
##    print('ajdust is', adjust)
#    fg_real = fg*adjust
#    return fg_real
#
#

#
#plt.legend()
#plt.grid(True)
#plt.xlim([0,80])
##plt.xlim([20,25])
##plt.ylim([10,15])
#
#tl.Gr(1682)/R/1682