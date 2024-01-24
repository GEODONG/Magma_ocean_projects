#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:39:34 2019

Designer: Jie Deng
Modified by: Donghao Zheng
"""

import src.IW_buffer as iw
from src.geotherm_MgSiO3_melt import geotherm
import src.tools as tl
import src.tools2 as tl2
import numpy as np
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt
import pandas as pd

def extrapolate(P,material,pivot = 25):
    f_fo2 = interp1d(P,material)
    f_iw  = interp1d(P,np.log10(fiw))
    Pm = np.linspace(1e-4,pivot,80)
    return Pm,f_fo2(Pm) - f_iw(Pm),f_fo2(Pm)

R = 8.314
A3_to_cm3   = 1e-30*1e6*6.022e23
cm3_to_A3   = 1e30*1e-6/6.022e23

Pearth = 55
# should be 55
Pmars = 14
Pmoon = 5


Tm0  = 2100
P    = np.linspace(1e5,60e9,100)/1e9  ## should have more than 50 bins
Tgeo = geotherm(P*1e9,Tm0=Tm0) 
#Tgeo = geotherm(55*1e9,Tm0=Tm0)

_, _,_, _, dv_fe_25   = tl.cal_dV_this_study(P,Tgeo,name='Fe_25',flag=False)
_, _,_, _, dv_fe_12p5 = tl.cal_dV_this_study(P,Tgeo,name='Fe_12p5',flag=False)

#fiw = iw.f3_vector(P,Tgeo)
fiw = iw.f2_vector(P,Tgeo)
fiw_cold = fiw

# Earth: IW-2, 25 GPa
# Mars: IW-1.5, 15 GPa
# Moon: IW-2, 5 GPa
# Mercury: IW-3,5 GPa


## for Fe 25% with new G
x_fe_25_cold_with_new_G = tl2.X_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,flag=False,fit='fit5',method='earth', PV_cal='12p5_cold')
f_25_cold_with_new_G    = interp1d(P,x_fe_25_cold_with_new_G);
r_25_cold_new_G    = f_25_cold_with_new_G(Pearth)
logf_25_cold_new_G = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_25_cold_new_G,flag=False,fit='fit5',method='earth', PV_cal='12p5_cold')
p_earth_cold_new_G,f_earth_cold_new_G,fo_earth_cold_new_G = extrapolate(P,logf_25_cold_new_G,Pearth);

x_fe_25_cold_with_new_G_moon = tl2.X_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,flag=False,fit='fit5',method='moon', PV_cal='12p5_cold')
f_25_cold_with_new_G_moon    = interp1d(P,x_fe_25_cold_with_new_G_moon);
r_25_cold_new_G_moon    = f_25_cold_with_new_G_moon(5)
logf_25_cold_new_G_moon = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_25_cold_new_G_moon,flag=False,fit='fit5',method='moon', PV_cal='12p5_cold')
p_cold_new_G_moon,f_cold_new_G_moon,fo_cold_new_G_moon = extrapolate(P,logf_25_cold_new_G_moon,5);


x_fe_25_cold_with_new_G_mars = tl2.X_cal_vector(P,Tgeo,dv_fe_25,np.log10(fiw)-1.5,flag=False,fit='fit5',method='mars', PV_cal='25_cold')
f_25_cold_with_new_G_mars    = interp1d(P,x_fe_25_cold_with_new_G_mars);
r_25_cold_new_G_mars    = f_25_cold_with_new_G_mars(14)
logf_25_cold_new_G_mars = tl2.fo2_cal_vector(P,Tgeo,dv_fe_25,r = r_25_cold_new_G_mars,flag=False,fit='fit5',method='mars', PV_cal='25_cold')
p_cold_new_G_mars,f_cold_new_G_mars,fo_cold_new_G_mars = extrapolate(P,logf_25_cold_new_G_mars,14);

## for Fe 25%
x_fe_25_cold = tl2.X_cal_vector(P,Tgeo,dv_fe_25,np.log10(fiw)-1.5,flag=False,method='mars', PV_cal='25_cold')
f_25_cold    = interp1d(P,x_fe_25_cold);
r_25_cold    = f_25_cold(14)
logf_25_cold = tl2.fo2_cal_vector(P,Tgeo,dv_fe_25,r = r_25_cold,flag=False,method='mars', PV_cal='25_cold')

# for Fe 12.5%
x_12p5_earth_cold = tl2.X_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,flag=False,method='earth', PV_cal='12p5_cold')
x_12p5_moon_cold  = tl2.X_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,flag=False,method='moon',  PV_cal='12p5_cold')

f_12p5_earth_cold = interp1d(P,x_12p5_earth_cold);r_12p5_earth_cold = f_12p5_earth_cold(Pearth);
f_12p5_moon_cold  = interp1d(P,x_12p5_moon_cold);r_12p5_moon_cold = f_12p5_moon_cold(5);

logf_12p5_earth_cold = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_12p5_earth_cold,flag=False,method='earth', PV_cal='12p5_cold')
logf_12p5_moon_cold  = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_12p5_moon_cold,flag=False,method='moon', PV_cal='12p5_cold')

p_mars_cold,f_mars_cold,fo_mars_cold   = extrapolate(P,logf_25_cold,14);
p_earth_cold,f_earth_cold,fo_earth_cold = extrapolate(P,logf_12p5_earth_cold,Pearth);
p_moon_cold,f_moon_cold,fo_moon_cold   = extrapolate(P,logf_12p5_moon_cold,5);

### Z17
x_Z17_fe_25_cold = tl2.X_Z17(P,Tgeo,np.log10(fiw)-1.5)
f_Z17_25_cold    = interp1d(P,x_Z17_fe_25_cold);
r_Z17_25_cold    = f_Z17_25_cold(14)
logf_Z17_25_cold = tl2.fo2_Z17(P,Tgeo,r = r_Z17_25_cold)

x_Z17_12p5_cold = tl2.X_Z17(P,Tgeo,np.log10(fiw)-2)

f_Z17_12p5_cold = interp1d(P,x_Z17_12p5_cold);
r_Z17_12p5_earth_cold = f_Z17_12p5_cold(25);
r_Z17_12p5_moon_cold = f_Z17_12p5_cold(5);

logf_Z17_12p5_earth_cold = tl2.fo2_Z17(P,Tgeo,r = r_Z17_12p5_earth_cold)
logf_Z17_12p5_moon_cold = tl2.fo2_Z17(P,Tgeo,r = r_Z17_12p5_moon_cold)

p_Z17_mars_cold,f_Z17_mars_cold,fo_Z17_mars_cold   = extrapolate(P,logf_Z17_25_cold,14);
p_Z17_earth_cold,f_Z17_earth_cold,fo_Z17_earth_cold = extrapolate(P,logf_Z17_12p5_earth_cold,25);
p_Z17_moon_cold,f_Z17_moon_cold,fo_Z17_moon_cold   = extrapolate(P,logf_Z17_12p5_moon_cold,5);

### O16
x_O06_fe_25_cold = tl2.X_O06(P,Tgeo,np.log10(fiw)-1.5,method='mars')
f_O06_25_cold    = interp1d(P,x_O06_fe_25_cold);
r_O06_25_cold    = f_O06_25_cold(14)
logf_O06_25_cold = tl2.fo2_O06(P,Tgeo,r = r_O06_25_cold,method='mars')

x_O06_12p5_earth_cold = tl2.X_O06(P,Tgeo,np.log10(fiw)-2,method='earth')
x_O06_12p5_moon_cold = tl2.X_O06(P,Tgeo,np.log10(fiw)-2,method='moon')

f_O06_12p5_earth_cold = interp1d(P,x_O06_12p5_earth_cold);
f_O06_12p5_moon_cold  = interp1d(P,x_O06_12p5_moon_cold);

r_O06_12p5_earth_cold = f_O06_12p5_earth_cold(Pearth);
r_O06_12p5_moon_cold  = f_O06_12p5_moon_cold(5);

logf_O06_12p5_earth_cold = tl2.fo2_O06(P,Tgeo,r = r_O06_12p5_earth_cold,method='earth')
logf_O06_12p5_moon_cold = tl2.fo2_O06(P,Tgeo,r = r_O06_12p5_moon_cold,method='moon')

p_O06_mars_cold,f_O06_mars_cold,fo_O06_mars_cold   = extrapolate(P,logf_O06_25_cold,14);
p_O06_earth_cold,f_O06_earth_cold,fo_O06_earth_cold = extrapolate(P,logf_O06_12p5_earth_cold,Pearth);
p_O06_moon_cold,f_O06_moon_cold,fo_O06_moon_cold   = extrapolate(P,logf_O06_12p5_moon_cold,5);

### A18
x_A18_fe_25_cold = tl2.X_A18(P,Tgeo,np.log10(fiw)-1.5,method='mars')
f_A18_25_cold    = interp1d(P,x_A18_fe_25_cold);
r_A18_25_cold    = f_A18_25_cold(14)
logf_A18_25_cold = tl2.fo2_A18(P,Tgeo,r = r_A18_25_cold,method='mars')

x_A18_12p5_earth_cold = tl2.X_A18(P,Tgeo,np.log10(fiw)-2,method='earth')
x_A18_12p5_moon_cold  = tl2.X_A18(P,Tgeo,np.log10(fiw)-2,method='moon')

f_A18_12p5_earth_cold = interp1d(P,x_A18_12p5_earth_cold);
f_A18_12p5_moon_cold  = interp1d(P,x_A18_12p5_moon_cold);

r_A18_12p5_earth_cold = f_A18_12p5_earth_cold(25);
r_A18_12p5_moon_cold  = f_A18_12p5_moon_cold(5);

logf_A18_12p5_earth_cold = tl2.fo2_A18(P,Tgeo,r = r_A18_12p5_earth_cold,method='earth')
logf_A18_12p5_moon_cold = tl2.fo2_A18(P,Tgeo,r = r_A18_12p5_moon_cold,method='moon')

p_A18_mars_cold,f_A18_mars_cold,fo_A18_mars_cold   = extrapolate(P,logf_A18_25_cold,14);
p_A18_earth_cold,f_A18_earth_cold,fo_A18_earth_cold = extrapolate(P,logf_A18_12p5_earth_cold,25);
p_A18_moon_cold,f_A18_moon_cold,fo_A18_moon_cold   = extrapolate(P,logf_A18_12p5_moon_cold,5);


## Hirschmann 2022
x_fe_25_cold_H22 = tl2.X_H22_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,method='earth', PV_cal='12p5_cold')
f_25_cold_H22    = interp1d(P,x_fe_25_cold_H22);
r_25_cold_H22   = f_25_cold_H22(Pearth)
logf_25_cold_H22 = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_25_cold_H22,method='earth', PV_cal='12p5_cold')
p_earth_cold_H22,f_earth_cold_H22,fo_earth_cold_H22= extrapolate(P,logf_25_cold_H22,Pearth);

x_fe_25_mars_H22 = tl2.X_H22_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-1.5,method='mars', PV_cal='12p5_cold')
f_25_mars_H22    = interp1d(P,x_fe_25_mars_H22);
r_25_mars_H22   = f_25_mars_H22(Pmars)
logf_25_mars_H22 = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_25_cold_H22,method='mars', PV_cal='12p5_cold')
p_mars_cold_H22,f_mars_cold_H22,fo_mars_cold_H22= extrapolate(P,logf_25_mars_H22,Pmars);


x_fe_25_moon_H22 = tl2.X_H22_cal_vector(P,Tgeo,dv_fe_12p5,np.log10(fiw)-2,method='moon', PV_cal='12p5_cold')
f_25_moon_H22    = interp1d(P,x_fe_25_moon_H22);
r_25_moon_H22   = f_25_moon_H22(Pmoon)
logf_25_moon_H22 = tl2.fo2_cal_vector(P,Tgeo,dv_fe_12p5,r = r_25_moon_H22,method='moon', PV_cal='12p5_cold')
p_moon_cold_H22,f_moon_cold_H22,fo_moon_cold_H22= extrapolate(P,logf_25_moon_H22,Pmoon);



fig = plt.figure(figsize=(14,7))
ax = fig.add_subplot(1,3,2)
ax2 = fig.add_subplot(1,3,3)
ax3 = fig.add_subplot(1,3,1)


import plot_tool
plot_tool.load_default_setting()
ax.plot(f_mars_cold,p_mars_cold,'black',linestyle = '-',linewidth = 1,label='')
ax.plot(f_A18_mars_cold,p_A18_mars_cold,'C0-.',alpha=0.7,label='')
ax.plot(f_Z17_mars_cold,p_Z17_mars_cold,'C1',linestyle= '--',alpha=0.7,label='')
ax.plot(f_O06_mars_cold,p_O06_mars_cold,'c:',alpha=0.7,label='')
ax.plot(f_cold_new_G_mars,p_cold_new_G_mars,'r-',linewidth = 6,label='',alpha = 0.8)


ax2.plot(f_earth_cold,p_earth_cold,'black',linestyle = '-',linewidth = 1,label='')
ax2.plot(f_A18_earth_cold,p_A18_earth_cold,'C0-.',alpha=0.7,label='')
ax2.plot(f_Z17_earth_cold,p_Z17_earth_cold,'C1--',alpha=0.7,label='')
ax2.plot(f_O06_earth_cold,p_O06_earth_cold,'c:',alpha=0.7,label='')
ax2.plot(f_earth_cold_new_G,p_earth_cold_new_G,'r-',linewidth = 6,label='',alpha = 0.8)



ax3.plot(f_moon_cold,p_moon_cold,'black',linestyle =  '-',linewidth = 1,label='')
ax3.plot(f_A18_moon_cold,p_A18_moon_cold,'C0-.',alpha=0.7,label='')
ax3.plot(f_Z17_moon_cold,p_Z17_moon_cold,'C1--',alpha=0.7,label='')
ax3.plot(f_O06_moon_cold,p_O06_moon_cold,'c:',alpha=0.7,label='')
ax3.plot(f_cold_new_G_moon,p_cold_new_G_moon,'r-',linewidth = 6,label='',alpha = 0.8)


ax.plot([0,1],[-1,-2],'k-',linewidth = 5)
ax.plot([0,1],[-1,-2],'k-',linewidth = 1.5)
ax.plot([0,1],[-1,-2],'k-.',alpha=0.7)
ax.plot([0,1],[-1,-2],'k--',alpha=0.7)
ax.plot([0,1],[-1,-2],'k:',alpha=0.7)


ax2.plot([0,1],[-1,-2],'C1--',alpha=0.7,label='Z17')#'Zhang et al., 2017')
ax2.plot([0,1],[-1,-2],'c:',alpha=0.7,label='O06')#'O\'Neill et al., 2006')
ax2.plot([0,1],[-1,-2],'C0-.',alpha=0.7,label='A19')#'Armstrong et al., 2018')
# Load the dataset from the Excel file.
data = pd.read_excel("k23.xlsx",sheet_name="fO2_MO_adiabat")

# Extract the columns
X_02 = data["  Fe3+/ΣFe=0.2"].values.reshape(-1, 1)
y = data["P (GPa)"].values

# Plot the data and the fitted curve
ax2.plot(X_02[:-1], y[:-1], color = "C5",linestyle = "--", label='K23')
ax2.plot([0,1],[-1,-2],'k-',linewidth = 1.5,label='D20')#'Deng et al., 2020')
ax2.plot([0,1],[-1,-2],'r-',linewidth = 5,label='This study')


ax3.plot([0,1],[-1,-2],'k-',linewidth = 5)
ax3.plot([0,1],[-1,-2],'k-',linewidth = 1.5)
ax3.plot([0,1],[-1,-2],'k-.',alpha=0.7)
ax3.plot([0,1],[-1,-2],'k--',alpha=0.7)
ax3.plot([0,1],[-1,-2],'k:',alpha=0.7)


ax3.set_ylabel('Pressure (GPa)',fontsize=22, fontname = "Arial")
ax.set_ylim([0,Pmars])
ax.set_xlim([-5,3])
ax.invert_yaxis()
yticks1 = np.linspace(0,Pmars,6)
yticks1 = np.append(yticks1,[Pmars])
ax.set_yticks(yticks1)
ax.set_yticklabels([str(int(i)) for i in yticks1],fontsize=15, fontname = "Arial")
ax.text(0.6,0.95*Pmars,'(b) Mars',fontsize=16,fontname= "Arial")



ax2.set_ylim([0,Pearth])
ax2.set_xlim([-6,6])
ax2.invert_yaxis()
yticks = np.linspace(0,50,6)
yticks = np.append(yticks,[55])
ax2.set_yticks(yticks)
ax2.set_yticklabels([str(int(i)) for i in yticks],fontsize=15, fontname = "Arial")
ax2.text(2.4,0.95*Pearth,'(c) Earth',fontsize=16,fontname= "Arial")

ax3.set_ylim([0,Pmoon])
ax3.set_xlim([-4,0])
ax3.invert_yaxis()
yticks3 = np.linspace(0,Pmoon,6)
yticks3 = np.append(yticks3,[Pmoon])
ax3.set_yticks(yticks3)
ax3.set_yticklabels([str(int(i)) for i in yticks3],fontsize=15, fontname = "Arial")
ax3.text(-1.2,0.95*Pmoon,'(a) Moon',fontsize=16,fontname= "Arial")

xticks = np.arange(-5,3.1,2)
ax.set_xticks(xticks)
ax.set_xticklabels([str(int(i)) for i in xticks],fontsize=15, fontname = "Arial")

xticks = np.arange(-6,6.1,2)
ax2.set_xticks(xticks)
ax2.set_xticklabels([str(int(i)) for i in xticks],fontsize=15, fontname = "Arial")

xticks = np.arange(-4,0.1,1)
ax3.set_xticks(xticks)
ax3.set_xticklabels([str(int(i)) for i in xticks],fontsize=15, fontname = "Arial")


ax_extra_y = ax.twinx()
ax_extra_y2 = ax2.twinx()
ax_extra_y3 = ax3.twinx()
ax_extra_y2.set_ylabel('Depth (km)',fontsize=22, fontname = "Arial")

import burnman
prem = burnman.seismic.PREM()
y_depth1 = prem.depth(yticks1*1e9)   ### need revise the burnman code from xx to xx.any()
y_depth1 = np.round(y_depth1/1e3)
ax_extra_y.set_yticks(yticks1)
ax_extra_y.set_yticklabels([str(int(i)) for i in y_depth1],fontsize=15, fontname = "Arial")
ax_extra_y.invert_yaxis()
ax_extra_y.minorticks_on()


plot_tool.show_minor_ticks(ax)


prem = burnman.seismic.PREM()
y_depth2 = prem.depth(yticks*1e9)   ### need revise the burnman code from xx to xx.any()
y_depth2 = np.round(y_depth2/1e3)
ax_extra_y2.set_yticks(yticks)
ax_extra_y2.set_yticklabels([str(int(i)) for i in y_depth2],fontsize=15, fontname = "Arial")
ax_extra_y2.invert_yaxis()
ax_extra_y2.minorticks_on()


plot_tool.show_minor_ticks(ax2)



prem = burnman.seismic.PREM()
y_depth3 = prem.depth(yticks3*1e9)   ### need revise the burnman code from xx to xx.any()
y_depth3 = np.round(y_depth3/1e3)
ax_extra_y3.set_yticks(yticks3)
ax_extra_y3.set_yticklabels([str(int(i)) for i in y_depth3],fontsize=15, fontname = "Arial")
ax_extra_y3.invert_yaxis()
ax_extra_y3.minorticks_on()

ax.set_xlabel('Redox within MO (' + r'$\mathrm{log } f_{\mathrm{O_{2}}}$'+  ' - IW)',fontsize=18, fontname = "Arial")
ax2.set_xlabel('Redox within MO (' + r'$\mathrm{log } f_{\mathrm{O_{2}}}$'+  ' - IW)',fontsize=18, fontname = "Arial")
ax3.set_xlabel('Redox within MO (' + r'$\mathrm{log } f_{\mathrm{O_{2}}}$'+  ' - IW)',fontsize=18, fontname = "Arial")

plot_tool.show_minor_ticks(ax3)


plt.tight_layout()
fig.legend(fancybox = False, loc ='upper center', ncol = 3, edgecolor = "black", bbox_to_anchor=(0.5, 0.1), bbox_transform=plt.gcf().transFigure, prop={'family': 'Arial', 'size': 15})
plt.subplots_adjust(bottom=0.2)
plt.savefig('./redox_profile.pdf',bbox_inches='tight')
plt.savefig('./redox_profile.png',bbox_inches='tight',dpi=600)
plt.show()




