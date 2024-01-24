#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:38:47 2019

@author: jiedeng
"""
import numpy as np
from src.geotherm_MgSiO3_melt import geotherm   
import matplotlib.pyplot as plt
from burnman import Mineral



P    = np.linspace(1e5,140e9,100)/1e9
#T25 = geotherm(P*1e9,Tm0=2500)
T21 = geotherm(P*1e9,Tm0=2100)


fig,ax = plt.subplots(1,1,figsize=(4,3))
 

#from vatic.plots import plot_tool
#plot_tool.load_default_setting()

ax.plot(T21,P,'k-')
#ax.plot(P,T25,'k--',label='hot')
#ax.legend(fontsize=13)
yticks = np.linspace(0,140,8)
plt.yticks(yticks,fontsize=13,family='Arial')
xticks = np.linspace(2000,5000,7)
plt.xticks(xticks,fontsize=13,family='Arial')
#plt.xtickslabels([str(int(i)) for i in xticks])
plt.ylim([0,140])
plt.xlim([2100,5000])
plt.ylabel("P (GPa)",fontsize=12,fontdict={'family':'Arial'})
ax.invert_yaxis()
plt.xlabel("T (K)",fontsize=12,fontdict={'family':'Arial'})
plt.minorticks_on()
#plot_tool.set_major_axis_font(ax,12)

fig.savefig('out/geotherm.pdf',bbox_inches='tight')
fig.savefig('out/geotherm.png',bbox_inches='tight',dpi=600)