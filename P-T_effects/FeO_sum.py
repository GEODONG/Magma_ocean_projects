#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 17:00:37 2019
FeO l: K14 Karki mine Frost 2010
FeO s: K14 SLB2011, campbell 2009 which likely already included in SLB2011, Frost 2010
@author: jiedeng
"""

from burnman import minerals
import numpy as np
import matplotlib.pyplot as plt

P  = np.linspace(1e5,100e9,100)
P2 = np.linspace(17e9,100e9,100)
T  = np.ones(P.shape)*4000

### solid 
feo_1 = minerals.Komabayashi_2014.feo()
feo_2 = minerals.SLB_2011.wuestite()
feo_3 = minerals.HHPH_2013.fper()


v_feo_1 = feo_1.evaluate(['molar_volume'],P,T)
g_feo_1 = feo_1.method.gibbs_free_energy_vector(P,T[0],feo_1.params)
v_feo_2,g_feo_2 = feo_2.evaluate(['molar_volume', 'gibbs'],P2,T) #P<2e9 numerical error raised
v_feo_3,g_feo_3 = feo_3.evaluate(['molar_volume','gibbs'],P,T)

### liquid
feo_1_l   = minerals.Komabayashi_2014.feo_l()
v_feo_1_l = feo_1_l.evaluate(['molar_volume'],P,T)
g_feo_1_l = feo_1_l.method.gibbs_free_energy_vector(P,T[0],feo_1_l.params)



################ cal fo2 ##################
fe_4  = minerals.HP_2011_ds62.iron()
feo_4 = minerals.HP_2011_ds62.fper()
o2_4  = minerals.HP_2011_fluids.O2()

v_fe_4,g_fe_4 = fe_4.evaluate(['molar_volume', 'gibbs'],P,T) #P<2e9 numerical error raised
v_feo_4,g_feo_4 = feo_4.evaluate(['molar_volume','gibbs'],P,T)
v_o2_4,g_o2_4 = o2_4.evaluate(['molar_volume','gibbs'],P,T)


plt.figure(figsize=(14,6))
plt.subplot(121)
plt.plot(P/1e9,v_feo_1.T*1e6,label='FeO Komabayashi 2014')
plt.plot(P2/1e9,v_feo_2.T*1e6,label='FeO SLB_2011')
#plt.plot(P/1e9,v_feo_3.T*1e6,label='FeO HHPH_2013')
plt.plot(P/1e9,v_feo_4.T*1e6,label='FeO HP_2011_ds62')
plt.plot(P/1e9,v_feo_1_l.T*1e6,'--',label='FeO melt Komabayashi 2014')

plt.legend()
plt.grid(True)
plt.subplot(122)
plt.plot(P/1e9,g_feo_1/1e3,'-',label='FeO Komabayashi 2014')
plt.plot(P2/1e9,g_feo_2/1e3,'-',label='FeO SLB_2011')
#plt.plot(P/1e9,g_feo_3/1e3,'-',label='FeO HHPH_2013')
plt.plot(P/1e9,g_feo_4.T/1e3,'-',label='FeO HP_2011_ds62')
plt.plot(P/1e9,g_feo_1_l/1e3,'--',label='FeO melt Komabayashi 2014')
plt.legend()
plt.grid(True)

dG_4 = g_feo_4 - g_fe_4 - 1/2*g_o2_4
dG_4_x = g_feo_4 - g_fe_4 

### O'Neil and Eggins 2002 Table 6 ####
dG_OG02 = -244118 + 115.559*T[0]-8.474*T[0]*np.log(T[0])
dG_OG02 = -244118 + 115.559*T[0]-8.474*T[0]*np.log(T[0])

T = Tgeo
f= -115997 + 27.036*T + 3.124*T*np.log(T)

#plt.plot(T,-f/T/8.31)