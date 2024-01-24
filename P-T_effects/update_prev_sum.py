#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 10:43:12 2019

include some routines to calculate the PV term, which is not performed in
all the codes, but may be useful to recall how we get this term.


@author: jiedeng
"""

import IW_buffer as iw
import numpy as np
import tools as tl
import tools2 as tl2

R = 8.314

np.log10(iw.Ru_RuO2_single(p=23,t=2473))
np.log10(iw.Ru_RuO2_single(p=10,t=2173))

np.log10(iw.Ru_RuO2_single(p=4,t=1873))
np.log10(iw.Ru_RuO2_single(p=6.5,t=2023))

np.log10(iw.Ru_RuO2_2_single(P=6.5,T=2023))

P = 10
T = 2023
np.log10(iw.Ru_RuO2_single(p=10,t=2173))
np.log10(iw.Ru_RuO2_single(p=10,t=2023))
np.log10(iw.Ru_RuO2_2_single(10,2023))
np.log10(iw.Ru_RuO2_single(p=15,t=2373))
np.log10(iw.Ru_RuO2_single(p=15,t=2323))
np.log10(iw.Ru_RuO2_2_single(15,2323))


_, _,_, _, dv_fe_25 = tl.cal_dV_this_study(23,2473,name='Fe_12p5',flag=False)

#P = 6.5
#T = 2023
#
P = 23
T = 2473

#P = 90
#T = 3300
Parray = np.linspace(1e-4,P)
Tarray = np.ones(Parray.shape)*T
dV     = tl2.dv_12p5(Parray,Tarray) # all high P are at 12.5%
PV     = tl.cal_PV(dV,Parray,T,P,minP=0)
print(PV)

dV     = tl2.dv_12p5(Parray,Tarray-50) # all high P are at 12.5%
PV     = tl.cal_PV(dV,Parray,T-50,P,minP=0)
print(PV)


dV     = tl2.dv_25(Parray,Tarray+50) # all high P are at 12.5%
PV     = tl.cal_PV(dV,Parray,T+50,P,minP=0)
print(PV)


dV     = tl2.dv_12p5(Parray,Tarray-500) # all high P are at 12.5%
PV     = tl.cal_PV(dV,Parray,T-500,P,minP=0)
print(PV)



#### combined the G model !!! => combine Gr_janaf and Gr_janaf old
import matplotlib.pyplot as plt
tin = np.linspace(200,10000,100)

plt.plot(tin,tl.Gr_janaf_old(tin),label='old')
plt.plot(tin,tl.Gr_janaf(tin),label='new')
plt.legend()
plt.grid(True)

### before 2100, use the new, after 2100, use the old one
new = tl.Gr_janaf(tin)
old = tl.Gr_janaf_old(tin)


Told = 3500
Tnew = 2500

dt = Told -Tnew

k = (tl.Gr_janaf_old(Told) - tl.Gr_janaf(Tnew))/dt

Gcombined  = np.zeros_like(new)
  
for i in range(len(tin)):
    if tin[i]<Tnew:
        Gcombined[i] = tl.Gr_janaf(tin[i])
    elif tin[i] > Told:
        Gcombined[i] = tl.Gr_janaf_old(tin[i])
    else:
        Gcombined[i] = (tin[i] - Tnew)*k + tl.Gr_janaf(Tnew)

plt.plot(tin,tl.Gr_janaf_old(tin),label='old')
plt.plot(tin,Gcombined,label='combined')
plt.plot(tin,tl.Gr_janaf(tin),label='new')
plt.legend()
plt.grid(True)


def cal_G(T,a,b,c,d,e,f):
    """
    """
    G = a+b*T+c*T*np.log(T) + d*(T**2) + e/T + f*(T**.5)
    return G


from lmfit import Parameters
from lmfit import report_fit,Model

params = Parameters()
params.add('a', value = -245000)
params.add('b', value = 231) 
params.add('c', value = -46) 
params.add('d', value = -5e-3) 
params.add('e', value = 0) 
params.add('f', value = 0) 
model = Model(cal_G)

#data = uct2pd(dGr)

#res = model.fit(data['value'].values,params,T=tin,weights=1./data['dev'].values)
res = model.fit(Gcombined,params,T=tin)

res.plot_fit()

plt.plot(tin,tl.Gr_janaf_old(tin),label='old')
#plt.plot(tin,Gcombined,label='combined')
plt.plot(tin,tl.Gr_janaf(tin),label='new')
#uplot.uplot(tin,G_fitted,'k')
#uplot.uplot(tin,dGr,'r')

#a =

from vatic import uplot as uplot

plt.figure(figsize=(5,4))
#plt.plot(tin,tl.Gr_janaf_old(tin),label='old')
#uplot.uplot(tin,tl.Gr_janaf_old(tin,True),label='old')
uplot.uplot(tin,tl.Gr_janaf_com(tin,True),label='combined')

#plt.plot(tin,res.best_fit,label='combined')
#plt.plot(tin,tl.Gr_janaf(tin),label='new')

#plt.plot(tin,tl.Gr5(tin),'b-',label='G03')
#plt.plot(tin,tl.Gr(tin),'g-',label='Model J04')
#plt.plot(tin,-16201*R + 8.031*tin*R,'k-',label='Fit J04')

plt.legend()
plt.xlim([1000,5000])
plt.ylim([-80000,80000])
plt.xlabel("T (K)",fontsize=13)
plt.ylabel(r'$\Delta$'+r'$G_{r}^{0}$'+ '(J/mol)',fontsize=13)



tl.Gr_janaf_old(tin[25])/2010/R
tl.Gr_janaf(tin[25])/2010/R
res.best_fit[25]/2010/R


TT = 2100
res.model.eval(res.params,T=TT)/R/TT  
tl.Gr_janaf_old(TT)/R/TT
tl.Gr_janaf(TT)/R/TT



########### input dV error ###########
Pin = 10
Tin = 2123
Parray = np.linspace(1e-4,Pin)
Tarray = np.ones(Parray.shape)*Tin

_,_,_,_,dV_test = tl.cal_dV_this_study(Parray,Tarray,name='Fe_12p5',flag=True)

tl.cal_PV_uct(dV_test,Parray,Tin,maxP=Pin,minP=0)


###



def cal_PV_error(Pin,Tin):
    Parray = np.linspace(1e-4,Pin)
    Tarray = np.ones(Parray.shape)*Tin
    
    _,_,_,_,dV_test = tl.cal_dV_this_study(Parray,Tarray,name='Fe_12p5',flag=True)
    
    PV_out  = tl.cal_PV_uct(dV_test,Parray,Tin,maxP=Pin,minP=0)
    std_dev = PV_out.std_dev
    return PV_out.nominal_value, std_dev

def cal_PV_error_vector(Pin,Tin):
    f_v = np.vectorize(cal_PV_error)
    return f_v(Pin,Tin)


### info of high P by running figS5 ###
#Pin = Z17.iloc[:-1]['P(GPa)'].values
#Tin = Z17.iloc[:-1]['T(K)'].values

Pin = O06['P(GPa)'].values
Tin = O06['T(K)'].values
val,std=cal_PV_error_vector(Pin,Tin)
