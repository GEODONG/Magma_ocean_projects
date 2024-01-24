#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 00:29:48 2020

@author: jiedeng
"""

import src.melt_eos as ms
import matplotlib.pyplot as plt
import numpy as np
from burnman import minerals
from scipy.interpolate import interp1d

import pandas as pd
from vatic.plots import plot_tool

plot_tool.load_default_setting()
#header  = ['P(GPa)','Depth(km)','T(K)','Vphi(km/s)','Vs(km/s)','Vp(km/s)','Vs_att(km/s)','Vp_att(km/s)',
#          'S(kJ/g)','S(J/g/K','alpha(1e-5K-1)','Cp(J/g/K)','KT(GPa)','Qs','Qp','rho(g/cm^3)','phase']
#elastic = pd.read_table('/Users/jiedeng/GD/papers/paper14_partition_ff/codes8/fort.56',sep=r"\s+",names=header,engine='python',index_col=2)

header=            ['P(GPa)' ,    'depth'    ,    'T(K)'      ,     'Vol'     ,       'KS',  
            'KT(GPa)'   ,      'alpha'     ,   'Heat C'   , 
            'thet'      ,       'g'     ,        'q'   ,
            'Velocity'   ,     'Ppart1'    ,    'Ppart2']
elastic = pd.read_table('/Users/jiedeng/GD/papers/paper14_partition_ff/codes8/fort.59',sep=r"\s+",names=header,engine='python',index_col=2)


def interp(T,P,param = 'KT(GPa)'):
    f = interp1d(elastic.loc[T]['P(GPa)'],elastic.loc[T][param])
    out = f(P)
    return out



def KT_K18(P,T):
    """
    P : array in GPa
    T : scalar in K
    """
    V_std,rho_std,V_A3 = ms.MgSiO3_pt2v(P,T)        
    KTs = -(P[1:] - P[:-1])/(V_std[1:] - V_std[:-1])*V_std[1:]
    return KTs

    
def KT_LF2012(P,T):
    """
    P : array in GPa
    T : scalar in K
    """
    pv_l  = minerals.LF_2012.mgsio3()
    V_std = (pv_l.evaluate(['V'],P*1e9,T*np.ones(P.shape))).reshape(P.shape)
    KTs   = -(P[1:] - P[:-1])/(V_std[1:] - V_std[:-1])*V_std[1:]
    return KTs

def KT_lars(P,T):
    """
    P : array in GPa
    T : scalar in K
    #use Lars' result

    """
    return interp(T,P)
    

P = np.linspace(20,140,141)

# bulk modulus at target pressure from 
P_target = 78


plt.figure()
#plt.plot(P[1:],KT_K18(P,3000),'g-',label='Karki et al 18')
##plt.plot(P[1:],KT_LF2012(P,3000),'g--',label='Liebske & Frost 2012')
##plt.plot(P,KT_lars(P,3000),'g:',label='HeFESTo')
#
#plt.plot(P[1:],KT_K18(P,4000),'b-')
#plt.plot(P[1:],KT_LF2012(P,4000),'b--')
#plt.plot(P,KT_lars(P,4000),'b:')

plt.plot(P[1:],KT_K18(P,5000),'k-',label='Karki et al., 2018')
plt.plot(78,279.8,'ro',label='This study')
#plt.plot(P[1:],KT_LF2012(P,5000),'r--')
#plt.plot(P,KT_lars(P,5000),'r:')

plt.legend()
plt.grid()
plt.xlabel('P (GPa)')
plt.ylabel(r'$K_{T}$ (GPa)')
plt.savefig('/Users/jiedeng/GD/ppt/2020/fig10_kt.pdf',bbox_inches = 'tight')


pv_l  = minerals.LF_2012.mgsio3()
V_std = (pv_l.evaluate(['V'],P*1e9,4000*np.ones(P.shape))).reshape(P.shape)

KTs   = -(P[1:] - P[:-1])/(V_std[1:] - V_std[:-1])*V_std[1:]