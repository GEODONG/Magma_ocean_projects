#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 00:18:02 2019

@author: jiedeng
"""
import burnman
from burnman import minerals
import numpy as np
from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt

def geotherm(P,Tm0=2500):
    pv_l = minerals.LF_2012.mgsio3()
    V0  = pv_l.evaluate(['molar_volume'],[0],[Tm0])[0];
    V   = np.linspace(V0/3,V0,100);
    tmp = np.log(V/V0)*(pv_l.params['grueneisen_prime'] - pv_l.params['grueneisen_0']) - \
          pv_l.params['grueneisen_prime']/V0*(V-V0)
    Tm  = Tm0*np.exp(tmp);
    Pm  = np.zeros(V.shape);
    for i in range(len(V)):
        Pm[i] = pv_l.method.pressure(Tm[i], V[i], pv_l.params)   
#    print("Pm is",Pm)
#    print("Tm is",Tm)

#    rho   = pv_l.params['molar_mass']/V;
    f   = interp1d(Pm[:,0],Tm[:,0])
#    f   = interp1d(Pm,Tm)

    Tin = f(P)
    return Tin

#pv_l = minerals.LF_2012.mgsio3()
#Tm0 = 2500;
#P   = np.linspace(1e5,100e9,100)
#V0  = pv_l.evaluate(['molar_volume'],[0],[Tm0])[0];
#V   = np.linspace(V0/2,V0,100);
#tmp = np.log(V/V0)*(pv_l.params['grueneisen_prime'] - pv_l.params['grueneisen_0']) - \
#      pv_l.params['grueneisen_prime']/V0*(V-V0)
#Tm  = Tm0*np.exp(tmp);
#Pm  = np.zeros(V.shape);
#for i in range(len(V)):
#    Pm[i] = pv_l.method.pressure(Tm[i], V[i], pv_l.params)   
#rho   = pv_l.params['molar_mass']/V;

#plt.figure(figsize=(12,8))
#plt.title('Geotherm of magma ocean')
#plt.subplot(221)
#plt.plot(Pm/1e9,Tm)
#plt.xlabel('Pressure (GPa)')
#plt.ylabel('Tm (K)')
#plt.grid(True)




