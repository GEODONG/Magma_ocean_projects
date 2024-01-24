#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:57:42 2020

constants

@author: jiedeng
"""
import pandas as pd

R = 8.314
A3_to_cm3   = 1e-30*1e6*6.022e23
cm3_to_A3   = 1e30*1e-6/6.022e23


### zhang et al., 2017 model
z17    = {'a':[-6.376,-6.627],'b':[107257,110729],'c':[15095,15243],
          'd':[0.0827,0.1137],'K0':[36.61, 27.11],'K0_old':[30.3, 21.1]}
df_z17 = pd.DataFrame(z17,index=[2.92,3.69])


A18_FeO_re = {'V0':13.65,'T0':1673,'dVdT':2.92e-3,'K0':37,'Kp':8}
A18_FeO_ox = {'V0':21.07,'T0':1673,'dVdT':4.54e-3,'K0':12.6,'Kp':1.3}

"""
xlxs


"""

### for detailed info and for errors, see fitted sheet in oxidation_4
par_12p5_re = {'a':35.79397483, 'b':71.10313668,'c':36.59545225,'V0':1180.114014,'T0':3000,
               'P0':0,'K0':26.94713861, 'Kp':2.802531871, 'Kdp':0.012313472}

par_12p5_ox = {'a':34.52616394, 'b':68.64429623,'c':35.27069116,'V0':1204.763652,'T0':3000,
               'P0':0,'K0':23.19530062, 'Kp':3.216089358, 'Kdp':0.009340183}

par_25_re = {'a':31.34712676, 'b':62.48520005,'c':32.4675829,'V0':1192.011066,'T0':3000,
               'P0':0,'K0':23.95435759, 'Kp':	3.32104996, 'Kdp':-0.008912497}

par_25_ox = {'a':30.38414264, 'b':59.10950152,'c':29.64971394,'V0':1256.727179,'T0':3000,
               'P0':0,'K0':16.1261390, 'Kp':4.584011905, 'Kdp':-0.177152954}

