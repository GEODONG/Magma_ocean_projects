#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:09:04 2019

@author: jiedeng
"""



import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals


fe_fcc  = minerals.Komabayashi_2014.fe_fcc()
fe_hcp  = minerals.Komabayashi_2014.fe_hcp()
fe_l    = minerals.Komabayashi_2014.fe_l()
feo     = minerals.Komabayashi_2014.feo()
feo_l   = minerals.Komabayashi_2014.feo_l()

p = np.linspace(1e5,100e9,3)
t = np.ones(p.shape)*5000

v_fe_fcc = fe_fcc.evaluate(['molar_volume'], p , t)
fe_fcc.method.volume(1e9 , 5000, fe_fcc.params)
v = v_fe_fcc.T[:,0]

#fe_fcc.evaluate(['gibbs'], p , t)
fe_fcc.method.gibbs_free_energy(1e9, 3000, fe_fcc.params)

feo.method.gibbs_free_energy(100e9, 5000, feo.params)
feo_l.method.gibbs_free_energy(100e9, 5000, feo_l.params)

feo.method.gibbs_free_energy(100e9, 4000, feo.params)
feo_l.method.gibbs_free_energy(100e9, 4000, feo_l.params)

feo.evaluate(['molar_volume'], [100e9], [5000])
feo_l.evaluate(['molar_volume'], [100e9], [5000])


#xxx = np.array([0.97883707, 0.98137553, 0.98226149])
fe_fcc.method.pressure(t[0], v, fe_fcc.params)

#burnman.eos.modified_vinet()
#p_fe_fcc = fe_fcc.evaluate(['pressure'], v , t)
#fe_fcc.evaluate(['pressure'], v , t)
fe_fcc.evaluate(['shear_modulus'], v_fe_fcc.T , t)

v_fe_hcp = fe_hcp.evaluate(['molar_volume'], p , t)
v_fe_l   = fe_l.evaluate(['molar_volume'], p , t)
v_feo    = feo.evaluate(['molar_volume'], p , t)
v_feo_l  = feo_l.evaluate(['molar_volume'], p , t)



plt.figure()
#plt.plot(P/1e9,k14_p2v(P,T,name = 'Fe_l')*1e6,label='Fe_l')
plt.plot(p/1e9,v_fe_fcc.T*1e6,label='fe fcc')
plt.plot(p/1e9,v_fe_hcp.T*1e6,label='fe hcp')
plt.plot(p/1e9,v_fe_l.T*1e6,label='fe l')
plt.plot(p/1e9,v_feo.T*1e6,label='feo')
plt.plot(p/1e9,v_feo_l.T*1e6,label='feo l')

plt.legend()
plt.grid(True)


### mu_O2
#Binding energy (kJ/mol/formula unit. Schimka et al (2011) J. Chem. Phys. Table VI PBEsol)
Ebind = -650.0 

#Rotational constant (kJ/mol)
Be = 0.017346

Rgas    = 8.3145  # J/mol/K
avn     = 6.022e23 # avogadro's number
hplanck = 6.62607015e-34
boltzk  = 1.3806e-23
wat     = 15.9994 # atomic weight  g/mol

J2eV      = 6.242e18  # joule to eV
kJ_mol2eV = 1e3*J2eV/avn

reduced_hplanck = hplanck/2/np.pi
wet_kg = 5.31352e-26

thermal_wavelength = (2*np.pi*(reduced_hplanck**2)/(wet_kg*boltzk*t))**(.5)
Vq = thermal_wavelength**3 #m^3
log10fo2 = -27489./t + 6.702 + 0.045*p*1e4/t  # here p in bar

w0 = 1556 # cm-1
w0 = w0*0.012  ## conversion factor 0.012  == hplanck*100*3e8*avn/1000 
vib =  boltzk*t*J2eV*np.log(w0*1e3/avn/boltzk/t/2)
# in eV
mu_P0_T = Ebind*kJ_mol2eV + boltzk*t*J2eV*np.log(1e5*Vq/boltzk/t) + \
            boltzk*t*J2eV*np.log(2*Be*1e3/avn/boltzk/t) + vib
# J/mol

mu_P0_T = mu_P0_T*96e3




