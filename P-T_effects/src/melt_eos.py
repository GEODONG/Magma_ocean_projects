#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 09:33:41 2018

@author: jiedeng
"""
from .eos_JD import *
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import scipy.constants as con

melt_info   = {25:{'Mg':14,'Fe':2,'Si':16,'O':49,'Ncell':16,'annot':'ox25','chem':'14 MgSiO3 +  2FeSiO3.5'}, 
             12.5:{'Mg':14,'Fe':2,'Si':16,'O':49,'Ncell':16,'annot':'ox12p5','chem':'14 MgSiO3 +  2FeSiO3.5'}, 
             6.25:{'Mg':30,'Fe':2,'Si':32,'O':97,'Ncell':32,'annot':'ox6p25','chem':'30 MgSiO3 +  2FeSiO3.5'}} 
melt_re_info   =   {25:{'Mg':14,'Fe':2,'Si':16,'O':49,'Ncell':16,'annot':'re25','chem':'14 MgSiO3 +  2FeSiO3.5'}, 
                  12.5:{'Mg':14,'Fe':2,'Si':16,'O':48,'Ncell':16,'annot':'re2p5','chem':'14 MgSiO3 +  2FeSiO3.5'}, 
                  6.25:{'Mg':30,'Fe':2,'Si':32,'O':96,'Ncell':32,'annot':'re6p25','chem':'30 MgSiO3 +  2FeSiO3.5'},
                  0:{'Mg':16,'Fe':0,'Si':16,'O':48,'Ncell':16,'annot':'re0','chem':'16 MgSiO3'}} 

pars = pd.read_excel('db/oxidation_4.xlsx',sheet_name='fitted',
                     usecols = list(range(10)),nrows = 5,index_col=0).dropna()
#print(pars)
def MgSiO3_vt2p(vol,T):
    """
    MgSiO3 volume + temperature -> pressure
    ----------------------------------------
    Input list
    ----------------------------------------
    vol    : one formula unit (A^3)
    T      : temperature (K)
    ----------------------------------------
    Output list
    ----------------------------------------    
    P       : pressure (GPa)
    ---------------------------------------- 
    Ref. Karki et al. 2018 
    """
    
    par = pars.loc['MgSiO3']
#    if vol<par['V0']: ##vol>par['V0'] does not work!!
#        P = 0
#    else:
    P = BM4_TH(vol,K0=par['K0'],Kp=par['Kp'],Kdp=par['Kdp'],V0=par['V0'],
               P0=par['P0'],T=T,T0=par['T0'],a=par['a'],b=par['b'],c=par['c'])        
    return P
#def MgSiO3_vt2p(vol,T):
#    par = pars.loc['MgSiO3']
#    P = bm(vol,K0=par['K0'],Kp=par['Kp'],Kdp=par['Kdp'],V0=par['V0'],
#               P0=par['P0'],T=T,T0=par['T0'],a=par['a'],b=par['b'],c=par['c'])
#    return P
def MgO_vt2p(vol,T):
    """
    MgO volume + temperature -> pressure
    ----------------------------------------
    Input list
    ----------------------------------------
    vol    : one formula unit (A^3)
    T      : temperature (K)
    ----------------------------------------
    Output list
    ----------------------------------------    
    P       : pressure (GPa)
    ---------------------------------------- 
    Ref. Ghosh and Karki 2016 
    """
    par  = pars.loc['MgFeO']
    Pc   = BM3(vol,K0=par['K0'],Kp=par['Kp'],V0=par['V0'],P0=par['P0'])
    temp = (vol-par['V0'])/par['V0']
    Pth  = 0.11/vol*(0.75-temp)
    return Pc+Pth

    
def MgO_pt2v(P,T):
    """
    MgO pressure + temperature -> temperature
    ----------------------------------------
    Input list
    ----------------------------------------
    P       : pressure (GPa)
    T       : temperature (K)
    ----------------------------------------
    Output list
    ----------------------------------------    
    V_std   : V in standard unit  m^3/mol
    rho_std : density in standard unit  kg/m^3
    V_A3    : V in A^3, 1 formula unit
    ---------------------------------------- 
    Ref.
    ----------------------------------------    
    Both Karki et al., 2018 and Karki et Ghosh and Karki 2016 have Fe_num dependecy
    Here I follow Kari et al., 2018
    """
    par         = pars.loc['MgFeO']
    alpha_guess = 1e-5
    v0_guess    = par['V0']*np.exp((T-par['T0'])*alpha_guess)
    V_in        = np.linspace(v0_guess*1.2,v0_guess*0.2,1000)
    P_out       = MgO_vt2p(V_in,T)
    f           = interp1d(P_out,V_in)
    V       = f(P)
    V_std   = V*1e-30*con.N_A
    V_A3    = V
    Mm      = 24.31 + 16
    rho_std = Mm*1e-3/V_std  
    return V_std,rho_std,V_A3

def MgFeO_pt2v(Fe_num,P,T,spin = 0):
    """
    MgO volume + temperature -> pressure
    ----------------------------------------
    Input list
    ----------------------------------------
    Fe_num : 0-100
    P      : pressure (GPa) list or scalar
    T      : temperature (K) scalar
    spin   : 0-1, 0 high spin, 1 low spin
    ----------------------------------------
    Output list
    ----------------------------------------    
    V_std   : V in standard unit  m^3/mol
    rho_std : density in standard unit  kg/m^3
    V_A3    : V in A^3, 1 formula unit
    ---------------------------------------- 
    Ref.  Karki et al., 2018 

    """
    x = Fe_num/100
    V_MgO_std, rho_MgO_std,_ = MgO_pt2v(P,T)
    Mm       = 24.31*(1-x) + 55.85*x + 16
    rho_std  = rho_MgO_std*(1+x*(0.66+0.12*spin))
    V_std    = Mm*1e-3/rho_std
    V_A3     = V_std/con.N_A*1e30
    return V_std,rho_std,V_A3

def MgFeO_pt2v_vector(Fe_num,P,T,spin = 0):
    """
    """
    f_v = np.vectorize(MgFeO_pt2v,excluded=[0,3])
    return f_v(Fe_num,P,T,spin)
    
def MgSiO3_pt2v(P,T):
    par = pars.loc['MgSiO3']
    #Kp = pars.loc['Kp']
    #v_guess = (Kp*(P-par['P0'])/pars.loc['K0']+1)**(-1/Kp)
    alpha_guess = 4e-5
    v0_guess = par['V0']*np.exp((T-par['T0'])*alpha_guess)
    V_in     =  np.linspace(v0_guess*1.5,v0_guess*0.2,1000)
    P_out    =  MgSiO3_vt2p(V_in,T)
    # if min of P_out is not first or last one, turnover occurs
    minPindex = np.argmin(P_out)
    P_useful  = P_out[minPindex:]
    V_useful  = V_in[minPindex:]
    f         = interp1d(P_useful,V_useful)
    # A^3
    ## if the target P is smaller than the 
#    V = np.zeros(np.shape(P))
#    if (P<P_out[minPindex]).any():
#         nearminP         = np.argmin(np.abs(P - P_out[minPindex]))
#         V[(nearminP+1):] = f(P[(nearminP+1):])
#         V[:(nearminP+1)] = np.ones((nearminP+1))*V_in[minPindex]
#    else:
    V = f(P)
    # m^3/mol, 16 formula unit
    V_std = V*1e-30/16*con.N_A
    V_A3  = V/16
    # g/mol
    Mm = 24.31 + 28.09*1 + 16*3
    # kg/m^3
    rho_std = Mm*1e-3/V_std  
    return V_std,rho_std,V_A3

#def MgSiO3_pt2v_single(P,T):
#    import scipy.optimize as opt
#    par = pars.loc['MgSiO3']
#    func = lambda vol: MgSiO3_vt2p(vol,T) - P
#    alpha_guess = 10e-5
#    v0_guess = par['V0']*np.exp((T-par['T0'])*alpha_guess)
#    V = opt.brentq(func, v0_guess*0.2, 1.5 * v0_guess)
#    return V
#
#def MgSiO3_pt2v(P,T):
#    f_v = np.vectorize(MgSiO3_pt2v_single)
#    V   = f_v(P,T)
#
#    V_std = V*1e-30/16*con.N_A
#    V_A3  = V/16
#    # g/mol
#    Mm = 24.31 + 28.09*1 + 16*3
#    # kg/m^3
#    rho_std = Mm*1e-3/V_std  
#    return V_std,rho_std,V_A3    

def MgFeSiO3_pt2v(Fe_num,P,T,spin = 0):
    """
    Fe2+
    spin = high spin 0, low spin 1 in MgFeSiO3
    ref. Karki 2018
    #### Jan2119 ###
    Fe_num = Fe_num/10 -> Fe_num = Fe_num/100
    #### Jan2119 ###
    """
    Fe_num = Fe_num/100
    V_MgSiO3_std, rho_MgSiO3_std,_ = MgSiO3_pt2v(P,T)
    Mm      = 24.31*(1-Fe_num) + 55.85*Fe_num + 28.09*1 + 16*3 
    rho_std = rho_MgSiO3_std*(1+Fe_num*(0.28+0.04*spin))
    # m^3/mol
    V_std   = Mm*1e-3/rho_std
    V_A3    = V_std/con.N_A*1e30
    return V_std,rho_std,V_A3

def MgFe3SiO3_vt2p(Fe_num,vol,T):
    """
    Fe3+
    """
    par   = pars.loc[melt_info[Fe_num]['annot']]
    P = BM4_TH(vol,K0=par['K0'],Kp=par['Kp'],Kdp=par['Kdp'],V0=par['V0'],
               P0=par['P0'],T=T,T0=par['T0'],a=par['a'],b=par['b'],c=par['c'])
    return P
  
def MgFe3SiO3_pt2v(Fe_num,P,T):
    """
    P, T -> V function for Fe3+ bearing MgFeSiO3 melt
    parameters:
    P GPa
    """
    par         = pars.loc[melt_info[Fe_num]['annot']]
#     print(par)
    alpha_guess = 1e-5
    v0_guess    = par['V0']*np.exp((T-par['T0'])*alpha_guess)
    V_in        = np.linspace(v0_guess*1.2,v0_guess*0.2,1000)
    P_out       = MgFe3SiO3_vt2p(Fe_num,V_in,T)
    f           = interp1d(P_out,V_in)
    # A^3
    V = f(P)
    # m^3/mol, cell_unit[Fe_num] formula unit
    V_std = V*1e-30/melt_info[Fe_num]['Ncell']*con.N_A
    V_A3  = V/melt_info[Fe_num]['Ncell']
    # g/mol
    Mm    =  (24.31*melt_info[Fe_num]['Mg'] + 55.85*melt_info[Fe_num]['Fe']+
          28.09*melt_info[Fe_num]['Si'] + 16*melt_info[Fe_num]['O'])/melt_info[Fe_num]['Ncell']
    # kg/m^3
    rho_std = Mm*1e-3/V_std  
    return V_std,rho_std,V_A3
    
    
    
##
#Px = np.linspace(-100,100,100)
#T = 4000
####MgFe3SiO3_pt2v(6.25,Px,T)
###MgO_pt2v(0,T)
#Fe_num =6.25
#V_A3=MgFeSiO3_pt2v(Fe_num,Px,T=4000)[-1]*melt_info[Fe_num]['Ncell']
#
##print(V_A3)
#
##v0 = pars.loc['MgSiO3']['V0']
##vol = np.linspace(v0*1.5,v0*0.9,100)
##p4 = MgSiO3_vt2p(vol,T)
#import matplotlib.pyplot as plt
#plt.figure()
#plt.plot(Px,V_A3)
#plt.grid(True)
#par = pars.loc['MgSiO3']
#p3 = BM3_TH(vol,K0=par['K0'],Kp=par['Kp'],V0=par['V0'],
#           P0=par['P0'],T=T,T0=par['T0'],a=par['a'],b=par['b'],c=par['c'])
#
#plt.figure()
#plt.plot(p3,vol*2,p4,vol*2)
#plt.grid(True)
