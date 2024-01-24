#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 16:16:53 2019

@author: jiedeng
"""
import pandas as pd
from lmfit import Parameters
import numpy as np
from src.eos_JD import BM4, BM4_TH,BM3
from scipy.interpolate import interp1d 
import scipy.optimize as opt

xlsx  = 'db/oxidation_4.xlsx'  
## note here should be 'db/oxidation_4.xlsx' rather than 'db/oxidation_4.xlsx'
## path relative to the file being executed 
sheet = 'fitted'

############################## plot ##############################

skiprow_dic = {'FeO_4k':9, 'Fe_25_4k':12, 'Fe_25_3k':15,  
               'Fe_12p5_4k':18, 'Fe_12p5_3k':21,
               'Fe_6p25_4k':24, 'Fe_6p25_3k':27,
               'previous':45}

def read_par(skiprows,usecol=7):
    """
    read par for xlsx
    """
    par_in  = pd.read_excel(xlsx,sheet_name='fitted',
                              usecols = list(range(usecol)), 
                              index_col = 0,
                              skiprows = list(range(skiprows)), 
                              nrows = 2).dropna()
    return par_in

def find_par(name):
    """
    """
    if name[-1] == 'k':
        print("No oxidized or reduced is specified")
        pass
    elif name[0] == 'P': ## P means previous
        name = name[1:]
        skiprow = skiprow_dic['previous']
        nrows   = 12
    elif name[-2:] == 're' or 'ox':
        skiprow = skiprow_dic[name[:-3]]
        nrows = 2
           
    par = pd.read_excel(xlsx,sheet_name=sheet, 
                      usecols = list(range(7)), 
                      index_col = 0,
                      skiprows = list(range(skiprow)), 
                      nrows = nrows,).dropna()
#    print("par is",par)
    return par.loc[name]

def vt2p(v,t,par):
    """
    """
    k0 = par['K0']
    kp = par['Kp']
    if 'P0' in par.index:
        p0 = par['P0']
    else:
        p0 = 0
        
    if 'Kdp' in par.index:
        kdp = par['Kdp']
    else:
        kdp = 0
        
    if 'dVdT' in par.index:
        dVdT = par['dVdT']
        v0 = par['V0'] + (t-par['T0'])*dVdT
    else:
        v0 = par['V0']
           
    p  = BM4(v,k0,kp, kdp, v0,p0)
    return p

def pt2v_single(p,t,name):
    """
    calculate dV at given P,T
    """
    par = find_par(name)
    if 'dVdT' in par.index:
        dVdT = par['dVdT']
        v0 = par['V0'] + (t-par['T0'])*dVdT
    else:
        v0 = par['V0']
    func = lambda vol: vt2p(vol,t,par) - p
    v = opt.brentq(func, v0, v0*0.2)
    return v


def pt2v(p,t,name):
    """
    calculate dV at given P,T
    """
    f_v = np.vectorize(pt2v_single, excluded=[2])
    return f_v(p,t,name)


def pt2dv(p,t,re,ox):
    """
    calculate dV at given P,T

    input list
    ----------
    P,T
    re -> reduced name
    ox -> oxidized name
    
    """
    re_par = find_par(re)
    ox_par = find_par(ox)
    re_v = pt2v(p,t,re_par)
    ox_v = pt2v(p,t,ox_par)
    dV = ox_v - re_v
    return dV

#######
    
def bm2_pt2v_single(p,t,v0,t0,dvdt,K0,Kp):
    vt0 = v0 + (t-t0)*dvdt
    V = vt0*(1+p*(Kp/K0))**(-1/Kp)
    return V

def bm2_vt2p_single(v,t,v0,t0,dvdt,K0,Kp):
    vt0 = v0 + (t-t0)*dvdt
    P = ((v/vt0)**(1/Kp) - 1)/(Kp/K0)
    return P
def pt2v_prev(p,t,name):
    """
    calculate dV at given P,T
    p : pressure
    t : temperature
    
    """
    par = find_par(name)
#    print("par is",par)
    f_v = np.vectorize(bm2_pt2v_single, excluded=[2,3,4,5,6])
    return f_v(p,t, par['V0'], par['T0'],par['dVdT'],par['K0'], par['Kp'])

 
def pt2v_this_study(p,t,name,flag=True):
    """
    calculate V at given P, T
    """
    skiprow_dic  = {'FeO':93, 'Fe_25':90,'Fe_12p5':87, 'Fe_6p25':84}
    skiprow_read = skiprow_dic[name[:-3]]
    par          = read_par(skiprow_read,usecol=17)
    if name[-2:] == 're':
        par_selected = par.iloc[0]
    else:
        par_selected = par.iloc[1]
#    print("par is")
#    print(par_selected)
    try:
        v = BM4_TH_pt2v_vector_uct(par_selected,p,t,flag=flag)
    except:
        print("Brentq cannot process it! Turn to extrapolation method")
        vx   = np.linspace(par_selected['V0']*0.2,par_selected['V0']*1.8,1000)
        v    = np.zeros(p.shape)
        if np.unique(t).size == 1:
            print("Temperatures are all equal to",t[0])
            px   = BM4_TH_vt2p_vector_uct(par_selected,vx,t[0],flag=flag)
            minP = np.argmin(px)
            if minP < len(px):
                print("The min P of {} at {:2} K is around {:2.4} GPa, corresponding V is {:2.2} A^3".format(name,t[0], px[minP], vx[minP]))
            f    = interp1d(px[:minP],vx[:minP])
            v    = f(p)        
        else:
            print("Temperatures are varying, SLOW!!")
            for i in range(len(p)):
                px   = BM4_TH_vt2p_vector_uct(par_selected,vx,t[i],flag=flag)
                minP = np.argmin(px)
    #        if minP > len(px):
    #            print("The min P of reduced at {:2} K is around {:2.4} GPa, corresponding V is {:2.2} A^3".format(tint, px[minP], vx[minP]))
                f = interp1d(px[:minP],vx[:minP])
                v[i] = f(p[i])
    return v

    
def BM4_TH_pt2v(V0,P0,T0,a,b,c,K0,Kp,Kdp,p,t):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    """

    func = lambda vol: BM4_TH(vol, K0, Kp, Kdp, V0, P0, t, T0, a, b, c) - p
#    V = opt.brentq(func, 0.1 * V0, 1.5 * V0)
#    print("V0 is",V0)
#    print("0.1*V yields",BM4_TH(0.1 * V0, K0, Kp, Kdp, V0, P0, t, T0, a, b, c) - p)
#    print("1.5*V yields",BM4_TH(1.5 * V0, K0, Kp, Kdp, V0, P0, t, T0, a, b, c) - p)
    try:
        V = opt.brentq(func, 0.1 * V0, 1.5 * V0)
    except:
        print("Brentq does not work")
    return V

import uncertainties as uct
from uncertainties import umath


def BM4_TH_pt2v_vector(V0,P0,T0,a,b,c,K0,Kp,Kdp,p,t):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[GPa]`. Returns molar volume in :math:`[m^3]`
    
    unit of p/v should be consistent with K0/V0 in par
    both p,t can be a scalar or vector/ ufloat or float
    """
    f_v = np.vectorize(uct.wrap(BM4_TH_pt2v), excluded=list(range(9)))
    return f_v(V0,P0,T0,a,b,c,K0,Kp,Kdp,p,t)



def BM4_TH_pt2v_vector_uct(par,p,t,flag=True):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    """
    V0,P0,T0,a,b,c,K0,Kp,Kdp = unpack_par(par,t,flag)
    return BM4_TH_pt2v_vector(V0,P0,T0,a,b,c,K0,Kp,Kdp,p,t)

def BM4_TH_vt2p_vector_uct(par,v,t,flag=True):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    """
    V0,P0,T0,a,b,c,K0,Kp,Kdp = unpack_par(par,t,flag)
    f_v = np.vectorize(uct.wrap(BM4_TH), excluded=[1,2,3,4,5,7,8,9,10])    
    return f_v(v, K0, Kp, Kdp, V0, P0, t, T0, a, b, c)


def unpack_par(par,t,flag):
    """
    unpack the par, 
    if t == T0, no uncertainties for a,b,c
    if t != T0, no uncertainties for K0,Kp,Kdp
    implicitly, if t is a array, I treat as t != T0
    """
    V0 = par['V0']
    P0 = par['P0']
    T0 = par['T0']
#    print(par)
    if flag:
        a = uct.ufloat(par['a'], par['aerr'])
        b = uct.ufloat(par['b'], par['berr'])
        c = uct.ufloat(par['c'], par['cerr'])
        K0  = uct.ufloat(par['K0'], par['K0err'])
        Kp  = uct.ufloat(par['Kp'], par['Kperr'])
        Kdp = uct.ufloat(par['Kdp'], par['Kdperr'])

    else:
        a = par['a']
        b = par['b']
        c = par['c']
        K0  = par['K0']
        Kp  = par['Kp']
        Kdp = par['Kdp']            
    return V0,P0,T0,a,b,c,K0,Kp,Kdp

###############################################################################    
A3_to_cm3   = 1e-30*1e6*6.022e23
cm3_to_A3   = 1e30*1e-6/6.022e23
def cal_dV_this_study(P,T,name='FeO',flag=False):
    """
    calculate dV = VFeO1.5 - VFeO
    
    params
    ----------
    P : pressure, list
    T : temperature, list
    name : FeO, Fe_25, Fe_12p5
    
    reutrn list
    -----------
    v_re : reduced volume in cm3/mol
    v_ox : oxidized volume in cm3/mol
    v_re_org : reduced volume in A3 with orginal formula unit
    v_re_org : oxidized volume in A3 with orginal formula unit
    dv: v_ox - v_re in in cm3/mol
    """
    if name == 'FeO':
        fu = 32
    elif name == 'Fe_25':
        fu = 4
    else:
        fu = 2
    v_re_org = pt2v_this_study(P,T,name+'_re',flag)
    v_re     = v_re_org/fu*A3_to_cm3
    v_ox_org = pt2v_this_study(P,T,name+'_ox',flag)
    v_ox     = v_ox_org/fu*A3_to_cm3
    dv       = v_ox - v_re
    return v_re, v_ox,v_re_org,v_ox_org, dv

def cal_PV(dV,P,T,maxP,minP=0):
    """
    integrate dV*P at cutoff pressure maxP and temperature T
    Parameters
    ----------
    dV: cm^3/mol
    P : GPa
    maxP : max P
    Return
    ------
    PV/R/T
    """
    R = 8.314
    min_ind = np.argmin(np.abs(P-minP))
    max_ind = np.argmin(np.abs(P-maxP))
#    print("min_ind is",min_ind)
#    print("max_ind is",max_ind)
#    if min_ind !=0:
#        print("Note the input P does not start from 0")    
    intg = np.trapz(dV[min_ind:(max_ind+1)]*1e-6,P[min_ind:(max_ind+1)]*1e9)
#    print("integration is", intg)
    return intg/R/T

def cal_PV_uct(dV,P,T,maxP,minP=0):
    cal_PV_uct = uct.wrap(cal_PV)
    return cal_PV_uct(dV,P,T,maxP,minP=0)


##############################FeO + 1/4O2 = Fe2O3###############################    
def Gr_J04_fit(T):
    """
    added by Jie -16201*R + 8.031*np.linspace(1278,1908,100)*R
    """
    R = 8.314
    return -16201*R + 8.031*T*R
    
def Gr(T):
    """
    calculate the free energy of reaction FeO + 1/4O2 = FeO1.5 for pure endmembers by Jayasuriya et al., 2004
    
    FeO (l) + 1/4O2 = FeO1.5 (l)
    return list
    -----------
    J/mol
    -----------
    Jayasuriya et al., 2004
    """
    return -115997 + 27.036*T + 3.124*T*np.log(T)

def Gr5(T):
    """
    calculate the free energy of reaction FeO + 1/4O2 = FeO1.5 for pure endmembers by Gaillard et al., 2003

    Fe + 1/2O2 = FeO(l) ---(1)
    2Fe + 3/2O2 = Fe2O3 (l)---(2)
    
    FeO (l) + 1/4O2  = FeO1.5 (l) 
    
    ----
    Gaillard does not directly gives the way of calcualting the Gr,
     but it gives the thermoydnamics of several intermediate reactions.
     Gr5 means we use the combination of (1) and (2). There are many other ways,
     Gr3,Gr4... which are deprecated because they are way too off
    ref Gaillard et al., 2003
    """
    G1 = -226244 + 42.29*T
    G2 = -0.056*(T**2) + 374.59*T - 846564
    return G2/2 -G1


def Gr_janaf(T,flag=False):
    """
    calculate the free energy of reaction FeO + 1/4O2 = FeO1.5 for pure endmembers by this study

    Params
    ----
    T :    list or scalar, temperature
    flag : whether or not include uncertainties
    
    Output
    ----
    G : energy list or scalar
    ----
    
    
    Note: It is a combination of Gr_janaf_new and Gr_janaf_old.
    Refer to the previous versions of `tools` of how to we derive the coefficients a,b,c,d based on these two old methods
    Several ways to improve it in the
    """
    a = uct.ufloat(-331035.9211346371,1.72e+02)
    b = uct.ufloat(-190.3795512883899,0.484)
    c = uct.ufloat(14.785873706952849,0.0859)
    d = uct.ufloat(-0.0016487959655627517,44e-6)
    e = uct.ufloat(9348044.389346942,1.22e+03)
    f = uct.ufloat(10773.299613088355,1.44)
    G = cal_G(T,a,b,c,d,e,f)

    if flag:
        return G
    else:
        try:
            return [uct.nominal_value(i) for i in G]
        except:
            return uct.nominal_value(G)
        
def cal_G(T,a,b,c,d,e,f):
    """
    """
    G = a+b*T+c*T*np.log(T) + d*(T**2) + e/T + f*(T**.5)
    return G


def cal_G_dz(T,a=-331035.9211346371,b=-190.3795512883899,c=14.785873706952849,d=-0.0016487959655627517,e=9348044.389346942,f=10773.299613088355):
    """
    """
    a = -331035.9211346371
    b = -190.3795512883899
    c = 14.785873706952849
    d = -0.0016487959655627517
    e = 9348044.389346942
    f = 10773.299613088355
    G = a+b*T+c*T*np.log(T) + d*(T**2) + e/T + f*(T**.5)
    return G