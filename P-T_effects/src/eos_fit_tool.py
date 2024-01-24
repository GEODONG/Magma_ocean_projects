#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 16:16:53 2019

This toolbox is originally from tools.py. It contains rich functions of EoS
fitting. I thus copy the original tools.py as this one.
It is not used in the oxidation paper, but can be very useful in the future.


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
def eos_fit(data,bmmodel,V0=0,P0=0,Kdp=0,free=False):
    """
    fit eos
    """
    params = Parameters()
    params.add('K0', value=12.9510408)
    params.add('Kp', value=4)
    if Kdp == 0:
        params.add('Kdp', value=-1.092579323,min=-1.4)
    else:
        params.add('Kdp', value=Kdp,vary=False)
        
    #hold as constant
    if V0 ==0 and P0==0:
        params.add('V0', value=max(data['V(A3)']), vary=False)
        params.add('P0', value=min(data['P(GPa)']), vary=False) 
    else:
        params.add('V0', value=V0, vary=False)
        params.add('P0', value=P0, vary=False)         
    result = bmmodel.fit(data['P(GPa)'],params,vol=data['V(A3)'])
    return result

#############
def adjustP(FeO_re,pivot=0.5):
    """
    adjust P, if it is within 0.8 GPa, 
    adjust since usually we can make it fit, 
    if not, make it approach to the best:
    """
    tempP = np.zeros(FeO_re['P(GPa)'].shape)
    for i in range(len(FeO_re['P(GPa)'])):
        diff_P = FeO_re.iloc[i]['P(GPa)'] - FeO_re.iloc[i]['best_P']
        if np.abs(diff_P) <= pivot:
            tempP[i]= FeO_re.iloc[i]['best_P']
        elif FeO_re.iloc[i]['P(GPa)'] > FeO_re.iloc[i]['best_P']:
            tempP[i] = FeO_re.iloc[i]['P(GPa)'] - pivot
        elif FeO_re.iloc[i]['P(GPa)'] < FeO_re.iloc[i]['best_P']:
            tempP[i] = FeO_re.iloc[i]['P(GPa)'] + pivot
    FeO_re['ajudst_P'] =  tempP
    return FeO_re

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

def set_par(data,bmmodel,ref = 3000,material='FeO_re',free=False,static=False,V0vary=True,Kdp = 0):
    """
    set parameters of fitting
    """
    params = Parameters()
    if static:
        T = float(input("input temperatures: "))
        data = data.loc[T]
        params.add('a', value=0,vary=False)
        params.add('b', value=0,vary=False)
        params.add('c', value=0,vary=False)
    else:
        params.add('a', value=20)
        params.add('b', value=30)
        params.add('c', value=6)
        
    if 'd' in bmmodel.param_names and not static:
        params.add('d', value=6)
    elif 'd' in bmmodel.param_names and static:
        params.add('d', value=0,vary=False)   
    if (static and T == 3000) or (not static and ref == 3000):
        skiprow = inp.loc[material]['skiprows_write_3k']
    else:
        skiprow = inp.loc[material]['skiprows_write_4k']
    par = read_par(skiprow)
    print("par is,", par)
    if material[-2:] == 're':
        par_single = par.iloc[0]
    else:
        par_single = par.iloc[1]
        
    if not free and not static:
    #hold as constant
        params.add('K0',  value=par_single['K0'], vary=False)
        params.add('Kp',  value=par_single['Kp'], vary=False)
        params.add('Kdp', value=par_single['Kdp'], vary=False)
        params.add('V0',  value=par_single['V0'], vary=False)
    elif free or static:
        params.add('K0',  value=par_single['K0'], min=1,vary=True)
        params.add('Kp',  value=par_single['Kp'], min = 0,vary=True)
        if Kdp == 0:
            params.add('Kdp', value=par_single['Kdp'], min=-4, vary=True)
        else:
            params.add('Kdp', value=par_single['Kdp'], min=Kdp, vary=True)            
        params.add('V0',  value=par_single['V0'], vary=V0vary)
    if static:
        params.add('P0',  value=0, vary=False)
        params.add('T0',  value=T, vary=False)  
    else:
        params.add('P0',  value=par_single['P0'], vary=False)
        params.add('T0',  value=par_single['T0'], vary=False) 
    return params
        
def eos_fit_th(data,bmmodel,ref = 3000,material = 'FeO_re',free=False,static=False,V0vary=True,Kdp = 0):
    """
    fit eos
    """
    params = set_par(data,bmmodel,ref = ref,material = material,free=free,static=static,V0vary=V0vary,Kdp=Kdp)   
    result = bmmodel.fit(data['P(GPa)'],params,vol=data['V(A3)'],T = data['T(K)'])
    return result

def eos_fit_th_adjust(data,bmmodel,ref = 3000,material = 'FeO_re',free=False,static=False,V0vary=True,Kdp = 0):
    """
    fit eos
    """
    print("here")
    print("material is",material)
    params = set_par(data,bmmodel,ref = ref,material = material,free=free,static=static,V0vary=V0vary,Kdp=Kdp)   
    result = bmmodel.fit(data['ajudst_P'],params,vol=data['V(A3)'],T = data['T(K)'])
    return result
############################## load data ##############################
def load_data(sheet_name = 'melt6_FeO',material = 'FeO_re'):
    FeO = pd.read_excel(xlsx,sheet_name=sheet_name,
                          usecols = list(range(inp.loc[material]['cols_start'],
                                               inp.loc[material]['cols_end'])), 
                          index_col = 0,
                          skiprows = list(range(inp.loc[material]['skiprows_read'])), 
                          nrows = inp.loc[material]['nrows'],).dropna()
    return FeO

temp  = pd.read_excel(xlsx,sheet_name='fitted',
                          usecols = list(range(8)), 
                          index_col = 0,
                          skiprows = list(range(70)), 
                          nrows = 9).dropna()
inp = temp.astype('int64')
############################## load data ##############################


def extract_th(result,index = ['FeO'],T0=4000):
    """
    extract the result from fitting to ['V0', 'P0', 'K0', 'Kp', 'Kdp', 'T0','redchi']
    """
    key_list = ['a','b','c','V0', 'P0', 'K0', 'Kp', 'Kdp', 'T0','redchi']
    temp     = result.params.valuesdict()
    temp2    = pd.DataFrame(temp,index=index)
#    temp2['T0']     = T0
    temp2['redchi'] = result.redchi 
    return temp2[key_list]

def extract_th_err(result,index = ['FeO'],T0=4000):
    """
    extract the result from fitting to ['V0', 'P0', 'K0', 'Kp', 'Kdp', 'T0','redchi']
    """
    key_list = ['a','b','c','V0', 'P0', 'K0', 'Kp', 'Kdp', 'T0',
                'aerr','berr','cerr','V0err','K0err','Kperr','Kdperr','redchi']
    temp     = result.params.valuesdict()
    temp2    = pd.DataFrame(temp,index=index)
#    temp2['T0']     = T0
    err = list(result.params.values())
    for ele in err:
        ind = ele.name+'err'
        temp2[ind] = ele.stderr
    temp2['redchi'] = result.redchi 
    return temp2[key_list]

def make_df(px,vx):
    temp = np.zeros((len(px),2))
    temp[:,0] = px
    temp[:,1] = vx
    df = pd.DataFrame(temp,columns=['P(GPa)','V(A3)'])
    return df


##################

############################## plot ##############################
import matplotlib.pyplot as plt


#def extrapolate(result_FeO_re,result_FeO_ox,vx,T):
#    p_FeO_re = result_FeO_re.model.eval(result_FeO_re.params,vol=vx,T=T)
#    p_FeO_ox = result_FeO_ox.model.eval(result_FeO_ox.params,vol=vx,T=T)
#    return p_FeO_re,p_FeO_ox
    
def plot_raw_fit(result_FeO_re,result_FeO_ox,FeO_re,FeO_ox,vx,xmin,xmax,index = 'org'):
    """
    plot raw and fitted data
    """
    fit_index = result_FeO_re.data.index.unique()
    org_index = FeO_re.index.unique()
    if index == 'org':
        index = org_index
    else:
        index = fit_index 
    
    from vatic import plot_tool as tool
    fig,ax = plt.subplots(2,1,figsize=(6,7),sharex=True,gridspec_kw={'height_ratios': (3,1)})
    for T in index:
        dfit_re = FeO_re['P(GPa)'].loc[T]-result_FeO_re.model.eval(result_FeO_re.params,vol=FeO_re['V(A3)'].loc[T],T=T)
        dfit_ox = FeO_ox['P(GPa)'].loc[T]-result_FeO_ox.model.eval(result_FeO_ox.params,vol=FeO_ox['V(A3)'].loc[T],T=T)
        re,ox,_,_   = extrapolate(result_FeO_re,result_FeO_ox,vx,T)
        re_line = ax[0].plot(re,vx,'-',label=str(T) + ' re')       
        
        ax[0].plot(ox,vx,'--',color=re_line[0].get_color(),label=str(T) + ' ox')
        ax[0].plot(FeO_re.loc[T]['P(GPa)'],FeO_re.loc[T]['V(A3)'],'o',color=re_line[0].get_color(),label='')
        ax[0].plot(FeO_ox.loc[T]['P(GPa)'],FeO_ox.loc[T]['V(A3)'],'^',color=re_line[0].get_color(),markerfacecolor='w',label='')

        ax[1].plot(FeO_re.loc[T]['P(GPa)'],dfit_re.loc[T],'o',color=re_line[0].get_color(),label= str(T) + ' re')
        ax[1].plot(FeO_ox.loc[T]['P(GPa)'],dfit_ox.loc[T],'^',color=re_line[0].get_color(),markerfacecolor='w',label= str(T) + ' ox')
    ax[0].set_ylabel('V (A^3)')
    ax[1].set_xlabel('P (GPa)'); ax[1].set_ylabel('Pexp - Ppre')
    ax[0].legend()
    ax[0].set_xlim([xmin,xmax]); ax[1].set_xlim([xmin,xmax])
    ax[0].grid(True); ax[1].grid(True)
    tool.autoscale_y(ax[0])
    tool.autoscale_y(ax[1]) 

################## Two different methods to calculate dV ######################
def extrapolate(result_FeO_re,result_FeO_ox,vx,T):
    p_FeO_re = result_FeO_re.model.eval(result_FeO_re.params,vol=vx,T=T)    
    p_FeO_ox = result_FeO_ox.model.eval(result_FeO_ox.params,vol=vx,T=T)
    minP_re = np.argmin(p_FeO_re)
    minP_ox = np.argmin(p_FeO_ox)
    if minP_re < len(p_FeO_re):
        print("The min P of reduced at {:2} K is around {:2.4} GPa, corresponding V is {:2.2} A^3".format(T, p_FeO_re[minP_re], vx[minP_re]))
    if minP_ox < len(p_FeO_ox):
        print("The min P of oxidized at {:2} K is around {:2.4} GPa, corresponding V is {:2.2} A^3".format(T, p_FeO_ox[minP_ox], vx[minP_ox]))
#    print("p_FeO_re is",p_FeO_re)
#    print("p_FeO_ox is",p_FeO_ox)
    return p_FeO_re,p_FeO_ox,minP_re,minP_ox

def cal_dV(result_FeO_re,result_FeO_ox,px,T,Fe_number):
    V0_re = result_FeO_re.best_values['V0']
    V0_ox = result_FeO_ox.best_values['V0']
    vx    = np.linspace(V0_re*0.2,V0_ox*1.8,1000)
    pre,pox,minP_re,minP_ox = extrapolate(result_FeO_re,result_FeO_ox,vx,T)
    if min(px)<pre[minP_re] or min(px)<pox[minP_ox]:
        print("min of input P is smaller the plausible input, ", min(pre[minP_re],pox[minP_ox]))

    f_FeO_re = interp1d(pre[:minP_re],vx[:minP_re])
    f_FeO_ox = interp1d(pox[:minP_ox],vx[:minP_ox])
    dV = (f_FeO_ox(px) - f_FeO_re(px))/Fe_number
    return dV


def plot_dV_fu(result_FeO_re,result_FeO_ox,px,Fe_number,xmin=0,xmax=80,index = 'org'):
    """
    plot dV per formula unit
    """
    fit_index = result_FeO_re.data.index.unique()
    
    if index == 'org':
        tarray = fit_index
    else:
        tmp = input("input temperatures (e.g., 2500 3000 4000): ")
        t = tmp.split()
        tarray = np.zeros(len(t))
        for i in range(len(t)):
            tarray[i] = float(t[i])
        
    from vatic import plot_tool as tool
    fig,ax = plt.subplots(1,2,sharey=False,figsize=(8,4))
    for T in tarray:
#        re,ox = extrapolate(result_FeO_re,result_FeO_ox,vx,T)
        dV    = cal_dV(result_FeO_re,result_FeO_ox,px,T,Fe_number)
        ax[0].plot(px,dV,label=str(T))
        ax[1].plot(px,dV*1e-30*1e6*6.022e23,label=str(T))
    ax[0].set_xlabel('P (GPa)')
    ax[0].set_ylabel('V (A^3)')
    ax[0].grid(True)
    ax[0].legend()
    ax[0].set_xlim([xmin,xmax])   
    ax[1].set_xlabel('P (GPa)')
    ax[1].set_ylabel('V (cm^3/mol)')
    ax[1].grid(True)
    ax[1].set_xlim([xmin,xmax])
    ax[1].legend()
    tool.autoscale_y(ax[0])
    tool.autoscale_y(ax[1])


def plot_raw_fit_k18(result_FeO_re,result_FeO_ox,FeO_re,FeO_ox,vx,Si_number,Fe_content,xmin,xmax,index = 'org',Fe_number=None):
    """
    plot raw and fitted data
    """
    fit_index = result_FeO_re.data.index.unique()
    org_index = FeO_re.index.unique()
    if index == 'org':
        index = org_index
    else:
        index = fit_index 
    from vatic import plot_tool as tool
    from melt_eos import MgFeSiO3_pt2v, MgFeO_pt2v
    fig,ax = plt.subplots(1,1,figsize=(8,6),sharex=False)

    for T in index:
        re,ox,_,_  = extrapolate(result_FeO_re,result_FeO_ox,vx,T)   
        p_re   = re[re>=0]  #MgFeSiO3_pt2v can only work for positive P
        if Si_number > 0:
            v_k18  = MgFeSiO3_pt2v(Fe_content,p_re,T=T,spin=0)[-1]*Si_number
        else:
            v_k18  = MgFeO_pt2v(Fe_content,p_re,T=T,spin=0)[-1]*Fe_number
        re_line = ax.plot(re,vx,'-',label=str(T) + ' re')        
        ax.plot(p_re,v_k18,':',color=re_line[0].get_color(),label=str(T) + ' re k18')        
        ax.plot(ox,vx,'--',color=re_line[0].get_color(),label=str(T) + ' ox')
        ax.plot(FeO_re.loc[T]['P(GPa)'],FeO_re.loc[T]['V(A3)'],'o',color=re_line[0].get_color(),label='')
        ax.plot(FeO_re.loc[T]['P(GPa)'],FeO_re.loc[T]['V(A3)'],'o',color=re_line[0].get_color(),label='')
        ax.plot(FeO_ox.loc[T]['P(GPa)'],FeO_ox.loc[T]['V(A3)'],'^',color=re_line[0].get_color(),markerfacecolor='w',label='')

    ax.grid(True); ax.grid(True)
    ax.set_xlabel('P (GPa)'); ax.set_ylabel('V (A^3)')
    ax.legend()
    ax.set_xlim([xmin,xmax]);
    tool.autoscale_y(ax)


def plot_raw_k18(FeO,FeO1p5,p_FeO,p_FeO1p5,vx,T,FeO_karki,px,formula_unit,xmin=-3,xmax=80):
    from vatic import plot_tool as tool
    fig,ax = plt.subplots()
    ax.plot(FeO.loc[T]['P(GPa)'],FeO.loc[T]['V(A3)'],'ro',label='reduced')
    ax.plot(FeO1p5.loc[T]['P(GPa)'],FeO1p5.loc[T]['V(A3)'],'bo',label='oxidized')
    ax.plot(p_FeO,vx,'r',label='reduced')
    ax.plot(p_FeO1p5,vx,'b',label='oxidized')
    ax.plot(px,FeO_karki*formula_unit,'r--',label='karki18')
    ax.set_xlim([xmin,xmax])
    ax.legend()
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Volume (' + r'$\AA^3$'+')')
    tool.autoscale_y(ax)
    plt.grid(True)
    plt.show()

def plot_dV_k18(dV,dV_karki,px,xmin=0,xmax=80):
    plt.figure()
    plt.subplot(121)
    plt.plot(px,dV,'-',label='This study')
    plt.plot(px,dV_karki,'--',label='Karki18')
    plt.xlabel('P (GPa)')
    plt.ylabel('V (A^3)')
#    plt.ylim([4,10])
    plt.xlim([xmin,xmax])
    plt.grid(True)
    plt.legend()
    plt.subplot(122)
    plt.plot(px,dV*1e-30*1e6*6.022e23,'-',label='This study')
    plt.plot(px,dV_karki*1e-30*1e6*6.022e23,'--',label='Karki18')
    plt.xlabel('P (GPa)')
    plt.ylabel('V (cm^3/mol)')
    plt.legend()
    plt.grid(True)

def plot_e(FeO_re,FeO_ox,xmin,xmax):
    """
    plot energy
    """
    index = FeO_re.index.unique()    
    from vatic import plot_tool as tool
    fig,ax = plt.subplots(1,1,figsize=(6,4))
    for T in index:        
        re_line = ax.plot(FeO_re.loc[T]['P(GPa)'],FeO_re.loc[T]['E(eV)'],'o',label='re '+str(T))
        ax.plot(FeO_ox.loc[T]['P(GPa)'],FeO_ox.loc[T]['E(eV)'],'^',color=re_line[0].get_color(),markerfacecolor='w',label='ox '+str(T))
    ax.set_ylabel('V (A^3)')
    ax.set_xlabel('P (GPa)'); ax.set_ylabel('E (eV)')
    ax.legend()
    ax.set_xlim([xmin,xmax]); 
    ax.grid(True);
    tool.autoscale_y(ax)
############################## plot ##############################

skiprow_dic = {'FeO_4k':9, 'Fe_25_4k':12, 'Fe_25_3k':15,  
               'Fe_12p5_4k':18, 'Fe_12p5_3k':21,
               'Fe_6p25_4k':24, 'Fe_6p25_3k':27,
               'previous':45}
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
    dv: v_ox - v_re
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
    Gaillard does not directly gives the way of calcualting the Gr, but it gives the thermoydnamics of several intermediate reactions. Gr5 means we use the combination of (1) and (2). There are many other ways, Gr3,Gr4... which are deprecated because they are way too off
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
    
    
    Note: It is a combination of Gr_janaf_new and Gr_janaf_old. Refer to the previous versions of `tools` of how to we derive the coefficients a,b,c,d based on these two old methods
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
