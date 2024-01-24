#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 28 11:29:42 2023

@author: jiedeng
@contributor: Donghao
"""
import numpy as np
import src.tools as tl
import src.tools2 as tl2
import pandas as pd
import matplotlib.pyplot as plt
import uncertainties as uct


##################################################################################################################################################################
##################################################   Functions       #############################################################################################
##################################################################################################################################################################

R = 8.314
a,b,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti = tl2.unpack_par(tl2.W.loc['fit3'],False)


def Gr(T, flag=False, composition= 'Fe25'):

    if composition == 'Fe25':
        a = uct.ufloat(-1278219.9542471424,24115.49293819573/100)
        b = uct.ufloat(-10905.74745614666,104.18153341404528/100)
        c = uct.ufloat(1177.9056464610862,10.10480094734614/100)
        d = uct.ufloat(-0.13961103642216463,0.00044179951813194326/100)
        e = uct.ufloat(103093738.42517175,2992805.064752209/100)
        f = uct.ufloat(126986.78970194193,1769.4906381762617/100)

        G = cal_G(T, a, b, c, d, e, f)
    else:
        a = uct.ufloat(-1225780.9400738582,24079.730032439846/100)
        b = uct.ufloat(-10421.923156287665,104.02515056964849/100)
        c = uct.ufloat(1125.1830774419097,10.089609394931353/100)
        d = uct.ufloat(-0.1340137144591397,0.00044112393023898927/100)
        e = uct.ufloat(98803019.7325353,2988396.1966464696/100)
        f = uct.ufloat(121815.10012724593,1766.8503852696128/100)

        G = cal_G(T, a, b, c, d, e, f)
    if flag:
        return G
    else:
        try:
            return [uct.nominal_value(i) for i in G]
        except:
            return uct.nominal_value(G)


def cal_G(T, a, b, c, d, e, f):
    """
    """
    G = a + b * T + c * T * np.log(T) + d * (T ** 2) + e / T + f * (T ** .5)
    return G

def X_cal(dat,a,b,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti, W_Mg_Al,
 W_Mg_Si, W_Mg_Ca, W_Al_Si, W_Al_Ca, W_Si_Ca, method = 'Deng'):
    """
    from
    /Users/jiedeng/Google Drive/papers/paper10_oxidation/revision_Oct1619/figS5_fit_all_2_JANAF.py
    -----
    input  list
    Fe - >tot iron, Fe
    r -> (Fe3/XFe)
    """
#    Fe3_Fe2 = []
#    print('enter X_cal')
    T = dat['T(K)'].values
    PV    = dat['PV'].values
    PV_uct = np.array([uct.ufloat( dat['PV'].values[i],dat['PV_err'].values[i]) for i in range(len(PV))])
    if method == 'Deng':
        Gterm = tl.Gr_janaf(T,flag=True)/R/T
        Wterm = (W_Fe * (dat['Fe2'] - dat['Fe3']) +
                 W_Mg * dat['MgO'] + W_Si * dat['SiO2'] +
                 W_Al * dat['Al2O3'] + W_Ca * dat['CaO'] +
                 W_Na * dat['Na2O'] + W_K * dat['K2O'] +
                 W_Ph * dat['P2O5'] + W_Ti * dat['TiO2']) / R / T
    else:
        Gterm = Gr(T, flag=True, composition= method)/R/T
        Wterm = (W_Fe * (dat['Fe2'] - dat['Fe3']) +
                 W_Mg * dat['MgO'] + W_Si * dat['SiO2'] +
                 W_Al * dat['Al2O3'] + W_Ca * dat['CaO'] +
                 W_Na * dat['Na2O'] + W_K * dat['K2O'] +
                 W_Ph * dat['P2O5'] + W_Ti * dat['TiO2']+
                 W_Mg_Al* dat['MgO']*dat['Al2O3']+
                 W_Mg_Si* dat['MgO']* dat['SiO2']+
                 W_Mg_Ca* dat['MgO'] * dat['CaO']+
                 W_Al_Si* dat['Al2O3']* dat['SiO2']+
                 W_Al_Ca* dat['Al2O3'] * dat['CaO']+
                 W_Si_Ca* dat['SiO2'] * dat['CaO']) / R / T

    tmp = -Gterm + np.log(10**dat['fo2'])/4 - Wterm - PV_uct
    print(tmp[:3])
    Fe3_Fe2_out = [uct.umath.exp(i) for i in tmp]
    
#    Fe3_Fe2.append(Fe3_Fe2)
    return Fe3_Fe2_out

def load_data(sheet_name = '',skiprow = list(range(81)), usecol = list(range(18)), index_col = 0):            
    FeO = pd.read_excel(tl2.xlsx,sheet_name=sheet_name,
                          usecols = usecol, 
                          index_col = index_col,
                          skiprows = skiprow, 
                          nrows = 100,).dropna()
    
    FeO['Fe3/Fe2'] = FeO['Fe3/Fe']/(1-FeO['Fe3/Fe'])
    FeO['Fe3'] = FeO['FeO']*FeO['Fe3/Fe']
    FeO['Fe2'] = FeO['FeO']*(1-FeO['Fe3/Fe'])
    return FeO

"""

This is for DelG_12p5_new
a =  -1225780.9400738582
b =  -10421.923156287665
c =  1125.1830774419097
d =  -0.1340137144591397
e =  98803019.7325353
f =  121815.10012724593
a_err =  24079.730032439846
b_err =  104.02515056964849
c_err =  10.089609394931353
d_err =  0.00044112393023898927
e_err =  2988396.1966464696
f_err =  1766.8503852696128
This is for DelG_25
a =  -1278219.9542471424
b =  -10905.74745614666
c =  1177.9056464610862
d =  -0.13961103642216463
e =  103093738.42517175
f =  126986.78970194193
a_err =  24115.49293819573
b_err =  104.18153341404528
c_err =  10.10480094734614
d_err =  0.00044179951813194326
e_err =  2992805.064752209
f_err =  1769.4906381762617


"""



######################################################################################################################################################################################################################
#################################################################################   ploting       ################################################################################################################
######################################################################################################################################################################################################################

"""
Fe25 
 W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti, 

 -3.25614156e+05, 1.54767751e+03,-5.09341345e+04,5.35953071e+05,
 2.72834359e+05,-3.22110805e+04 ,-9.52835576e+04,0, 9.77265699e+04, 
 3.21794220e+05,  3.91257558e+05, -1.66959988e+06, -8.65206480e+05,
 -3.37696151e+06,  4.41975302e+05
 
0.9571956668994056


Fe12p5
W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti, 
-330699.47877531, 29523.22324146 ,-615143114, 491979.99196477 ,
305681.86964438, -36643.54727279,  -101697.60588497, 0, 88606.53279936 ,     
432992.27444423,   340229.3436809,  -1846613.3796063 ,  -783519.55820583,
 -3660061.45921339,   467417.04514842



0.9518259768958255




"""


O06 = load_data(sheet_name='O06',skiprow=None)
Z17 = load_data(sheet_name='Z17',skiprow=None)
A18 = load_data(sheet_name='A18',skiprow=None);
A18_with_corr = load_data(sheet_name='A18',skiprow=None,usecol=list(range(19)));
Ku2023 = load_data(sheet_name='Ku2023',skiprow=None)
#Reported_A19 = load_data(sheet_name='Reported_A19',skiprow=None,usecol=list(range(14)))

org = pd.concat([O06,Z17,A18],sort=False).fillna(0)
#org = pd.concat([Thornber80,Kress88,Kress91,Sack81,Kilinc83,Moore95,J04,high_P],'sort=False').fillna(0)

all_dat = org.iloc[:-8]
#all_dat = org
Fe3_Fe  = all_dat['Fe3/Fe']
Fe2_Fe  = 1-Fe3_Fe
Fe3_Fe2 = Fe3_Fe/Fe2_Fe

Fe3_Fe2_model  = X_cal(org,a,b,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti,0,0,0,0,0,0)
Fe3_Fe_model = [i/(1+i) for i in Fe3_Fe2_model] 
#Fe3_Fe_model   = Fe3_Fe2_model/(1+Fe3_Fe2_model)
Fe3_Fe_model_uct = [uct.std_dev(i) for i in Fe3_Fe_model]
Fe3_Fe_model_val = [uct.nominal_value(i) for i in Fe3_Fe_model]



New_model = pd.concat([O06,Z17,A18,Ku2023],sort=False).fillna(0)

all_dat = New_model.iloc[:-8]

Fe3_Fe_new_model  = all_dat['Fe3/Fe']
Fe2_Fe_new_model  = 1-Fe3_Fe
Fe3_Fe2_new_model = Fe3_Fe/Fe2_Fe

"""
Fe3_Fe2_model_new  = X_cal(New_model,a,b,-3.25614156e+05, 1.54767751e+03,-5.09341345e+04,5.35953071e+05,
 2.72834359e+05,-3.22110805e+04 ,-9.52835576e+04,0, 9.77265699e+04,
3.21794220e+05,  3.91257558e+05, -1.66959988e+06, -8.65206480e+05,
 -3.37696151e+06,  4.41975302e+05,method='Fe25')#,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti
"""
Fe3_Fe2_model_new  = X_cal(New_model,a,b,W_Fe,W_Mg,W_Si,W_Al,W_Ca,W_Na,W_K,W_Ph,W_Ti,0,0,0,0,0,0,method='Fe125')

"""
Fe3_Fe2_model_new  = X_cal(New_model,a,b, -329265.56357094,111970.0839,-60507.06644171,
  490360.43466179 ,299021.21556838 , -42360.80231379 ,
   -108944.31026724   , 0,  41957.64883924 ,
    182733.23181499,   236505.02448532, -1971060.75950741 , -769302.85660747,
  -3495173.53104364,   462589.49532712,method='Fe125')"""
Fe3_Fe_model_new = [i/(1+i) for i in Fe3_Fe2_model_new]

Fe3_Fe_model_uct_new = [uct.std_dev(i) for i in Fe3_Fe_model_new]
Fe3_Fe_model_val_new = [uct.nominal_value(i) for i in Fe3_Fe_model_new]

fig,ax = plt.subplots(2,1,figsize=(8.5,5.5),sharey=True)



ax[0].errorbar(Z17['P(GPa)'],Z17['Fe3/Fe'],yerr = Z17['uct']*2,fmt='s',color='C1',label='Z17',alpha=0.6)
ax[0].errorbar(O06['P(GPa)'],O06['Fe3/Fe'],yerr = O06['uct']*4,fmt='o',color='c',label='O06',alpha=0.5)
ax[0].errorbar(A18['P(GPa)'],A18['Fe3/Fe'],yerr = A18['uct'],fmt='^',color='C0',label='A19')
ax[0].errorbar(Ku2023['P(GPa)'],Ku2023['Fe3/Fe'],yerr = Ku2023['uct'],fmt='h',color='brown',label='K23',alpha=0.5)# added by Donghao
ax[0].errorbar(org['P(GPa)'],Fe3_Fe_model_val,yerr = Fe3_Fe_model_uct,fmt='d',color='k',label='D20',alpha=0.4)
ax[0].errorbar(New_model['P(GPa)'],Fe3_Fe_model_val_new,yerr = Fe3_Fe_model_uct_new,fmt='D',color='r',alpha = 0.7,label='This study')#plt.plot(org['P(GPa)'],Fe3_Fe_model,'*',label='This')

ax[0].set_xlabel(r'$P$ (GPa)'+"\n\n",fontsize=16,fontname= "Arial")
ax[0].set_xlim([0,30])
ax[0].set_ylim([0,1])# modified by Donghao

# change the temperature of A18 shown in her figure to her table, those two are not the same.
# refer to prev_sum for the changes
T_A18_table = A18['T(K)'].values
T_A18_table[3] = T_A18_table[3]#+50
T_A18_table[7:9] = T_A18_table[7:9] #+ 100
T_A18_table[-2]    = T_A18_table[-2] #+ 200
T_org_corr         = org['T(K)'].values 
T_org_corr[-11:]   =  T_A18_table

T_new_model_corr = New_model['T(K)'].values
#T_new_model_corr[-11:]   =  T_A18_table
ax[1].errorbar(Z17['T(K)'],Z17['Fe3/Fe'],yerr = Z17['uct']*2,fmt='s',color='C1',label='Z17',alpha=0.4)
ax[1].errorbar(O06['T(K)'],O06['Fe3/Fe'],yerr = O06['uct']*4,fmt='o',color='c',label='O06',alpha=0.5)
ax[1].errorbar(T_A18_table,A18['Fe3/Fe'],yerr = A18['uct'],fmt='^',color='C0',label='A19')
ax[1].errorbar(Ku2023['T(K)'],Ku2023['Fe3/Fe'],yerr = Ku2023['uct'],fmt='h',color='brown',label='Ku2023',alpha=0.5)# added by Donghao

ax[1].errorbar(T_org_corr,Fe3_Fe_model_val,yerr = Fe3_Fe_model_uct,fmt='d',color='k',label='D20',alpha=0.4)
ax[1].errorbar(T_new_model_corr,Fe3_Fe_model_val_new,yerr = Fe3_Fe_model_uct_new,fmt='D',color='r',alpha = 0.7,label='This Study')
ax[1].legend(loc = 'lower center', ncol=3 ,fancybox = False, edgecolor = 'black',fontsize=14,bbox_to_anchor=(0.5,-1))


ax[1].set_xlabel(r'$T$ (K)',fontsize=16,fontname= "Arial")
ax[1].set_ylabel(r'$\mathrm{Fe}^{3+} /\Sigma \mathrm{Fe}$',fontsize=15,fontname= "Arial")
ax[0].set_ylabel(r'$\mathrm{Fe}^{3+} /\Sigma \mathrm{Fe}$',fontsize=15,fontname= "Arial")
ax[1].set_xlim([1600,2900])
ax[1].set_ylim([0,1])#modified by Donghao
ax[0].minorticks_on()
ax[1].minorticks_on()
ax[0].text(0.25,0.9,'(a)',fontsize=16,fontname= "Arial")
ax[1].text(1615,0.85,'(b)',fontsize=16,fontname= "Arial")

ax[0].tick_params(axis='both', which='major', labelsize=15, labelfontfamily = "Arial") 
ax[1].tick_params(axis='both', which='major', labelsize=15, labelfontfamily = "Arial") 
plt.subplots_adjust(wspace=0.5, hspace=0.6)

fig.savefig('./P-T_effect.pdf',bbox_inches='tight')
fig.savefig('./P-T_effect.png',bbox_inches='tight',dpi=600)

