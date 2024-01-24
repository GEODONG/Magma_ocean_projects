import numpy as np
import pandas as pd
import random
import  matplotlib.pyplot as plt
from scipy.optimize import  curve_fit

def func(T, a, b, c, d, e,f):
    G = a+b*T+c*T*np.log(T) + d*(T**2) + e/T + f*(T**.5)
    return G

def Cp(T, x,y):
    Cp_tot = x+y*T
    return Cp_tot

my_dataset = pd.read_excel('cv_cal.xlsx',#previous file is ''eos.xlsx
                           sheet_name='de_cp_with_5000K_new_data')
for fe_name in ['Fe_12p5','Fe_25']:
    data = my_dataset[my_dataset['Fe'] == fe_name]
    print(data['delta_Cp'])
    popt, pcov = curve_fit(Cp,data['Temperature'], data['delta_Cp'])
    #perr = np.sqrt(np.diag(pcov))
    print('This is for '+ fe_name)
    x = popt[0]
    y = popt[1]
    print('x = ',x)
    print('y = ',y)


width = 2

dataset2 = pd.read_excel('./copy_of_magma_ocean_redox_4_3 (copy).xlsx', sheet_name='Sheet2')
fig2 = plt.figure(figsize=(10,7))

subset1 = dataset2[dataset2['Temperature']<= 2300]

subset2 = dataset2[dataset2['Temperature']> 2300]

ax6 = fig2.add_subplot(1,1,1)
ax6.plot(dataset2['Temperature'], dataset2['del(G)']/1000, color='black',linestyle =':', linewidth=width, label='Deng et al., 2020')

ax6.plot(subset1['Temperature'], subset1['Jaysurya 2004']/1000, color='red',linestyle = '-', linewidth=width)

ax6.plot(subset1['Temperature'], subset1['Jaysurya 2004']/1000, color='red',linestyle = '-', linewidth=width*5,alpha =.6)
ax6.plot(subset2['Temperature'], subset2['Jaysurya 2004']/1000, color='red',linestyle = '--', linewidth=width, label='Jayasuriya et al., 2004')

ax6.plot(dataset2['Temperature'][:-300], dataset2['Gaillard 2003: Fe-Ir equilibrium'][:-300]/1000, color='#260b9e',linestyle='-.', linewidth=width,
         label='Gaillard et al., 2003')

for (name, clr,lab) in [('DelG_12p5_new', '#f58d42',r'$X_{\mathrm{FeO^T}}$ = 12.5 mol.%'),('DelG_25','#42b9f5',r'$X_{\mathrm{FeO^T}}$ = 25.0 mol.%')]:
    popt, pcov = curve_fit(func,dataset2['Temperature'], dataset2[name])
    perr = np.sqrt(np.diag(pcov))
    print('This is for '+ name)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    d = popt[3]
    e = popt[4]
    f = popt[5]

    G_vals = func(dataset2['Temperature'],a,b,c,d,e,f)

    print('a = ',a)
    print('b = ',b)
    print('c = ',c)
    print('d = ',d)
    print('e = ',e)
    print('f = ',f)

    a_err = perr[0]
    b_err = perr[1]
    c_err = perr[2]
    d_err = perr[3]
    e_err = perr[4]
    f_err = perr[5]

    print('a_err = ',perr[0])
    print('b_err = ',perr[1])
    print('c_err = ',perr[2])
    print('d_err = ',perr[3])
    print('e_err = ',perr[4])
    print('f_err = ',perr[5])

    ax6.plot(dataset2['Temperature'],dataset2[name]/1000,color = clr,linewidth = width,label = 'This study '+lab,linestyle = '-')
   
ax6.yaxis.set_ticks_position('left')
ax6.spines['left'].set_position(('data',1500))
ax6.set_xlabel(r'$T$ (K)',font = 'Arial',fontsize = 20)
ax6.set_ylabel(r'$\Delta$$G$$_r^0$ (kJ mol$^{-1})$',fontsize = 20)
ax6.legend(fancybox = False, edgecolor = 'black', prop={'family': 'Arial', 'size': 16})
ax_xlabel = ax6.get_xticklabels()
[x1_label_temp.set_fontname('Arial') for x1_label_temp in ax_xlabel]
ax_ylabel = ax6.get_yticklabels()
[y1_label_temp.set_fontname('Arial') for y1_label_temp in ax_xlabel]


plt.xlim(1500,5000)
plt.ylim(-50,1e2)
ax6.set_xticks(np.arange(1500,5100,100),minor=True)
ax6.set_yticks(np.arange(-40,100,5),minor=True)
plt.tick_params(axis='x',labelsize = 20,pad = 10)
plt.tick_params(axis='y',labelsize = 20,pad = 10)

plt.tight_layout()

plt.savefig("./Delta_G.pdf")
plt.savefig("./Delta_G_new.png",dpi = 600)

