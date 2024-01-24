import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize as op

A3_to_m3   = 1e-30
ev_to_J  =  1.602177*1e-19
N = 6.02*1e23
bar_width = 100

def model(x, b, c):
    return  c+b/(x**2)

my_dataset = pd.read_excel('cv_cal.xlsx',
                           sheet_name='de_cp_with_5000K_new_data')
fig,ax = plt.subplots(figsize = (10,9))

ax.errorbar(1452+(1854-1452)/2,(199.7/2-78.94), xerr = (1854-1452)/2,fmt= 'v', label = 'Richet and Bottinga, 1985', markersize = 10, color = 'black')#'R84')###new experiments
ax.errorbar(1200+(1850-1200)/2,(229/2-78.9), xerr =(1850-1200)/2,fmt= '^',label = 'Stebbins et. al., 1984' , markersize = 10, color = 'grey')#'S84')
ax.errorbar(1140+(1674-1140)/2,(240.9/2-78.8), xerr =(1674-1140)/2,fmt= 'o',label = 'Lange and Navrotsky, 1992', markersize = 10, color = 'black', alpha = 0.3)#'L91')

colors = ['b','g', 'r','c']
for (name, clr,mole,lab, mark) in [('Fe_12p5', 'red',8,r'$X_{\mathrm{FeO^T}}$ = 12.5 mol.$\%$','s'),('Fe_25','blue',4,r'$X_{\mathrm{FeO^T}}$ = 25.0 mol.$\%$','D')]:#,('Fe_25_new','c'),('Fe_12p5_new','orange'),('Fe_12p5_with_rechecked_ox','pink')]:
    my_data = my_dataset[my_dataset.Fe == name]

    b,c = op.curve_fit(model, my_data.Temperature, my_data.delta_Cp*mole)[0]
    print(b)
    print('\n'+name)
    print(c)
    x = np.arange(1400,5100,1)
    ax.plot(x,c+b/(x**2)
, color= clr, linestyle='-')
    ax.scatter(my_data.Temperature, my_data.delta_Cp*mole
     , color= clr, label = 'This study ' + lab, marker = mark, s=100, edgecolor = clr)
    ax.errorbar(my_data.Temperature, my_data.delta_Cp*mole, yerr = 3, fmt = 'none', color = clr)
   
ax.set_xlabel(r'$T$ (K)',fontsize = 22, fontname ='Arial')
ax.set_ylabel(r'$\Delta$'+r'$C_{P}$ (J '+'mol'+'$^{-1}$'+'K'+r'$^{-1}$'+')',fontsize = 22, fontname ='Arial')
ax.set_xlim(1000,5600)
ax.set_ylim(0,50)
ax.set_xticks(np.arange(1000,5600,1000))#, fontsize = 24)
ax.set_xticks(np.arange(1000,5600,100),minor = True)#, fontsize = 24)
ax.set_yticks(np.arange(0,50,2),minor = True)#, fontsize = 24)
plt.tick_params(axis='x',labelsize = 20,pad = 16)
plt.tick_params(axis='y',labelsize = 20,pad = 16)
ax_xlabel = ax.get_xticklabels()
[x1_label_temp.set_fontname('Arial') for x1_label_temp in ax_xlabel]
ax_ylabel = ax.get_yticklabels()
[y1_label_temp.set_fontname('Arial') for y1_label_temp in ax_xlabel]

ax.xaxis.labelpad = 10  # adjust the space to 10 points
ax.yaxis.labelpad = 10  # adjust the space to 10 points

ax.legend(ncol=1, loc="upper right", fontsize= 16,fancybox='false',edgecolor='k')
ax = plt.gca()


plt.savefig('./Delta_Cp.pdf')
plt.savefig('./Delta_Cp.png',dpi = 600)