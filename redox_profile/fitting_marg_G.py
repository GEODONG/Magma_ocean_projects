import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

# Define the gas constant
R = 8.314

# read in the data
df = pd.read_excel('Supplementary_Excel.xlsx', sheet_name='Table S1')

# Prepare the data
T = df['T (K)'].values
Fe2 = df['FeO (mol.%)'].values / 100
Fe3 = df['Fe2O3 (mol.%)'].values / 100 / 2
Mg = df['MgO (mol.%)'].values / 100
Al = df['Al2O3 (mol.%)'].values / 100 / 2
Si = df['SiO2 (mol.%)'].values / 100
Ca = df['CaO (mol.%)'].values / 100
Na = df['Na2O (mol.%)'].values / 100 / 2
K = df['K2O (mol.%)'].values / 100 / 2
Ph = df['P2O5 (mol.%)'].values / 100 / 2
Ti = df['TiO2 (mol.%)'].values / 100
ln_gamma_r = df['LNGAMMA'].values


#Plot the comparison
import matplotlib.pyplot as plt
plt.figure(figsize=(8,8))

# Prepare the data
T = df['T (K)'].values
Mg = df['MgO (mol.%)'].values / 100
Al = df['Al2O3 (mol.%)'].values / 100
Si = df['SiO2 (mol.%)'].values / 100
Ca = df['CaO (mol.%)'].values / 100
Na = df['Na2O (mol.%)'].values / 100
K = df['K2O (mol.%)'].values / 100 
Ph = df['P2O5 (mol.%)'].values / 100 
Ti = df['TiO2 (mol.%)'].values / 100
ln_gamma_r = df['LNGAMMA'].values

def G95_new(variables, *W):
    #T,Fe2,Fe3, Mg, Al, Si, Ca, Na, Ti ,K, Ph= variables
    return (W[0]*K + W[1]*Mg + W[2]*Al + W[3]*Si +
            W[4]*Ca + W[5]*Na  + W[6]*Ti + W[7]*Ph+ + W[8]*Fe2 + W[9]*Fe3
            +W[10]*Si*Ti+W[11]*Si*Al+W[12]*Si*Fe3+W[13]*Si*Fe2+W[14]*Si*Mg+W[15]*Si*Ca
            +W[16]*Si*Na+W[17]*Si*K+W[18]*Si*Ph
            +W[19]*Ti*Al+W[20]*Ti*Fe3+W[21]*Ti*Fe2+W[22]*Ti*Mg+W[23]*Ti*Ca
            +W[24]*Ti*Na+W[25]*Ti*K+W[26]*Ti*Ph
            +W[27]*Al*Fe3+W[28]*Al*Fe2+W[29]*Al*Mg+W[30]*Al*Ca
            +W[31]*Al*Na+W[32]*Al*K+W[33]*Al*Ph
            +W[34]*Fe3*Fe2+W[35]*Fe3*Mg+W[36]*Fe3*Ca
            +W[37]*Fe3*Na+W[38]*Fe3*K+W[39]*Fe3*Ph
            +W[40]*Fe2*Mg+W[41]*Fe2*Ca
            +W[42]*Fe2*Na+W[43]*Fe2*K+W[44]*Fe2*Ph
            +W[45]*Mg*Ca+W[46]*Mg*Na+W[47]*Mg*K+W[48]*Mg*Ph
            +W[49]*Ca*Na+W[50]*Ca*K+W[51]*Ca*Ph
            +W[52]*Na*K+W[53]*Na*Ph+W[54]*K*Ph
            ) / R / T
           
 # Combine the independent variables
variables = (T, K,Mg, Al, Si, Ca, Na, Ti, Ph,Fe2,Fe3,
             Si*Ti,Si*Al,Si*Fe3,Si*Fe2,Si*Mg,Si*Ca,Si*Na,Si*K,Si*Ph,
             Ti*Al,Ti*Fe3,Ti*Fe2,Ti*Mg,Ti*Ca,Ti*Na,Ti*K,Ti*Ph,
             Al*Fe3,Al*Fe2,Al*Mg,Al*Ca,Al*Na,Al*K,Al*Ph,
             Fe3*Fe2,Fe3*Mg,Fe3*Ca,Fe3*Na,Fe3*K,Fe3*Ph,
             Fe2*Mg,Fe2*Ca,Fe2*Na,Fe2*K,Fe2*Ph,
             Mg*Ca,Mg*Na,Mg*K,Mg*Ph,
             Ca*Na,Ca*K,Ca*Ph,
             Na*K,Na*Ph,K*Ph)
initial_guess = [1.0] * 55
# Perform the curve fitting
optimal_parameters, covariance = curve_fit(
    G95_new, 
    variables, 
    ln_gamma_r, initial_guess)

print("Optimal parameters (W):", optimal_parameters)
#print("Covariance:", covariance)

# Calculate the standard deviation of the parameters
stdev = np.sqrt(np.diag(covariance))
#print("Standard deviation:", stdev)

# Calculate the R^2 value
residuals = ln_gamma_r - G95_new(variables, *optimal_parameters)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ln_gamma_r - np.mean(ln_gamma_r))**2)
r_squared = 1 - (ss_res / ss_tot)
print("R^2:", r_squared)

# Calculate the RSME value
n = len(ln_gamma_r)
rmse = np.sqrt(ss_res / n)
print("RMSE:", rmse)
plt.scatter(ln_gamma_r, G95_new(variables, *optimal_parameters), marker='^',facecolor = 'C0',edgecolor='C0',s=50, alpha=0.5, label= 'No truncation\n'+'RMSE=%.3f, ' % (rmse)+r'$R^2$'+'=%.3f'%(r_squared))

     


plt.xlabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Reference value', fontdict={'family': 'Arial', 'size': 24})
plt.ylabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Predicted value', fontdict={'family': 'Arial', 'size': 24})
plt.xlim(-1, 4)
plt.ylim(-1, 4)
plt.tick_params(labelsize=28,labelfontfamily='Arial')
plt.minorticks_on()
plt.plot(np.linspace(-1,3,100),np.linspace(-1,3,100),color='black',linestyle='-')
plt.text(3, 3, '1:1 line', rotation=45,fontsize=16, fontdict={'family': 'Arial'})
plt.plot(np.linspace(3.5,4,10),np.linspace(3.5,4,10),color='black',linestyle='-')

plt.legend(loc='lower right', fontsize=16, fancybox=False, edgecolor='black')
plt.tight_layout()
plt.savefig('W_fittings_G.png', dpi=600)
plt.savefig('W_fittings_G.pdf')

plt.show()

