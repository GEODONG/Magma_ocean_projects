import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

# Define the gas constant
R = 8.314

#Plot the comparison
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
lb_size = 20
tick_size = 24
text_size = 14

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


def J(variables, *W):
    #T, Fe2, Fe3, Mg, Al, Si, Ca, Na, Ti = variables
    #T, Fe2, Fe3, Mg, Al, Si, Ca,Na, Ti = variables
    return (W[0]*(Fe2 - Fe3) + W[1]*Mg + W[2]*Al + W[3]*Ca +
          + W[4]*Na  + W[5]*K + W[6]*Ph) / R / T

# Combine the independent variables
variables = (T, Fe2-Fe3, Mg, Al, Ca, Na, K, Ph)
initial_guess = [1.0] * 7
# Perform the curve fitting
optimal_parameters, covariance = curve_fit(
    J,
    variables,
    ln_gamma_r, initial_guess)
print('J04 model')
print("Optimal parameters (W):", optimal_parameters)
#print("Covariance:", covariance)

# Calculate the standard deviation of the parameters
stdev = np.sqrt(np.diag(covariance))
#print("Standard deviation:", stdev)

# Calculate the R^2 value
residuals = ln_gamma_r - J(variables, *optimal_parameters)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ln_gamma_r - np.mean(ln_gamma_r))**2)
r_squared = 1 - (ss_res / ss_tot)
print("\n\nR^2:", r_squared)

# Calculate the RSME value
n = len(ln_gamma_r)
rmse = np.sqrt(ss_res / n)
print("RMSE:", rmse)


ax.scatter(ln_gamma_r, J(variables, *optimal_parameters), marker= 's', color='red', s=50, alpha=0.5, label='J04 refitted model \n'+'RMSE=%.3f\n' % (rmse)+r'$R^2$'+'=%.3f'%(r_squared))
ax.set_xlabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Reference value', fontdict={'family': 'Arial', 'size': lb_size})
ax.set_ylabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Predicted value', fontdict={'family': 'Arial', 'size': lb_size})
ax.set_xlim(-1, 4)
ax.set_ylim(-1, 4)
ax.tick_params(labelsize=tick_size,labelfontfamily='Arial')
ax.legend(loc='lower right',fontsize=text_size, fancybox=False, edgecolor='black')
ax.minorticks_on()
ax.plot(np.linspace(-1,3,100),np.linspace(-1,3,100),color='black',linestyle='-')
ax.text(-0.8, 3.7, '(a)',fontsize=text_size, fontdict={'family': 'Arial'})
ax.text(3, 3, '1:1 line', rotation=45,fontsize=text_size, fontdict={'family': 'Arial'})
ax.plot(np.linspace(3.55,4,10),np.linspace(3.55,4,10),color='black',linestyle='-')


def J_new(variables, *W):
    #T, Fe2, Fe3, Mg, Al, Si, Ca, Na, Ti = variables
    #T, Fe2, Fe3, Mg, Al, Si, Ca,Na, Ti = variables
    return (W[0]*(Fe2 - Fe3) + W[1]*Mg + W[2]*Al + W[3]*Si +
          + W[4]*Na  + W[5]*Ti +
            W[6]*Ca*Mg + W[7]*Si*Al) / R / T

# Combine the independent variables
variables = (T, Fe2-Fe3, Mg, Al, Si, Na, Ti, Ca*Mg, Si*Al)
initial_guess = [1.0] * 8
# Perform the curve fitting
optimal_parameters, covariance = curve_fit(
    J_new, 
    variables, 
    ln_gamma_r, initial_guess)
print('J04 model inferred from ML')
print("Optimal parameters (W):", optimal_parameters)
#print("Covariance:", covariance)

# Calculate the standard deviation of the parameters
stdev = np.sqrt(np.diag(covariance))
#print("Standard deviation:", stdev)

# Calculate the R^2 value
residuals = ln_gamma_r - J_new(variables, *optimal_parameters)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ln_gamma_r - np.mean(ln_gamma_r))**2)
r_squared = 1 - (ss_res / ss_tot)
print("\n\nR^2:", r_squared)

# Calculate the RSME value
n = len(ln_gamma_r)
rmse = np.sqrt(ss_res / n)
print("RMSE:", rmse)

ax2.scatter(ln_gamma_r, J_new(variables, *optimal_parameters), marker= 's', color='red', s=50, alpha=0.5, label='J04 Model inferred from ML\n'+'RMSE=%.3f\n' % (rmse)+r'$R^2$'+'=%.3f'%(r_squared))
ax2.set_xlabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Reference value', fontdict={'family': 'Arial', 'size': lb_size})
ax2.set_ylabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Predicted value', fontdict={'family': 'Arial', 'size': lb_size})
ax2.set_xlim(-1, 4)
ax2.set_ylim(-1, 4)
ax2.tick_params(labelsize=tick_size,labelfontfamily='Arial')
ax2.legend(loc='lower right',fontsize=text_size, fancybox=False, edgecolor='black')
ax2.minorticks_on()
ax2.plot(np.linspace(-1,3,100),np.linspace(-1,3,100),color='black',linestyle='-')
ax2.text(-0.8, 3.7, '(b)',fontsize=text_size, fontdict={'family': 'Arial'})
ax2.text(3, 3, '1:1 line', rotation=45,fontsize=text_size, fontdict={'family': 'Arial'})
ax2.plot(np.linspace(3.55,4,10),np.linspace(3.55,4,10),color='black',linestyle='-')


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

def H22(variables, *W):
    return (W[0]*K + W[1]*Mg  + W[2]*Si +
            W[3]*Ca + W[4]*Na  + W[5]*Ti +
            W[6]*Si*Mg + W[7]*Si*Al + W[8]*Ph) / R / T


# Combine the independent variables
variables = (T, K,Mg, Si, Ca, Na, Ti, Si*Mg,Si*Al, Ph)
initial_guess = [1.0] * 9
# Perform the curve fitting
optimal_parameters, covariance = curve_fit(
    H22, 
    variables, 
    ln_gamma_r, initial_guess)
print('H22 model')
print("Optimal parameters (W):", optimal_parameters)
#print("Covariance:", covariance)

# Calculate the standard deviation of the parameters
stdev = np.sqrt(np.diag(covariance))
#print("Standard deviation:", stdev)

# Calculate the R^2 value
residuals = ln_gamma_r - H22(variables, *optimal_parameters)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ln_gamma_r - np.mean(ln_gamma_r))**2)
r_squared = 1 - (ss_res / ss_tot)
print("R^2:", r_squared)

# Calculate the RSME value
n = len(ln_gamma_r)
rmse = np.sqrt(ss_res / n)
print("RMSE:", rmse)

ax3.scatter(ln_gamma_r, H22(variables, *optimal_parameters), facecolor = 'none',edgecolor='black',s=50, alpha=0.7, label='H22 refitted model \n'+'RMSE=%.3f\n' % (rmse)+r'$R^2$'+'=%.3f'%(r_squared))
ax3.set_xlabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Reference value', fontdict={'family': 'Arial', 'size': lb_size})
ax3.set_ylabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Predicted value', fontdict={'family': 'Arial', 'size': lb_size})
ax3.set_xlim(-1, 4)
ax3.set_ylim(-1, 4)
ax3.tick_params(labelsize=tick_size,labelfontfamily='Arial')
ax3.legend(loc='lower right',fontsize=text_size, fancybox=False, edgecolor='black')
ax3.minorticks_on()
ax3.plot(np.linspace(-1,3,100),np.linspace(-1,3,100),color='black',linestyle='-')
ax3.text(-0.8, 3.7, '(c)',fontsize=text_size, fontdict={'family': 'Arial'})
ax3.text(3, 3, '1:1 line', rotation=45,fontsize=text_size, fontdict={'family': 'Arial'})
ax3.plot(np.linspace(3.55,4,10),np.linspace(3.55,4,10),color='black',linestyle='-')

def H22_new(variables, *W):
    return ( W[0]*Mg  + W[1]*Si + W[2]*Al +
             + W[3]*Na  + W[4]*Ti +
            W[5]*Ca*Mg + W[6]*Si*Al +W[7]*Fe2 + W[8]*Fe3) / R / T


# Combine the independent variables
variables = (T, Mg, Si, Al, Na, Ti, Ca*Mg, Si*Al, Fe2, Fe3)
initial_guess = [1.0] * 9
# Perform the curve fitting
optimal_parameters, covariance = curve_fit(
    H22_new, 
    variables, 
    ln_gamma_r, initial_guess)
print('H22 model inferred from ML')
print("Optimal parameters (W):", optimal_parameters)
#print("Covariance:", covariance)

# Calculate the standard deviation of the parameters
stdev = np.sqrt(np.diag(covariance))
#print("Standard deviation:", stdev)

# Calculate the R^2 value
residuals = ln_gamma_r - H22_new(variables, *optimal_parameters)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ln_gamma_r - np.mean(ln_gamma_r))**2)
r_squared = 1 - (ss_res / ss_tot)
print("R^2:", r_squared)

# Calculate the RSME value
n = len(ln_gamma_r)
rmse = np.sqrt(ss_res / n)
print("RMSE:", rmse)

ax4.scatter(ln_gamma_r, H22_new(variables, *optimal_parameters), facecolor = 'none',edgecolor='black',s=50, alpha=0.7, label='H22 Model inferred from ML\n'+'RMSE=%.3f\n' % (rmse)+r'$R^2$'+'=%.3f'%(r_squared))
ax4.set_xlabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Reference value', fontdict={'family': 'Arial', 'size': lb_size})
ax4.set_ylabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Predicted value', fontdict={'family': 'Arial', 'size': lb_size})
ax4.set_xlim(-1, 4)
ax4.set_ylim(-1, 4)
ax4.tick_params(labelsize=tick_size,labelfontfamily='Arial')
ax4.legend(loc='lower right',fontsize=text_size, fancybox=False, edgecolor='black')
ax4.minorticks_on()
ax4.plot(np.linspace(-1,3,100),np.linspace(-1,3,100),color='black',linestyle='-')
ax4.text(-0.8, 3.7, '(d)',fontsize=text_size, fontdict={'family': 'Arial'})
ax4.text(3, 3, '1:1 line', rotation=45,fontsize=text_size, fontdict={'family': 'Arial'})
ax4.plot(np.linspace(3.55,4,10),np.linspace(3.55,4,10),color='black',linestyle='-')

"""

plt.xlabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Reference value', fontdict={'family': 'Arial', 'size': lb_size})
plt.ylabel(r'$\ln{\frac{\gamma_{\mathrm{FeO_{1.5}}}}{\gamma_\mathrm{FeO}}}$ Predicted value', fontdict={'family': 'Arial', 'size': lb_size})
plt.xlim(-1, 4)
plt.ylim(-1, 4)
plt.tick_params(labelsize=tick_size,labelfontfamily='Arial')
plt.minorticks_on()
plt.plot(np.linspace(-1,3,100),np.linspace(-1,3,100),color='black',linestyle='-')
plt.text(3, 3, '1:1 line', rotation=45,fontsize=text_size, fontdict={'family': 'Arial'})
plt.plot(np.linspace(3.5,4,10),np.linspace(3.5,4,10),color='black',linestyle='-')

plt.legend(loc='lower right',fontsize=text_size, fancybox=False, edgecolor='black')
"""
plt.tight_layout()
plt.savefig('W_fittings_new.png', dpi=600)
plt.savefig('W_fittings_new.pdf')

plt.show()

