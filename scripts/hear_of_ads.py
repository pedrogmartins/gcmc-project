import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#This script extracts adsorption isotherms stored in arrays, and
#   calculates the isoteric heat of adsorption versus CO2 occupancy.
#   We fit the isotherm with a Langmuit dual-site adsorption
#   model from 3 different isotherms temperatures simultaneously to
#   extract the isosteric heat of adsoprtion according to [1].

# [1] Mason, Jarad A., et al. "Evaluating metal–organic frameworks
#    for post-combustion carbon dioxide capture via temperature swing
#    adsorption." Energy & Environmental Science 4.8 (2011): 3030-3040.

#Physical Constants
R = 8.3145e-3
T = 298

#Import data and concatenate pressures and concentrations
y1 = conc_298_avg
y2 = conc_313_avg
y3 = conc_323_avg
comboY = np.append(y1, np.append(y2, y3))

h = ps_set
comboX = np.append(h, np.append(h, h))

#Define initial needed functions

# Model at 298
def mod_298(p, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B):

        T = 298
        R = 8.3145*1e-3
        b_A = b_A0*np.exp(E_A/R/T)
        b_B = b_B0*np.exp(E_B/R/T)
        term_1 = (q_sat_A*b_A*p)/(1+b_A*p)
        term_2 = (q_sat_B*b_B*p)/(1+b_B*p)

        return term_1 + term_2

# Model at 313
def mod_313(p, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B):

        T = 313
        R = 8.3145*1e-3
        b_A = b_A0*np.exp(E_A/R/T)
        b_B = b_B0*np.exp(E_B/R/T)
        term_1 = (q_sat_A*b_A*p)/(1+b_A*p)
        term_2 = (q_sat_B*b_B*p)/(1+b_B*p)

        return term_1 + term_2

# Model at 323
def mod_323(p, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B):

        T = 323
        R = 8.3145*1e-3
        b_A = b_A0*np.exp(E_A/R/T)
        b_B = b_B0*np.exp(E_B/R/T)
        term_1 = (q_sat_A*b_A*p)/(1+b_A*p)
        term_2 = (q_sat_B*b_B*p)/(1+b_B*p)

        return term_1 + term_2

#Combinding all 3 models into a function combo for simulatanous fit
def comboFunc(comboData, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B):
    # single data set passed in, extract separate data
    extract1 = comboData[:len(y1)] # first data
    extract2 = comboData[len(y1):2*len(y2)] # second data
    extract3 = comboData[2*len(y2):] # second data

    result1 = mod_298(extract1, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B)
    result2 = mod_313(extract2, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B)
    result3 = mod_323(extract3, q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B)

    return np.append(result1, np.append(result2, result3))

#Define initial numerical guesses
initialParameters = np.array([1, 1e-6, 10, 1, 1e-6, 10])

#Curve fit the combined data to the combined function
fittedParameters, pcov = curve_fit(comboFunc, comboX, comboY, initialParameters, maxfev = 1400)

#Error minimizing constant in the model
q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B = fittedParameters

#Defining arrays of model prediction
ps_support = range(0, 60, 1)
y_fit_1 = mod_298(np.asarray(ps_support), q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B) # first data set, first equation
y_fit_2 = mod_313(np.asarray(ps_support), q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B) # second data set, second equation
y_fit_3 = mod_323(np.asarray(ps_support), q_sat_A, b_A0, E_A, q_sat_B, b_B0, E_B) # second data set, second equation

#Plotting data and model, and printing our parameters
plt.plot(comboX, comboY, 'D') # plot the raw data
plt.plot(ps_support, y_fit_1) # plot the equation using the fitted parameters
plt.plot(ps_support, y_fit_2) # plot the equation using the fitted parameters
plt.plot(ps_support, y_fit_3) # plot the equation using the fitted parameters
plt.rcParams['axes.facecolor'] = 'white'

#plt.xlim(0, 1)
plt.show()
print(fittedParameters)



#Isosteric heat of adsoprtion according to [1]

#Physical constant and algebraic operations
q = np.linspace(0.000001, 5)
T = 298
R = 8.3145*1e-3
b_A = b_A0*np.exp(E_A/R/T)
b_B = b_B0*np.exp(E_B/R/T)
alpha = (q_sat_A + q_sat_B - q)*b_A*b_B
beta = (q_sat_A - q)*b_A + (q_sat_B - q)*b_B
dbabbdt = -b_A0*b_B0*np.exp((E_A + E_B)/(R*T))*((E_A + E_B)/(R*T**2))
dadt = (q_sat_A + q_sat_B - q)*dbabbdt
dbadt = -b_A0*np.exp(E_A/R/T)*(E_A/R/T**2)
dbbdt = -b_B0*np.exp(E_B/R/T)*(E_B/R/T**2)
dbdt = (q_sat_A - q)*dbadt +  (q_sat_B - q)*dbbdt
der = 1/(np.sqrt(beta**2 + 4*alpha*q) - beta) * ( ( 1 / ( 2*np.sqrt(beta**2 + 4*alpha*q) ) ) * ( 2*beta*dbdt + 4*q*dadt) - (dbdt) ) - 1/alpha*dadt
Q = R*T**2*der



#Plotting isosteric heat of adsorption and DFT predicted binding energies
figure(figsize=(15, 12), dpi=80)
plt.rcParams['axes.facecolor'] = 'white'
plt.plot(q, Q, '-',  mfc='none', linewidth = 3)
#plt.errorbar(ps_setIV, ads_isotherm_mc[5][0], yerr = ads_isotherm_mc[10][1], fmt=',k', elinewidth = 0.8, capsize = 0.5);
plt.xlabel('CO$_2$ molecule / Mg metal site', size = 24)
plt.ylabel('|$Q_{st}$| (kJ/mol)', size = 24)
#plt.title('Isosteric Heat of Adsoprtion: CO$_2$ in MOF-274', size = 33)
plt.axvline(x = 1, c = "r", linestyle='dashed')
plt.axhline(y = 53.89, c = "k", linestyle='dotted')
#plt.axhline(y = 19.9, c = "k", linestyle='dotted')
#plt.axhline(y = 24, c = "#ff7f0e", linestyle='dotted')
#plt.axhline(y = 42, c = "#ff7f0e", linestyle='dotted')
plt.tick_params(axis = 'x', labelsize = 22.5)
plt.tick_params(axis = 'y', labelsize = 22.5)
plt.ylim(0, 65)
plt.xlim([0, 3])
plt.rcParams['axes.facecolor'] = 'white'
#plt.xscale('log')
plt.show()
