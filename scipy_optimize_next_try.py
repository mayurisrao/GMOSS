import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

#The docstring lines below are working piece of code.

frequency = np.array([22,45,150,408,1420,23000]) #x_values
b_temp = np.array([5.49282937e+04, 1.04875149e+04, 5.57289208e+02, 4.39957717e+01, 4.37429422e+00, 1.07386277e-03]) #y_values

#Defining the function that I want to fit
def concave_func(x, C_1, C_2, alpha_one, alpha_two, I_x, nu_t, T_e): 
    one = x**(-alpha_one)
    two = (C_2/C_1)*(x**(-alpha_two))
    three = I_x*(x**-2.1)
    expo = np.exp(-1*((nu_t/x)**2.1))
    eqn_one = C_1*(one + two + three)*expo
    eqn_two = T_e*(1 - expo)
    return eqn_one + eqn_two

#Defining chi_square function
def chisq(params, xobs, yobs):
    ynew = concave_func(xobs, *params)
    #yerr = np.sum((ynew- yobs)**2)
    yerr = np.sum(((yobs- ynew)/ynew)**2)
    print(yerr)
    return yerr

result = minimize(chisq, [1,2,1,1,1,1,1], args = (frequency,b_temp),  method = 'Nelder-Mead', options = {'disp' : True, 'maxiter': 10000})
#ynw3 = concave_func(frequency, *result.x)
x = np.linspace(-300,24000,1000)
plt.plot(x,concave_func(x, *result.x) )
print(result.x)
print(result)
#plt.plot(frequency, ynw3)
plt.plot(frequency, b_temp, 'r*')
plt.grid()
plt.show()

'''
frequency = np.array([22,45,150,408,1420,23000]) #x_values
b_temp = np.array([5.49282937e+04, 1.04875149e+04, 5.57289208e+02, 4.39957717e+01, 4.37429422e+00, 1.07386277e-03]) #y_values
log_frequency = np.log10(frequency)
log_b_temp = np.log10(b_temp)
#Defining the function that I want to fit
def concave_func(x, C_1, C_2, alpha_one, alpha_two, I_x, nu_t, T_e): 
    one = x**(-alpha_one)
    two = (C_2/C_1)*(x**(-alpha_two))
    three = I_x*(x**-2.1)
    expo = np.exp(-1*((nu_t/x)**2.1))
    eqn_one = C_1*(one + two + three)*expo
    eqn_two = T_e*(1 - expo)
    return eqn_one + eqn_two

#Defining chi_square function
def chisq(params, xobs, yobs):
    ynew = concave_func(xobs, *params)
    yerr = np.sum(((yobs- ynew)/ynew)**2)
    print(yerr)
    return yerr

result = minimize(chisq, [2,2,1,1,1,1,1], args = (log_frequency,log_b_temp), method = 'Nelder-Mead', options = {'disp' : True, 'maxiter': 10000})
#ynw3 = concave_func(log_frequency, *result.x)
print(result.x)
x = np.linspace(min(log_frequency),max(log_frequency),1000)
plt.plot(x,concave_func(x, *result.x) )
#plt.plot(log_frequency, ynw3)
plt.plot(log_frequency, log_b_temp, 'r*')
print(result)
plt.show()
'''