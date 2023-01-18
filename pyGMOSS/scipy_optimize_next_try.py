import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

#The docstring lines below are working piece of code.

frequency = np.array([22,45,150,408,1420,23000]) #x_values
#b_temp = np.array([5.49282937e+04, 1.04875149e+04, 5.57289208e+02, 4.39957717e+01, 4.37429422e+00, 1.07386277e-03]) #y_values
b_temp = [2.55080863e+04, 4.90777800e+03, 2.28984753e+02, 2.10842949e+01, 3.58631166e+00, 5.68716056e-04]

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

<<<<<<< HEAD
result = minimize(chisq, [1.0,2.0,1.0,1.0,1.0,1.0,1.0], args = (frequency,b_temp),  method = 'Nelder-Mead', options = {'disp' : True, 'maxiter': 10000})
=======
result = minimize(chisq, [1,2,2,2,1,1e6,8000], args = (frequency,b_temp),  method = 'Nelder-Mead', options = {'disp' : True, 'maxiter': 10000})
#result = minimize(chisq, [1,2,1,1,1,1,1], args = (frequency,b_temp),  method = 'Nelder-Mead', options = {'disp' : True, 'maxiter': 10000})
>>>>>>> ee6ea3822d4f8d856f86eae8972f3d33c0bf8e33
#ynw3 = concave_func(frequency, *result.x)
x = np.linspace(-300,24000,1000)
plt.yscale("log")
plt.xscale("log")
plt.plot(x,concave_func(x, *result.x) )
print(result.x)
print(result)
#plt.plot(frequency, ynw3)
<<<<<<< HEAD
plt.xscale("log")
plt.yscale("log")
plt.plot(frequency, b_temp, 'r*')
=======
plt.plot(frequency, b_temp, 'r*' )
plt.xlabel("log Frequency[MHz]")
plt.ylabel("log Temp[K]")
plt.title('log Temparature vs log Frequency')

>>>>>>> ee6ea3822d4f8d856f86eae8972f3d33c0bf8e33
plt.grid()
plt.savefig('the_plot_2060.png')


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

result = minimize(chisq, [2,2,1,1,1,1,1], args = (frequency,b_temp), method = 'Nelder-Mead', options = {'disp' : True, 'maxiter': 10000})
#ynw3 = concave_func(log_frequency, *result.x)
print(result.x)
x = np.linspace(10,24000,10)
plt.plot(x,concave_func(x, *result.x))
#plt.plot(log_frequency, ynw3)
plt.plot(frequency, b_temp, 'r*')
print(result)
#plt.yscale("log")
#plt.xscale("log")
#print(x)
#print(concave_func(frequency, *result.x))
plt.show()

'''
