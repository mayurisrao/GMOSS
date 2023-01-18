#apply corrections to 150,408,1420
#corrections_applied_to = np.array([150,408,1420])

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gamma
import scipy.special as special
import scipy.integrate as integrate
import math
import scipy as sp

#convex shape at pixel 36 and we will use that
frequency = np.array([22,45,150,408,1420,23000]) #x_values
frequency = np.array([np.float32(f*10**-3) for f in frequency])
#b_temp = [9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04] #y_values
b_temp = np.array([9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04])

cvel = 2.99792458e+08 # m s^-1
m_e = 9.1e-31
q_e = 1.6e-19
sin_alph = 1.0
Bmag = 1e-9 # Tesla == 10 micro-Gauss
#PI = 3.14159265358979
# FTOL = 1.0e-6
T_e = 8000.0 # Thermal electron temperature
# TCMB = 2.72548
# NHPIX = 3072
GSPAN = 100.0



#some important stuff
scale_gam_nu = (3.0*q_e*Bmag*sin_alph)/(4.0*np.pi*m_e*cvel)

Te = 8000.0
nu_t = 0.001
nu_break = np.sqrt(0.150*0.408)*1e9
extn = np.exp(-1.0*np.power((nu_t/frequency[4]),2.1))
alpha1, alpha2 = 2.6728667075093107, 2.7477254162083455
GSPAN = 100
nu = 1.420
nu_min = nu*1e9/GSPAN
nu_max = nu*1e9*GSPAN
gama_min = np.sqrt((nu_min)/scale_gam_nu)
gama_max = np.sqrt((nu_max)/scale_gam_nu)
gama_break = np.sqrt((nu_break)/scale_gam_nu)

xb = gama_break
xl = gama_min
xu = gama_max

def integrand_for_param(gama, alpha):
    nu_c = scale_gam_nu * (gama**2)
    x = nu/nu_c
    integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
    return integrand_


if xl < xb:
    C1 = alpha2
    I, _ = integrate.quad(integrand_for_param, gama_min, gama_max, args = (alpha1))
    I *= np.power(gama_break, 2*C1-3)

elif xu < xb:
    C1 = alpha2
    I, _ = integrate.quad(integrand_for_param, gama_min, gama_max, args = (alpha2))
    I *= np.power(gama_break, 2*C1-3)

else:
    xu = xb
    C1 = alpha2
    I1, _ = integrate.quad(integrand_for_param, gama_min, gama_max, args = (alpha1))
    I1 *= np.power(gama_break, 2*C1-3)




def F(x):
    if x<3: 
        one = (np.pi*x)/3
        two =  (9*(x**(11/3))*gamma(-2/3))/(160*2**(2/3)) 
        three = ((x**(1/3))*(16+(3*x**2))*gamma(-1/3))/(24*2**(1/3))
        return -one + two -three  
    else:
        expo = np.exp(-x)/(967458816*np.sqrt(2)*x**(5/3))
        const = 13*np.sqrt(np.pi)
        quad_term = 2429625 + 2*x*(-1922325 + (5418382*x) + 83221732*(x**2))
        error_function_term =  119630621+6*np.exp(x)*np.pi*(x**(7/2))+sp.special.erfc(np.sqrt(x))
        return expo*((const*quad_term) - error_function_term) 

def integrand(gamma, x, alpha):
    result = F(x)*x*gamma**(-(2*alpha) - 3)
    return result

def convex_func(nus, C_1, alpha1, alpha2, nu_break, I_x, T_e, nu_t):
    b_temps = []
    global scale_gam_nu, GSPAN
    print(f"alpha1 is {alpha1}")

    # nu_min = 1420*1e6/GSPAN
    # nu_max = 1420*1e6*GSPAN
    # gama_min = np.sqrt((nu_min)/scale_gam_nu)
    # gama_max = np.sqrt((nu_max)/scale_gam_nu)
    # gama_break = np.sqrt((nu_break)/scale_gam_nu)
    # #print(f"gama_break is {gama_break}")

    gam_alpha1_term  = (gama_break**((2*alpha1) - 3))
    #print(f"gam_alpha1_term is {gam_alpha1_term}")
    gam_alpha2_term  = (gama_break**((2*alpha2) - 3))
    #print(f"gam_alpha2_term is {gam_alpha2_term}")
    
    #nu = 1420.
    for nu in nus:
        gamma = np.sqrt(nu/scale_gam_nu)
        gamma_t = np.sqrt(nu_t/scale_gam_nu)
        x = gamma/gamma_t
        integ1, _ = integrate.quad(integrand, gama_min, gama_max, args = (x, alpha1))
        integ2, _ = integrate.quad(integrand, gama_min, gama_max, args = (x, alpha2))
        expo = np.exp(-1*((gamma_t/gamma)**2.1))
        three = I_x*(gamma**-2.1)
        result = C_1*((gamma**-2)*(gam_alpha1_term*integ1 + gam_alpha2_term*integ2) + three)* expo + T_e*(1 - expo)
        b_temps.append(result)
        
    return result, b_temps

#Defining chi_square function
def chisq(params, xobs, yobs):
    ynew, _ = convex_func(xobs, *params)
    #yerr = np.sum((ynew- yobs)**2)
    yerr = np.sum(((yobs- ynew)/ynew)**2)
    print(f"y error is {yerr}")
    return yerr
#bounds = ([0,100], [2,3], [2,3], [0, 1e12], [0, 1e-15], [0,5000], [0,1e7])

result = minimize(chisq,args = (frequency, b_temp), x0 = [25000, 2.5, 2.5, 0.36e6,8.39e-10, 2060, 0.3e6] , method='Nelder-Mead', options={'verbose': 1, 'maxiter': 100000})


a0= [75, 2.5, 2.5, 0.36e9,8.39e-10, 2060, 0.3e6]
xs = np.linspace(1,24000,100)
_, yinitial = convex_func(xs, *np.array(a0))
_, ys = convex_func(xs, *result.x)
plt.yscale("log")
plt.xscale("log")
plt.plot(xs, ys, label='best fit')
plt.plot(xs, yinitial,label='guess')
print(result.x)
print(result)
plt.plot(frequency, b_temp, 'r*' )
plt.xlabel("log Frequency[MHz]")
plt.ylabel("log Temp[K]")
plt.legend()
plt.title('log Temparature vs log Frequency')

plt.grid()
plt.show()