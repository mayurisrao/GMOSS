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
frequency = np.array([22.,45.,150.,408.,1420.,23000.]) #x_values
frequency = np.array([np.float32(f*10**-3) for f in frequency])
#b_temp = [9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04] #y_values
print(frequency)
b_temp = np.array([9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04])

cvel = 2.99792458e+08 # m s^-1
m_e = 9.1e-31
q_e = 1.6e-19
sin_alph = 1.0
Bmag = 1e-9 # Tesla == 10 micro-Gauss
#PI = 3.14159265358979
# FTOL = 1.0e-6
#T_e = 8000.0 # Thermal electron temperature
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

def integrand_for_param(gama, alpha):
    nu_c = (scale_gam_nu * (gama**2))/1e9
    x = nu/nu_c
    integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
    return integrand_

def integrand_for_convex(gama, alpha, nu):
    nu_c = (scale_gam_nu * (gama**2))/1e9
    x = nu/nu_c
    integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
    return integrand_


if xl > xb:
    C1 = alpha2
    I, _ = integrate.quad(integrand_for_param, xl, xb, args = (alpha1))
    I *= np.power(gama_break, 2*C1-3)

elif xu < xb:
    C1 = alpha1
    I, _ = integrate.quad(integrand_for_param, xl, xb, args = (alpha2))
    I *= np.power(gama_break, 2*C1-3)

else:
    xu = xb
    C1 = alpha1
    I1, _ = integrate.quad(integrand_for_param, xl, xb, args = (alpha1))
    I1 *= np.power(gama_break, 2*C1-3)
    xl = xb
    xu = gama_max
    C1 = alpha2
    I2, _ = integrate.quad(integrand_for_param, xl, xb, args = (alpha2))
    I2 *= np.power(gama_break, 2*C1-3)
    I = I1 + I2

fnorm = (b_temp[4] - (Te*(1.0- extn)))/((np.power(nu,-2.0)*I)*extn)
print(fnorm)
print(I)
print(extn)

nu = 22.690
nu_min = nu*1e9/GSPAN
print(nu_min)
nu_max = nu*1e9*GSPAN
print(nu_max)
gama_min = np.sqrt((nu_min)/scale_gam_nu)
gama_max = np.sqrt((nu_max)/scale_gam_nu)
gama_break = np.sqrt((nu_break)/scale_gam_nu)

xb = gama_break
xl = gama_min
xu = gama_max

if xl > xb:
    C1 = alpha2
    I, _ = integrate.quad(integrand_for_param, xl,xb, args = (alpha1))
    I *= np.power(gama_break, 2*C1-3)

elif xu < xb:
    C1 = alpha1
    I, _ = integrate.quad(integrand_for_param, xl,xb, args = (alpha2))
    I *= np.power(gama_break, 2*C1-3)

else:
    xu = xb
    C1 = alpha1
    I1, _ = integrate.quad(integrand_for_param, xl,xb, args = (alpha1))
    I1 *= np.power(gama_break, 2*C1-3)
    xl = xb
    xu = gama_max
    C1 = alpha2
    I2, _ = integrate.quad(integrand_for_param, xl,xb, args = (alpha2))
    I2 *= np.power(gama_break, 2*C1-3)
    I = I1 + I2

extn = np.exp(-1.0*np.power((nu_t/nu),2.1))
temp1 = fnorm*extn

Tx = (((b_temp[5] - Te*(1.0-extn))/temp1)-(np.power(nu,-2.0)*I))/np.power(frequency[4],-2.1)
if Tx <= 0:
    Tx = 1.0e-10

print(f"the parameters are fnorm = {fnorm}, alpha1 = {alpha1}, alpha2 = {alpha2}, nu_break = {nu_break/1e6}, Tx = {Tx}, Te = {Te}, nu_t = {nu_t}")

x_ini = []
#x_ini.extend([np.log10(fnorm), np.log(alpha1), np.log10(alpha2), np.log10(nu_break), np.log10(Tx), np.log10(Te), np.log10(nu_t)])
x_ini.extend([fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t])
print(x_ini)
print(x_ini)

print(nu_min)
print(nu_max)
def convex_func(nus, C_1, alpha1, alpha2, nu_break, I_x, Te, nu_t):
    b_temps = []
    global scale_gam_nu, GSPAN
    #print(f"alpha1 is {alpha1}")

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
    
    
    for fre in nus:
        # gamma = np.sqrt(nu/scale_gam_nu)
        # gamma_t = np.sqrt(nu_t/scale_gam_nu)
        # x = nu/nu_c
        integ1, _ = integrate.quad(integrand_for_convex, gama_min, gama_break, args = (alpha1, fre))
        print(f"integ1 = {integ1}")
        integ2, _ = integrate.quad(integrand_for_convex, gama_break, gama_max, args = (alpha2, fre))
        print(f"integ2 = {integ2}")
        print(f"nu_t = {nu_t}")
        print(f"fre = {fre}")
        expo = np.exp(-1*((nu_t/fre)**2.1))
        print(f"expo = {expo}")
        three = I_x*np.power(fre, -2.1)
        print(f"three = {three}")
        result = C_1*((fre**-2)*(gam_alpha1_term*integ1 + gam_alpha2_term*integ2) + three)* expo + Te*(1 - expo)
        print(f"result = {result}")
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
bounds = ([-np.inf, np.inf], [-np.inf, np.inf], [2,3], [2,3], [-np.inf, np.inf], [0,10000],[-np.inf, np.inf])
result = minimize(chisq,args = (frequency, b_temp), x0 = x_ini , method='Nelder-Mead',bounds = bounds, options={'verbose': 1, 'maxiter': 100000})


#a0= [75, 2.5, 2.5, 0.36e9,8.39e-10, 2060, 0.3e6]
xs = np.linspace(22e-3,24,100)
#_, yinitial = convex_func(xs, *np.array(a0))
_, ys = convex_func(xs, *result.x)
# plt.yscale("log")
# plt.xscale("log")
plt.plot(xs, ys, label='best fit')
#plt.plot(xs, yinitial,label='guess')
print(result.x)
print(result)
plt.plot(frequency, b_temp, 'r*' )
plt.xlabel("log Frequency[MHz]")
plt.ylabel("log Temp[K]")
plt.legend()
plt.title('log Temparature vs log Frequency')
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.show()