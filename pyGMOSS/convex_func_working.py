import pdb
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gamma
import scipy.special as special
import scipy.integrate as integrate
import math
import scipy as sp

frequency = np.array([22.,45.,150.,408.,1420.,23000.]) #x_values
frequency = np.array([np.float32(f*10**-3) for f in frequency])
b_temp = np.array([9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04])#y_values
GSPAN = 100
#adopting Dr. Rao's method and writing a function called func that has to be minimized
#therefore essentially a code that computes chi_square
NHPIX = 3072

cvel = 2.99792458e+08 # m s^-1
m_e = 9.1e-31
q_e = 1.6e-19
sin_alph = 1.0
Bmag = 1e-9 # Tesla == 10 micro-Gauss
scale_gam_nu = (3.0*q_e*Bmag*sin_alph)/(4.0*np.pi*m_e*cvel)

#x_ini = np.array([ 9.13659851e-07,  9.66942286e+00,  2.00531364e+00,  8.86653501e+01,
#        -3.66142917e-10,  9.99999996e+03,  1.19717340e-02])
#x_ini = np.array([ -1.976754550046829e-07, 2.6728667075093107, 2.7477254162083455, 247386337.5370596/1e9, 1e-10, 8000.0, 0.001])
#x_ini = np.array([-5.761185054330454e-08, 2.6728667075093107, 2.7477254162083455, 247386337.5370596/1e9, 1e-10, 8000.0, 0.001]) #working
x_ini = np.array([-5.766426064650115e-08, 2.6728667075093107, 2.7477254162083455, 247386337.5370596/1e9, 1e-10, 8000.0, 0.001]) # broken leg
#x_ini = np.array([ -6.928880, 2.6728667075093107, 2.7477254162083455, 247.3863375370596, 1e-10, 8000.0, 0.001])

def F(x):
    if x<3: 
        one = (np.pi*x)/np.sqrt(3)
        two =  (9*(x**(11/3))*gamma(-2/3))/(160*2**(2/3)) 
        three = ((x**(1/3))*(16+(3*x**2))*gamma(-1/3))/(24*2**(1/3))
        return -one + two -three  
    else:
        exponential_term = np.exp(-x)/(967458816*np.sqrt(2)*x**(5/2))
        const = 13*np.sqrt(np.pi)
        quad_term = 2429625 + 2*x*(-1922325 + (5418382*x) + 83221732*(x**2))
        error_function_term =  1196306216*np.exp(x)*np.pi*(x**(7/2))*sp.special.erfc(np.sqrt(x))
        return exponential_term*((const*quad_term) - error_function_term)

def integrand_for_convex(gama, alpha, nu):
    nu_c = (scale_gam_nu * (gama**2))/1e9
    x = nu/nu_c
    integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
    return integrand_

def func(pp: np.ndarray) -> float: #pp is an array
    #unpacking pp
    fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = pp
    #blowing up the  
    if alpha1 < 2.0 or alpha1 > 3.0: alpha1 = 1000000000
    if alpha2 < 2.0 or alpha2 > 3.0: alpha2  = 1000000000
    if Te < 0.0 or Te > 10000: Te = 1000000000
    
    #computing chi square
    chisq = 0.0 
    for i in range(len(frequency)):
        nu = frequency[i]
        nu_min = nu*1e9/GSPAN
        nu_max = nu*1e9*GSPAN
        gama_min = np.sqrt(nu_min/scale_gam_nu) 
        gama_max = np.sqrt(nu_max/scale_gam_nu)
        gama_break = np.sqrt(nu_break/scale_gam_nu)
        xl = gama_min
        xu = gama_max
        xb = gama_break

        if xl > xb:
            C1 = alpha2
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args = (alpha1,nu))
            I *= np.power(gama_break, 2*C1-3)
        
        elif xu < xb:
            C1 = alpha1
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args = (alpha1,nu))
            I *= np.power(gama_break, 2*C1-3)
        else:
            xu = xb
            C1 = alpha1
            I1, _ = integrate.quad(integrand_for_convex, xl,xu, args = (alpha1, nu))
            I1 *= np.power(gama_break, 2*C1-3)
            xl = xb
            xu = gama_max
            C1 = alpha2
            I2, _ = integrate.quad(integrand_for_convex, xl,xu, args = (alpha2, nu))
            I2 *= np.power(gama_break, 2*C1-3)
            I = I1 + I2
        extn = np.exp(-1.0*np.power((nu_t/nu),2.1))
        FFIT = fnorm * ((np.power(nu,-2)*I) + Tx *np.power(nu,-2.1))*extn + Te *(1.0 - extn)
        DDIF = (b_temp[i] - FFIT) / (b_temp[i])
        if i<=5: chisq += (DDIF*DDIF)
    chisq /= 6.0
    print(chisq)
    return chisq
bounds = ([-np.inf, np.inf], [2,3], [2,3], [0.001, np.inf],[0.001, np.inf], [0.001,10000],[0.001, np.inf])
result = minimize(func, x0 = x_ini, method='Nelder-Mead',bounds = bounds, options={'verbose': 1, 'maxiter': 10000})
print(x_ini)
print(type(result.x))
print(result)
xs = np.linspace(22e-3,24,100)

def convex_func(frequency, pp):
    #unpacking pp
    fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = pp
    #blowing up the  
    # if alpha1 < 2.0 or alpha1 > 2.0: alpha1 = 1000000000000
    # if alpha2 < 2.0 or alpha2 > 3.0: alpha2  = 1000000000000
    # if Te < 0.0 or Te > 10000: Te = 1000000000000
    b_temps = []
    
    for i in range(len(frequency)):
        nu = frequency[i]
        nu_min = nu*1e9/GSPAN
        nu_max = nu*1e9*GSPAN
        gama_min = np.sqrt(nu_min/scale_gam_nu) 
        gama_max = np.sqrt(nu_max/scale_gam_nu)
        gama_break = np.sqrt(nu_break/scale_gam_nu)
        xl = gama_min
        xu = gama_max
        xb = gama_break

        if xl > xb:
            C1 = alpha2
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args = (alpha1,nu))
            I *= np.power(gama_break, 2*C1-3)
        
        elif xu < xb:
            C1 = alpha1
            I, _ = integrate.quad(integrand_for_convex, xl, xu, args = (alpha1,nu))
            I *= np.power(gama_break, 2*C1-3)
        else:
            xu = xb
            C1 = alpha1
            I1, _ = integrate.quad(integrand_for_convex, xl,xu, args = (alpha1, nu))
            I1 *= np.power(gama_break, 2*C1-3)
            xl = xb
            xu = gama_max
            C1 = alpha2
            I2, _ = integrate.quad(integrand_for_convex, xl,xu, args = (alpha2, nu))
            I2 *= np.power(gama_break, 2*C1-3)
            I = I1 + I2
        extn = np.exp(-1.0*np.power((nu_t/nu),2.1))
        FFIT = fnorm * ((np.power(nu,-2) * I) + Tx * np.power(nu,-2.1)) * extn + Te * (1.0 - extn)
        b_temps.append(FFIT)
    return b_temps


xs = np.linspace(22e-3,24,100)
plt.plot(xs, convex_func(xs, np.array(x_ini)), label = 'this is the initial guess')
plt.plot(frequency, b_temp, 'r*')
plt.plot(frequency, np.array(convex_func(frequency, result.x)), label = 'the curve after fitting')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("log Frequency[MHz]")
plt.ylabel("log Temp[K]")
plt.legend()
plt.grid()
plt.show()
# plt.plot(result.x)



