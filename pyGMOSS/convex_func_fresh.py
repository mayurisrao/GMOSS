import numpy as np
from scipy.optimize import minimize
from scipy.special import gamma
import scipy as sp
from scipy import integrate

frequency = np.array([22.,45.,150.,408.,1420.,23000.]) #x_values
frequency = np.array([np.float32(f*10**-3) for f in frequency])

b_temp = np.array([9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04])

def F(x):
    vals = np.array([])
    for i in x:
        if i<3: 
            one = (np.pi*i)/3
            two =  (9*(i**(11/3))*gamma(-2/3))/(160*2**(2/3)) 
            three = ((i**(1/3))*(16+(3*i**2))*gamma(-1/3))/(24*2**(1/3))
            vals = np.append(vals, -one + two -three)
            #return -one + two -three
        else:
            expo = np.exp(-i)/(967458816*np.sqrt(2)*i**(5/2))
            const = 13*np.sqrt(np.pi)
            quad_term = 2429625 + 2*i*(-1922325 + (5418382*i) + 83221732*(i**2))
            error_function_term =  1196306216*np.exp(i)*np.pi*(i**(7/2))*sp.special.erfc(np.sqrt(i))
            #return expo*((const*quad_term) - error_function_term)
            vals = np.append(vals,expo*((const*quad_term) - error_function_term))
    return vals

def convex_func(nus, C_1, alpha1, alpha2, nu_break, I_x, Te, nu_t):
    b_temps = []
    global scale_gam_nu, GSPAN
    gam_alpha1_term  = (gama_break**((2*alpha1) - 3))
    gam_alpha2_term  = (gama_break**((2*alpha2) - 3))

    integ1, _ = integrate.quad(integrand_for_convex, gama_min, gama_break, args = (alpha1, nus))
    integ2, _ = integrate.quad(integrand_for_convex, gama_break, gama_max, args = (alpha2, nus))

    expo = np.exp(-1*((nu_t/nus)**2.1))

    three = I_x*np.power(nus, -2.1)
    result = C_1*((nus**-2)*(gam_alpha1_term*integ1 + gam_alpha2_term*integ2) + three)* expo + Te*(1 - expo)
    print(f"result = {result}")
    b_temps.append(result)
        
    return result, b_temps


