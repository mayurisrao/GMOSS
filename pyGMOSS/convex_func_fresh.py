import numpy as np
from scipy.optimize import minimize
from scipy.special import gamma


frequency = np.array([22.,45.,150.,408.,1420.,23000.]) #x_values
frequency = np.array([np.float32(f*10**-3) for f in frequency])

b_temp = np.array([9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04])

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
s = F(np.array([5,7,9]))
print(s)
# def convex_model(frequency, params):
