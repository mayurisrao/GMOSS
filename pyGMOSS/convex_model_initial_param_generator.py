import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import gamma
import scipy.special as special
import scipy.integrate as integrate
import math
import scipy as sp
import consts


frequency = np.array([22.,45.,150.,408.,1420.,23000.]) #x_values
frequency = np.array([np.float32(f*10**-3) for f in frequency])

# file_data = np.genfromtxt('convexity.txt', delimiter = ' ')
# print(file_data)

# arrays = [np.array(map(int, line.split())) for line in open('convexity.txt')]
# print(arrays)


f = open("convexity.txt", "r")
array = []
line = f.readline()
index = 0
while line:
    line = line.strip("\n")
    line = line.split()
    array.append([])
    flag = 0
    for item in line:
        if flag == 0 or flag == 1 or flag == 2:
            array[index].append(item)
        elif flag == 3:
            array[index].append(array(item))
        else:
            array[index].append(str(item))
    line = f.readline()
    index += 1
f.close()

print(array)

# cvel = 2.99792458e+08 # m s^-1
# m_e = 9.1e-31
# q_e = 1.6e-19
# sin_alph = 1.0
# Bmag = 1e-9 # Tesla == 10 micro-Gauss
# GSPAN = 100.0

# scale_gam_nu = (3.0*q_e*Bmag*sin_alph)/(4.0*np.pi*m_e*cvel)

# Te = 8000.0
# nu_t = 0.001
# nu_break = np.sqrt(0.150*0.408)*1e9
# extn = np.exp(-1.0*np.power((nu_t/frequency[4]),2.1))
# alpha1, alpha2 = 2.6728667075093107, 2.7477254162083455

# nu = 1.420
# nu_min = nu*1e9/GSPAN
# nu_max = nu*1e9*GSPAN
# gama_min = np.sqrt((nu_min)/scale_gam_nu)
# gama_max = np.sqrt((nu_max)/scale_gam_nu)
# gama_break = np.sqrt((nu_break)/scale_gam_nu)


# class Initial_param_convex:

#     def __init__(self):
#         pass
    
#     def f(x):
#         if x<3: 
#             one = (np.pi*x)/np.sqrt(3)
#             two =  (9*(x**(11/3))*gamma(-2/3))/(160*2**(2/3)) 
#             three = ((x**(1/3))*(16+(3*x**2))*gamma(-1/3))/(24*2**(1/3))
#             return -one + two -three  
#         else:
#             exponential_term = np.exp(-x)/(967458816*np.sqrt(2)*x**(5/2))
#             const = 13*np.sqrt(np.pi)
#             quad_term = 2429625 + 2*x*(-1922325 + (5418382*x) + 83221732*(x**2))
#             error_function_term =  1196306216*np.exp(x)*np.pi*(x**(7/2))*sp.special.erfc(np.sqrt(x))
#             return exponential_term*((const*quad_term) - error_function_term) 

#     def integrand_for_param(gama, alpha):
#         nu_c = (scale_gam_nu * (gama**2))/1e9
#         x = nu/nu_c
#         integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
#         return integrand_

#     def integrand_for_convex(gama, alpha, nu):
#         nu_c = (scale_gam_nu * (gama**2))/1e9
#         x = nu/nu_c
#         integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
#         return integrand_
    
#     def I_generator(self):
#         if xl > xb:
#             C1 = alpha2
#             I, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha1))
#             I *= np.power(gama_break, 2*C1-3)

#         elif xu < xb:
#             C1 = alpha1
#             I, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha2))
#             I *= np.power(gama_break, 2*C1-3)

#         else:
#             xu = xb
#             C1 = alpha1
#             I1, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha1))
#             I1 *= np.power(gama_break, 2*C1-3)
#             xl = xb
#             xu = gama_max
#             C1 = alpha2
#             I2, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha2))
#             I2 *= np.power(gama_break, 2*C1-3)
#             I = I1 + I2


#________________________________________________________________################################________________________________________________________

# if __name__ == "__main__":



# frequency = np.array([22.,45.,150.,408.,1420.,23000.]) #x_values
# frequency = np.array([np.float32(f*10**-3) for f in frequency])
# b_temp = np.array([9.18573758e+04, 1.77507604e+04, 7.10610657e+02, 6.49989393e+01, 2.11183872e+00, 9.89014738e-04])#y_values

# cvel = 2.99792458e+08 # m s^-1
# m_e = 9.1e-31
# q_e = 1.6e-19
# sin_alph = 1.0
# Bmag = 1e-9 # Tesla == 10 micro-Gauss
# GSPAN = 100.0

# scale_gam_nu = (3.0*q_e*Bmag*sin_alph)/(4.0*np.pi*m_e*cvel)

# Te = 8000.0
# nu_t = 0.001
# nu_break = np.sqrt(0.150*0.408)*1e9
# extn = np.exp(-1.0*np.power((nu_t/frequency[4]),2.1))
# alpha1, alpha2 = 2.6728667075093107, 2.7477254162083455

# nu = 1.420
# nu_min = nu*1e9/GSPAN
# nu_max = nu*1e9*GSPAN
# gama_min = np.sqrt((nu_min)/scale_gam_nu)
# gama_max = np.sqrt((nu_max)/scale_gam_nu)
# gama_break = np.sqrt((nu_break)/scale_gam_nu)

# xb = gama_break
# xl = gama_min
# xu = gama_max

# def F(x):
#     if x<3: 
#         one = (np.pi*x)/np.sqrt(3)
#         two =  (9*(x**(11/3))*gamma(-2/3))/(160*2**(2/3)) 
#         three = ((x**(1/3))*(16+(3*x**2))*gamma(-1/3))/(24*2**(1/3))
#         return -one + two -three  
#     else:
#         exponential_term = np.exp(-x)/(967458816*np.sqrt(2)*x**(5/2))
#         const = 13*np.sqrt(np.pi)
#         quad_term = 2429625 + 2*x*(-1922325 + (5418382*x) + 83221732*(x**2))
#         error_function_term =  1196306216*np.exp(x)*np.pi*(x**(7/2))*sp.special.erfc(np.sqrt(x))
#         return exponential_term*((const*quad_term) - error_function_term) 

# def integrand_for_param(gama, alpha):
#     nu_c = (scale_gam_nu * (gama**2))/1e9
#     x = nu/nu_c
#     integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
#     return integrand_

# def integrand_for_convex(gama, alpha, nu):
#     nu_c = (scale_gam_nu * (gama**2))/1e9
#     x = nu/nu_c
#     integrand_ = F(x)*x*np.power(gama, -1*(2*alpha - 3)) 
#     return integrand_


# if xl > xb:
#     C1 = alpha2
#     I, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha1))
#     I *= np.power(gama_break, 2*C1-3)

# elif xu < xb:
#     C1 = alpha1
#     I, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha2))
#     I *= np.power(gama_break, 2*C1-3)

# else:
#     xu = xb
#     C1 = alpha1
#     I1, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha1))
#     I1 *= np.power(gama_break, 2*C1-3)
#     xl = xb
#     xu = gama_max
#     C1 = alpha2
#     I2, _ = integrate.quad(integrand_for_param, xl, xu, args = (alpha2))
#     I2 *= np.power(gama_break, 2*C1-3)
#     I = I1 + I2

# fnorm = (b_temp[4] - (Te*(1.0- extn)))/((np.power(nu,-2.0)*I)*extn)
# print(fnorm)
# print(I)
# print(extn)

# nu = 22.690
# nu_min = nu*1e9/GSPAN
# print(nu_min)
# nu_max = nu*1e9*GSPAN
# print(nu_max)
# gama_min = np.sqrt((nu_min)/scale_gam_nu)
# gama_max = np.sqrt((nu_max)/scale_gam_nu)
# gama_break = np.sqrt((nu_break)/scale_gam_nu)

# xb = gama_break
# xl = gama_min
# xu = gama_max

# if xl > xb:
#     C1 = alpha2
#     I, _ = integrate.quad(integrand_for_param, xl,xu, args = (alpha1))
#     I *= np.power(gama_break, 2*C1-3)

# elif xu < xb:
#     C1 = alpha1
#     I, _ = integrate.quad(integrand_for_param, xl,xu, args = (alpha2))
#     I *= np.power(gama_break, 2*C1-3)

# else:
#     xu = xb
#     C1 = alpha1
#     I1, _ = integrate.quad(integrand_for_param, xl,xu, args = (alpha1))
#     I1 *= np.power(gama_break, 2*C1-3)
#     xl = xb
#     xu = gama_max
#     C1 = alpha2
#     I2, _ = integrate.quad(integrand_for_param, xl,xu, args = (alpha2))
#     I2 *= np.power(gama_break, 2*C1-3)
#     I = I1 + I2

# extn = np.exp(-1.0*np.power((nu_t/nu),2.1))
# temp1 = fnorm*extn

# Tx = (((b_temp[5] - Te*(1.0-extn))/temp1)-(np.power(nu,-2.0)*I))/np.power(frequency[4],-2.1)
# if Tx <= 0:
#     Tx = 1.0e-10

# print(f"the parameters are fnorm = {fnorm}, alpha1 = {alpha1}, alpha2 = {alpha2}, nu_break = {nu_break/1e6}, Tx = {Tx}, Te = {Te}, nu_t = {nu_t}")

# x_ini = []
# #x_ini.extend([np.log10(fnorm), np.log(alpha1), np.log10(alpha2), np.log10(nu_break), np.log10(Tx), np.log10(Te), np.log10(nu_t)])
# x_ini.extend([fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t])
# print(x_ini)
