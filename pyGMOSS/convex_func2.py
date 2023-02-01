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
GSPAN = 100
#adopting Dr. Rao's method and writing a function that has to be minimized
#therefore essentially a code that computes chi_square
NHPIX = 3072
def func(pp): #pp is an array
    #unpacking pp
    fnorm, alpha1, alpha2, nu_break, Tx, Te, nu_t = pp
    #blowing up the  
    if alpha1 < 2.0 and alpha2 > 2.0: alpha1 1000000000
    if alpha2 < 2.0 and alpha2 > 3.0: alpha2 1000000000
    for i in range(len(frequency)):
        nu = frequency[i]
        nu_min = nu*1e9/GSPAN
        nu_max = 




