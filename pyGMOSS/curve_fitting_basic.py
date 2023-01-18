:import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import scipy as sp
import seaborn as sns
from scipy.optimize import curve_fit
sns.set()


x = np.array([22,45,150,408,1420,23000])
y = np.array([5.49282937e+04, 1.04875149e+04, 5.57289208e+02, 4.39957717e+01, 4.37429422e+00, 1.07386277e-03])
#y = np.array([11, 10, 12, 14, 15, 15, 14, 13, 14, 11, 10, 11, 7, 9, 8, 6, 5])  
#x = np.array(range(len(y)))  

#plt.grid()
plt.scatter(x,y)
plt.plot(x,y, c = "red",linestyle='dashed')
#plt.show()

#def model_f(x,a,b,c):
#    return a*(x-b)**2 + c

# defining objective functions  
def mapping1(x, a, b, c):  
    return a * x**2 + b * x + c  
  
def mapping2(x, a, b, c):  
    return a * x**3 + b * x + c  
  
def mapping3(x, a, b, c):  
    return a * x**3 + b * x**2 + c  
  
def mapping4(x, a, b, c):  
    return a * np.exp(b * x) + c  

#using cuve_fit fuciton from scipy
args, covar = curve_fit(mapping1, x, y, maxfev = 5000)  
print("Arguments: ", args)  
print("Co-Variance: ", covar)

args, _ = curve_fit(mapping1, x, y,maxfev = 5000)  
a, b, c = args[0], args[1], args[2]  
y_fit1 = a * x**2 + b * x + c  
  
args, _  = curve_fit(mapping2, x, y, maxfev = 5000)  
a, b, c = args[0], args[1], args[2]  
y_fit2 = a * x**3 + b * x + c  
  
args, _  = curve_fit(mapping3, x, y, maxfev = 5000)  
a, b, c = args[0], args[1], args[2]  
y_fit3 = a * x**3 + b * x**2 + c  
  
#args, _  = curve_fit(mapping4, x, y, maxfev = 5000)  
#a, b, c = args[0], args[1], args[2]  
#y_fit4 = a * np.exp(x * b) + c  


# plotting the graph  
plt.plot(x, y, 'bo', label="y - original")  
plt.plot(x, y_fit1, label="y = a * x^2 + b * x + c")  
plt.plot(x, y_fit2, label="y = a * x^3 + b * x + c")  
plt.plot(x, y_fit3, label="y = a * x^3 + b * x^2 * c")  
#plt.plot(x, y_fit4, label="y = a * exp(b * x) + c")  
plt.xlabel('x')  
plt.ylabel('y')  
#plt.xlim([min(x)-100,max(x)+100])
plt.legend(loc = 'best', fancybox = True, shadow = True)  
plt.grid(True)  
plt.show()

'''The problem is cause by the exponential function- Try to tackle it'''