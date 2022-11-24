from scipy.optimize import minimize
import numpy as np
import pandas as pd

def chi_square(temperature_data):
    global temperature_model
    chi_sq = np.sum(np.power(((temperature_data - temperature_model) / temperature_data), 2))
    chi_sq /= len(chi_sq)
    return chi_sq

df = pd.read_csv("pixel_fits_degree_1.csv", converters = {"Log_Frequency": eval, "Log_Frequency": eval, "Fit": eval})
pixel_number = df.loc[:, "Pixel_Number"][0]
log_frequency = df.loc[:, "Log_Frequency"][0]
log_brightness_temperature = df.loc[:, "Log_Brightness_Temperature"][0]
fit = df.loc[:, "Fit"][0]

print(type(log_brightness_temperature))
print(type(log_frequency))

#temperature_model = (fit[0] * log_frequency) + fit[1]

#res = minimize(method = "Nelder-Mead", fun = chi_square, initial_simplex = fit.T)

#print(res.x)