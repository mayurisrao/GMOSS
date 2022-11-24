import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import seaborn as sns
sns.set()
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import glob
import time

pixels = 3072 # Total Number of Pixels
DEGREE = 1

start = time.time()

text_files = glob.glob("*.txt")
map_files = []
for file in text_files:
    if file.startswith("map"):
        map_files.append(file)

frequency = np.array([])

for file in map_files:
    frequency = np.append(frequency, float(file.split("_")[1]))

log_frequency = np.log10(frequency)

pixel_fits = []

print(f"Fitting a {DEGREE} Degree Polynomial")

for pixel in range(pixels):
    print(f"Pixel Number: {pixel}")
    brightness_temperature = np.array([])
    pixel_fit = {}
    for file in map_files:
        data = np.genfromtxt(file)
        brightness_temperature = np.append(brightness_temperature, data[pixel])

    log_brightness_temperature = np.log10(brightness_temperature)
    
    mean_second_derivative = np.mean(np.diff(np.diff(log_brightness_temperature)))

    z = np.polyfit(log_frequency, log_brightness_temperature, deg = DEGREE)
    f = np.poly1d(z)

    pixel_fit["Pixel_Number"] = pixel
    pixel_fit["Log_Frequency"] = log_frequency
    pixel_fit["Log_Brightness_Temperature"] = log_brightness_temperature
    pixel_fit["Fit"] = z
    if mean_second_derivative < 0:
        pixel_fit["Concave/Convex"] = "Concave"
    elif mean_second_derivative > 0:
        pixel_fit["Concave/Convex"] = "Convex"
    else:
        pixel_fit["Concave/Convex"] = "Straight"

    pixel_fits.append(pixel_fit)

df = pd.DataFrame(pixel_fits)
df.to_csv(f"pixel_fits_degree_{DEGREE}.csv", quoting = csv.QUOTE_ALL)

end = time.time()

print(f"Time Taken to fit {DEGREE} polynomials for {pixels} pixels = {end - start} s")
