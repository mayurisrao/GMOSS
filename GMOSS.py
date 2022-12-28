import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time as time
import glob

start = time.time()
pixels = 3072

#some important constants
kk = 1.3806488e-23 # m^2 kg s^-1
hh = 6.62606957e-34 # m^2 kg s^-2 K^-1
cvel = 2.99792458e+08 # m s^-1
m_e = 9.1e-31
q_e = 1.6e-19
sin_alph = 1.0
Bmag = 1e-9 # Tesla == 10 micro-Gauss
#PI = 3.14159265358979
FTOL = 1.0e-6
T_e = 8000.0 # Thermal electron temperature
TCMB = 2.72548
NHPIX = 3072
GSPAN = 100.0

correction_150offset = 21.4
correction_150scaling = 1.05
#apply corrections to 150,408,1420
#corrections_applied_to = np.array([150,408,1420])

#some important stuff
scale_gam_nu = (3.0*q_e*Bmag*sin_alph)/(4.0*np.pi*m_e*cvel)

#appending the maps to the list
text_files = glob.glob("*.txt")
map_files = []
for file in text_files:
    if file.startswith("map"):
        map_files.append(file)

frequency = np.array([])

#obtaining the sorted list of maps and frequencies
for file in map_files:
    frequency = np.append(frequency, float(file.split("_")[1]))

frequency = np.sort(frequency)
frequency = [int(i) for i in frequency]

map_files = [f'map_{i}_r4_5deg_nested_galactic_Kelvin.txt' for i in frequency]

#obtaining frequency in GHz and the log of the frequencyies
frequency_in_GHz = [np.float32(f*10**-3) for f in frequency] 
#print(frequency)

f = open("convexity.txt", "w")
length_ = len(map_files)

######################Starting our main program#####################
for pixel in range(pixels):
    print(f"Pixel Number: {pixel+1}")
    b_temp = np.array([])
    
    for i in range(length_):
        data = np.genfromtxt(map_files[i])
        if frequency[i] ==  22: b_temp= np.append(b_temp, data[pixel])
        if frequency[i] == 45: b_temp= np.append(b_temp, data[pixel])
        if frequency[i] == 150:
            data[pixel] = (data[pixel] - correction_150offset)*correction_150scaling
            data[pixel] = data[pixel] - TCMB
            b_temp= np.append(b_temp, data[pixel])
        if frequency[i] == 408: b_temp= np.append(b_temp, data[pixel] - TCMB)
        if frequency[i] == 1420: b_temp= np.append(b_temp, data[pixel] - TCMB)
        if frequency[i] == 23000: b_temp= np.append(b_temp, data[pixel])

    #to obtain concavity or convexity
    #Spectral index between 45 MHz and 150 MHz
    alpha1 = (np.log10(b_temp[1]) - np.log10(b_temp[2]))/(np.log10(frequency_in_GHz[2]) - np.log10(frequency_in_GHz[1]))
    if alpha1 < 2.0: alpha1 = 2.
    if alpha1 > 3.0: alpha1 = 3.
    
    #to obtain concave or convexity
    #Spectral index between 408 MHz and 1420 MHz
    alpha2 = (np.log10(b_temp[3]) - np.log10(b_temp[4]))/(np.log10(frequency_in_GHz[4]) - np.log10(frequency_in_GHz[3]))
    if alpha2 < 2.0: alpha2 = 2.
    if alpha2 > 3.0: alpha2 = 3.

    if alpha1 < alpha2: f.write(f"{pixel+1} convex\n")
    elif alpha1 > alpha2: f.write(f"{pixel+1} concave\n")

f.close()
