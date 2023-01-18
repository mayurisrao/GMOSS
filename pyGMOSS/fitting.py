import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#sns.set()

data1 = np.genfromtxt("map_22_r4_5deg_nested_galactic_Kelvin.txt")
data2 = np.genfromtxt("map_45_r4_5deg_nested_galactic_Kelvin.txt")
data3 = np.genfromtxt("map_150_r4_5deg_nested_galactic_Kelvin.txt")
data4 = np.genfromtxt("map_408_r4_5deg_nested_galactic_Kelvin.txt")
data5 = np.genfromtxt("map_1420_r4_5deg_nested_galactic_Kelvin.txt")
data6 = np.genfromtxt("map_23000_r4_5deg_nested_galactic_Kelvin.txt")

pixel_number = int(input("Enter the pixel number from 1 to 3072: ")) - 1
print(pixel_number)

for_graphing = np.zeros(6)
len_of_data = len(data1)

for_graphing[0] = data1[pixel_number]
for_graphing[1] = data2[pixel_number]
for_graphing[2] = data3[pixel_number]
for_graphing[3] = data4[pixel_number]
for_graphing[4] = data5[pixel_number]
for_graphing[5] = data6[pixel_number]
plt.yscale("log")
plt.xscale("log")
plt.plot(for_graphing)
#plt.yscale("log")
plt.savefig("roughplot.png")
print(for_graphing)

