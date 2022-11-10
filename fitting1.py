import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
sns.set()

x = [10,10**2,10**3,10**4]
plt.yscale("log")
plt.xscale("log")
plt.plot(x)


data1 = np.genfromtxt("map_22_r4_5deg_nested_galactic_Kelvin.txt")
data2 = np.genfromtxt("map_45_r4_5deg_nested_galactic_Kelvin.txt")
data3 = np.genfromtxt("map_150_r4_5deg_nested_galactic_Kelvin.txt")
data4 = np.genfromtxt("map_408_r4_5deg_nested_galactic_Kelvin.txt")
data5 = np.genfromtxt("map_1420_r4_5deg_nested_galactic_Kelvin.txt")
data6 = np.genfromtxt("map_23000_r4_5deg_nested_galactic_Kelvin.txt")
coods = np.genfromtxt("coods.txt", delimiter=" ", skip_header = 1)

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
plt.scatter([22,45,150,408,1420,23000],for_graphing)
#plt.yscale("log")
#plt.savefig("roughplot.png")
print(for_graphing)

coordinates = np.genfromtxt("coods.txt", delimiter=" ", skip_header = 1)
coordinates[0]

coordinates[:,0]

type(coordinates[:,1][0])

gal = SkyCoord(coordinates[:,0], coordinates[:,1], frame='galactic', unit=u.deg)

plt.subplot(111, projection='aitoff')
plt.grid(True)
plt.scatter(gal.l.wrap_at('180d').radian, gal.b.radian)

# Use the heatmap function from the seaborn package
sns.heatmap(coordinates,fmt="",cmap='RdYlGn',linewidths=0.30)

# Display the Pharma Sector Heatmap
plt.show()