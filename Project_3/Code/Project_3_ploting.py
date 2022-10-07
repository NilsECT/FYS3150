import matplotlib.pyplot as plt
import numpy as np

"""
All the independent variables are stored as individual matrices in the shape of columns individual particles ans rows as time step in the simulation. 
Shape of array [time step, particle number]
"""

## Data extraction:
x = np.loadtxt("x.txt")
y = np.loadtxt("y.txt")
y = np.loadtxt("z.txt")

vx = np.loadtxt("vx.txt")
vy = np.loadtxt("vy.txt")
vz = np.loadtxt("vz.txt")


## Plotting:
plt.title("Z-movement over time")

plt.title("Two particles in xy-plane with interaction")

plt.title("Two particles in xy-plane without interaction")

plt.title("Phase space for two particles with interaction x and vx")

plt.title("Phase space for two particles with interaction y and vy")

plt.title("Phase space for two particles with interaction z and vz")

plt.title("Phase space for two particles without interaction x and vx")

plt.title("Phase space for two particles without interaction y and vy")

plt.title("Phase space for two particles without interaction z and vz")

plt.title("3D-plot of 2 particles with interaction")

plt.title("3D-plot of 2 particles without interaction")

plt.title("Single for different h-values")

plt.title("Single for different h-values with Forward Euler")

## 100 particle 
plt.title("")

