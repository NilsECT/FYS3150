import matplotlib.pyplot as plt
import numpy as np

"""
All the independent variables are stored as individual matrices in the shape of columns individual particles ans rows as time step in the simulation. 
Shape of array [time step, particle number]
"""

## Data extraction with interaction:
x = np.loadtxt("x.txt")
y = np.loadtxt("y.txt")
y = np.loadtxt("z.txt")

vx = np.loadtxt("vx.txt")
vy = np.loadtxt("vy.txt")
vz = np.loadtxt("vz.txt")



## Two particles with interactions:
plt.title("Two particles in xy-plane with interaction")

plt.title("Phase space for two particles with interaction x and vx")

plt.title("Phase space for two particles with interaction y and vy")

plt.title("Phase space for two particles with interaction z and vz")

plt.title("3D-plot of 2 particles with interaction")



## Data extraction without interaction:
x = np.loadtxt("x.txt")
y = np.loadtxt("y.txt")
y = np.loadtxt("z.txt")

vx = np.loadtxt("vx.txt")
vy = np.loadtxt("vy.txt")
vz = np.loadtxt("vz.txt")

## Two particles without interactions:
plt.title("Two particles in xy-plane without interaction")

plt.title("Phase space for two particles without interaction x and vx")

plt.title("Phase space for two particles without interaction y and vy")

plt.title("Phase space for two particles without interaction z and vz")

plt.title("3D-plot of 2 particles without interaction")

## Data extraction single particle:
x = np.loadtxt("x.txt")
y = np.loadtxt("y.txt")
y = np.loadtxt("z.txt")

vx = np.loadtxt("vx.txt")
vy = np.loadtxt("vy.txt")
vz = np.loadtxt("vz.txt")

## Single particle:
plt.title("Z-movement over time")
plt.plot(z)
plt.xlabel("time in s")
plt.ylabel("Movement on z-axis")

plt.title("Single for different h-values for Runge-Kutta-4")

plt.title("Single for different h-values with Forward Euler")

## 100 particle 
plt.title("")

