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
