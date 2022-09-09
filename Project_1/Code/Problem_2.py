import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Problem_2.txt") #numpy's txt-reading

x,y = data[:,0], data[:,1] #putting x and y into arrays

plt.plot(x,y,label="u(x)") # plotting y over x

#Plotting stuff:
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("Problem_2_plot.svg")