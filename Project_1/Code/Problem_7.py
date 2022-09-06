import numpy as np
import matplotlib.pyplot as plt

x,g,v = np.loadtxt("Problem_7.txt", unpack=True) #numpy's txt-reading

# x,g,v = data[:,0], data[:,1], data[:,2] #putting x and y into arrays

plt.plot(x,g,label="g") # plotting y over x
plt.plot(x,v,label="v") 
#Plotting stuff:
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("Problem_7_plot.pdf")