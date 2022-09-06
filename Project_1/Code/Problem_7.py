import numpy as np
import matplotlib.pyplot as plt

x,u = np.loadtxt("Problem_2.txt", unpack=True) #numpy's txt-reading
x1, v1 = np.loadtxt("out10.txt", unpack=True)
x2, v2 = np.loadtxt("out100.txt", unpack=True)
x3, v3 = np.loadtxt("out1000.txt", unpack=True)
# x,g,v = data[:,0], data[:,1], data[:,2] #putting x and y into arrays

plt.plot(x,u,label="u") # plotting y over x
plt.plot(x1,v1,label="n=10")
plt.plot(x2,v2,label="n=100")
plt.plot(x3,v3,label="n=1000")
# plt.plot(x,v,label="v") 
#Plotting stuff:
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("Problem_7_plot.pdf")