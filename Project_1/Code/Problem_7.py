import numpy as np
import matplotlib.pyplot as plt

x,u = np.loadtxt("Problem_2.txt", unpack=True) #numpy's txt-reading
x1, v1 = np.loadtxt("out_e1.txt", unpack=True)
x2, v2 = np.loadtxt("out_e2.txt", unpack=True)
x3, v3 = np.loadtxt("out_e3.txt", unpack=True)

plt.plot(x,u,label="u") # plotting y over x
plt.plot(x1,v1,"x",linewidth=0.4,label="n=10")
plt.plot(x2,v2,"x",linewidth=0.4,label="n=100")
plt.plot(x3,v3,"--",linewidth=0.4,label="n=1000")

#Plotting stuff:
plt.title("Exact solution $u(x)$ and approximate $v(x)$\n with increasing number of steps $n$")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.savefig("Problem_7_plot.pdf")