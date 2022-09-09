import numpy as np
import matplotlib.pyplot as plt

x,u = np.loadtxt("Problem_2.txt", unpack=True) #numpy's txt-reading


def error(d2u_approx, d2u_exact):
    n = int((len(d2u_exact))/(len(d2u_approx)-1))
    d2u_exact_1 = np.empty(0)
    d2u_exact_1 = np.append(d2u_exact[::n],d2u_exact[-1])
    
    abs_err = np.absolute(d2u_approx - d2u_exact_1)
    rel_err = np.absolute(abs_err / d2u_exact_1)
    # print(d2u_exact_1[-1])
    x_1 = np.linspace(0,1,len(d2u_exact_1))
    
    return x_1, rel_err


liste = ["out10.txt","out100.txt","out1000.txt","out1e4.txt","out1e5.txt","out1e6.txt","out1e7.txt"]
max_rel = []
# print(max_rel)
for filename in liste:
    x,v = np.loadtxt(filename,unpack=True)
    x_1, rel_err = error(v,u)
    a = np.nanmax(np.abs(rel_err[:-1]))
    max_rel.append(a)

    # print(a)

k = int(1e7)
x = np.logspace(1,7,7)

plt.loglog(x,max_rel)
plt.xlabel("n_step")
plt.ylabel("relativ_error")
plt.savefig("Problem_8_c_plot.svg")


