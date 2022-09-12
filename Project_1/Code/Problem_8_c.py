import numpy as np
import matplotlib.pyplot as plt

x,u = np.loadtxt("Problem_2.txt", unpack=True) #numpy's txt-reading

def error(d2u_approx, d2u_exact):
    n = int((len(d2u_exact))/(len(d2u_approx)-1))
    d2u_exact_1 = np.empty(0)
    d2u_exact_1 = np.append(d2u_exact[::n],d2u_exact[-1])
    
    abs_err = np.absolute(d2u_approx - d2u_exact_1)
    rel_err = np.absolute(abs_err / d2u_exact_1)
    x_1 = np.linspace(0,1,len(d2u_exact_1))
    
    return x_1, rel_err

max_rel = []
for i in range(1,8):
    x,v = np.loadtxt(f"out_e{i}.txt",unpack=True)
    x_1, rel_err = error(v,u)
    a = np.nanmax(np.abs(rel_err[:-1]))
    max_rel.append(a)

x = np.logspace(1,7,7)

plt.loglog(x,max_rel)
plt.title("Plot of the maximum relative error $\\max(\\epsilon_i)$ for each choice of $n_{steps}$")
plt.xlabel("$n_{steps}$")
plt.ylabel("relative_error")
plt.savefig("Problem_8_c_plot.pdf")



