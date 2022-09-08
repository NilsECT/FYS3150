import numpy as np
import matplotlib.pyplot as plt

x,u = np.loadtxt("Problem_2.txt", unpack=True) #numpy's txt-reading
x1, v1 = np.loadtxt("out10.txt", unpack=True)
x2, v2 = np.loadtxt("out100.txt", unpack=True)
x3, v3 = np.loadtxt("out1000.txt", unpack=True)

# Compute numpy arrays with the absolute error, relative error and log10(relative error)
def error(d2u_approx, d2u_exact):
    n = int((len(d2u_exact))/(len(d2u_approx)-1))
    print(n)
    d2u_exact_1 = np.append(d2u_exact[::n],d2u_exact[-1])
    print(len(d2u_exact_1))
    abs_err = np.absolute(d2u_approx - d2u_exact_1)
    rel_err = np.absolute(abs_err / d2u_exact_1)
    log10_abs_err = np.log10(abs_err)
    log10_rel_err = np.log10(rel_err)

    x_1 = np.linspace(0,1,len(d2u_exact_1))
    
    return x_1, log10_abs_err, log10_rel_err

#Def error arrays:
x_1, log10_abs_err_10, log10_rel_err_10 = error(v1,u)
x_2, log10_abs_err_100, log10_rel_err_100 = error(v2,u)
x_3, log10_abs_err_1000, log10_rel_err_1000 = error(v3,u)

#Plot abs erro:
plt.figure(1)
plt.plot(x_1, log10_abs_err_10,label="n=10")
plt.plot(x_2, log10_abs_err_100,label="n=100")
plt.plot(x_3, log10_abs_err_1000,label="n=1000")

plt.xlabel("x")
plt.ylabel("log(abs(v-u))")
plt.legend
plt.savefig("log10_abs_err.pdf")

#Plot rel error:
plt.figure(2)
plt.plot(x_1, log10_rel_err_10,label="n=10")
plt.plot(x_2, log10_rel_err_100,label="n=100")
plt.plot(x_3, log10_rel_err_1000,label="n=1000")

plt.xlabel("x")
plt.ylabel("log(abs(v-u)/u))")
plt.legend
plt.savefig("log10_rel_err.pdf")


