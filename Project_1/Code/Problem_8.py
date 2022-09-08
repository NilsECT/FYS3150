import numpy as np
import matplotlib.pyplot as plt

x,u = np.loadtxt("Problem_2.txt", unpack=True) #numpy's txt-reading
x1, v1 = np.loadtxt("out10.txt", unpack=True)
x2, v2 = np.loadtxt("out100.txt", unpack=True)
x3, v3 = np.loadtxt("out1000.txt", unpack=True)
# Compute numpy arrays with the absolute error, relative error and log10(relative error)
def error(d2u_approx, d2u_exact, name):
    d2u_exact_1 = np.resize(d2u_exact,int(len(d2u_approx)))
    abs_err = np.abs(d2u_approx - d2u_exact_1)
    rel_err = np.abs(abs_err / d2u_exact_1)
    log10_rel_err = np.log10(rel_err)

    x_1 = np.resize(x,int(len(d2u_approx)))
    plt.plot(x_1, log10_rel_err, '--', c="0.8", linewidth=1.5)
    plt.plot(x_1, log10_rel_err, '.', c="black", markersize=10)
    # plt.title("Log10(relative error) " + common_title_string, fontsize=10)
    plt.ylabel("log10(relative error)")
    plt.xlabel("x")
    plt.savefig(f"log10_rel_err_vs_x_{name}.pdf")
    

error(v1,u,"n=10")


error(v2,u,"n=100")


error(v3,u,"n=1000")
