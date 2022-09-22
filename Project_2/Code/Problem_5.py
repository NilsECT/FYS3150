import matplotlib.pyplot as plt
import numpy as np


N_val, num_trans = np.loadtxt("Problem_5.txt", unpack=True)

plt.title("Number of transformations needed for an NxN tridiagonal matrix")
plt.plot(N_val, num_trans, label="Number of transformations")
plt.xlabel("N")
plt.ylabel("Number of transformations")
plt.legend()
plt.savefig("Problem_5_plot.pdf")







