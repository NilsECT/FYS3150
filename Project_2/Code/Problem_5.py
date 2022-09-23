import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress


N_val, num_trans = np.loadtxt("Problem_5.txt", unpack=True) # Getting the data

# Slicing fata:
slice = len(N_val)//2
linreg_model = linregress(x=np.log(N_val[slice:]), y=np.log(num_trans[slice:]))

# Plotting:
plt.title("Number of similarity transformations needed for an NxN tridiagonal matrix")
plt.axvline(N_val[slice], color='green', label='Point from which we slice the data')
plt.plot(N_val, num_trans, label="Number of transformations")
plt.loglog(N_val, np.exp(linreg_model.intercept)*N_val**linreg_model.slope, '--', label=f'Linear regression model, slope = {linreg_model.slope:.2f}')
plt.xlabel('$N$')
plt.title(f'Number of similarity transformations \nneeded for an NxN tridiagonal matrix')
plt.ylabel('Number of transformations')
plt.savefig('Problem_5_plot.pdf')
plt.legend()
plt.show()
