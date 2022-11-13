import numpy as np
import matplotlib.pyplot as plt
import sys

# Plot problem 6:
data_ordered = np.loadtxt(f'histogram_ordered.txt', skiprows=1, delimiter=',').T
data_unordered = np.loadtxt(f'histogram_unordered.txt', skiprows=1, delimiter=',').T
seed, temperature_o, energy_o, freq_o = data_ordered
seed, temperature_u, energy_u, freq_u = data_unordered

for T in [1, 2.4]:
    filter_o = temperature_o == T
    filter_u = temperature_u == T

    lw = 0.8
    plt.plot(energy_o[filter_o], freq_o[filter_o], label='Ordered', alpha=.5, lw=lw)
    plt.plot(energy_u[filter_u], freq_u[filter_u], label='Unordered', alpha=.5, lw=lw)
    plt.title(f'Histogram of $\epsilon$ for T = {T} $J/k_B$')
    plt.xlabel('$\epsilon  [1/J] $')
    plt.ylabel('Frequency, normalized')
    plt.show()
