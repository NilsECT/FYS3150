import numpy as np
import matplotlib.pyplot as plt
import sys

# Plot problem 6:
data = np.loadtxt(f'histogram.txt', skiprows=1, delimiter=',').T
seed, temperature, energy, freq = data

for T in [1, 2.4]:
    filter = temperature == T
    energies = energy[filter]
    freqs = freq[filter]

    plt.plot(energies, freqs, lw=.4)
    plt.title(f'Histogram of $\epsilon$ for T = {T} $J/k_B$')
    plt.xlabel('$\epsilon  [1/J] $')
    plt.ylabel('Frequency, normalized')
    plt.show()
