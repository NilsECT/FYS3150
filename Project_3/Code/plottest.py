import sys
import numpy as np
import matplotlib.pyplot as plt

N, inter, h = sys.argv[1:]

filename_prefix = f'{N}_{inter}_{h}_'

## Data extraction with interaction:
x = np.loadtxt(filename_prefix + "x.txt")
y = np.loadtxt(filename_prefix + "y.txt")
y = np.loadtxt(filename_prefix + "z.txt")

vx = np.loadtxt(filename_prefix + "vx.txt")
vy = np.loadtxt(filename_prefix + "vy.txt")
vz = np.loadtxt(filename_prefix + "vz.txt")

plt.plot(x[:,:-1], y[:,:-1],"blue")
plt.plot(x[:,-1], y[:,-1],"red")
plt.ylabel('$y$')
plt.xlabel('$x$')
plt.title('Test!')
plt.savefig(f"{N}_particle_{inter}_interacting.pdf")


