from pyarma import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import seaborn as sns
import pandas as pd

A = cube()
A.load("Experiment_3.bin", arma_binary)
print(size(A))
P = np.array(A)

u = cube()
u.load('Experiment_3_Re.bin', arma_binary)
u_real = np.array(u)

u.load('Experiment_3_Im.bin', arma_binary)
u_imag = np.array(u)

# B = mat()
# B.load("pot.bin", arma_binary)
# V = np.array(B)

sns.set_theme()

########### Plot sum of probabilities over time:

A = cube()
A.load("Experiment_1.bin", arma_binary)
P_zero = np.array(A)

A.load("Experiment_2.bin", arma_binary)
P_double = np.array(A)

# Array of time points
dt = 2.5e-5
T = 0.008
t_points = np.arange(0, T, dt)

prob = np.abs(np.sum(P_zero, axis=(1, 2)) - 1)
sns.lineplot(y = prob, x = t_points, label='No slits')

prob = np.abs(np.sum(P_double, axis=(1, 2)) - 1)
sns.lineplot(y = prob, x = t_points, label='2 slits')

plt.xlabel('Time')
plt.ylabel(f'$| |p(x, y; t)|^2 - 1 |$')

plt.show()


########### Plot colour maps:
h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)

# Array of time points
T = 0.002
t_points = np.arange(0, T, dt)

# Some settings
fontsize = 12
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]


time_points = np.array([0, 0.001, 0.0015])

for t in time_points:

    ### Plot colourmap of probability:

    t_idx = int(t / dt)

    # Create figure
    fig = plt.figure(figsize=(8, 7))

    ax = plt.gca()
    ax.grid(False)

    norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(P[t_idx, :, :]))

    # Plot the first frame
    img = ax.imshow(P[t_idx, :, :], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("hot"), norm=norm)
    img.set_norm(norm)
    # Axis labels
    plt.xlabel("$x$", fontsize=fontsize)
    plt.ylabel("$y$", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(f"$p(x,y; t = ${t})", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    plt.show()

    ### Plot real part:

    # Create figure
    fig = plt.figure(figsize=(8, 7))

    ax = plt.gca()
    ax.grid(False)

    norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(P[t_idx, :, :]))

    # Plot the first frame
    img = ax.imshow(u_real[t_idx, :, :], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("hot"), norm=norm)
    img.set_norm(norm)
    # Axis labels
    plt.xlabel("$x$", fontsize=fontsize)
    plt.ylabel("$y$", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(f"Real part of $u(x,y; t = {t})$", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    plt.show()

    ### Plot imaginary part:

    # Create figure
    fig = plt.figure(figsize=(8, 7))

    ax = plt.gca()
    ax.grid(False)

    norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(P[t_idx, :, :]))

    # Plot the first frame
    img = ax.imshow(u_imag[t_idx, :, :], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("hot"), norm=norm)
    img.set_norm(norm)
    # Axis labels
    plt.xlabel("$x$", fontsize=fontsize)
    plt.ylabel("$y$", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(f"Imaginary part of $u(x,y; t = {t})$", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    plt.show()


########### Plot "detector screen" for different slit numbers:

A = cube()
A.load("Experiment_3.bin", arma_binary)
P_double = np.array(A)

A.load("Experiment_4.bin", arma_binary)
P_single = np.array(A)

A.load("Experiment_5.bin", arma_binary)
P_triple = np.array(A)

# Find indices corresponding to t = 0.002, x = 0.8
x_coord = 0.8
t = 0.001    # change this!
t_idx = -1
x_idx = int(x_coord / h)

# Set plotting parameters
palette = sns.color_palette("colorblind")
sns.set_palette(palette=palette, desat=0.9, color_codes=True)

# Create figure
plt.figure(figsize=(12, 7))

x_points = np.linspace(y_min, y_max, len(P[t_idx, x_idx, :]))

sns.lineplot(y=P_double[t_idx, x_idx, :], x=x_points, label='Two slits')
sns.lineplot(y=P_single[t_idx, x_idx, :], x=x_points, label='One slit')
sns.lineplot(y=P_triple[t_idx, x_idx, :], x=x_points, label='Three slits')

plt.xlabel('$y$')
plt.ylabel(f'$p(y | x = 0.8; t = {t})$')
plt.show()

# fig = sns.heatmap(data = P[0, :, :])
# fig.set_xticks(np.linspace(0, 1, len(x_points)))

# for idx, (label, label2) in enumerate(zip(fig.get_xticklabels(), fig.get_yticklabels())):
#     if idx % 2 == 0:
#         label.set_visible(True)
#         label2.set_visible(True)
#     else:
#         label.set_visible(False)
#         label2.set_visible(False)

# plt.show()
