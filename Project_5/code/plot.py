from pyarma import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import seaborn as sns
import pandas as pd

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams['font.size'] = '16'

sns.set_theme("notebook", "whitegrid")

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
t_points = np.arange(0, T + dt, dt)

plt.figure(figsize=(8, 5))
fontsize = 14

prob = np.abs(np.sum(P_zero, axis=(1, 2)) - 1)
sns.lineplot(y = prob, x = t_points, label='No slits')

prob = np.abs(np.sum(P_double, axis=(1, 2)) - 1)
sns.lineplot(y = prob, x = t_points, label='2 slits')

plt.legend()
plt.xlabel('Time', fontsize=fontsize)
plt.ylabel(f'$| |p(x, y; t)|^2 - 1 |$', fontsize=fontsize)
plt.tick_params(axis='both', which='major', labelsize=fontsize-1)
plt.savefig('figs/probability.pdf')


########### Plot colour maps:
h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)

# Array of time points
T = 0.002
t_points = np.arange(0, T, dt)

# Some settings
fontsize = 15
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

##### Plot colourmaps:

time_points = [0, 0.001, 0.002]
filenames = ['3', '3_Re', '3_Im']

for filename, label in zip(filenames, ['full', 'real', 'imag']):

    A = cube()
    A.load(f'Experiment_{filename}.bin', arma_binary)
    probability = np.array(A)

    fig, axs = plt.subplots(1, 3, figsize=(16, 7))
    fontsize = 15

    for i, t in enumerate(time_points):

        t_idx = int(t / dt)

        axs[i].grid(False)

        norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(probability[t_idx, :, :]))
        img = axs[i].imshow(probability[t_idx, :, :], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("hot"))#, norm=norm)
        img.set_norm(norm)

        cbar = fig.colorbar(img, ax=axs[i],orientation="horizontal", pad=.2, fraction = 0.05)

        cbar.set_label(f'$p(x,y; t)$', fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)
        axs[i].tick_params(axis='both', which='major', labelsize=fontsize)

        axs[i].set_xlabel("$x$", fontsize=fontsize)
        axs[i].set_ylabel("$y$", fontsize=fontsize)


        axs[i].title.set_text(f'$t = {t}$')

    plt.tight_layout()
    # Make space for title
    plt.subplots_adjust(top=0.99)

    plt.savefig(f'figs/colourmap_{label}.pdf')


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
t = 0.002
<<<<<<< HEAD
t_idx = int(t / dt)
print(t_idx)
=======
>>>>>>> 14e28c4e9e26bd232cfa79306b56c88795ddb3a3
x_idx = int(x_coord / h)


# Set plotting parameters
palette = sns.color_palette("colorblind")
sns.set_palette(palette=palette, desat=0.9, color_codes=True)

# Create figure
plt.figure(figsize=(12, 7))

x_points = np.linspace(y_min, y_max, len(P_double[t_idx, x_idx, :]))

<<<<<<< HEAD
sns.lineplot(y=P_double[t_idx, x_idx, :] / np.sum(P_double[t_idx, x_idx, :]), x=x_points, label='Two slits')
sns.lineplot(y=P_single[t_idx, x_idx, :] / np.sum(P_single[t_idx, x_idx, :]), x=x_points, label='One slit')
sns.lineplot(y=P_triple[t_idx, x_idx, :] / np.sum(P_triple[t_idx, x_idx, :]), x=x_points, label='Three slits')
=======
sns.lineplot(y=P_double[t_idx, :, x_idx] / np.sum(P_double[t_idx, :, x_idx]), x=x_points, label='Two slits')
sns.lineplot(y=P_single[t_idx, :, x_idx] / np.sum(P_single[t_idx, :, x_idx]), x=x_points, label='One slit')
sns.lineplot(y=P_triple[t_idx, :, x_idx] / np.sum(P_triple[t_idx, :, x_idx]), x=x_points, label='Three slits')
>>>>>>> 14e28c4e9e26bd232cfa79306b56c88795ddb3a3

fontsize = 22
plt.xlabel('$y$', fontsize=fontsize)
plt.ylabel(f'$p(y | x = 0.8; t = {t})$', fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(f'figs/detectionscreen.pdf')

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
