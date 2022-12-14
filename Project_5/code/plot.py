from pyarma import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams['font.size'] = '16'

sns.set_theme("notebook", "whitegrid")

def turn_the_lights_down_low():
    """
    Function for setting the plotting environment equal for all plots
    """

    context = {'font.size': 22.0,
                'axes.labelsize': 22.0,
                'axes.titlesize': 22.0,
                'xtick.labelsize': 22.0,
                'ytick.labelsize': 22.0,
                'legend.fontsize': 22.0,
                'legend.title_fontsize': None,
                'axes.linewidth': 0.8,
                'grid.linewidth': 0.8,
                'lines.linewidth': 1.5,
                'lines.markersize': 6.0,
                'patch.linewidth': 1.0,
                'xtick.major.width': 0.8,
                'ytick.major.width': 0.8,
                'xtick.minor.width': 0.6,
                'ytick.minor.width': 0.6,
                'xtick.major.size': 3.5,
                'ytick.major.size': 3.5,
                'xtick.minor.size': 2.0,
                'ytick.minor.size': 2.0}

    plt.rcParams['text.usetex'] = True
    sns.set_theme(context=context, palette="colorblind", font="sans-serif", font_scale=1)
    sns.set_style("whitegrid", {'axes.linewidth': 2, 'axes.edgecolor':'black'})

turn_the_lights_down_low()

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
plt.xlabel('Time')
plt.ylabel(f'$| |p(x, y; t)|^2 - 1 |$')
plt.tick_params(axis='both', which='major', labelsize=fontsize-1)
plt.savefig('probability.pdf')


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

    for i, t in enumerate(time_points):

        t_idx = int(t / dt) # index for current time-step

        axs[i].grid(False)

        norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(probability[t_idx, :, :]))
        img = axs[i].imshow(probability[t_idx, :, :], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("hot"))
        img.set_norm(norm)

        cbar = fig.colorbar(img, ax=axs[i],orientation="horizontal", pad=.2, fraction = 0.05)

        cbar.set_label(f'$p(x,y; t)$')
        cbar.ax.tick_params(labelsize=fontsize)
        axs[i].tick_params(axis='both', which='major', labelsize=fontsize)

        axs[i].set_xlabel("$x$")
        axs[i].set_ylabel("$y$")

        axs[i].title.set_text(f'$t = {t}$')

    plt.tight_layout()
    # Make space for title
    plt.subplots_adjust(top=0.99)

    plt.savefig(f'colourmap_{label}.pdf')


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
t_idx = int(t / dt)
x_idx = int(x_coord / h)


# Set plotting parameters
palette = sns.color_palette("colorblind")
sns.set_palette(palette=palette, desat=0.9, color_codes=True)

# Create figure
plt.figure(figsize=(12, 7))

x_points = np.linspace(y_min, y_max, len(P_double[t_idx, :, x_idx]))

sns.lineplot(y=P_double[t_idx, :, x_idx] / np.sum(P_double[t_idx, :, x_idx]), x=x_points, label='Two slits')
sns.lineplot(y=P_single[t_idx, :, x_idx] / np.sum(P_single[t_idx, :, x_idx]), x=x_points, label='One slit')
sns.lineplot(y=P_triple[t_idx, :, x_idx] / np.sum(P_triple[t_idx, :, x_idx]), x=x_points, label='Three slits')

plt.xlabel('$y$')
plt.ylabel(f'$p(y | x = 0.8; t = {t})$')
plt.legend()
plt.savefig(f'detectionscreen.pdf')
