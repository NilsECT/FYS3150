import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def set_params():
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

plotcolors = ['red', 'orange', 'black']
sns.set_theme("notebook", "whitegrid")

# system parameters:
q = 1
V0 = 2.41e6
B0 = 9.65*10
m = 40.078
d = 500

omega_z = np.sqrt(2*q*V0/(m*d*d))
omega_0 = q*B0/m

fig, ax = plt.subplots(nrows=3, figsize=(8, 7))
set_params()

# time-varying applied potential for amplitudes f = 0.1, 0.4, 0.7:
for i, f in enumerate([0.1, 0.4, 0.7]):
    N_particles = 25
    coulomb_force = 0

    # load data containing fraction of particles for each amplitude:
    data = np.loadtxt(f"Perturbation_{f}00000_{N_particles}_{coulomb_force}.txt")
    particles_trapped = (N_particles-data[:, 1])/N_particles

    ax[i].plot(data[:, 0], particles_trapped, color='black', lw=.8)#label=f'$f$ = {f}', lw=.8)

    # add vertical lines corresponding to resonance frequencies:
    for j in range(3):
        ax[j].axvline(omega_z*(i+1), linestyle='dashed', color=plotcolors[i-1],
            label=f'{i+1} $\cdotp \omega_z$')

    ax[i].set_ylabel(f"$f=$ {f}")
    plt.xlabel("Angular frequency $\omega_V$ [MHz]")
    plt.legend()

plt.savefig(f"Perturbation_{N_particles}_{coulomb_force}.pdf")

# plot the "zoom-in" on angular frequencies, f = 0.1:
for i, f in enumerate([0.1]):
    N_particles = 50

    plt.figure()
    set_params()

    # with particle interaction:
    coulomb_force = 1
    data = np.loadtxt(f"Perturbation_{f}00000_{N_particles}_{coulomb_force}.txt")
    particles_trapped = (N_particles-data[:, 1])/N_particles
    plt.plot(data[:, 0], particles_trapped, color='red', label=f'With Coulomb interaction', lw=.9)

    # without particle interaction:
    coulomb_force = 0
    data = np.loadtxt(f"Perturbation_{f}00000_{N_particles}_{coulomb_force}.txt")
    particles_trapped = (N_particles-data[:, 1])/N_particles
    plt.plot(data[:, 0], particles_trapped, color='black', label=f'No Coulomb interaction', lw=.9)

    # add vertical line for resonance frequency:
    plt.axvline(omega_z*2, linestyle='dashed', color=plotcolors[1],
        label=f'$\omega_V$ = 2 $\cdotp \omega_z$')

    plt.ylabel("Fraction of particles still trapped")
    plt.xlabel("Angular frequency $\omega_V$ [MHz]")
    plt.legend()
    plt.savefig(f"Perturbation_{f}_{N_particles}_Zoom.pdf")
