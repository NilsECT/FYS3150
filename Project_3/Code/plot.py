import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# cm = 1 / 2.54
# plt.rcParams["figure.figsize"] = (12 * cm, 8 * cm)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

sns.set_theme("notebook", "whitegrid")

def Perturbation_plot(f,coulomb_force,zoom):
    N_particles = 50
    data = np.loadtxt(f"Perturbation_{f}00000_{N_particles}_{coulomb_force}.txt")
    particles_trapped = (N_particles-data[:, 1])/N_particles
    plt.figure()
    plt.plot(data[:, 0], particles_trapped)
    plt.ylabel("Fraction of particles still trapped")
    plt.xlabel("Angular frequency $\omega_V$ [MHz]")
    plt.savefig(f"Perturbation_{f}_{N_particles}_{coulomb_force}{zoom}.pdf")

# Plot of frequency scan:
f = [0.1,0.4,0.7] # Amplitudes
coulomb_force = 0
zoom = ""
for i in f:
    Perturbation_plot(f,coulomb_force,zoom)
# Plot of more detailed scan.
## Without coulomb force:
f = 0.1
zoom = "zoom"
Perturbation_plot(f,coulomb_force,zoom)
## With coulomb force:
coulomb_force = 1
zoom = "zoom"
Perturbation_plot(f,coulomb_force,zoom)

# Forward Euler, has coulomb force
data_2p_FE_x_hasC = np.loadtxt("Forward_Euler2_y_0.010000_x.txt")
data_2p_FE_y_hasC = np.loadtxt("Forward_Euler2_y_0.010000_y.txt")
data_2p_FE_z_hasC = np.loadtxt("Forward_Euler2_y_0.010000_z.txt")

# Forward Euler, no coulomb force
data_2p_FE_x = np.loadtxt("Forward_Euler2_n_0.010000_x.txt")
data_2p_FE_y = np.loadtxt("Forward_Euler2_n_0.010000_y.txt")
data_2p_FE_z = np.loadtxt("Forward_Euler2_n_0.010000_z.txt")

# Runge Kutta, has coulomb force
data_2p_RK_x_hasC = np.loadtxt("RK42_y_0.010000_x.txt")
data_2p_RK_y_hasC = np.loadtxt("RK42_y_0.010000_y.txt")
data_2p_RK_z_hasC = np.loadtxt("RK42_y_0.010000_z.txt")

# Runge Kutta, no coulomb force
data_2p_RK_x = np.loadtxt("RK42_n_0.010000_x.txt")
data_2p_RK_y = np.loadtxt("RK42_n_0.010000_y.txt")
data_2p_RK_z = np.loadtxt("RK42_n_0.010000_z.txt")

# One particle with and without interaction in 3D:
plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], data_2p_RK_z_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
ax.plot3D(data_2p_RK_x_hasC[0,0], data_2p_RK_y_hasC[0,0], data_2p_RK_z_hasC[0,0], "ro", label="Particle 1 initial position", markersize=7)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_rk4_with_and_without_coulomb_3d.pdf")

# One particles Runge-Kutta with and without interaction xy-plane:
plt.figure()
plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x_hasC[0,0], data_2p_RK_y_hasC[0,0], "ro", label="Particle 1 initial position", markersize=7)
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.axis('equal')
plt.legend()
plt.savefig("p1_rk4_with_and_without_coulomb_2d.pdf")

# Two particles without interaction 3D:
plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
ax.plot3D(data_2p_RK_x[:,1], data_2p_RK_y[:,1], data_2p_RK_z[:,1], label="Particle 2: RK4, no coulomb", linewidth=0.8)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_and_p2_rk4_without_coulomb_3d.pdf")

# Two particles with interaction 3D:
plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], data_2p_RK_z_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
ax.plot3D(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], data_2p_RK_z_hasC[:,1], label="Particle 2: RK4, has coulomb", linewidth=0.8)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_and_p2_rk4_with_coulomb_3d.pdf")

# Two particles without interaction for both forward Euler and Runge-Kutta 3D:
plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
ax.plot3D(data_2p_FE_x[:,0], data_2p_FE_y[:,0], data_2p_FE_z[:,0], label="Particle 1: FE, no coulomb", linewidth=0.8)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_rk4_and_fe_without_coulomb_3d.pdf")

# Two particles without interaction for Runge-Kutta xy-plane:
plt.figure()
plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x[:,1], data_2p_RK_y[:,1], label="Particle 2: RK4, no coulomb", linewidth=0.8)
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.axis('equal')
plt.legend()
plt.savefig("p1_and_p2_rk4_without_coulomb_2d.pdf")

# Two particles with interaction for Runge-Kutta xy-plane:
plt.figure()
plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], label="Particle 2: RK4, has coulomb", linewidth=0.8)
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.axis('equal')
plt.legend()
plt.savefig("p1_and_p2_rk4_with_coulomb_2d.pdf")

########################################################################
#PHASE SPACE PLOTS
########################################################################

p1 = "Particle 1"
p2 = "Particle 2"

# without interaction
RK4_x = np.loadtxt("RK42_n_0.010000_x.txt")
RK4_z = np.loadtxt("RK42_n_0.010000_z.txt")
N = len(RK4_x[:, 0])

RK4_vx = np.loadtxt("RK42_n_0.010000_vx.txt")
RK4_vz = np.loadtxt("RK42_n_0.010000_vz.txt")

# Runge-Kutta phase plot for x-direction Two particles:
plt.figure()
plt.plot(RK4_x[:, 0], RK4_vx[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_x[0, 0], RK4_vx[0, 0], 30, 'red', label=p1)
plt.plot(RK4_x[:, 1], RK4_vx[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_x[0, 1], RK4_vx[0, 1], 30, 'blue', label=p2)
plt.xlabel("position x")
plt.ylabel("velocity in x-direction")
plt.axis('equal')
plt.legend()
plt.savefig("RK4_phase_space_x_n.pdf")

# Runge-Kutta phase plot for z-direction:
plt.figure()
plt.plot(RK4_z[:, 0], RK4_vz[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_z[0, 0], RK4_vz[0, 0], 30, 'red', label=p1)
plt.plot(RK4_z[:, 1], RK4_vz[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_z[0, 1], RK4_vz[0, 1], 30, 'blue', label=p2)
plt.xlabel("position z")
plt.ylabel("velocity in z-direction")
plt.axis('equal')
plt.legend()
plt.savefig("RK4_phase_space_z_n.pdf")

# With Coulomb

RK4_x = np.loadtxt("RK42_y_0.010000_x.txt")
RK4_z = np.loadtxt("RK42_y_0.010000_z.txt")
N = len(RK4_x[:, 0])
RK4_vx = np.loadtxt("RK42_y_0.010000_vx.txt")
RK4_vz = np.loadtxt("RK42_y_0.010000_vz.txt")
cut = -int(np.ceil(N/2))

# Runge-Kutta phase plot with interaction x-direction:
plt.figure()
plt.plot(RK4_x[:, 0], RK4_vx[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_x[0, 0], RK4_vx[0, 0], 30, 'red', label=p1)
plt.plot(RK4_x[:, 1], RK4_vx[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_x[0, 1], RK4_vx[0, 1], 30, 'blue', label=p2)
plt.xlabel("position x")
plt.ylabel("velocity in x-direction")
plt.axis('equal')
plt.legend()
plt.savefig("RK4_phase_space_x_y.pdf")

# Runge-Kutta phase plot for z-direction with interaction:
plt.figure()
plt.plot(RK4_z[:, 0], RK4_vz[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_z[0, 0], RK4_vz[0, 0], 30, 'red', label=p1)
plt.plot(RK4_z[:, 1], RK4_vz[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_z[0, 1], RK4_vz[0, 1], 30, 'blue', label=p2)
plt.xlabel("position z")
plt.ylabel("velocity in z-direction")
plt.axis('equal')
plt.legend()
plt.savefig("RK4_phase_space_z_y.pdf")
