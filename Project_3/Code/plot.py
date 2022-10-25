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
    
# Plot of fist frequency scan:
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
"""

# 100 particles, RK4, has coulomb force
data_100p_RK_x_hasC = np.loadtxt("RK4100_y_0.010000_x.txt")
data_100p_RK_y_hasC = np.loadtxt("RK4100_y_0.010000_y.txt")
data_100p_RK_z_hasC = np.loadtxt("RK4100_y_0.010000_z.txt")

plt.figure()
ax = plt.axes(projection = "3d")
N = len(data_100p_RK_x_hasC[0,:])
for i in range(N): # , label="Particle 1: RK4, has coulomb"
    ax.plot3D(data_100p_RK_x_hasC[:,i], data_100p_RK_y_hasC[:,i], data_100p_RK_z_hasC[:,i], linewidth=0.4)
    # ax.plot3D(data_100p_RK_x_hasC[0,i], data_100p_RK_y_hasC[0,i], data_100p_RK_z_hasC[0,i], label="Particle 1 initial position", markersize=7)

ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_rk4_with_and_without_coulomb_3d.pdf")
plt.show()

"""
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

plt.figure()
plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x_hasC[0,0], data_2p_RK_y_hasC[0,0], "ro", label="Particle 1 initial position", markersize=7)
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.axis('equal')
plt.legend()
plt.savefig("p1_rk4_with_and_without_coulomb_2d.pdf")

plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
ax.plot3D(data_2p_RK_x[:,1], data_2p_RK_y[:,1], data_2p_RK_z[:,1], label="Particle 2: RK4, no coulomb", linewidth=0.8)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_and_p2_rk4_without_coulomb_3d.pdf")

plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], data_2p_RK_z_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
ax.plot3D(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], data_2p_RK_z_hasC[:,1], label="Particle 2: RK4, has coulomb", linewidth=0.8)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_and_p2_rk4_with_coulomb_3d.pdf")

plt.figure()
ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
ax.plot3D(data_2p_FE_x[:,0], data_2p_FE_y[:,0], data_2p_FE_z[:,0], label="Particle 1: FE, no coulomb", linewidth=0.8)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.savefig("p1_rk4_and_fe_without_coulomb_3d.pdf")

plt.figure()
plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x[:,1], data_2p_RK_y[:,1], label="Particle 2: RK4, no coulomb", linewidth=0.8)
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.axis('equal')
plt.legend()
plt.savefig("p1_and_p2_rk4_without_coulomb_2d.pdf")

plt.figure()
plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb", linewidth=0.8)
plt.plot(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], label="Particle 2: RK4, has coulomb", linewidth=0.8)
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.axis('equal')
plt.legend()
plt.savefig("p1_and_p2_rk4_with_coulomb_2d.pdf")


# plt.plot(data_2p_FE_x_hasC[:,0], data_2p_FE_y_hasC[:,0], label="Particle 1: FE, has coulomb")
# plt.plot(data_2p_FE_x[:,0], data_2p_FE_y[:,0], label="Particle 1: FE, no coulomb")
# plt.plot(data_2p_FE_x[0,0], data_2p_FE_y[0,0], "ro", label="Particle 1 initial position", markersize=10)
# plt.plot(data_2p_FE_x[1,1], data_2p_FE_y[1,1], "ro", label="Particle 1 initial position", markersize=10)
# plt.legend()
# plt.show()


# plt.plot(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], label="Particle 2: RK4, has coulomb")
# plt.plot(data_2p_FE_x[:,1], data_2p_RK_y[:,1], label="Particle 2: RK4, no coulomb")
# plt.plot(data_2p_FE_x_hasC[:,1], data_2p_FE_y_hasC[:,1], label="Particle 2: FE, has coulomb")
# plt.plot(data_2p_FE_x[:,1], data_2p_FE_y[:,1], label="Particle 2: FE, no coulomb")

# plt.legend()
# plt.show()

########################################################################
#PHASE SPACE PLOTS
########################################################################

p1 = "Particle 1"
p2 = "Particle 2"

# without interaction

## Data extraction with interaction:
#fe_x = np.loadtxt("Forward_Euler2_n_0.010000_x.txt")
#e_y = np.loadtxt("Forward_Euler2_n_0.010000_y.txt")
#fe_y = np.loadtxt("Forward_Euler2_n_0.010000_z.txt")

#fe_vx = np.loadtxt("Forward_Euler2_n_0.010000_vx.txt")
#fe_vy = np.loadtxt("Forward_Euler2_n_0.010000_vy.txt")
#fe_vz = np.loadtxt("Forward_Euler2_n_0.010000_vz.txt")

RK4_x = np.loadtxt("RK42_n_0.010000_x.txt")
#RK4_y = np.loadtxt("RK42_n_0.010000_y.txt")
RK4_z = np.loadtxt("RK42_n_0.010000_z.txt")
N = len(RK4_x[:, 0])

RK4_vx = np.loadtxt("RK42_n_0.010000_vx.txt")
#RK4_vy = np.loadtxt("RK42_n_0.010000_vy.txt")
RK4_vz = np.loadtxt("RK42_n_0.010000_vz.txt")

#cut = -int(np.ceil(N/2))

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
# plt.show()

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
# plt.show()

# With Coulomb

plt.figure()
RK4_x = np.loadtxt("RK42_y_0.010000_x.txt")
#RK4_y = np.loadtxt("RK42_y_0.010000_y.txt")
RK4_z = np.loadtxt("RK42_y_0.010000_z.txt")
N = len(RK4_x[:, 0])

RK4_vx = np.loadtxt("RK42_y_0.010000_vx.txt")
#RK4_vy = np.loadtxt("RK42_y_0.010000_vy.txt")
RK4_vz = np.loadtxt("RK42_y_0.010000_vz.txt")

cut = -int(np.ceil(N/2))

plt.plot(RK4_x[:, 0], RK4_vx[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_x[0, 0], RK4_vx[0, 0], 30, 'red', label=p1)

plt.plot(RK4_x[:, 1], RK4_vx[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_x[0, 1], RK4_vx[0, 1], 30, 'blue', label=p2)

plt.xlabel("position x")
plt.ylabel("velocity in x-direction")

plt.axis('equal')
plt.legend()
plt.savefig("RK4_phase_space_x_y.pdf")
# plt.show()

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
# plt.show()
