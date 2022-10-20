import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cm = 1 / 2.54
plt.rcParams["figure.figsize"] = (12 * cm, 8 * cm)
sns.set_theme("notebook", "whitegrid")



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

ax = plt.axes(projection = "3d")

ax.plot3D(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], data_2p_RK_z_hasC[:,0], label="Particle 1: RK4, has coulomb")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb")
ax.plot3D(data_2p_RK_x_hasC[0,0], data_2p_RK_y_hasC[0,0], data_2p_RK_z_hasC[0,0], "ro", label="Particle 1 initial position", markersize=10)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.show()

plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb")
plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb")
plt.plot(data_2p_RK_x_hasC[0,0], data_2p_RK_y_hasC[0,0], "ro", label="Particle 1 initial position", markersize=10)
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.show()

ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb")
ax.plot3D(data_2p_RK_x[:,1], data_2p_RK_y[:,1], data_2p_RK_z[:,1], label="Particle 2: RK4, no coulomb")
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.show()

ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], data_2p_RK_z_hasC[:,0], label="Particle 1: RK4, has coulomb")
ax.plot3D(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], data_2p_RK_z_hasC[:,1], label="Particle 2: RK4, has coulomb")
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.show()

ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb")
ax.plot3D(data_2p_FE_x[:,0], data_2p_FE_y[:,0], data_2p_FE_z[:,0], label="Particle 1: FE, no coulomb")
ax.set_xlabel("x [$\\mu$m]")
ax.set_ylabel("y [$\\mu$m]")
ax.set_zlabel("z [$\\mu$m]")
plt.legend()
plt.show()

plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb")
plt.plot(data_2p_RK_x[:,1], data_2p_RK_y[:,1], label="Particle 2: RK4, no coulomb")
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.legend()
plt.show()

plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb")
plt.plot(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], label="Particle 2: RK4, has coulomb")
plt.xlabel("x [$\\mu$m]")
plt.ylabel("y [$\\mu$m]")
plt.legend()
plt.show()

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

plt.plot(RK4_x[:, 0], RK4_vx[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_x[0, 0], RK4_vx[0, 0], 30, 'red', label=p1)

plt.plot(RK4_x[:, 1], RK4_vx[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_x[0, 1], RK4_vx[0, 1], 30, 'blue', label=p2)

plt.xlabel("position x")
plt.ylabel("velocity in x-direction")

plt.legend()
plt.savefig("RK4_phase_space_x_n.pdf")
plt.show()

plt.plot(RK4_z[:, 0], RK4_vz[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_z[0, 0], RK4_vz[0, 0], 30, 'red', label=p1)

plt.plot(RK4_z[:, 1], RK4_vz[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_z[0, 1], RK4_vz[0, 1], 30, 'blue', label=p2)

plt.xlabel("position z")
plt.ylabel("velocity in z-direction")

plt.legend()
plt.savefig("RK4_phase_space_z_n.pdf")
plt.show()

# With Coulomb

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

plt.legend()
plt.savefig("RK4_phase_space_x_y.pdf")
plt.show()

plt.plot(RK4_z[:, 0], RK4_vz[:, 0], 'red', linewidth=0.8)
plt.scatter(RK4_z[0, 0], RK4_vz[0, 0], 30, 'red', label=p1)

plt.plot(RK4_z[:, 1], RK4_vz[:, 1], 'blue', linewidth=0.8)
plt.scatter(RK4_z[0, 1], RK4_vz[0, 1], 30, 'blue', label=p2)

plt.xlabel("position z")
plt.ylabel("velocity in z-direction")

plt.legend()
plt.savefig("RK4_phase_space_z_y.pdf")
plt.show()