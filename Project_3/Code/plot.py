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
plt.legend()
plt.show()

plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb")
plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb")
plt.plot(data_2p_RK_x_hasC[0,0], data_2p_RK_y_hasC[0,0], "ro", label="Particle 1 initial position", markersize=10)
plt.legend()
plt.show()

ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb")
ax.plot3D(data_2p_RK_x[:,1], data_2p_RK_y[:,1], data_2p_RK_z[:,1], label="Particle 2: RK4, no coulomb")
plt.legend()
plt.show()

ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], data_2p_RK_z_hasC[:,0], label="Particle 1: RK4, has coulomb")
ax.plot3D(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], data_2p_RK_z_hasC[:,1], label="Particle 2: RK4, has coulomb")
plt.legend()
plt.show()

ax = plt.axes(projection = "3d")
ax.plot3D(data_2p_RK_x[:,0], data_2p_RK_y[:,0], data_2p_RK_z[:,0], label="Particle 1: RK4, no coulomb")
ax.plot3D(data_2p_FE_x[:,0], data_2p_FE_y[:,0], data_2p_FE_z[:,0], label="Particle 1: FE, no coulomb")
plt.legend()
plt.show()

plt.plot(data_2p_RK_x[:,0], data_2p_RK_y[:,0], label="Particle 1: RK4, no coulomb")
plt.plot(data_2p_RK_x[:,1], data_2p_RK_y[:,1], label="Particle 2: RK4, no coulomb")
plt.legend()
plt.show()

plt.plot(data_2p_RK_x_hasC[:,0], data_2p_RK_y_hasC[:,0], label="Particle 1: RK4, has coulomb")
plt.plot(data_2p_RK_x_hasC[:,1], data_2p_RK_y_hasC[:,1], label="Particle 2: RK4, has coulomb")
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