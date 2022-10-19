import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cm = 1 / 2.54
plt.rcParams["figure.figsize"] = (12 * cm, 8 * cm)
sns.set_theme("notebook", "whitegrid")


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