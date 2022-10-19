import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
cm = 1/2.54
plt.rcParams["figure.figsize"] = (12*cm,8*cm)
sns.set_theme("notebook","whitegrid")

## Data extraction with interaction for Forward Euler:
data_1p_FE_x = np.loadtxt("Forward_Euler1_n_0.010000_x.txt")
data_1p_FE_y = np.loadtxt("Forward_Euler1_n_0.010000_y.txt")
data_1p_FE_z = np.loadtxt("Forward_Euler1_n_0.010000_z.txt")

data_1p_RK4_x = np.loadtxt("RK41_n_0.010000_x.txt")
data_1p_RK4_y = np.loadtxt("RK41_n_0.010000_y.txt")
data_1p_RK4_z = np.loadtxt("RK41_n_0.010000_z.txt")

data_1p_a_x, data_1p_a_y, data_1p_a_z = np.loadtxt("analytical.txt",unpack=True)

plt.figure(1)
plt.plot(data_1p_FE_z,label="Forward Euler")
plt.plot(data_1p_RK4_z,label="Runge Kutta 4")
plt.plot(data_1p_a_z,label="Analytic")

plt.ylabel('$z$')
plt.xlabel('$t$')
plt.legend()
plt.savefig(f"1_particle_z_direction.pdf")




n = np.array([4000,8000,16000,32000])

def make_r(n):
    dt = 50/n
    data_1p_FE_x = np.loadtxt(f"FE1_n_{dt:.6f}_x.txt",unpack=True)
    data_1p_FE_y = np.loadtxt(f"FE1_n_{dt:.6f}_y.txt",unpack=True)
    data_1p_FE_z = np.loadtxt(f"FE1_n_{dt:.6f}_z.txt",unpack=True)

    data_1p_RK_y = np.loadtxt(f"RK1_n_{dt:.6f}_y.txt",unpack=True)
    data_1p_RK_z = np.loadtxt(f"RK1_n_{dt:.6f}_z.txt",unpack=True)
    data_1p_RK_x = np.loadtxt(f"RK1_n_{dt:.6f}_x.txt",unpack=True)
    r_F = np.zeros((n,3))
    for i in range(n):
        r_F[i] = np.array([data_1p_FE_x[i], data_1p_FE_y[i], data_1p_FE_z[i]])
    r_R = np.zeros((n,3))
    for i in range(n):
        r_R[i] = np.array([data_1p_RK_x[i], data_1p_RK_y[i], data_1p_RK_z[i]])
    return r_F, r_R

r_a = np.zeros(((len(data_1p_a_x)),3))
for i in range(len(data_1p_a_x)):
    r_a[i] = np.array([data_1p_a_x[i],data_1p_a_y[i],data_1p_a_z[i]])

r_f_1 , r_r_1 = make_r(n[0])
r_f_2 , r_r_2 = make_r(n[1])
r_f_3 , r_r_3 = make_r(n[2])
r_f_4 , r_r_4 = make_r(n[3])

max_error= r_a[::2.5]-r_f_1

    

