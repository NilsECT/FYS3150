import numpy as np
import matplotlib.pyplot as plt
from matplotlib.markers import *
import seaborn as sns
# cm = 1/2.54
# plt.rcParams["figure.figsize"] = (12*cm,8*cm)
sns.set_theme("notebook","whitegrid")
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

## Data extraction with interaction for Forward Euler:
dt = 100/10000
data_1p_Forward_Euler_x = np.loadtxt(f"Forward_Euler1_n_{dt:.6f}_x.txt")
data_1p_Forward_Euler_y = np.loadtxt(f"Forward_Euler1_n_{dt:.6f}_y.txt")
data_1p_Forward_Euler_z = np.loadtxt(f"Forward_Euler1_n_{dt:.6f}_z.txt")

data_1p_RK4_x = np.loadtxt(f"RK41_n_{dt:.6f}_x.txt")
data_1p_RK4_y = np.loadtxt(f"RK41_n_{dt:.6f}_y.txt")
data_1p_RK4_z = np.loadtxt(f"RK41_n_{dt:.6f}_z.txt")

data_4000_1p_a_x, data_4000_1p_a_y, data_4000_1p_a_z = np.loadtxt("4000_analytical.txt",unpack=True)
data_8000_1p_a_x, data_8000_1p_a_y, data_8000_1p_a_z = np.loadtxt("8000_analytical.txt",unpack=True)
data_10000_1p_a_x, data_10000_1p_a_y, data_10000_1p_a_z = np.loadtxt("10000_analytical.txt",unpack=True)
data_16000_1p_a_x, data_16000_1p_a_y, data_16000_1p_a_z = np.loadtxt("16000_analytical.txt",unpack=True)
data_32000_1p_a_x, data_32000_1p_a_y, data_32000_1p_a_z = np.loadtxt("32000_analytical.txt",unpack=True)
t = np.linspace(0,100,10000)

plt.figure(1)
plt.plot(t,data_1p_Forward_Euler_z,label="Forward Euler", c = "r")
plt.plot(t,data_1p_RK4_z,label="Runge Kutta 4", c = "b")
plt.plot(t,data_10000_1p_a_z,label="Analytic", c ="g")

plt.ylabel('z [$\\mu$m]')
plt.xlabel('t [$\\mu$s]')
plt.legend()
plt.savefig(f"1_particle_z_direction.pdf")



def error_plot(num_time_steps, method):

    # take out analytical datapoints
    data_1p_a_x, data_1p_a_y, data_1p_a_z = np.loadtxt(f"{num_time_steps}_analytical.txt",unpack=True)

    dt = 50/num_time_steps

    # take out either FE or RK4 datapoints
    data_1p_x = np.loadtxt(f"{method}1_n_{dt:.6f}_x.txt",unpack=True)
    data_1p_y = np.loadtxt(f"{method}1_n_{dt:.6f}_y.txt",unpack=True)
    data_1p_z = np.loadtxt(f"{method}1_n_{dt:.6f}_z.txt",unpack=True)

    # storage for the fe/rk4 positions
    r_F = np.zeros((num_time_steps,3))
    for i in range(num_time_steps):
        r_F[i] = np.array([data_1p_x[i], data_1p_y[i], data_1p_z[i]])


    r_a = np.zeros((num_time_steps,3))
    for i in range(num_time_steps):
        r_a[i] = np.array([data_1p_a_x[i],data_1p_a_y[i],data_1p_a_z[i]])

    max_error = r_a-r_F

    error = np.zeros(num_time_steps)
    t = np.linspace(0,50,num_time_steps)
    error = np.zeros(num_time_steps)
    Delta_max_error = np.zeros(num_time_steps)

    for i in range(num_time_steps):
        r_F = np.array([data_1p_x[i], data_1p_y[i], data_1p_z[i]])
        r_a = np.array([data_1p_a_x[i],data_1p_a_y[i],data_1p_a_z[i]])
        max_error = np.linalg.norm(r_a-r_F)
        error[i] = max_error/np.linalg.norm(r_a)
        Delta_max_error[i] = np.linalg.norm(r_a - r_F)

    plt.plot( t, error,".",label=f"n_step ={num_time_steps}")
    return np.amax(Delta_max_error)

## Error Plot:
n = np.array([4000,8000,16000,32000])
delta_max_error_Forward_Euler = np.zeros(4)
delta_max_error_RK4 = np.zeros(4)

plt.figure(2)
plt.subplot(1,2,1)
plt.title("Error over time for Forward Euler ")
for i in range(len(n)):
    delta_max_error_Forward_Euler[i] = error_plot(n[i],"Forward_Euler")

plt.xlabel("t [$\\mu$s]")
plt.ylabel("y [$\\mu$m]")
plt.yscale('log')
plt.legend()

plt.subplot(1,2,2)
plt.title("and Runge Kutta 4")
for i in range(len(n)):
    delta_max_error_RK4[i] = error_plot(n[i], "RK4")

plt.xlabel("t [$\\mu$s]")
plt.ylabel("y [$\\mu$m]")
plt.yscale('log')
plt.legend()

plt.savefig("error_plot.pdf")

## Error convergence rate for Forward Euler
r_err_Forward_Euler = 0
for i in range(1,4):
    r_err_Forward_Euler += np.log(delta_max_error_Forward_Euler[i]/delta_max_error_Forward_Euler[i-1])/np.log((50/n[i])/(50/n[i-1]))

print(f"Error convergence rate for Forward Euler: {(1/3)*r_err_Forward_Euler}")

## Error convergence rate for Runge Kutta 4
r_err_RK4 = 0
for i in range(1,4):
    r_err_RK4 += np.log(delta_max_error_RK4[i]/delta_max_error_RK4[i-1])/np.log((50/n[i])/(50/n[i-1]))
print(f"Error convergence rate for Runge Kutta 4: {(1/3)*r_err_RK4}")
