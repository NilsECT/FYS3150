import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
# cm = 1/2.54
# plt.rcParams["figure.figsize"] = (12*cm,8*cm)
# sns.set_theme("notebook","whitegrid")

## Data extraction with interaction for Forward Euler:
data_1p_Forward_Euler_x = np.loadtxt("Forward_Euler1_n_0.010000_x.txt")
data_1p_Forward_Euler_y = np.loadtxt("Forward_Euler1_n_0.010000_y.txt")
data_1p_Forward_Euler_z = np.loadtxt("Forward_Euler1_n_0.010000_z.txt")

data_1p_RK4_x = np.loadtxt("RK41_n_0.010000_x.txt")
data_1p_RK4_y = np.loadtxt("RK41_n_0.010000_y.txt")
data_1p_RK4_z = np.loadtxt("RK41_n_0.010000_z.txt")

data_4000_1p_a_x, data_4000_1p_a_y, data_4000_1p_a_z = np.loadtxt("4000_analytical.txt",unpack=True)
data_8000_1p_a_x, data_8000_1p_a_y, data_8000_1p_a_z = np.loadtxt("8000_analytical.txt",unpack=True)
data_16000_1p_a_x, data_16000_1p_a_y, data_16000_1p_a_z = np.loadtxt("16000_analytical.txt",unpack=True)
data_32000_1p_a_x, data_32000_1p_a_y, data_32000_1p_a_z = np.loadtxt("32000_analytical.txt",unpack=True)

plt.figure(1)
plt.plot(data_1p_Forward_Euler_z,label="Forward Euler")
plt.plot(data_1p_RK4_z,label="Runge Kutta 4")
plt.plot(data_4000_1p_a_z,label="Analytic")

plt.ylabel('$z$')
plt.xlabel('$t$')
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

    for i in range(num_time_steps):
        error[i] = np.linalg.norm(max_error[i]) / (np.linalg.norm(r_a[i]))

    plt.plot( t, error, marker = "o", markersize = 0.05,label=f"n_step ={num_time_steps}")

n = np.array([4000,8000,16000,32000])

plt.figure(2)
plt.subplot(2,1,1) 
plt.title("Error over time for Forward Euler")
for i in n:
    error_plot(i,"Forward_Euler")
plt.xlabel("t [$\\mu$s]")
plt.ylabel("y [$\\mu$m]")
plt.yscale('log')
plt.legend()

plt.subplot(2,1,2)
plt.title("Error over time for Runge Kutta 4")
for i in n:
    error_plot(i,"RK4")
plt.xlabel("t [$\\mu$s]")
plt.ylabel("y [$\\mu$m]")
plt.yscale('log')
plt.legend()

plt.savefig("error_plot.pdf")