import matplotlib.pyplot as plt
import numpy as np

def plot_stuf(name,N):
    data = np.loadtxt(f"{name}.txt")
    data_shape = data.shape
    data_new = np.zeros((data_shape[0]+2, data_shape[1])) #declearing new array with boundry conditions
    for i in range(data_new.shape[1]):
        data_new[1:-1,i] = data[:,i] #adding the vectors to the new arrays
    
    #Plotting:
    plt.figure(N)
    plt.title(f"3 lowest eigen values for {N+1} steps:")
    for i in range(0,3):
        plt.plot(data_new[:,i], label=f"lambda_{i+1}")
        plt.plot(data_new[:,i+3], label=f"anal_{i+1}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.savefig(f"{name}_plot.pdf")

plot_stuf("Problem_6_a",9)
plot_stuf("Problem_6_b",99)


