import matplotlib.pyplot as plt
import numpy as np

def plot_stuf(name,N):
    data = np.loadtxt(f"{name}.txt")
    for i in data.T:
        np.insert(i, 0, 0)
        np.insert(i, -1, 0)
    plt.figure(N)
    plt.title(f"3 lowest eigen values for {N+1} steps:")
    for i in range(0,3):
        plt.plot(data.T[i], label=f"lambda_{i+1}")
        plt.plot(data.T[i+3], label=f"anal_{i+1}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.savefig(f"{name}_plot.pdf")

plot_stuf("Problem_6_a",9)
plot_stuf("Problem_6_b",99)


