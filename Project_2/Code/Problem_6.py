import matplotlib.pyplot as plt
import numpy as np

def plot_stuff(name,N):
    '''Read three eigenvectors (numerical and analytical) from {name}.txt,
    and plot them against corresponding x-values. N is the dimension
    of the matrix, such that there are N+2 points in each eigenvectors
    when we include the boundary points (which are zero).'''
    data = np.loadtxt(f"{name}.txt")
    xaxis = np.linspace(0, 1, N+2)

    # declaring new array that includes boundary conditions:
    data_new = np.zeros((data.shape[0]+2, data.shape[1]))
    for i in range(data_new.shape[1]):
        data_new[1:-1,i] = data[:,i] # add the vectors to the new array

    # plotting the eigenvectors:
    plt.figure()
    plt.title(f"Eigenvectors v corresponding to the 3 lowest eigenvalues for {N+1} steps")
    for i in range(0,3):
        plt.plot(xaxis, data_new[:,i+3], alpha=.4, label=f"$v_{i+1}$, analytic")
        plt.plot(xaxis, data_new[:,i], '--', label=f"$v_{i+1}$, numerical")
    plt.xlabel("$x$")
    plt.ylabel("$v$")
    plt.legend()
    plt.savefig(f"{name}_plot.pdf")
    plt.show()

plot_stuff("Problem_6_a",9)
plot_stuff("Problem_6_b",99)
