import matplotlib.pyplot as plt
import numpy as np

def plot_stuf(name,N):
    lambda_1, lambda_2, lambda_3, anal_1, anal_2 ,anal_3 = np.loadtxt(f"{name}.txt", unpack=True)

    np.insert(lambda_1,0,0)
    np.insert(lambda_2,0,0)
    np.insert(lambda_3,0,0)
    np.insert(anal_1,0,0)
    np.insert(anal_1,0,0)
    np.insert(anal_1,0,0)

    np.insert(lambda_1,-1,0)
    np.insert(lambda_2,-1,0)
    np.insert(lambda_3,-1,0)
    np.insert(anal_1,-1,0)
    np.insert(anal_1,-1,0)
    np.insert(anal_1,-1,0)
    plt.figure(N)
    plt.title(f"3 lowest eigen values for {N+1} steps:")
    plt.plot(lambda_1, label="lambda_1")
    plt.plot(anal_1, label="anal_1")
    plt.plot(lambda_2, label="lambda_2")
    plt.plot(anal_2, label="anal_2")
    plt.plot(lambda_3, label="lambda_3")
    plt.plot(anal_3, label="anal_3")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.savefig(f"{name}_plot.pdf")

plot_stuf("Problem_6_a",9)
plot_stuf("Problem_6_b",99)


