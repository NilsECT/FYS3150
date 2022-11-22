# Project 4
Welcome to our world of code. In this code we run Monte Carlo simulation to study properties of the discrete 2D Ising model. We do this by running Monte Carlo cycles with random walks for our important sampling. We run various investigations for exploring the system.

For compiling:
```
g++ main.cpp src/* -I include/ -o main.exe -larmadillo -O3 -fopenmp
```
(The keys are for a run in ubuntu. For mac users, use what you need to compile with armadillo and openmp.)

For running: (This is where we run the simulation, and generate all our data. You
should expect this part to take very long.)
```
time ./main.exe
```
("time" is optional)

All results are written to txt files so that they could be used for plotting in python.

For plotting in python

```
python3 plot.py
```