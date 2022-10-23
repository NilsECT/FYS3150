# Mappe for project 3 code:

Hi and welcome to our code for project 3. It is centered around object oriented programming for simulating a penning trap containing Calcium ions. The simulation is done with a few different parameters that has all been specified in the simulation main functions. Of these there are running a simulation for 1,2 and 25 particles, with and without particle-particle interactions with coulomb force. These simulations are also done with two different numerical integration methods, Forward Euler and Runge Kutta to the 4'th order. There is also an analytical solution for 1 particle in the penning trap that can be compeard with to check our results. For the las part we wil be alternating the elecrtical field of the penningtrep with different frequencies and angular frequencies. This will be done for 25 particles so the simulation time not to be too long. 

## Simulation Generator.
For the first part we run the Simulation_generator.cpp to generate all the different ways we wand to simulate the system. This includes; 1 particle analyticaly, with Runge Kutta 4 and with Forward Euler, 2 particles with and without particle interaction with both Forward Euler and Runge Kutta 4. 

### To run:
To run this we must include the Particle and Penningtrap class where we have all the functionality for running these simulation. What we are doing is that for each time step we wright out the position coordinate and the velocities for each particle in the simulation. There is one txt file for each of the coordinates with each row in the file being one particles trajectory in that coordinat form for each time step. We wright this as:

####    -   g++ Simulation_generator.cpp src/Particle.cpp src/Penningtrap.cpp -I include/ -o Simulation_generator.exe -larmadillo && ./Simulation_generator.exe

(Here there have been difficulties with armadillo not working for some people running the code, so we include the -larmadillo for safety.)(The prints in the terminal is to let you know how far the code have gotten in the simulations. Not in time steps but what kind of simulation it is running.)

## Simulation perturbation.
This simulation is a little longer so we include a optimization key -O2 this will help the runtime and still keep our calculations ok. Here there is still a cout to let you know how far the simulation has come.

####    -   g++ Simulation_perturbation.cpp src/Particle.cpp src/Penningtrap.cpp -I include/ -o Simulation_perturbation.exe -O2 -larmadillo && ./Simulation_perturbation.exe

## Plotting.
There are two plot codes for these for simulations. plot_one_particle.py is meant for interpreting the one particle result of the simulation. It generates one plot for one particle in the z-direction for the  analytical, Forward Euler and Runge Kutta 4 simulation. It also generates a plot for the relative error in the simulation. This is done in two subplots for both Runge Kutta 4 and Forward Euler, and compering them to the analytical for the different number of time steps. This is done For 4000, 8000, 16000, 32000 number of time steps in a 50 micro second simulation. It also uses the max of these error to return the values the convergence rate for the error for both the Runge Kutta and the Forward Euler.
####    -   python3 plot_one_particle.py 

Now opening these plots straight forward:
For one particle in z-direction:
#### - xgd-open 1_particle_t_direction
For the relative error:
#### - xdg-open error_plot.pdf

For plot.py there is a lot of plots that are meant for checking that the results from the simulation make sense from visualization. The run is straight forward:
#### - python3 plot.py

Opening the plots:


