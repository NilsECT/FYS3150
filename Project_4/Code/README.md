# Project 4
Welcome to our world of code. In this code we run Monte Carlo integration to study properties of the discrete 2D Ising model. We do this by running Monte Carlo cycles with random walks for our important sampling. We run various investigations for exploring the system. 

For compiling:
```
g++ main.cpp src/* -I include/  -o main.exe -larmadillo -O2 -fopenmp
```
(The keys are for a run in ubuntu. For mac users, use what you need to compile with armadillo and open mp.)

For running:
```
time ./main.exe
```
("time" is optional)

All results are written to txt files so that they could be used for plotting in python.

For running 


For running:
```
./main <L> <number of threads> <number of MCMC cycles>
```
Resultatene skrives til en tekstfil

Viktig: Dette er bare "grunnstrukturen" til hvordan vi kan lage resultatene. Det mangler å faktisk lage funksjonene/idk som svarer på oppgavene.

### To do:
* Endre løkka slik at vi ikke regner ut exp() hver bidige runde
* Brage snakket om å bytte ut double (jeg "backer")
* Mer
