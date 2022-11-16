# Project 4

For compiling:
```
g++ main.cpp src/Grid.cpp -I include/  -o main -larmadillo -O2
```
Ta med flagget du trenger for å få OMP

'''
For Sigurd:
Comp:   g++ main.cpp src/Grid.cpp src/Ising.cpp -I include/  -o main.exe -larmadillo -O2 -fopenmp
exe:    ./main.exe 10 10 4000000
plot:   python3 plotting.py
vis:    xdg-open Histogram_of_epsilon_for_T_1.pdf && xdg-open Histogram_of_epsilon_for_T_2.4.pdf
'''


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
