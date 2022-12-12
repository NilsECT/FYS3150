# Project 5
Welcome to the fifth and final project in this repo. In this project we explore the world of partial differential equations(PEDs). More specifically we explore the double slit experiment numerically with the Crank-Nicolson(CN) scheme for evolving/solving PDEs. We conduct simulations of various versions of the experiment and produce plots and an animation to visually confirm that the method works. We also investigate the stability for the CN method by plotting the total probability over time. 

## C++
main.cpp generates all the data that is needed for the python programs. It run 5 different experiments for different conditions of the experiment such as wall, one slit, double slit, triple slit and various widths of the wave in the simulation. The code saves all the experiments data as binary files in form the Armadillo cube type.

To compile and run:
```
g++ main.cpp src/* -I include/ -o main.exe -larmadillo && ./main.exe
```
(The compiling includes -larmadillo to make use of the Armadillo modula in c++. The code takes about 20 min to execute.)

## Python
plot.py produces all visual data for the report and essentially all of results. It converts Armadillo cubes to numpy arrays and generates various plots. 

To run:
```
python3 plot.py
``` 

animate.py produces an animation of the experiment for visually checking the reliability of the method and the results it produces.

To run:
```
python3 animate.py
```

Thank you for following our projects, we hope you have enjoyed our code and results. :)