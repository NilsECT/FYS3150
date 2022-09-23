# Code for project 2

Hola and welcome to our project code. This is an instruction on how to run our code to reproduce the results in our report.

Running c++ programs that uses armadillo might need the "-larmadillo" included in the compiling of the code.

## Problem 2:
The program Problem_2.cpp is checking wether two matrices have the same eigen vectors.
Running the executable Problem_2.exe cout's a statment telling us if the two matrices of eigenvectors are the same.  

How to run Problem_2:
#### 1. g++ Problem_2.cpp -o Problem_2.exe -larmadillo
#### 2. ./Problem_2.exe


## Problem 3:
The constructed class Jacobi.cpp contains a test function for the matrix in Problem 7 b. It raises an assertion error if it can not access the k and l element of the matrix.

How to run Problem_3:
#### 1. g++ Problem_3.cpp Jacobi.cpp -I ./ -o Problem_3.exe -larmadillo
#### 2. ./Problem_3.exe

## Problem 4:
The code in Problem 4 generates the matrix from the problem text, puts it into the Jacobi class and gives us a matrix of eigenvectors for the matrix. The matrix is then arranged in order of raising eigenvalue and compared with armadillos method for finding eigenvectors. It returns a statement telling us if the matrix of eigenvectors from the jacobi rotation method is consistent with Armadillos matrix.

How to run Problem_4:
#### 1. g++ Problem_4.cpp Jacobi.cpp -I ./ -o Problem_4.exe -larmadillo
#### 2. ./Problem_4.exe

## Problem 5:
The code for problem 5 keeps track of how many rotations the jacobi rotation method has to do to solve the eigen problem. It runs over the same matrix as problem 4 but for N ranging from 2 to 99. it returns the number of transformations for each N to a txt file. Problem_5.py takes the txt file and plots the number of transformations for each N. It also makes a linear fit to the plot helping us estimate the number of transformations needed for an N times N matrix.

How to run Problem_5:
#### 1. g++ Problem_5.cpp Jacobi.cpp -I ./ -o Problem_5.exe -larmadillo
#### 2. ./Problem_5.exe
#### 3. python3 Problem_5.py
#### 4. xdg-open Problem_5_plot.pdf

## Problem 6:
For problem 6 a the code uses the jacobi class on an 9 times 9 matrix, and returns the eigenvectors of the 3 lowest eigenvalues to a txt file. The python script then takes the eigenvectors and adds the boundary conditions and plots the eigenvectors with the analytical eigenvectors from Armadillos eig_sym method.
The code for b does the same only for a 99 times 99 matrix.

How to run Problem_6:
#### 1. g++ Problem_6.cpp Jacobi.cpp -I ./ -o Problem_6.exe -larmadillo
#### 2. ./Problem_6.exe
#### 1. g++ Problem_6_b.cpp Jacobi.cpp -I ./ -o Problem_6_b.exe -larmadillo
#### 2. ./Problem_6_b.exe
#### 3. python3 Problem_6.py
#### 4. xdg-open Problem_6_a_plot.pdf
#### 4. xdg-open Problem_6_b_plot.pdf


fin.