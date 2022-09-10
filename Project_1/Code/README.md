# Code

Hello and welcome to our getting to learn c++ through a series of interesting challenges closely related to linear algebra. This is our first encounter to c++ so the code might be messy now but in the next projects we hope our flaws are rectified.

## The compiling and running of the program is here the ones used for the terminal in Ubuntu 20.04.1, for MacOS you might have to add "-larmadillo" as such: g++ -larmadillo foo.cpp -o foo.exe

A few notes, the packages needed to run these programs are listed in the course-page of FYS3150/4150. The page contains a link to the github repository of the course which has so much nice info!
The commands starting with "xdg-open" are to see the plots created.
There will be few .txt files created through the running of these programs so that you're aware of that.

So lets get down to business! To run these programs!

## Problem 2

The executable 'Problem_2.exe' will produce a .txt file named 'Problem_2.txt' which contains the data used for plotting by 'Problem_2.py'. The name can be changed in the end of the respective .cpp. Here are the steps for the command line to run the programs and scrips related to Problem 2:

#### 1. g++ Problem_2.cpp -o Problem_2.exe
#### 2. ./Problem_2.exe
#### 3. python3 Problem_2.py
#### 4. xdg-open Problem_2_plot.svg

## Problem 7

Here we will create quite a few files, some are needed for this problem and others will be for Problem 8. All the files are created automagically by the Problem_7.exe program and are called "out_e*.txt" where the star is a number representing the order of the number of steps. Anyhow, let§s get to it:

#### 1. g++ Problem_7.cpp -o Problem_7.exe
#### 2. ./Problem_7.exe
#### 3. python3 Problem_7.py
#### 4. xdg-open Problem_7_plot.svg

## Problem 8

Here it's just to see the plots from the data created in the previous problem.

#### 1. python3 Problem_8.py
#### 2. xdg-open log10_abs_err.svg
#### 3. xdg-open log10_rel_err.svg
#### 4. python3 Problem_8_c.py
#### 5. xdg-open Problem_8_c_plot.svg

## Problem 9

Here it is just if you need to, the code is written again in Problem_10_master.cpp. The program we will create here prints out the result to the terminal so if you want to add it to a file just add "> name.txt" in the command line when you run the program (you choose the name).

#### 1. g++ Problem_9.cpp -o Problem_9.exe
#### 2. ./Problem_9.exe

## Problem 10

Here we will create the tables wanted through a script that takes both problem 7 and 9. The output is printed to the terminal and commented in the .pdf we handed in.

#### 1. g++ Problem_10_master.cpp -o Problem_10.exe
#### 2. g++ Problem_10.exe
