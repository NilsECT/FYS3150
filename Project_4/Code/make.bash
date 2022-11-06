g++omp main.cpp src/Grid.cpp -I include/ -o main -larmadillo


g++dillo -Xpreprocessor -fopenmp main.cpp src/Grid.cpp -I include/ -o main -larmadillo

CPATH=/opt/homebrew/Cellar/libomp/15.0.4/include LIBRARY_PATH=/opt/homebrew/Cellar/libomp/15.0.4/lib g++ -std=c++11

g++dillo -std=c++11 -I /opt/homebrew/opt/libomp/include -L /opt/homebrew/opt/libomp/lib
