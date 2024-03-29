# Landscape_Smoothing_ILS_TSP
# Landscape Smoothing Iterated Local Search (LSILS) for solving the TSP
# Jialong Shi (jialong.shi@xjtu.edu.cn)



## Introduction

In this program, Landscape Smoothing Iterated Local Search (LSILS) is implemented to optimize the Traveling Salesman Problem (TSP). 

The code is distributed for research use. The author reserves all rights to the code.

Relevant literature:

[1] J. Shi, J. Sun, Q. Zhang and K. Ye, Homotopic Convex Transformation: A New Landscape Smoothing Method for the Traveling Salesman Problem, in IEEE Transactions on Cybernetics, vol. 52, no. 1, pp. 495-507, Jan. 2022, doi: 10.1109/TCYB.2020.2981385.

## File list (the following files should be in the same directory)

- Source code: main.cpp  tspProblem.h  tspProblem.cpp  tspSolution.h  tspSolution.cpp  solver.h  solver.cpp

- Problem file examples: rd400.tsp

- cmake file: CMakeLists.txt

## How to compile it

CMakeLists.txt has been provided for compiling the program by cmake on linux (such as Ubuntu), so make sure that cmake has been installed. Then you can compile the project by running

```
cd <directory-of-this-program>
cmake CMakeLists.txt
make
```

Then an executable named 'lsils_tsp_3opt' will appear in the directory. 

!!! If the cmake command or the make command is failed, you can directly compile the program by the following command. 

```
g++  -o  lsils_tsp_3opt  -O2  main.cpp  tspProblem.cpp  tspSolution.cpp  solver.cpp
```

## How to use it

Before running, please create a subfolder named 'results'.

```
mkdir results
```

You can run the program by

```
./lsils_tsp_3opt
```

The above command runs ILS, LSILS with constant lambda 0.06, LSILS with dynamic lambda (0 to 0.09) on the TSPLIB instance rd400. The function evaluation budget is 1000000.



## Outputs

The outputs are located in the 'results' folder. The 'log_....txt' files contain crucial information regarding the execution process, while the 'fit_eval_....txt' files record the fitness (first column) versus the function evaluation (second column).

```



