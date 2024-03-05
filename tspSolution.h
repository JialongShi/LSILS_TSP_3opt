/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem.
IEEE Transactions on Cybernetics, 52(1), 495-507.*/

#pragma once

#include<vector>
#include <stdlib.h>
#include"tspProblem.h"

class tspSolution
{
public:
	std::vector<int> x;
	int fitness;

public:
	void clear();
	void random(const tspProblem& instance);
	bool calFitness(const tspProblem& inst);
	tspSolution& operator = (const tspSolution &another);
};

