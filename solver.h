/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem. 
IEEE Transactions on Cybernetics, 52(1), 495-507.*/

#pragma once

#include <vector>
#include <climits>
#include <list>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include "tspProblem.h"
#include "tspSolution.h"

class flipper_class
{
public: //protected:
	int reversed;
	int cycle_size;
	int short_size;
	int *cyc;
	int *cyc_inv;

public:
	flipper_class();
	~flipper_class();

	int		flipper_init(const std::vector<int> &incyc, int n);
	void	flipper_cycle(std::vector<int> &x);
	int		flipper_next(int x);
	int		flipper_prev(int x);
	void	flipper_flip(int x, int y);
	void	flipper_reverse();
	int		flipper_sequence(int x, int y, int z);
};

bool move2Opt(const int& P, const int& Q, tspSolution& sol);
bool doubleBridgeKickBits(tspSolution &sol, const tspProblem &inst, std::vector<bool> &bits);
bool genNearestList(std::vector< std::vector<int> > &nearestList, int nearestNum, const tspProblem &inst);
bool homotopicConvexTransfer(tspProblem &pseudoTSP, const tspProblem originalTSP, const tspSolution &localOpt, const double lambda);

int ILS_3Opt(const char* tspName, const int iniSolSeed, const int seed, const int globalOptFit, const double maxFunEval,
	const int nearestNum, const char* outputLabel, const int printFunEvalIntvl);

int LSILS_3Opt_ConstantLambda(const char* tspName, const int iniSolSeed, const int seed, const int globalOptFit, const double maxFunEval,
	const int nearestNum, const double lambda, const char* outputLabel, const int printFunEvalIntvl);

int LSILS_3Opt_DynamicLambda(const char* tspName, const int iniSolSeed, const int seed, const int globalOptFit, const double maxFunEval,
	const int nearestNum, const double maxLambda, const char* outputLabel, const int printFunEvalIntvl);