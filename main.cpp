/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem.
IEEE Transactions on Cybernetics, 52(1), 495-507.*/


#include "tspProblem.h"
#include "tspSolution.h"
#include "solver.h"


using namespace std;

int main(int argc, char *argv[]) {

	int iniSeed = 1;
	int seed = 1;
	int globalOptFit = 15281;
	double maxFunEval = 1e6;
	int nearestNum = 20;
	int printFunEvalIntvl = 1000;

	ILS_3Opt("rd400.tsp", iniSeed, seed, globalOptFit, maxFunEval, nearestNum, "ILS_rd400_test", printFunEvalIntvl);

	double lambda = 0.06;
	LSILS_3Opt_ConstantLambda("rd400.tsp", iniSeed, seed, globalOptFit, maxFunEval, nearestNum, lambda, "LSILS-C_rd400_test", printFunEvalIntvl);
	
	double maxLambda = 0.09;
	LSILS_3Opt_DynamicLambda("rd400.tsp", iniSeed, seed, globalOptFit, maxFunEval, nearestNum, maxLambda, "LSILS-D_rd400_test", printFunEvalIntvl);
	return 0;
}












