/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem.
IEEE Transactions on Cybernetics, 52(1), 495-507.*/

#include "tspSolution.h"

using namespace std;

//clear the solution x.
void tspSolution::clear() {
	x.clear();
}

//generate a random solution
void tspSolution::random(const tspProblem& inst) {
	clear();
	int n = inst.getN();
	x.resize(n);
	for (int i = 0; i < n; i++) {
		x[i] = i;
	}
	int temp;
	int randPos;
	for (int i = 0; i < n; ++i) {
		randPos = rand() % (n - i) + i;
		temp = x[i];
		x[i] = x[randPos];
		x[randPos] = temp;
	}
	calFitness(inst);
}

//calculate the fitness
bool tspSolution::calFitness(const tspProblem& inst) {
	int n = inst.getN();
	if (x.size() != n) {
		return false;
	}

	int P, Q;
	fitness = 0;
	for (int i = 0; i < n; i++) {
		P = x[i];
		Q = x[(i == (n - 1)) ? 0 : (i + 1)];
		fitness += inst.getDist(P, Q);
	}
	return true;
}

//assign a tspSolution object to another
tspSolution& tspSolution::operator = (const tspSolution &another) {
	int i;
	int n = another.x.size();

	x.clear(); 
	x.resize(n);
	for (i = 0; i < n; i++) {
		x[i] = another.x[i];
	}

	fitness = another.fitness;


	return *this;
}
