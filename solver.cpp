/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem.
IEEE Transactions on Cybernetics, 52(1), 495-507.*/

#include "solver.h"

#if defined(__sun) || defined(linux)
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#else
#include <windows.h>
#endif


using namespace std;

#if defined(__sun) || defined(linux)
static struct tms  Old_Time;
static struct tms  New_Time;
static int ClkTck;
static clock_t clock_t_start;
static clock_t clock_t_end;
void start_timer()
{
	ClkTck = sysconf(_SC_CLK_TCK);
	clock_t_start = times(&Old_Time);
}

double return_CPU_time()
{
	double cpu_time;
	clock_t_end = times(&New_Time);
	cpu_time = (New_Time.tms_utime - Old_Time.tms_utime) / ((double)ClkTck)
		+ (New_Time.tms_stime - Old_Time.tms_stime) / ((double)ClkTck)
		+ (New_Time.tms_cutime - Old_Time.tms_cutime) / ((double)ClkTck)
		+ (New_Time.tms_cstime - Old_Time.tms_cstime) / ((double)ClkTck);
	return cpu_time;
}
double return_elapsed_time()
{
	double real_time;
	clock_t_end = times(&New_Time);
	real_time = (clock_t_end - clock_t_start) / ((double)ClkTck);
	return real_time;
}
#else
static time_t Old_Time;
static time_t New_Time;
void start_timer()
{
	//time_t Old_Time = time(NULL);//
	time(&Old_Time);
}

double return_elapsed_time()
{
	time(&New_Time);
	//return (double)(New_Time.tms_utime - Old_Time.tms_utime)/sysconf(_SC_CLK_TCK);
	return (double)difftime(New_Time, Old_Time);
}
#endif

flipper_class::flipper_class()
{
}

//clean up
flipper_class::~flipper_class()
{
	if (cyc) delete[] cyc;
	if (cyc_inv) delete[] cyc_inv;
	cycle_size = 0;
	short_size = 0;
	reversed = 0;
}

//initialize
int flipper_class::flipper_init(const std::vector<int> &incyc, int n)
{
	int i;

	cyc = new int[n];
	cyc_inv = new int[n];;

	for (i = 0; i < n; i++) {
		cyc[i] = incyc[i];
		cyc_inv[incyc[i]] = i;
	}

	cycle_size = n;
	short_size = n / 2;
	reversed = 0;
	return 1;
}

//reverse the cycle based on the 'reversed' flag
void flipper_class::flipper_cycle(std::vector<int> &x)
{
	int *p;

	if (reversed) {
		p = cyc + cycle_size;
		int index = 0;
		while (p > cyc) {
			x[index] = *--p;
			index++;
		
		}
	}
	else {
		p = cyc + cycle_size;
		int index = cycle_size - 1;
		while (p > cyc) {
			x[index] = *--p;
			index--;
		}
	}
}

//find the next element of x in the cycle 
int flipper_class::flipper_next(int x)
{
	int y;

	if (reversed) {
		y = cyc_inv[x] - 1;
		return (y >= 0) ? cyc[y] : cyc[cycle_size - 1];
	}
	else {
		y = cyc_inv[x] + 1;
		return (y < cycle_size) ? cyc[y] : cyc[0];
	}
}

//find the previous element of x in the cycle 
int flipper_class::flipper_prev(int x)
{
	int y;

	if (reversed) {
		y = cyc_inv[x] + 1;
		return (y < cycle_size) ? cyc[y] : cyc[0];
	}
	else {
		y = cyc_inv[x] - 1;
		return (y >= 0) ? cyc[y] : cyc[cycle_size - 1];
	}
}

//perform a flip between x and y
void flipper_class::flipper_flip(int x, int y)
{
	int xloc = cyc_inv[x];
	int yloc = cyc_inv[y];
	int temp;
	int gap;

	if (reversed) {
		CC_SWAP(xloc, yloc, temp);
	}
	gap = yloc - xloc;
	if (gap < 0)
		gap += cycle_size;
	if (gap > short_size) {
		CC_SWAP(xloc, yloc, temp);
		reversed ^= 1;
		xloc++;
		if (xloc >= cycle_size)
			xloc = 0;
		yloc--;
		if (yloc < 0)
			yloc = cycle_size - 1;
		gap = cycle_size - gap - 2;
	}

	if (xloc > yloc) {
		gap++;
		gap /= 2;
		for (; gap; gap--) {
			x = cyc[xloc];
			y = cyc[yloc];
			cyc[xloc] = y;
			cyc[yloc] = x;
			cyc_inv[x] = yloc--;
			cyc_inv[y] = xloc++;
			if (xloc >= cycle_size)
				xloc = 0;
			if (yloc < 0)
				yloc = cycle_size - 1;
		}
	}
	else {
		gap++;
		gap /= 2;
		for (; gap; gap--) {
			x = cyc[xloc];
			y = cyc[yloc];
			cyc[xloc] = y;
			cyc[yloc] = x;
			cyc_inv[x] = yloc--;
			cyc_inv[y] = xloc++;
		}
	}
}

//change the direction of the cycle
void flipper_class::flipper_reverse()
{
	reversed ^= 1;
}

//determine whether the x, y, and z occur in sequence
int flipper_class::flipper_sequence(int x, int y, int z)
{
	int xloc = cyc_inv[x];
	int yloc = cyc_inv[y];
	int zloc = cyc_inv[z];

	if (reversed) {
		if (xloc >= yloc)
			return yloc >= zloc || zloc >= xloc;
		else
			return yloc >= zloc && zloc >= xloc;
	}
	else {
		if (xloc <= yloc)
			return yloc <= zloc || zloc <= xloc;
		else
			return yloc <= zloc && zloc <= xloc;
	}
}

//2-Opt move
bool move2Opt(const int& P, const int& Q, tspSolution& sol) {
	int n = sol.x.size();
	int cityP = (P <= Q) ? P : Q;
	int cityQ = (P >= Q) ? P : Q;

	int right = (cityP == (n - 1)) ? 0 : (cityP + 1);
	int left = cityQ;
	int orderR, orderL;

	while (right < left) {
		if (right == n)
			right = 0;
		if (left == -1)
			left = n - 1;
		orderR = sol.x[right];
		orderL = sol.x[left];
		sol.x[left--] = orderR;
		sol.x[right++] = orderL;
	}
	return true;
}

//perturbation
bool doubleBridgeKickBits(tspSolution &sol, const tspProblem &inst, vector<bool> &bits) {
	 int n = sol.x.size();
	 if (n < 8) {
		 return false;
	 }
	 int t1, t2, t3, t4, t5, t6, t7, t8;
	 vector<int> sortSpace(4);

	 t1 = rand() % n;
	 t2 = t1 + 1;
	 if (t2 >= n) t2 = 0;

	 do {
		 t3 = rand() % n;
		 t4 = t3 + 1;
		 if (t4 >= n) t4 = 0;
	 } while (t3 == t1 || t3 == t2 || t4 == t1);

	 do {
		 t5 = rand() % n;
		 t6 = t5 + 1;
		 if (t6 >= n) t6 = 0;
	 } while (t5 == t1 || t5 == t2 || t5 == t3 || t5 == t4 || t6 == t1 || t6 == t3);


	 do {
		 t7 = rand() % n;
		 t8 = t7 + 1;
		 if (t8 >= n) t8 = 0;
	 } while (t7 == t1 || t7 == t2 || t7 == t3 || t7 == t4 || t7 == t5 || t7 == t6 || t8 == t1 || t8 == t3 || t8 == t5);

	 sortSpace[0] = t1;
	 sortSpace[1] = t3;
	 sortSpace[2] = t5;
	 sortSpace[3] = t7;

	 sort(sortSpace.begin(), sortSpace.end());
	 t1 = sortSpace[0];
	 t2 = t1 + 1;
	 t3 = sortSpace[1];
	 t4 = t3 + 1;
	 t5 = sortSpace[2];
	 t6 = t5 + 1;
	 t7 = sortSpace[3];
	 t8 = t7 + 1;
	 if (t8 >= n) t8 = 0;

	 sol.fitness +=
		 inst.getDist(sol.x[t1], sol.x[t6]) + inst.getDist(sol.x[t2], sol.x[t5]) +
		 inst.getDist(sol.x[t3], sol.x[t8]) + inst.getDist(sol.x[t4], sol.x[t7]) -
		 inst.getDist(sol.x[t1], sol.x[t2]) - inst.getDist(sol.x[t3], sol.x[t4]) -
		 inst.getDist(sol.x[t5], sol.x[t6]) - inst.getDist(sol.x[t7], sol.x[t8]);

	 bits[sol.x[t1]] = bits[sol.x[t2]] = bits[sol.x[t3]] = bits[sol.x[t4]] = bits[sol.x[t5]] = bits[sol.x[t6]] = bits[sol.x[t7]] = bits[sol.x[t8]] = true;

	 move2Opt(t1, t3, sol);
	 move2Opt(t3, t5, sol);
	 move2Opt(t5, t7, sol);
	 move2Opt(t1, t7, sol);

	 return true;
 }

//update a sorted list of nearest cities
void intsertNearlist(int m, int thisDist, vector<int> &nearList, vector<int> &nearDists) {
	 int i;
	 if (thisDist < nearDists[0])
	 {
		 for (i = 0; nearDists[i + 1] > thisDist; i++)
		 {
			 nearList[i] = nearList[i + 1];
			 nearDists[i] = nearDists[i + 1];
		 }
		 nearDists[i] = thisDist;
		 nearList[i] = m;
	 }
 }

//generate a list of nearest cities for each city 
bool genNearestList(vector< vector<int> > &nearestList, int nearestNum, const tspProblem &inst) {
	int n = inst.getN();
	nearestList.resize(n);
	if (nearestNum > n - 1) {
		nearestNum = n - 1;
	}
	for (int i = 0; i < n; i++)
		nearestList[i].resize(nearestNum);

	vector<int> nearDists(nearestNum + 1);

	for (int nn = 0; nn < n; nn++) {
		for (int i = 0; i < nearestNum; i++)
			nearDists[i] = INT_MAX;
		nearDists[nearestNum] = -INT_MAX;

		for (int m = nn - 1; m >= 0; m--) {
			int thisDist = inst.getDist(nn, m);
			intsertNearlist(m, thisDist, nearestList[nn], nearDists);
		}
		for (int m = nn + 1; m < n; m++) {
			int thisDist = inst.getDist(nn, m);
			intsertNearlist(m, thisDist, nearestList[nn], nearDists);
		}
	}

	return true;
}

//generate a pseudoTSP by combining the originalTSP and the localCircleTSP based on localOpt. 
bool homotopicConvexTransfer(tspProblem &pseudoTSP, const tspProblem originalTSP, const tspSolution &localOpt, const double lambda) {
	if (originalTSP.n != localOpt.x.size() || lambda > 1) {
		return false;
	}

	int n = originalTSP.n;
	int i, j;

	//generate the circle TSP based on the local optimum
	tspProblem localCircleTSP;
	localCircleTSP.edgeWeightType = CC_EUCLIDEAN;
	localCircleTSP.X.resize(n);
	localCircleTSP.Y.resize(n);
	localCircleTSP.dist.resize(n);
	for (i = 0; i < n; i++) {
		localCircleTSP.dist[i].resize(n);
	}

	double meanDist = 0;
	for (i = 0; i < n; i++){
		int minDistTemp = INT_MAX;
		for (j = 0; j < n; j++){
			if (originalTSP.getDist(i, j) < minDistTemp && j != i) {
				minDistTemp = originalTSP.getDist(i, j);
			}
		}
		meanDist += ((double)minDistTemp / n);
	}
	double r = meanDist / 2.0 / sin(PI / n);
	for (i = 0; i < n; i++) {
		double theta = (i - 1) * 2.0 * PI / n;
		localCircleTSP.X[localOpt.x[i]] = r*cos(theta) + r;
		localCircleTSP.Y[localOpt.x[i]] = r*sin(theta) + r;
	}
	for (i = 0; i < n; i++) {
		localCircleTSP.dist[i][i] = 0;
		for (j = 0; j < i; j++) {
			int distTemp = (int)floor(sqrt((localCircleTSP.X[i] - localCircleTSP.X[j])*(localCircleTSP.X[i] - localCircleTSP.X[j]) 
				+ (localCircleTSP.Y[i] - localCircleTSP.Y[j])*(localCircleTSP.Y[i] - localCircleTSP.Y[j])) + 0.5);
			localCircleTSP.dist[i][j] = distTemp;
			localCircleTSP.dist[j][i] = distTemp;
		}
	}

	//generate the pseudo TSP
	pseudoTSP.edgeWeightType = CC_MATRIXNORM;
	pseudoTSP.X.clear();
	pseudoTSP.Y.clear();
	pseudoTSP.dist.resize(n);
	for (i = 0; i < n; i++) {
		pseudoTSP.dist[i].resize(n);
	}
	for (i = 0; i < n; i++) {
		pseudoTSP.dist[i][i] = 0;
		for (j = 0; j < i; j++) {
			pseudoTSP.dist[i][j] = (int)floor((1 - lambda)*((double)originalTSP.dist[i][j]) + lambda*((double)localCircleTSP.dist[i][j]) + 0.5);
			pseudoTSP.dist[j][i] = pseudoTSP.dist[i][j];
		}
	}

	return true;
}

//Iterated Local Search with 3-Opt local search
int ILS_3Opt(const char* tspName, const int iniSolSeed, const int seed, const int globalOptFit, const double maxFunEval,
	const int nearestNum, const char* outputLabel, const int printFunEvalIntvl) {
	tspProblem tspInst;
	if (!tspInst.load(tspName)) {
		return 1;
	}

	int n = tspInst.getN();
	tspSolution sol, bestSol;
	double currentTime;
	vector<vector<int> > nearestList;
	int i, j, k;
	int A, B, C, D, E, F;
	int AB, CD, EF, AC, BF, DE, BD, CE, BC;//distances
	vector<int> A_neighbor;
	vector<int> D_neighbor;
	int delta;
	flipper_class fliptour;
	int flippercost;
	int search_index;
	bool stopFlag;
	double funEvalNum;
	double lsNum;
	vector<bool> bits;
	int improveType;


	char logFileName[200];
	sprintf(logFileName, "results/log_%s.txt", outputLabel);
	FILE* logFile = fopen(logFileName, "w");
	if (logFile == NULL) {
		return 1;
	}

	vector<int> funEvalPoint;
	funEvalPoint.resize(maxFunEval / printFunEvalIntvl - 1);
	funEvalPoint[0] = printFunEvalIntvl;
	for (i = 1; i < funEvalPoint.size(); i++) {
		funEvalPoint[i] = funEvalPoint[i - 1] + printFunEvalIntvl;
	}

	char fitFunEvalFileName[200]; //record how
	sprintf(fitFunEvalFileName, "results/fit_funeval_%s.txt", outputLabel);
	FILE* fitFunEvalFile = fopen(fitFunEvalFileName, "w");
	if (fitFunEvalFile == NULL) {
		return 1;
	}

	//initialization
	srand(iniSolSeed);
	sol.random(tspInst);
	fprintf(logFile, "%s\nMax fun eval=%.0f, seed=%d, Ini fitness=%d\n", tspName, maxFunEval, seed, sol.fitness);
	printf("ILS started on %s\nMax fun eval=%.0f, seed=%d, Ini fitness=%d\n", tspName, maxFunEval, seed, sol.fitness);
	srand(seed);
	genNearestList(nearestList, nearestNum, tspInst);
	bits.resize(n);
	for (i = 0; i < n; i++) bits[i] = true;

	fprintf(fitFunEvalFile, "%d 0\n", sol.fitness);

	bestSol = sol;
	stopFlag = false;
	start_timer();
	funEvalNum = 0;
	lsNum = 0;
	int outputFEIndex = 0;

	while (1) {
		//3-Opt local search
		fliptour.flipper_init(sol.x, n);
		flippercost = sol.fitness;

		vector<int> fixedStart(n);
		for (i = 0; i < n; i++) fixedStart[i] = i;
		while (1) {
			delta = 0;
			for (i = 0; i < n; i++) {
				A = fixedStart[i];
				if (bits[A]) {
					A_neighbor = nearestList[A];

					search_index = 0;
					while (search_index < 2) {
						if (search_index == 1)
							fliptour.flipper_reverse();
						B = fliptour.flipper_next(A);
						if (B != A_neighbor[nearestNum - 1]) {
							AB = tspInst.getDist(A, B);
							for (j = nearestNum - 1; j >= 0; j--) {//search in the neighbor list of A, from nearest to farthest
								C = A_neighbor[j];
								AC = tspInst.getDist(A, C);
								if (AC < AB) {
									D = fliptour.flipper_next(C);
									if (D != A) {

										CD = tspInst.getDist(C, D);
										BD = tspInst.getDist(B, D);

										//test 2-opt move
										delta = AC + BD - AB - CD;
										funEvalNum++;

										if (outputFEIndex < funEvalPoint.size()) {
											if (funEvalNum >= funEvalPoint[outputFEIndex]) {
												fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
												outputFEIndex++;
											}
										}

										if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
											stopFlag = true;//stopping criterion is met, stop the entire procedure
											break;
										}
										if (delta < 0) {
											fliptour.flipper_flip(D, A);
											improveType = 0;
											break;
										}

										//test 2.5 opt move
										E = fliptour.flipper_prev(C);
										if (E != B) {
											CE = tspInst.getDist(C, E);
											BC = tspInst.getDist(B, C);
											DE = tspInst.getDist(D, E);
											delta = AC + BC + DE - AB - CD - CE;
											funEvalNum++;

											if (outputFEIndex < funEvalPoint.size()) {
												if (funEvalNum >= funEvalPoint[outputFEIndex]) {
													fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
													outputFEIndex++;
												}
											}

											if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
												stopFlag = true;//stopping criterion is met, stop the entire procedure
												break;
											}
											if (delta < 0)
											{
												fliptour.flipper_flip(B, E);
												fliptour.flipper_flip(E, C);
												improveType = 1;
												break;
											}
										}

										//test 3 opt move
										D_neighbor = nearestList[D];
										for (k = nearestNum - 1; k >= 0; k--) {//search in the neighbor list of D, from nearest to farthest
											E = D_neighbor[k];
											if (E != C) {
												DE = tspInst.getDist(D, E);
												if (DE < AB + CD - AC) {
													if (fliptour.flipper_sequence(B, E, C))
														F = fliptour.flipper_next(E);
													else
														F = fliptour.flipper_prev(E);
													EF = tspInst.getDist(E, F);
													BF = tspInst.getDist(B, F);
													delta = AC + BF + DE - AB - CD - EF;
													funEvalNum++;

													if (outputFEIndex < funEvalPoint.size()) {
														if (funEvalNum >= funEvalPoint[outputFEIndex]) {
															fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
															outputFEIndex++;
														}
													}


													if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
														stopFlag = true;//stopping criterion is met, stop the entire procedure
														break;
													}
													if (delta < 0) {
														fliptour.flipper_flip(D, A);
														fliptour.flipper_flip(B, E);
														improveType = 2;
														break;
													}
												}
											}
											if (stopFlag == true) break;
										}
									}
								}
								if (stopFlag == true) break;
								if (delta < 0) break;
							}
						}
						if (stopFlag == true) break;
						if (delta < 0) break;
						search_index++;
					}
				}

				if (stopFlag == true) break;

				
				if (delta < 0) {
					flippercost += delta;
					if (flippercost < bestSol.fitness) {
						fliptour.flipper_cycle(bestSol.x);
						bestSol.fitness = flippercost;
					}

					if (improveType == 0) {//2opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = true;
					}
					else if (improveType == 1) {//2.5opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = bits[E] = true;
					}
					else if (improveType == 2) {//3opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = bits[E] = bits[F] = true;
					}
					else {
						fprintf(logFile, "ERROR! unclear improve type!!!\n");
					}

					break;
				}
				else {
					bits[A] = false;
				}
			}

			if (stopFlag == true) break;
			if (delta >= 0) break; //local optimum is found, stop current LS
		}
		fliptour.flipper_cycle(sol.x);
		sol.fitness = flippercost;

		lsNum++;

		if (stopFlag == true) break;

		for (i = 0; i < n; i++) {
			if (bits[i] == true) {
				fprintf(logFile, "ERROR! don't=look-bits is not all-zero after LS!!!\n");
				break;
			}
		}

		//perturbation 
		doubleBridgeKickBits(sol, tspInst, bits);		
	}

#if defined(__sun) || defined(linux)
	currentTime = return_CPU_time();
#else
	currentTime = return_elapsed_time();
#endif

	int the_fitness = sol.fitness;
	sol.calFitness(tspInst);
	int the_fitness2 = bestSol.fitness;
	bestSol.calFitness(tspInst);
	if (the_fitness != sol.fitness || the_fitness2 != bestSol.fitness) {
		fprintf(logFile, "ERROR! fitness miss match!!!\n");
	}

	fprintf(logFile, "Optimal fitness = %d, Best fitness=%d\n", globalOptFit, bestSol.fitness);
	fprintf(logFile, "Finish time=%.2f, LS num=%.0f, Fun eval num= %.0f\n\n", currentTime, lsNum, funEvalNum);
	fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);

	printf("ILS finished, best fitness=%d\n", bestSol.fitness);
	printf("Finish time=%.2fs, LS num=%.0f, Fun eval num= %.0f\n\n", currentTime, lsNum, funEvalNum);

	fclose(logFile);
	fclose(fitFunEvalFile);

	return 0;
}

//Landscape Smoothing Iterated Local Search with 3-Opt local search and a constant lambda value
int LSILS_3Opt_ConstantLambda(const char* tspName, const int iniSolSeed, const int seed, const int globalOptFit, const double maxFunEval,
	const int nearestNum, const double lambda, const char* outputLabel, const int printFunEvalIntvl) {
	tspProblem tspInst;
	if (!tspInst.load(tspName)) {
		return 1;
	}

	int n = tspInst.getN();
	tspProblem pseudoInst;
	tspSolution sol, bestSol;
	bool bestUpdated;
	double currentTime;
	vector<vector<int> > nearestList;
	int i, j, k;
	int A, B, C, D, E, F;
	int AB, CD, EF, AC, BF, DE, BD, CE, BC;//distances
	vector<int> A_neighbor;
	vector<int> D_neighbor;
	int delta, realDelta, solRealFit;
	flipper_class fliptour;
	int flippercost;
	int search_index;
	bool stopFlag;
	double funEvalNum;
	double lsNum;
	vector<bool> bits;
	int improveType;


	char logFileName[200];
	sprintf(logFileName, "results/log_%s.txt", outputLabel);
	FILE* logFile = fopen(logFileName, "w");
	if (logFile == NULL) {
		return 1;
	}

	vector<int> funEvalPoint;
	funEvalPoint.resize(maxFunEval / printFunEvalIntvl - 1);
	funEvalPoint[0] = printFunEvalIntvl;
	for (i = 1; i < funEvalPoint.size(); i++) {
		funEvalPoint[i] = funEvalPoint[i - 1] + printFunEvalIntvl;
	}

	char fitFunEvalFileName[200]; //record how
	sprintf(fitFunEvalFileName, "results/fit_funeval_%s.txt", outputLabel);
	FILE* fitFunEvalFile = fopen(fitFunEvalFileName, "w");
	if (fitFunEvalFile == NULL) {
		return 1;
	}

	//initialization
	srand(iniSolSeed);
	sol.random(tspInst);
	fprintf(logFile, "%s lambda=%.4f\nMax fun eval=%.0fs, seed=%d, Ini fitness=%d\n", tspName, lambda, maxFunEval, seed, sol.fitness);
	printf("LSILS-C started on %s, lambda=%.4f\nMax fun eval=%.0f, seed=%d, Ini fitness=%d\n", tspName, lambda, maxFunEval, seed, sol.fitness);
	srand(seed);
	pseudoInst = tspInst;
	genNearestList(nearestList, nearestNum, tspInst);
	bits.resize(n);
	for (i = 0; i < n; i++) bits[i] = true;

	fprintf(fitFunEvalFile, "%d 0\n", sol.fitness);

	bestSol = sol;
	bestUpdated = false;
	solRealFit = sol.fitness;
	stopFlag = false;
	start_timer();
	funEvalNum = 0;
	lsNum = 0;
	int outputFEIndex = 0;

	while (1) {
		//3-Opt local search on the pseudo funtion
		sol.calFitness(pseudoInst);
		fliptour.flipper_init(sol.x, n);
		flippercost = sol.fitness;

		vector<int> fixedStart(n);
		for (i = 0; i < n; i++) fixedStart[i] = i;//GenRandomTour(fixedstart);
		while (1) {
			delta = 0;
			for (i = 0; i < n; i++) {
				A = fixedStart[i];
				if (bits[A]) {
					A_neighbor = nearestList[A];

					search_index = 0;
					while (search_index < 2) {
						if (search_index == 1)
							fliptour.flipper_reverse();
						B = fliptour.flipper_next(A);
						if (B != A_neighbor[nearestNum - 1]) {
							AB = pseudoInst.getDist(A, B);
							for (j = nearestNum - 1; j >= 0; j--) {//search in the neighbor list of A, from nearest to farthest
								C = A_neighbor[j];
								AC = pseudoInst.getDist(A, C);
								if (AC < AB) {
									D = fliptour.flipper_next(C);
									if (D != A) {

										CD = pseudoInst.getDist(C, D);
										BD = pseudoInst.getDist(B, D);

										//test 2-opt move
										delta = AC + BD - AB - CD;
										funEvalNum++;

										if (outputFEIndex < funEvalPoint.size()) {
											if (funEvalNum >= funEvalPoint[outputFEIndex]) {
												fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
												outputFEIndex++;
											}
										}

										if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
											stopFlag = true;//stopping criterion is met, stop the entire procedure
											break;
										}
										if (delta < 0) {
											fliptour.flipper_flip(D, A);
											realDelta = tspInst.getDist(A, C) + tspInst.getDist(B, D) - tspInst.getDist(A, B) - tspInst.getDist(C, D);
											improveType = 0;
											break;
										}

										//test 2.5 opt move
										E = fliptour.flipper_prev(C);
										if (E != B) {
											CE = pseudoInst.getDist(C, E);
											BC = pseudoInst.getDist(B, C);
											DE = pseudoInst.getDist(D, E);
											delta = AC + BC + DE - AB - CD - CE;
											funEvalNum++;

											if (outputFEIndex < funEvalPoint.size()) {
												if (funEvalNum >= funEvalPoint[outputFEIndex]) {
													fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
													outputFEIndex++;
												}
											}

											if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
												stopFlag = true;//stopping criterion is met, stop the entire procedure
												break;
											}
											if (delta < 0)
											{
												fliptour.flipper_flip(B, E);
												fliptour.flipper_flip(E, C);
												realDelta = tspInst.getDist(A, C) + tspInst.getDist(B, C) + tspInst.getDist(D, E)
													- tspInst.getDist(A, B) - tspInst.getDist(C, D) - tspInst.getDist(C, E);
												improveType = 1;
												break;
											}
										}

										//test 3 opt move
										D_neighbor = nearestList[D];
										for (k = nearestNum - 1; k >= 0; k--) {//search in the neighbor list of D, from nearest to farthest
											E = D_neighbor[k];
											if (E != C) {
												DE = pseudoInst.getDist(D, E);
												if (DE < AB + CD - AC) {
													if (fliptour.flipper_sequence(B, E, C))
														F = fliptour.flipper_next(E);
													else
														F = fliptour.flipper_prev(E);
													EF = pseudoInst.getDist(E, F);
													BF = pseudoInst.getDist(B, F);
													delta = AC + BF + DE - AB - CD - EF;
													funEvalNum++;

													if (outputFEIndex < funEvalPoint.size()) {
														if (funEvalNum >= funEvalPoint[outputFEIndex]) {
															fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
															outputFEIndex++;
														}
													}

													if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
														stopFlag = true;//stopping criterion is met, stop the entire procedure
														break;
													}
													if (delta < 0) {
														fliptour.flipper_flip(D, A);
														fliptour.flipper_flip(B, E);
														realDelta = tspInst.getDist(A, C) + tspInst.getDist(B, F) + tspInst.getDist(D, E)
															- tspInst.getDist(A, B) - tspInst.getDist(C, D) - tspInst.getDist(E, F);
														improveType = 2;
														break;
													}
												}
											}
											if (stopFlag == true) break;
										}
									}
								}
								if (stopFlag == true) break;
								if (delta < 0) break;
							}
						}
						if (stopFlag == true) break;
						if (delta < 0) break;
						search_index++;
					}
				}

				if (stopFlag == true) break;


				if (delta < 0) {
					flippercost += delta;
					solRealFit += realDelta;
					if (solRealFit < bestSol.fitness) {
						fliptour.flipper_cycle(bestSol.x);
						bestSol.fitness = solRealFit;
						bestUpdated = true;
					}

					if (improveType == 0) {//2opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = true;
					}
					else if (improveType == 1) {//2.5opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = bits[E] = true;
					}
					else if (improveType == 2) {//3opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = bits[E] = bits[F] = true;
					}
					else {
						fprintf(logFile, "ERROR! unclear improve type!!!\n");
					}

					break;
				}
				else {
					bits[A] = false;
				}
			}

			if (stopFlag == true) break;
			if (delta >= 0) break; //local optimum is found, stop current LS
		}
		fliptour.flipper_cycle(sol.x);
		sol.fitness = solRealFit;

		lsNum++;

		if (stopFlag == true) break;

		for (i = 0; i < n; i++) {
			if (bits[i] == true) {
				fprintf(logFile, "ERROR! don't=look-bits is not all-zero after LS!!!\n");
				break;
			}
		}

		//perturbation 
		doubleBridgeKickBits(sol, tspInst, bits);
		solRealFit = sol.fitness;

		if (bestUpdated) {
			//generate new pseudoInst based on new best solution
			homotopicConvexTransfer(pseudoInst, tspInst, bestSol, lambda);
			bestUpdated = false;

			//do NOT update the nearest list, because it can refleact the original information 
		}

	}

#if defined(__sun) || defined(linux)
	currentTime = return_CPU_time();
#else
	currentTime = return_elapsed_time();
#endif

	int the_fitness = sol.fitness;
	sol.calFitness(tspInst);
	int the_fitness2 = bestSol.fitness;
	bestSol.calFitness(tspInst);
	if (the_fitness != sol.fitness || the_fitness2 != bestSol.fitness) {
		fprintf(logFile, "ERROR! fitness miss match!!!\n");
	}

	fprintf(logFile, "Optimal fitness = %d, Best fitness=%d\n", globalOptFit, bestSol.fitness);
	fprintf(logFile, "Finish time=%.2f, LS num=%.0f, Fun eval num= %.0f\n\n", currentTime, lsNum, funEvalNum);
	fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);

	printf("LSILS-C finished, best fitness=%d\n", bestSol.fitness);
	printf("Finish time=%.2fs, LS num=%.0f, Fun eval num= %.0f\n\n", currentTime, lsNum, funEvalNum);

	fclose(logFile);
	fclose(fitFunEvalFile);

	return 0;
}

//Landscape Smoothing Iterated Local Search with 3-Opt local search and dynamically changing lambda values
int LSILS_3Opt_DynamicLambda(const char* tspName, const int iniSolSeed, const int seed, const int globalOptFit, const double maxFunEval,
	const int nearestNum, const double maxLambda, const char* outputLabel, const int printFunEvalIntvl) { // with homotopic transfer
	tspProblem tspInst;
	if (!tspInst.load(tspName)) {
		return 1;
	}

	int n = tspInst.getN();
	tspProblem pseudoInst;
	tspSolution sol, bestSol;
	bool bestUpdated;
	double currentTime;
	vector<vector<int> > nearestList;
	int i, j, k;
	int A, B, C, D, E, F;
	int AB, CD, EF, AC, BF, DE, BD, CE, BC;//distances
	vector<int> A_neighbor;
	vector<int> D_neighbor;
	int delta, realDelta, solRealFit;
	flipper_class fliptour;
	int flippercost;
	int search_index;
	bool stopFlag;
	double funEvalNum;
	double lsNum;
	vector<bool> bits;
	int improveType;

	int lambdaNum = round(maxLambda * 100) + 1;
	vector<double> allLambda(lambdaNum);
	for (i = 0; i < lambdaNum; i++) allLambda[i] = 0.01 * i;
	double changePeriod = ceil(maxFunEval / lambdaNum);

	char logFileName[200];
	sprintf(logFileName, "results/log_%s.txt", outputLabel);
	FILE* logFile = fopen(logFileName, "w");
	if (logFile == NULL) {
		return 1;
	}

	vector<int> funEvalPoint;
	funEvalPoint.resize(maxFunEval / printFunEvalIntvl - 1);
	funEvalPoint[0] = printFunEvalIntvl;
	for (i = 1; i < funEvalPoint.size(); i++) {
		funEvalPoint[i] = funEvalPoint[i - 1] + printFunEvalIntvl;
	}

	char fitFunEvalFileName[200]; //record how
	sprintf(fitFunEvalFileName, "results/fit_funeval_%s.txt", outputLabel);
	FILE* fitFunEvalFile = fopen(fitFunEvalFileName, "w");
	if (fitFunEvalFile == NULL) {
		return 1;
	}

	//initialization
	srand(iniSolSeed);
	sol.random(tspInst);
	fprintf(logFile, "%s all lambda = [", tspName);
	for (i = 0; i < allLambda.size(); i++) fprintf(logFile, " %.2f", allLambda[i]);
	fprintf(logFile, " ], change period = %.0f\n", changePeriod);
	fprintf(logFile, "Max fun eval = %.0f, seed = %d, Ini fitness = %d\n", maxFunEval, seed, sol.fitness);

	printf("LSILS-D started on %s, all lambda = [", tspName);
	for (i = 0; i < allLambda.size(); i++) printf(" %.2f", allLambda[i]);
	printf(" ], change period = %.0f\n", changePeriod);
	printf("Max fun eval = %.0f, seed = %d, Ini fitness = %d\n", maxFunEval, seed, sol.fitness);

	srand(seed);
	pseudoInst = tspInst;
	genNearestList(nearestList, nearestNum, tspInst);
	bits.resize(n);
	for (i = 0; i < n; i++) bits[i] = true;

	fprintf(fitFunEvalFile, "%d 0\n", sol.fitness);

	bestSol = sol;
	bestUpdated = false;
	solRealFit = sol.fitness;
	stopFlag = false;
	start_timer();
	funEvalNum = 0;
	lsNum = 0;
	int outputFEIndex = 0;

	double lambda = allLambda[0];

	while (1) {
		//3-Opt local search on the pseudo funtion
		sol.calFitness(pseudoInst);
		fliptour.flipper_init(sol.x, n);
		flippercost = sol.fitness;

		vector<int> fixedStart(n);
		for (i = 0; i < n; i++) fixedStart[i] = i;//GenRandomTour(fixedstart);
		while (1) {
			delta = 0;
			for (i = 0; i < n; i++) {
				A = fixedStart[i];
				if (bits[A]) {
					A_neighbor = nearestList[A];

					search_index = 0;
					while (search_index < 2) {
						if (search_index == 1)
							fliptour.flipper_reverse();
						B = fliptour.flipper_next(A);
						if (B != A_neighbor[nearestNum - 1]) {
							AB = pseudoInst.getDist(A, B);
							for (j = nearestNum - 1; j >= 0; j--) {//search in the neighbor list of A, from nearest to farthest
								C = A_neighbor[j];
								AC = pseudoInst.getDist(A, C);
								if (AC < AB) {
									D = fliptour.flipper_next(C);
									if (D != A) {

										CD = pseudoInst.getDist(C, D);
										BD = pseudoInst.getDist(B, D);

										//test 2-opt move
										delta = AC + BD - AB - CD;
										funEvalNum++;

										if (outputFEIndex < funEvalPoint.size()) {
											if (funEvalNum >= funEvalPoint[outputFEIndex]) {
												fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
												outputFEIndex++;
											}
										}

										if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
											stopFlag = true;//stopping criterion is met, stop the entire procedure
											break;
										}
										if (delta < 0) {
											fliptour.flipper_flip(D, A);
											realDelta = tspInst.getDist(A, C) + tspInst.getDist(B, D) - tspInst.getDist(A, B) - tspInst.getDist(C, D);
											improveType = 0;
											break;
										}

										//test 2.5 opt move
										E = fliptour.flipper_prev(C);
										if (E != B) {
											CE = pseudoInst.getDist(C, E);
											BC = pseudoInst.getDist(B, C);
											DE = pseudoInst.getDist(D, E);
											delta = AC + BC + DE - AB - CD - CE;
											funEvalNum++;

											if (outputFEIndex < funEvalPoint.size()) {
												if (funEvalNum >= funEvalPoint[outputFEIndex]) {
													fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
													outputFEIndex++;
												}
											}

											if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
												stopFlag = true;//stopping criterion is met, stop the entire procedure
												break;
											}
											if (delta < 0)
											{
												fliptour.flipper_flip(B, E);
												fliptour.flipper_flip(E, C);
												realDelta = tspInst.getDist(A, C) + tspInst.getDist(B, C) + tspInst.getDist(D, E)
													- tspInst.getDist(A, B) - tspInst.getDist(C, D) - tspInst.getDist(C, E);
												improveType = 1;
												break;
											}
										}

										//test 3 opt move
										D_neighbor = nearestList[D];
										for (k = nearestNum - 1; k >= 0; k--) {//search in the neighbor list of D, from nearest to farthest
											E = D_neighbor[k];
											if (E != C) {
												DE = pseudoInst.getDist(D, E);
												if (DE < AB + CD - AC) {
													if (fliptour.flipper_sequence(B, E, C))
														F = fliptour.flipper_next(E);
													else
														F = fliptour.flipper_prev(E);
													EF = pseudoInst.getDist(E, F);
													BF = pseudoInst.getDist(B, F);
													delta = AC + BF + DE - AB - CD - EF;
													funEvalNum++;

													if (outputFEIndex < funEvalPoint.size()) {
														if (funEvalNum >= funEvalPoint[outputFEIndex]) {
															fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);
															outputFEIndex++;
														}
													}

													if (funEvalNum >= maxFunEval || bestSol.fitness <= globalOptFit) {
														stopFlag = true;//stopping criterion is met, stop the entire procedure
														break;
													}
													if (delta < 0) {
														fliptour.flipper_flip(D, A);
														fliptour.flipper_flip(B, E);
														realDelta = tspInst.getDist(A, C) + tspInst.getDist(B, F) + tspInst.getDist(D, E)
															- tspInst.getDist(A, B) - tspInst.getDist(C, D) - tspInst.getDist(E, F);
														improveType = 2;
														break;
													}
												}
											}
											if (stopFlag == true) break;
										}
									}
								}
								if (stopFlag == true) break;
								if (delta < 0) break;
							}
						}
						if (stopFlag == true) break;
						if (delta < 0) break;
						search_index++;
					}
				}

				if (stopFlag == true) break;

				if (delta < 0) {
					flippercost += delta;
					solRealFit += realDelta;
					if (solRealFit < bestSol.fitness) {
						fliptour.flipper_cycle(bestSol.x);
						bestSol.fitness = solRealFit;
						bestUpdated = true;
					}

					if (improveType == 0) {//2opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = true;
					}
					else if (improveType == 1) {//2.5opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = bits[E] = true;
					}
					else if (improveType == 2) {//3opt improve
						bits[A] = bits[B] = bits[C] = bits[D] = bits[E] = bits[F] = true;
					}
					else {
						fprintf(logFile, "ERROR! unclear improve type!!!\n");
					}

					break;
				}
				else {
					bits[A] = false;
				}
			}

			if (stopFlag == true) break;
			if (delta >= 0) break; //local optimum is found, stop current LS
		}
		fliptour.flipper_cycle(sol.x);
		sol.fitness = solRealFit;

		lsNum++;

		if (stopFlag == true) break;

		for (i = 0; i < n; i++) {
			if (bits[i] == true) {
				fprintf(logFile, "ERROR! don't=look-bits is not all-zero after LS!!!\n");
				break;
			}
		}

		//perturbation 
		doubleBridgeKickBits(sol, tspInst, bits);
		solRealFit = sol.fitness;

		if (bestUpdated) {
			//generate new pseudoInst based on new best solution
			lambda = allLambda[floor(funEvalNum / changePeriod)];
			homotopicConvexTransfer(pseudoInst, tspInst, bestSol, lambda);
			bestUpdated = false;

			//do NOT update the nearest list, because it can refleact the original information 
		}
		else if (lambda != allLambda[floor(funEvalNum / changePeriod)]) {
			lambda = allLambda[floor(funEvalNum / changePeriod)];
			homotopicConvexTransfer(pseudoInst, tspInst, bestSol, lambda);
		}

	}

#if defined(__sun) || defined(linux)
	currentTime = return_CPU_time();
#else
	currentTime = return_elapsed_time();
#endif

	int the_fitness = sol.fitness;
	sol.calFitness(tspInst);
	int the_fitness2 = bestSol.fitness;
	bestSol.calFitness(tspInst);
	if (the_fitness != sol.fitness || the_fitness2 != bestSol.fitness) {
		fprintf(logFile, "ERROR! fitness miss match!!!\n");
	}

	fprintf(logFile, "Optimal fitness = %d, Best fitness=%d\n", globalOptFit, bestSol.fitness);
	fprintf(logFile, "Finish time=%.2f, LS num=%.0f, Fun eval num= %.0f\n\n", currentTime, lsNum, funEvalNum);
	fprintf(fitFunEvalFile, "%d %.0f\n", bestSol.fitness, funEvalNum);


	printf("LSILS-D finished, best fitness=%d\n", bestSol.fitness);
	printf("Finish time=%.2fs, LS num=%.0f, Fun eval num= %.0f\n\n", currentTime, lsNum, funEvalNum);

	fclose(logFile);
	fclose(fitFunEvalFile);

	return 0;
}

