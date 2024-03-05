/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem.
IEEE Transactions on Cybernetics, 52(1), 495-507.*/

#pragma once

#include <vector>
#include <math.h>
#include<stdlib.h>

typedef short int sint;
typedef unsigned int uint;

#define PI 3.14159265
#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CC_EUCLIDEAN      0
#define CC_EUCLIDEAN_CEIL 1
#define CC_GEOGRAPHIC     2
#define CC_ATT            3
#define CC_MATRIXNORM	  4

#define MATRIX_FULL_MATRIX     0
#define MATRIX_UPPER_ROW       1
#define MATRIX_LOWER_ROW       2
#define MATRIX_UPPER_DIAG_ROW  3
#define MATRIX_LOWER_DIAG_ROW  4
#define MATRIX_UPPER_COL       5
#define MATRIX_LOWER_COL       6
#define MATRIX_UPPER_DIAG_COL  7
#define MATRIX_LOWER_DIAG_COL  8


class tspProblem
{
public:
	int n;
	int	edgeWeightType;  //edge weight type

	std::vector<double>	X; //x coord
	std::vector<double>	Y; //y coord
	std::vector< std::vector<int> > dist; //distance matrix

public:
	bool load(const char* filename);
	int getN() const;
	sint getDist(const int i, const int j) const;
};

