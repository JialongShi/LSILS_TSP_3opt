/*Shi, J., Sun, J., Zhang, Q., & Ye, K. (2020). Homotopic convex transformation: A new landscape smoothing method for the traveling salesman problem.
IEEE Transactions on Cybernetics, 52(1), 495-507.*/

#include "tspProblem.h"
#include <cstring>
#include <fstream>

using namespace std;

//truncate the decimal part of a double value x 
double dtrunc(double x)
{
	int k;

	k = (int)x;
	x = (double)k;
	return x;
}

//load a TSP instance from a file specified by the filename
bool tspProblem::load(const char* filename)
{
	ifstream tspStream(((string)filename).c_str());
	if (!tspStream) {
		return false;
	}
	int null, distTemp;
	int i, j;
	char key[256];

	while (strncmp(key, "TYPE", 4))
	{
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1))
	{
		tspStream >> key;
	}
	if (strncmp(key, "TSP", 3))
	{
		return false;
	}
	while (strncmp(key, "DIMENSION", 9))
	{
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1))
	{
		tspStream >> key;
	}
	n = atoi(key);

	while (strncmp(key, "EDGE_WEIGHT_TYPE", 16))
	{
		tspStream >> key;
	}
	memset(key, 0, sizeof(key));
	tspStream >> key;
	if (!strncmp(key, ":", 1))
	{
		tspStream >> key;
	}

	if (!strncmp(key, "EXPLICIT", 8))
	{
		edgeWeightType = CC_MATRIXNORM;
	}
	else if (!strncmp(key, "EUC_2D", 6))
	{
		edgeWeightType = CC_EUCLIDEAN;
	}
	else if (!strncmp(key, "CEIL_2D", 7))
	{
		edgeWeightType = CC_EUCLIDEAN_CEIL;
	}
	else if (!strncmp(key, "GEO", 3))
	{
		edgeWeightType = CC_GEOGRAPHIC;
	}
	else if (!strncmp(key, "ATT", 3))
	{
		edgeWeightType = CC_ATT;
	}
	else
	{
		return false; //just support 5 case now
	}

	dist.resize(n);
	for (i = 0; i < n; ++i) {
		dist[i].resize(n);
	}

	if (edgeWeightType != CC_MATRIXNORM)
	{
		while (strncmp(key, "NODE_COORD_SECTION", 18))
		{
			tspStream >> key;
		}

		X.resize(n);
		Y.resize(n);

		for (i = 0; i < n; i++)
		{
			double Xi, Yi;
			tspStream >> null >> Xi >> Yi;
			X[i] = Xi;
			Y[i] = Yi;
		}


		switch (edgeWeightType)
		{
		case CC_EUCLIDEAN:
		{
			for (i = 0; i < n; i++)
			{
				dist[i][i] = 0;
				for (j = 0; j < i; j++)
				{
					distTemp = (int)floor(sqrt((X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])) + 0.5);
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		case CC_EUCLIDEAN_CEIL:
		{
			for (int i = 0; i < n; i++)
			{
				dist[i][i] = 0;
				for (int j = 0; j < i; j++)
				{
					distTemp = (int)ceil(sqrt((X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])));
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		case CC_GEOGRAPHIC:
		{
			for (int i = 0; i < n; i++)
			{
				dist[i][i] = 0;
				for (int j = 0; j < i; j++)
				{
					double deg, min;
					double lati, latj, longi, longj;
					double q1, q2, q3;
					double x1 = X[i], x2 = X[j], yy1 = Y[i], yy2 = Y[j];

					deg = dtrunc(x1);
					min = x1 - deg;
					lati = PI * (deg + 5.0 * min / 3.0) / 180.0;
					deg = dtrunc(x2);
					min = x2 - deg;
					latj = PI * (deg + 5.0 * min / 3.0) / 180.0;

					deg = dtrunc(yy1);
					min = yy1 - deg;
					longi = PI * (deg + 5.0 * min / 3.0) / 180.0;
					deg = dtrunc(yy2);
					min = yy2 - deg;
					longj = PI * (deg + 5.0 * min / 3.0) / 180.0;

					q1 = cos(longi - longj);
					q2 = cos(lati - latj);
					q3 = cos(lati + latj);
					distTemp = (int)(6378.388 * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
						+ 1.0);
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		case CC_ATT:
		{
			for (int i = 0; i < n; i++)
			{
				dist[i][i] = 0;
				for (int j = 0; j < i; j++)
				{
					double xd = X[i] - X[j];
					double yd = Y[i] - Y[j];
					double rij = sqrt((xd * xd + yd * yd) / 10.0);
					double tij = dtrunc(rij);

					if (tij < rij)
						distTemp = (int)tij + 1;
					else
						distTemp = (int)tij;
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		}
	}
	else if (edgeWeightType == CC_MATRIXNORM)
	{
		int matForm;
		while (strncmp(key, "EDGE_WEIGHT_FORMAT", 17))
		{
			tspStream >> key;
		}
		memset(key, 0, sizeof(key));
		tspStream >> key;
		if (!strncmp(key, ":", 1))
		{
			tspStream >> key;
		}
		if (!strncmp(key, "FULL_MATRIX", 11))
		{
			matForm = MATRIX_FULL_MATRIX;
		}
		else if (!strncmp(key, "UPPER_ROW", 9))
		{
			matForm = MATRIX_UPPER_ROW;
		}
		else if (!strncmp(key, "LOWER_ROW", 9))
		{
			matForm = MATRIX_LOWER_ROW;
		}
		else if (!strncmp(key, "UPPER_DIAG_ROW", 14))
		{
			matForm = MATRIX_UPPER_DIAG_ROW;
		}
		else if (!strncmp(key, "LOWER_DIAG_ROW", 14))
		{
			matForm = MATRIX_LOWER_DIAG_ROW;
		}
		else if (!strncmp(key, "UPPER_COL", 9))
		{
			matForm = MATRIX_UPPER_COL;
		}
		else if (!strncmp(key, "LOWER_COL", 9))
		{
			matForm = MATRIX_LOWER_COL;
		}
		else if (!strncmp(key, "UPPER_DIAG_COL", 14))
		{
			matForm = MATRIX_UPPER_DIAG_COL;
		}
		else if (!strncmp(key, "LOWER_DIAG_COL", 14))
		{
			matForm = MATRIX_LOWER_DIAG_COL;
		}


		while (strncmp(key, "EDGE_WEIGHT_SECTION", 17))
		{
			tspStream >> key;
		}
		switch (matForm)
		{
		case MATRIX_FULL_MATRIX: // only for STSP!!!!
		{
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					tspStream >> distTemp;
					if (j <= i)
					{
						dist[i][j] = distTemp;
						dist[j][i] = distTemp;
					}
				}
			}
			break;
		}
		case MATRIX_UPPER_ROW:
		{
			for (i = 0; i < n; i++)
			{
				dist[i][i] = 0;
			}
			for (i = 0; i < n - 1; i++)
			{
				for (j = i + 1; j < n; j++)
				{
					tspStream >> distTemp;
					dist[j][i] = distTemp;
					dist[i][j] = distTemp;
				}
			}
			break;
		}
		case MATRIX_LOWER_ROW:
		{
			for (int i = 0; i < n; i++)
			{
				dist[i][i] = 0;
			}
			for (int i = 1; i < n; i++)
			{
				for (int j = 0; j < i; j++)
				{
					tspStream >> distTemp;
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		case MATRIX_UPPER_DIAG_ROW:
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = i; j < n; j++)
				{
					tspStream >> distTemp;
					dist[j][i] = distTemp;
					dist[i][j] = distTemp;
				}
			}
			break;
		}
		case MATRIX_LOWER_DIAG_ROW:
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					tspStream >> distTemp;
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		case MATRIX_UPPER_COL:
		{
			for (int i = 0; i < n; i++)
			{
				dist[i][i] = 0;
			}
			for (int j = 1; j < n; j++)
			{
				for (int i = 0; i < j; i++)
				{
					tspStream >> distTemp;
					dist[j][i] = distTemp;
					dist[i][j] = distTemp;
				}
			}
			break;
		}
		case MATRIX_LOWER_COL:
		{
			for (int i = 0; i < n; i++)
			{
				dist[i][i] = 0;
			}
			for (int j = 0; j < n - 1; j++)
			{
				for (int i = j + 1; i < n; i++)
				{
					tspStream >> distTemp;
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		case MATRIX_UPPER_DIAG_COL:
		{
			for (int j = 0; j < n; j++)
			{
				for (int i = 0; i <= j; i++)
				{
					tspStream >> distTemp;
					dist[j][i] = distTemp;
					dist[i][j] = distTemp;
				}
			}
			break;
		}
		case MATRIX_LOWER_DIAG_COL:
		{
			for (int j = 0; j < n; j++)
			{
				for (int i = j; i < n; i++)
				{
					tspStream >> distTemp;
					dist[i][j] = distTemp;
					dist[j][i] = distTemp;
				}
			}
			break;
		}
		default:
			return false;
		}
	}
	else
	{
		return false;
	}

	return true;
}

//get the the number of cities
int tspProblem::getN() const {
	return n;
}

//get the distance between city i and city j
sint tspProblem::getDist(const int i, const int j) const {
	return dist[i][j];
}



