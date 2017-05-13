#include "Monte_Carlo.h"


template <typename Gen>
std::vector<double> MCPROJ::monte_carlo(int n, Gen X)
{
	std::vector<double> result(3, 0.);
	double x;
	for (int j = 0; j < n; j++) {
		x = X(); // restriction sur X...
		result[0] += x;
		result[1] += x*x;
	}
	result[0] /= (double)n;
	result[1] = (result[1] - n*result[0] * result[0]) / (double)(n - 1);
	result[2] = 1.96*sqrt(result[1] / (double)n);
	return result;
};

template <typename Gen>
double MCPROJ::percentage_default(Gen X, double	q, double corr, double R, int Nb_CDS, double K1, double	K2) {
	double C = normal_CDF_inverse(q);
	double counter = 0;
	double x;
	
	double M = X();
	for (int i = 0; i < Nb_CDS; i++) {
		x = X();
		if (x < (C - corr*M) / sqrt(1 - corr*corr)) { counter++; }
	}
	counter = std::max((counter*(1 - R) / Nb_CDS) - K1, 0.0) / (K2 - K1);		//Normalization of counter in K1 and K2
	return counter;
};
