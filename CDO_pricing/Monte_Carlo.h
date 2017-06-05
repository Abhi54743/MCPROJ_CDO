/*!
Monte carlo Routines, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 11/03/2017

Version 1.0
*/
#pragma once
#ifndef _MONTE_CARLO_H_
#define _MONTE_CARLO_H_

#include <vector>

namespace MCPROJ {

	template <typename Gen>
	std::vector<double> monte_carlo(int n, Gen X)
	{
		std::vector<double> result(3, 0.0);
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
	double percentage_default(Gen M, Gen X, double	C, double corr, double R, int Nb_CDS, double K1, double	K2) {
		double counter = 0;
		double x;

		double Mvalue = M();

		for (int i = 0; i < Nb_CDS; i++) {
			x = X();
			if (x < (C - corr*Mvalue) / sqrt(1 - corr*corr)) { counter++; }
		}
		counter = std::max((counter*(1 - R) / Nb_CDS) - K1, 0.0) / (K2 - K1);		//Normalization of counter in K1 and K2
		//counter *= (1 - R) / Nb_CDS;
		return counter;
	};


};

#endif // !_MONTE_CARLO_H_

