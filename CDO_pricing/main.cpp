#define _USE_MATH_DEFINES

#include<cstdio>
#include<cmath>
#include<sstream>
#include<iostream>
#include<vector>
#include "GaussianPricer.h"
#include "NIGPricer.h"
#include "zero_search.h"


using namespace MCPROJ;

/*
template <typename Gen>
class percent_default_vr {
public:
	percent_default_vr(Gen X, double K1, double K2, std::vector<double> theta) :X(X), K1(K1), K2(K2), theta(theta) {};
	~percent_default_vr() {};
	double operator()() {
		double C = normal_CDF_inverse(q);
		double counter = 0;
		double x;

		double M = X();
		double theta_prod_X = theta[0]*M;
		double scalar_pdt_theta = theta[0] * theta[0];

		//if(Gen.name== normal)
		
		for (int i = 0; i < Nb_CDS; i++) {
			x = X();
			if (x+theta[i+1] < (C - corr*(M+theta[0])) / sqrt(1 - corr*corr)) { counter ++; }
			theta_prod_X += x*theta[1+i];
			scalar_pdt_theta += theta[1+i] * theta[1+i];
		}

		counter *= exp(-theta_prod_X-(scalar_pdt_theta / 2));
		counter = std::max((counter*(1 - R) / Nb_CDS) - K1, 0.0) / (K2 - K1);			//Normalization of counter in K1 and K2
		return counter;
	}
private:
	Gen X;
	double K1, K2;
	std::vector<double> theta;
};




template <typename Gen>
std::vector<double> grad_f(int i, std::vector<double> theta, std::vector<std::vector<double>> stored,
							int Nb_CDS, Gen X, double q, double corr, double K1, double K2) {//ALL these will be attributes 
	std::vector<double> res(Nb_CDS + 1, 0);
	double C = normal_CDF_inverse(q);
	for (int j = 0; j < Nb_CDS + 1; j++) {
		for (int k = 0; k < 1000; k++) {
			int counter = 0;
			double scalar_product = theta[0] * stored[k][0];
			double square_theta = theta[0] * theta[0];
			for (int r = 1; r < Nb_CDS+1; r++) {
				if (stored[k][r] < (C - corr*stored[k][0]) / sqrt(1 - corr*corr)) { counter++; }
				scalar_product += theta[r] * stored[k][r];
				square_theta += theta[r] * theta[r];
			}
			double A = std::max((counter*(1 - R) / Nb_CDS) - K1, 0.0) / (K2 - K1);
			A *= exp(-scalar_product + square_theta / 2);
			res[j] += A*((theta[j] - stored[k][j])*(theta[i] - stored[k][i]) + (i == j ? 1 : 0));
			res[j] /= 1000;
		}
	}
	return res;
}

template <typename Gen>
std::vector<double> test(std::vector<double> theta, std::vector<std::vector<double>> stored,
	int Nb_CDS, Gen X, double q, double corr, double K1, double K2) {//ALL these will be attributes 
	std::vector<double> res(Nb_CDS + 1, 0);
	double C = normal_CDF_inverse(q);
	for (int j = 0; j < Nb_CDS + 1; j++) {
		for (int k = 0; k < 1000; k++) {
			int counter = 0;
			double scalar_product = theta[0] * stored[k][0];
			double square_theta = theta[0] * theta[0];
			for (int r = 1; r < Nb_CDS + 1; r++) {
				if (stored[k][r] < (C - corr*stored[k][0]) / sqrt(1 - corr*corr)) { counter++; }
				scalar_product += theta[r] * stored[k][r];
				square_theta += theta[r] * theta[r];
			}
			double A = std::max((counter*(1 - R) / Nb_CDS) - K1, 0.0) / (K2 - K1);
			A *= exp(-scalar_product + square_theta / 2);
			res[j] += A*(theta[j] - stored[k][j]);
			res[j] /= 1000;
		}
	}
	return res;
}

template <typename Gen>
std::vector<double> optimizer_Theta(int Nb_CDS, Gen X, double q, double corr, double K1, double K2) //ALL these will be attributes
{


	double C = normal_CDF_inverse(q);
	std::vector<std::vector<double>> stored (1000, std::vector<double>(Nb_CDS+1));
	for (int i = 0; i < 1000; i++) {
		for (int j = 0; j < Nb_CDS + 1; j++)
		{
			stored[i][j] = X();
		}
	}



	//algorithm
	std::vector<double> theta0(Nb_CDS + 1, 0);

	generator gen;
	boost::random::uniform_int_distribution<> dist(0, 100);

	double step = 1e-1;

	double t = 0;
	std::vector<double> eps(Nb_CDS + 1, 1e-5);
	while (test(theta0, stored, Nb_CDS, X, q, corr, K1, K2)> eps) { //ALL these will be attributes ) > eps) {
		t += 1;
		std::cout << "iteration " << t << std::endl;
		int i = dist(gen);
		for(int y=0; y < Nb_CDS + 1; y++){
		theta0[y] -= step / (sqrt(t)) * grad_f(i, theta0, stored,
											Nb_CDS, X, q, corr, K1, K2)[y]; //ALL these will be attributes );
		}
	} 

	return theta0;

};
*/

int main() {
	
	double q = 0.5;   //default probability
	double corr = 0.3;	//correlation between CDS
	double R = 0;		//recovery rate
	int Nb_CDS = 100;
	double	K1 = 0.2;
	double	K2 = 0.7;
	double alpha = 0.5;
	double beta = 1;
	



	GaussianPricer Gauss(q, corr, R, Nb_CDS, K1, K2);

	std::cout << Gauss.expected_Loss() << std::endl ;

	std::vector<double> vect = Gauss.expected_LossMC(10000);

	std::cout << vect[0] <<"    " << vect[1] << "    " << vect[2] << std::endl;

	
	NIGPricer NIG(q, corr, R, Nb_CDS, K1, K2, alpha, beta);

	std::vector<double> vect2 = NIG.expected_LossMC(10000);

	std::cout << vect2[0] << "    " << vect2[1] << "    " << vect2[2] << std::endl;
	
	int lol;
	std::cin >> lol;

	return 0;

	/*generator gen;
	normal_rv G(gen, normal_dist(0, 1));
	gen.seed(static_cast<unsigned int>(std::time(0)));
	percent_default<normal_rv> pd1(G, 0.2, 0.7);                        //percentage default generator with normal random

	std::vector<double>theta = optimizer_Theta(Nb_CDS, G, q, corr, 0.2, 0.7);
	
	percent_default_vr<normal_rv> pd2(G, 0.2, 0.7, theta);   //percentage default generator with normal random

	
	std::vector<double> result = monte_carlo(1000, pd1);

	std::cout << "without vr" << result[0] <<'\t'<< result[1] << '\t' << result[2] << std::endl;

	result = monte_carlo(1000, pd2);
	std::cout << "vr" <<result[0] << '\t' << result[1] << '\t' << result[2] << std::endl;
	
	std::cout << expected_loss(0.2, 0.7) << std::endl;
	*/


}