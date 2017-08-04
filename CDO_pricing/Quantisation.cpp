#include "Quantisation.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

std::vector<double> MCPROJ::Quantisator::price(int N)
{
	if (m_dist_type == "Gaussian") {

		Quantisator Quant(m_dist_type, m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2);

		return MCPROJ::monte_carlo(N, Quant);
	}
	else if (m_dist_type == "NIG") {
		Quantisator Quant(m_dist_type, m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2, m_alpha, m_beta, m_gridSize);

		return MCPROJ::monte_carlo(N, Quant);
	}
	else throw std::exception("only Gaussian or NIG are supported");
}

double MCPROJ::Quantisator::operator()()
{
	if (m_dist_type == "Gaussian") {
		double u = m_normal_variable->operator()();
		double expectedPrice = conditionallyExpectedPrice(u, m_K1) - conditionallyExpectedPrice(u, m_K2);
		return expectedPrice / (m_K2 - m_K1);
	}
	else if (m_dist_type == "NIG") {
		double u = m_NIG_variable->operator()();
		double expectedPrice = conditionallyExpectedPrice(u, m_K1) - conditionallyExpectedPrice(u, m_K2);
		return expectedPrice / (m_K2 - m_K1);
	}
	else throw std::exception("only Gaussian or NIG are supported");
}

double MCPROJ::Quantisator::conditionallyExpectedPrice(double u, double K)
{
	double mu;
	double sigma;
	if (m_dist_type == "Gaussian") {

		boost::math::normal normal;

		double toto = (boost::math::cdf(normal, (m_C - m_corr*u) / sqrt(1 - (m_corr*m_corr))));
		mu = (1 - m_R)*toto / m_Nb_CDS;
		sigma = (1 - m_R)*(1 - m_R)*(toto - toto*toto) / (m_Nb_CDS*m_Nb_CDS);
	}

	double alpha = (1 - m_R) / m_Nb_CDS;

	std::vector<double> grid_old(1);
	std::vector<double> grid_new;
	std::vector<double> weights_old(1, 1.0);
	std::vector<double> weights_new;

	for (int k = 0; k < m_Nb_CDS; k++)
	{

		double mu_k = mu * (k+1);
		double sqrt_sigma_k = sqrt(sigma * (k + 1) * (k + 1));

		grid_new.resize(m_Grids[k].size() + 1);

		std::transform(m_Grids[k].begin(), m_Grids[k].end(), grid_new.begin() + 1, std::bind1st(std::multiplies<double>(), sqrt_sigma_k));
		std::transform(grid_new.begin() + 1, grid_new.end(), grid_new.begin() + 1, std::bind1st(std::plus<double>(), mu_k));

		grid_new.push_back(alpha * (k + 1));
		
		weights_new.resize(m_Grids[k].size() + 2);
		
		for (j = 0; j < grid_old.size(); j++)
		{
		}

		grid_old = grid_new;
		weights_old = weights_new;

	}
	
	
	
	
	
	std::vector<double> finalWeights;




	double res=0.0;
	for (int i = 0; i < m_Grids.back().size(); i++)
	{
		res += std::max((m_Grids.back()[i] - K), 0.0) * finalWeights[i];
	}
	return res;
}
