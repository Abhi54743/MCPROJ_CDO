#include "SteinPricer.h"

std::vector<double> GaussianSteinPricer::price(int N)
{
	if (m_dist_type == "Gaussian") {

		GaussianSteinPricer gaussStein(m_dist_type, m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2);

		return MCPROJ::monte_carlo(N, gaussStein);
	}
	else if (m_dist_type == "NIG") {
		GaussianSteinPricer gaussStein(m_dist_type, m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2, m_alpha, m_beta, m_gridSize);

		return MCPROJ::monte_carlo(N, gaussStein);
	}
	else throw std::exception("only Gaussian or NIG are supported");
}

double GaussianSteinPricer::operator()()
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

double GaussianSteinPricer::conditionallyExpectedPrice(double u, double K)
{
	double mu, sigma, third_moment;

	if (m_dist_type == "Gaussian") { 
		
		boost::math::normal normal;
		
		double toto = (boost::math::cdf(normal, (m_C-m_corr*u)/sqrt(1-(m_corr*m_corr))));
		mu = (1-m_R)*toto/m_Nb_CDS;
		sigma = (1 - m_R)*(1 - m_R)*(toto - toto*toto) / (m_Nb_CDS*m_Nb_CDS);

		third_moment = pow(1 - m_R, 3)* (toto - 3 * toto*toto + 2 * pow(toto, 3))/pow(m_Nb_CDS, 3);	
	}
	else if (m_dist_type == "NIG") { 
				
		double toto =  m_NIG_X->CDF((m_C - m_corr*u) / sqrt(1 - (m_corr*m_corr)));
		mu = (1 - m_R)*toto / m_Nb_CDS;
		sigma = (1 - m_R)*(1 - m_R)*(toto - toto*toto) / (m_Nb_CDS*m_Nb_CDS);

		third_moment = pow(1 - m_R, 3)* (toto - 3 * toto*toto + 2 * pow(toto, 3)) / pow(m_Nb_CDS, 3);


	}
	

	double displaced_K = K - m_Nb_CDS*mu;

	boost::math::normal normal(0, m_Nb_CDS*sigma);
	
	double gaussianapprox = (m_Nb_CDS*sigma*boost::math::pdf(normal, displaced_K)) -displaced_K*boost::math::cdf(normal, -displaced_K);

	double gaussianerrorcorrection = third_moment*displaced_K*boost::math::pdf(normal, displaced_K) / (6 * sigma);

	return gaussianapprox + gaussianerrorcorrection;

}
