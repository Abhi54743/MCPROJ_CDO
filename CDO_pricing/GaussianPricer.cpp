#include "GaussianPricer.h"


MCPROJ::GaussianPricer::GaussianPricer(double q, double corr, double R, int	Nb_CDS, double K1, double K2) 
	: m_q(q), m_corr(corr), m_R(R), m_Nb_CDS(Nb_CDS), m_K1(K1), m_K2(K2){
	

	m_G = new normal_rv(m_gen, normal_dist(0, 1));
	m_gen.seed(static_cast<unsigned int>(0));//to remove later
	

	m_C = normal_CDF_inverse(q);

};


double MCPROJ::GaussianPricer::expected_Loss() {

	if (m_R == 1 || m_K1 / (1 - m_R) > 1 || m_K2 / (1 - m_R) > 1) { throw std::exception("Unreasonable K and R parameters: R=1 or K/(1-R)>1"); }
	else {
		double K11 = normal_CDF_inverse(m_K1 / (1 - m_R));
		double K21 = normal_CDF_inverse(m_K2 / (1 - m_R));
		return (BivariateNormalCDF(-K11, m_C, -sqrt(1 - m_corr*m_corr)) - BivariateNormalCDF(-K21, m_C, -sqrt(1 - m_corr*m_corr))) * (1 - m_R) / (m_K2 - m_K1);
	}



};


double MCPROJ::GaussianPricer::operator()() {

	return MCPROJ::percentage_default(*m_G, *m_G, m_C, m_corr, m_R, m_Nb_CDS, m_K1, m_K2);

};

std::vector<double> MCPROJ::GaussianPricer::expected_LossMC(int N) {

	MCPROJ::GaussianPricer gauss(m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2);

	return MCPROJ::monte_carlo(N, gauss);

};
