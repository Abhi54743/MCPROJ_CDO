#include "NIGPricer.h"



MCPROJ::NIGPricer::NIGPricer(double q, double corr, double R, int	Nb_CDS, double K1, double K2, double alpha, double beta)
	: m_q(q), m_corr(corr), m_R(R), m_Nb_CDS(Nb_CDS), m_K1(K1), m_K2(K2), m_alpha(alpha), m_beta(beta) {

	double gamma = sqrt(alpha*alpha - beta*beta);

	m_NIG_M = NIG_rv(alpha, beta, -beta*gamma*gamma/(alpha*alpha), gamma*gamma*gamma/(alpha*alpha));

	double lol = sqrt(1 - (corr*corr)) / corr;
	m_NIG_X = NIG_rv(lol*alpha, lol*beta, -lol*beta*gamma*gamma / (alpha*alpha), lol*gamma*gamma*gamma / (alpha*alpha));

};





double MCPROJ::NIGPricer::operator()() {

	return MCPROJ::percentage_default(m_NIG_M, m_NIG_X, m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2);

};

std::vector<double> MCPROJ::NIGPricer::expected_LossMC(int N) {

	MCPROJ::NIGPricer NIG(m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2, m_alpha, m_beta);

	return MCPROJ::monte_carlo(N, NIG);

};
