#include "NIGPricer.h"



MCPROJ::NIGPricer::NIGPricer(double q, double corr, double R, int	Nb_CDS, double K1, double K2, double alpha, double beta, int gridSize)
	: m_q(q), m_corr(corr), m_R(R), m_Nb_CDS(Nb_CDS), m_K1(K1), m_K2(K2), m_alpha(alpha), m_beta(beta) {

	
	double gamma = sqrt(alpha*alpha - beta*beta);

	if (abs(beta) > alpha) { throw std::exception(" value of abs(beta) is greater than alpha "); }

	m_NIG_M = new NIG_rv(alpha, beta, -beta*gamma*gamma/(alpha*alpha), gamma*gamma*gamma/(alpha*alpha), gridSize);

	double factor = sqrt(1 - (corr*corr)) / corr;
	m_NIG_X = new NIG_rv(factor*alpha, factor*beta, -factor*beta*gamma*gamma / (alpha*alpha), factor*gamma*gamma*gamma / (alpha*alpha), gridSize);
	
	NIG_rv  NIG_A(alpha/corr, beta/corr, -beta*gamma*gamma / (alpha*alpha*corr), gamma*gamma*gamma / (alpha*alpha*corr), gridSize);

	m_C = NIG_A.inverseCDF(q); //
};





double MCPROJ::NIGPricer::operator()() {

	return MCPROJ::percentage_default(*m_NIG_M, *m_NIG_X, m_C, m_corr, m_R, m_Nb_CDS, m_K1, m_K2);

}
double MCPROJ::NIGPricer::expected_Loss()
{
	return 0.0;
}
;

std::vector<double> MCPROJ::NIGPricer::expected_LossMC(int N) {

	MCPROJ::NIGPricer NIG(m_q, m_corr, m_R, m_Nb_CDS, m_K1, m_K2, m_alpha, m_beta, m_NIG_M->getSizeListOfValues());

	return MCPROJ::monte_carlo(N, NIG);

};
