#include "Generators.h"

MCPROJ::NIG_rv::NIG_rv(double alpha, double beta, double mu, double delta):m_mu(mu), m_alpha(alpha), m_beta(beta), m_delta(delta){
	m_Gaussian = new normal_rv(m_Gen, normal_dist(0, 1));

	m_Gen.seed(static_cast<unsigned int>(0));//to remove later

	double gamma = sqrt(alpha*alpha - beta*beta);

	m_IG = inverse_gaussian_dist(m_delta*gamma, gamma*gamma);

}

double MCPROJ::NIG_rv::operator()()
{	

	normal_rv gauss = *m_Gaussian;

	double y = boost::math::quantile(m_IG, m_Gen());

	return (sqrt(y)*gauss()) + m_mu + (m_beta*y);
}
