#include "Generators.h"

MCPROJ::NIG_rv::NIG_rv(double alpha, double beta, double mu, double delta, int gridSize) :m_mu(mu), m_alpha(alpha), m_beta(beta), m_delta(delta) {
	// creating a normal random variable with mean 0 and variance 1
	m_Gaussian = new normal_rv(m_Gen, normal_dist(0, 1));
	// creating uniform random variable between 0 and 1
	m_U_dist = new uni_01_dist(0.0, 1.0);

	m_Gen.seed(static_cast<unsigned int>(0));//to remove later

	if (abs(beta) > alpha) { throw std::exception(" value of abs(beta) is greater than alpha "); }

	else if (delta < 1e-30) { throw std::exception(" value of delta is negative "); }
	double gamma = sqrt(alpha*alpha - beta*beta);

	m_IG = inverse_gaussian_dist(m_delta*gamma, gamma*gamma);

	for(int i=0; i<gridSize; i++)
	{
		m_Listofvalues.push_back(boost::math::quantile(m_IG, ((double)i)/gridSize));
	}

}

double MCPROJ::NIG_rv::operator()()
{	
	// getting uniform random variable value between 0 and 1
	double x = m_U_dist->operator()(m_Gen);
	// getting value of IG distribution

	//double y = boost::math::quantile(m_IG, x); //slower version
	double y = inverseIG(x, m_Listofvalues);

	double lol = m_Gaussian -> operator()();

	//returning value of NIG distribution
	return (lol*sqrt(y)) + m_mu + (m_beta*y);
}


int MCPROJ::NIG_rv::getSizeListOfValues() {
	return m_Listofvalues.size();

}
double MCPROJ::NIG_rv::CDF(double x)
{	
	std::vector<double> abs(m_Listofvalues.size()), weight(m_Listofvalues.size());
	gauleg(0, 1, abs, weight);
	
	boost::math::normal normal;

	double value = 0.0;

	for (int i = 0; i < m_Listofvalues.size(); i++) {
		
		double toto = log(abs[i]);
		value += (weight[i] *  boost::math::cdf(normal, (x-(m_mu - m_beta*toto)) / sqrt(-toto) ) * boost::math::pdf(m_IG, -toto) / abs[i]);
	
	}

	return value;
}



double MCPROJ::NIG_rv::inverseCDF(double q)
{
	double x0 = 0.5;
	double x1;

	double tol = 1.0e-12;
	double err = 1 + tol;

	double dx = 1.0e-3;

	while (err > tol) {
		
		double deriv = (CDF(x0 + dx) - CDF(x0 - dx)) / (2 * dx); 
		x1 = x0 - (CDF(x0) - q) / deriv;
		err = CDF(x1) - q;
		x0 = x1;
	
	}
	return x1;

}
