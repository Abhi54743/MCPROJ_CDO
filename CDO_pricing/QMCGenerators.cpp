#include "QMCGenerators.h"

double MCPROJ::Halton1D::operator()()
{
	size_t i = m_index;
	double f = 1.0;
	double r = 0.0;

	while (i>0){
		f /= m_base;
		r += f*(i%m_base);
		i /= m_base;
	}

	m_index += 1;
	return r;
}

std::vector<double> MCPROJ::Halton2D::operator()()
{
	std::vector<size_t> i(2, m_index);
	std::vector<double> f(2, 1.0);
	std::vector<double> r(2);

	while (i[0]>0) {
		f[0] /= m_base_1;
		r[0] += f[0]*(i[0]%m_base_1);
		i[0] /= m_base_1;
	}


	while (i[1]>0) {
		f[1] /= m_base_2;
		r[1] += f[1] * (i[1] % m_base_2);
		i[1] /= m_base_2;
	}

	m_index += 1;
	return r;
}

std::vector<double> MCPROJ::Halton2DGauss::operator()()
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;
	std::vector<double> z;
	double u1, u2;
	do
	{
		std::vector<double> temp = m_H2D ->operator()();
		u1 = temp[0];
		u2 = temp[1];
	} while (u1 <= epsilon);

	z.push_back((sqrt(-2.0 * log(u1)) * cos(two_pi * u2))*m_sigma + m_mu);
	z.push_back((sqrt(-2.0 * log(u1)) * sin(two_pi * u2))*m_sigma + m_mu);
	return z;

}

std::vector<double> MCPROJ::Halton2DNIG::operator()()
{
	// getting uniform random variable value between 0 and 1
	std::vector<double> x = m_H2D->operator()();
	// getting value of IG distribution

	std::vector<double> y;
	y.push_back(inverseIG(x[0], m_Listofvalues_1));
	y.push_back(inverseIG(x[1], m_Listofvalues_2));

	std::vector<double> lol = m_HG2D -> operator()();

	//returning value of NIG distribution
	std::vector<double> r;
	r.push_back((lol[0]*sqrt(y[0])) + m_mu_1 + (m_beta_1*y[0]));
	r.push_back((lol[1] * sqrt(y[1])) + m_mu_2 + (m_beta_2*y[1]));
	return r;
}

std::vector<double> MCPROJ::Kakutani2D::double2piadic(double x, size_t base)
{
	std::vector<int> res(m_decimals);
	div_t divresult;
	for (int i = 0; i < m_decimals; i++)
	{
		res[i] = (int)floor(x * base);
		x -= res[i];
	}
	return res;
}

double MCPROJ::Kakutani2D::piadic2double(std::vector<int> p, site_t base)
{
	double res = 0.0;

	for (int i = 0; i < m_decimals; i++)
	{
		res += p[i] * pow(base, -(i + 1));
	}

	return res;
}


std::vector<double> MCPROJ::Kakutani2D::operator()()
{
	std::vector<double> r(2);

	m_previous[0][0] += 1;
	m_previous[1][0] += 1;

	for (int i = 0; i < m_decimals - 1; i++)
	{
		if (m_previous[0][i] == m_base_1)
		{
			m_previous[0][i] = 0;
			m_previous[0][i + 1] + 1;
		}
		if (m_previous[1][i] == m_base_2)
		{
			m_previous[1][i] = 0;
			m_previous[1][i + 1] + 1;
		}
	}
	if (m_previous[0][m_decimals - 1] == m_base_1)
		m_previous[0][m_decimals - 1] = 0;
	if (m_previous[1][m_decimals - 1] == m_base_2)
		m_previous[1][m_decimals - 1] = 0;

	r[0] = padic2double(m_previous[0], m_base_1);
	r[1] = padic2double(m_previous[1], m_base_2);
	return r;
}

std::vector<double> MCPROJ::Kakutani2DGauss::operator()()
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;
	std::vector<double> z;
	double u1, u2;
	do
	{
		std::vector<double> temp = m_K2D ->operator()();
		u1 = temp[0];
		u2 = temp[1];
	} while (u1 <= epsilon);

	z.push_back((sqrt(-2.0 * log(u1)) * cos(two_pi * u2))*m_sigma + m_mu);
	z.push_back((sqrt(-2.0 * log(u1)) * sin(two_pi * u2))*m_sigma + m_mu);
	return z;
}

std::vector<double> MCPROJ::Kakutani2DNIG::operator()()
{
	// getting uniform random variable value between 0 and 1
	std::vector<double> x = m_K2D->operator()();
	// getting value of IG distribution

	std::vector<double> y;
	y.push_back(inverseIG(x[0], m_Listofvalues_1));
	y.push_back(inverseIG(x[1], m_Listofvalues_2));

	std::vector<double> lol = m_KG2D -> operator()();

	//returning value of NIG distribution
	std::vector<double> r;
	r.push_back((lol[0] * sqrt(y[0])) + m_mu_1 + (m_beta_1*y[0]));
	r.push_back((lol[1] * sqrt(y[1])) + m_mu_2 + (m_beta_2*y[1]));
	return r;
}