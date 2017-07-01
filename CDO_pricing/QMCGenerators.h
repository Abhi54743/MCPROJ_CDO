#pragma once
#ifndef _QMCGENERATOR_H_
#define _QMCGENERATOR_H_

#include "Generators.h"

namespace MCPROJ {

	class Halton1D {

		size_t		m_index;
		size_t		m_base;

	public:
		//constructor
		Halton1D(size_t base):m_base(base), m_index(1) {};
		
		double operator()();
		
	};

	class Halton2D {

		size_t		m_index;
		size_t		m_base_1;
		size_t		m_base_2;

	public:
		//constructor
		Halton2D(size_t base_1, size_t base_2) :m_base_1(base_1), m_base_2(base_2), m_index(1) {};

		std::vector<double> operator()();

	};

	class Halton2DGauss {

		size_t			m_index;
		Halton2D	*	m_H2D;
		double			m_mu;
		double			m_sigma;

	public:
		//constructor
		Halton2DGauss(double mu = 0.0, double sigma = 1.0, size_t base_1 = 2, size_t base_2 = 3) :m_index(1), m_mu(mu), m_sigma(sigma) {
			m_H2D = new Halton2D(base_1, base_2);
		};

		std::vector<double> operator()();

	};

	class Halton2DNIG {

		size_t			m_index;
		Halton2D	*	m_H2D;
		Halton2DGauss * m_HG2D;
		double			m_alpha_1, m_alpha_2;
		double			m_beta_1, m_beta_2;
		double			m_mu_1, m_mu_2;
		double			m_delta_1, m_delta_2;
		int				m_gridSize;
		std::vector<double> m_Listofvalues_1;
		std::vector<double> m_Listofvalues_2;
		inverse_gaussian_dist		m_IG_1;
		inverse_gaussian_dist		m_IG_2;

	public:
		//constructor
		Halton2DNIG(double			alpha_1, double alpha_2,
					double			beta_1, double beta_2,
					double			mu_1, double mu_2,
					double			delta_1, double delta_2, int grid_size, size_t base_1 = 2, size_t base_2 = 3, size_t base_3 = 5, size_t base_4 = 7)
			:m_index(1), m_alpha_1(alpha_1), m_alpha_2(alpha_2), m_beta_1(beta_1), m_beta_2(beta_2),
			m_mu_1(mu_1), m_mu_2(mu_2), 	m_delta_1(delta_1), m_delta_2(delta_2), m_gridSize(grid_size)
		{	
			if (abs(beta_1) > alpha_1 || abs(beta_2) > alpha_2) { throw std::exception(" value of abs(beta) is greater than alpha "); }

			else if (delta_1 < 1e-30|| delta_2 < 1e-30) { throw std::exception(" value of delta is negative "); }
			m_H2D = new Halton2D(base_1, base_2);
			m_HG2D = new Halton2DGauss(0.0, 1.0, base_3, base_4);

			double gamma_1 = sqrt(alpha_1*alpha_1 - beta_1*beta_1);
			double gamma_2 = sqrt(alpha_2*alpha_2 - beta_2*beta_2);

			m_IG_1 = inverse_gaussian_dist(m_delta_1*gamma_1, gamma_1*gamma_1);
			m_IG_2 = inverse_gaussian_dist(m_delta_2*gamma_2, gamma_2*gamma_2);

			for (int i = 0; i<m_gridSize; i++)
			{
				m_Listofvalues_1.push_back(boost::math::quantile(m_IG_1, ((double)i) / m_gridSize));
				m_Listofvalues_2.push_back(boost::math::quantile(m_IG_2, ((double)i) / m_gridSize));
			}
		};

		std::vector<double> operator()();

	};

}// end of namespace MCPROJ



#endif // !_QMCGENERATOR_H_
