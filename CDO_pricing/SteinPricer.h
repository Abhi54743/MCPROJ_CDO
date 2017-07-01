/*!
Class implimenting method of Stein: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 10/06/2017

Version 1.0
*/
#define _USE_MATH_DEFINES
#include<cstdio>
#include<cmath>
#include<sstream>
#include<iostream>
#include<vector>
#include "Generators.h"
#include "Monte_Carlo.h"

#pragma once
#ifndef _STEIN_PRICER_H_
#define _STEIN_PRICER_H_

class GaussianSteinPricer {
public:
	// default Constructor
	GaussianSteinPricer() {};

	// default destructor
	~GaussianSteinPricer() {};

	//constructor
	GaussianSteinPricer(std::string dist_type,
		double	q,
		double	corr,
		double	R,
		int		Nb_CDS,
		double	K1,
		double	K2,
		double	alpha = 1.0,
		double	beta = 0.5,
		int		gridSize = 100) :m_dist_type(dist_type), m_q(q), m_corr(corr), m_R(R), m_Nb_CDS(Nb_CDS), m_K1(K1), m_K2(K2), m_alpha(alpha), m_beta(beta), m_gridSize(gridSize) {
		m_Gen.seed(static_cast<unsigned int>(0));//to remove later
		m_normal_variable = new normal_rv(m_Gen, normal_dist(0, 1));

		double gamma = sqrt(alpha*alpha - beta*beta);
		if (abs(beta) > alpha) { throw std::exception(" value of abs(beta) is greater than alpha "); }
		m_NIG_variable = new MCPROJ::NIG_rv(alpha, beta, -beta*gamma*gamma / (alpha*alpha), gamma*gamma*gamma / (alpha*alpha), gridSize);
		
		if (m_dist_type == "Gaussian") {
			boost::math::normal normal;
			m_C = boost::math::quantile(normal, m_q);
		}
		else if (m_dist_type == "NIG") {
			double gamma = sqrt(m_alpha*m_alpha - m_beta*m_beta);
			MCPROJ::NIG_rv NIG_A(m_alpha / m_corr, m_beta / m_corr, -m_beta*gamma*gamma / (m_corr*m_alpha*m_alpha), pow(gamma, 3) / (m_corr*m_alpha*m_alpha), m_gridSize);
			m_C = NIG_A.inverseCDF(m_q);

		}

		else throw std::exception("only Gaussian or NIG are supported");
		double factor = sqrt(1 - (m_corr*m_corr)) / m_corr;
		m_NIG_X = new MCPROJ::NIG_rv(factor*m_alpha, factor*m_beta, -factor*m_beta*gamma*gamma / (m_alpha*m_alpha), factor*gamma*gamma*gamma / (m_alpha*m_alpha), m_gridSize);

	};

	std::vector<double> price(int N);

	double operator()();

	

private:

	double conditionallyExpectedPrice( double u, double K);

	std::string			m_dist_type;	// distribution type of M and X_i
	double				m_q;			// default probability
	double				m_C;
	double				m_corr;			// correlation between CDS
	double				m_R;			// recovery rate
	double				m_alpha;
	double				m_beta;
	int					m_gridSize;
	int					m_Nb_CDS;		// Number of CDS
	double				m_K1;			// lowest default probability of the selected tranche
	double				m_K2;			// highest default probability of the selected tranche
	normal_rv	*		m_normal_variable;
	MCPROJ::NIG_rv	*	m_NIG_variable;
	generator			m_Gen;
	MCPROJ::NIG_rv  *   m_NIG_X;

};


class PoissonSteinPricer {
public:
	// default Constructor
	PoissonSteinPricer() {};

	// default destructor
	~PoissonSteinPricer() {};

	//constructor
	PoissonSteinPricer(std::string dist_type,
		double	q,
		double	corr,
		double	R,
		int		Nb_CDS,
		double	K1,
		double	K2,
		double	alpha = 1.0,
		double	beta = 0.5,
		int		gridSize = 100) :m_dist_type(dist_type), m_q(q), m_corr(corr), m_R(R), m_Nb_CDS(Nb_CDS), m_K1(K1), m_K2(K2), m_alpha(alpha), m_beta(beta), m_gridSize(gridSize) {
		m_Gen.seed(static_cast<unsigned int>(0));//to remove later
		m_normal_variable = new normal_rv(m_Gen, normal_dist(0, 1));

		double gamma = sqrt(alpha*alpha - beta*beta);
		if (abs(beta) > alpha) { throw std::exception(" value of abs(beta) is greater than alpha "); }
		m_NIG_variable = new MCPROJ::NIG_rv(alpha, beta, -beta*gamma*gamma / (alpha*alpha), gamma*gamma*gamma / (alpha*alpha), gridSize);

		if (m_dist_type == "Gaussian") {
			boost::math::normal normal;
			m_C = boost::math::quantile(normal, m_q);
		}
		else if (m_dist_type == "NIG") {
			double gamma = sqrt(m_alpha*m_alpha - m_beta*m_beta);
			MCPROJ::NIG_rv NIG_A(m_alpha / m_corr, m_beta / m_corr, -m_beta*gamma*gamma / (m_corr*m_alpha*m_alpha), pow(gamma, 3) / (m_corr*m_alpha*m_alpha), m_gridSize);
			m_C = NIG_A.inverseCDF(m_q);

		}

		else throw std::exception("only Gaussian or NIG are supported");
		double factor = sqrt(1 - (m_corr*m_corr)) / m_corr;
		m_NIG_X = new MCPROJ::NIG_rv(factor*m_alpha, factor*m_beta, -factor*m_beta*gamma*gamma / (m_alpha*m_alpha), factor*gamma*gamma*gamma / (m_alpha*m_alpha), m_gridSize);

	};

	std::vector<double> price(int N);

	double operator()();



private:

	double conditionallyExpectedPrice(double u, double K);

	std::string			m_dist_type;	// distribution type of M and X_i
	double				m_q;			// default probability
	double				m_C;
	double				m_corr;			// correlation between CDS
	double				m_R;			// recovery rate
	double				m_alpha;
	double				m_beta;
	int					m_gridSize;
	int					m_Nb_CDS;		// Number of CDS
	double				m_K1;			// lowest default probability of the selected tranche
	double				m_K2;			// highest default probability of the selected tranche
	normal_rv	*		m_normal_variable;
	MCPROJ::NIG_rv	*	m_NIG_variable;
	generator			m_Gen;
	MCPROJ::NIG_rv  *   m_NIG_X;

};


#endif // !_STEIN_PRICER_H_
