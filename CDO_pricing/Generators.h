/*!
Tools used in porject: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 11/03/2017

Version 1.0
*/


#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <memory>

typedef boost::mt19937 generator;
typedef boost::normal_distribution<double> normal_dist;
typedef boost::variate_generator< generator&, normal_dist> normal_rv;
typedef boost::math::inverse_gaussian_distribution<double> inverse_gaussian_dist;

#pragma once
#ifndef _GENERATORS_H
#define _GENERATORS_H

namespace MCPROJ {

	class NIG_rv {
		
		// members

		double						m_mu, m_alpha, m_beta, m_delta;
		generator					m_Gen;
		normal_rv			*		m_Gaussian;
		inverse_gaussian_dist		m_IG;

	public:

		// default constructor

		NIG_rv() {};

		// constructor

		NIG_rv( double alpha, double beta, double mu, double delta);
		
		// destructor

		~ NIG_rv() {};

		// operator

		double operator()();

	};

}//end namespace MCPROJ

#endif // !_GENERATORS_H


