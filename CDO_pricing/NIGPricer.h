/*!
Gaussian Pricer: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 11/03/2017

Version 1.0
*/


#include "IPricer.h"
#include "Monte_Carlo.h"
#include "MCPROJ_Tools.h"

#pragma once
#ifndef _GAUSSIAN_PRICER_H
#define _GAUSSIAN_PRICER_H

namespace MCPROJ {

	class NIGPricer : public IPricer {
	public:
		// default Constructor
		NIGPricer() {};

		// default destructor
		~NIGPricer() {};

		//constructor
		NIGPricer(double	q,
				  double	corr,
				  double	R,
				  int		Nb_CDS,
				  double	K1,
				  double	K2,
				  double alpha,
				  double beta);

		double operator()();

		double expected_Loss();

		std::vector<double> expected_LossMC(int N);

	private:
		double				m_q;			// default probability
		double				m_corr;			// correlation between CDS
		double				m_R;			// recovery rate
		int					m_Nb_CDS;		// Number of CDS
		double				m_K1;			// lowest default probability of the selected tranche
		double				m_K2;			// highest default probability of the selected tranche
		double				m_alpha;		// NIG parameter
		double				m_beta;			// NIG parameter
		NIG_rv				m_NIG_X;		// normal inverse gaussian random variable generator for X
		NIG_rv				m_NIG_M;		// normal inverse gaussian random variable generator for M


	};


}// end namespace MCPROJ


#endif // !_GAUSSIAN_PRICER_H


