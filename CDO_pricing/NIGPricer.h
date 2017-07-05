/*!
Gaussian Pricer: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 11/03/2017

Version 1.0
*/


#include "IPricer.h"
#include "Monte_Carlo.h"
#include "MCPROJ_Tools.h"
#include "QMCGenerators.h"

#pragma once
#ifndef _NIG_PRICER_H
#define _NIG_PRICER_H

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
				  double	alpha,
				  double	beta,
				  int		gridSize);

		double operator()();

		double expected_Loss();

		std::vector<double> expected_LossMC(int N);

		double expected_LossQMC(int N, std::string Qtype);

		double percentage_defaultQMCH(double C, double corr, double R, int Nb_CDS, double K1, double K2);

		double percentage_defaultQMCK(double C, double corr, double R, int Nb_CDS, double K1, double K2);

	private:
		double				m_q;			// default probability
		double			    m_C;
		double				m_corr;			// correlation between CDS
		double				m_R;			// recovery rate
		int					m_Nb_CDS;		// Number of CDS
		double				m_K1;			// lowest default probability of the selected tranche
		double				m_K2;			// highest default probability of the selected tranche
		double				m_alpha;		// NIG parameter
		double				m_beta;			// NIG parameter
		NIG_rv		*		m_NIG_X;		// normal inverse gaussian random variable generator for X
		NIG_rv		*		m_NIG_M;		// normal inverse gaussian random variable generator for M
		Halton2DNIG *	    m_Halton2DNIG;  // 2D NIG QMC 
		Kakutani2DNIG  *       m_Kakutani2DNIG;// 2D NIG QMC 


	};


}// end namespace MCPROJ


#endif // !_NIG_PRICER_H


