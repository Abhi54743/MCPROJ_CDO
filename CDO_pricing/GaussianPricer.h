/*!
Gaussian Pricer: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo PULCINI

Date: 11/03/2017

Version 1.0
*/


#include "Generators.h"
#include "IPricer.h"
#include "Monte_Carlo.h"
#include "MCPROJ_Tools.h"
#include "QMCGenerators.h"

#pragma once
#ifndef _GAUSSIAN_PRICER_H
#define _GAUSSIAN_PRICER_H

namespace MCPROJ{

	class GaussianPricer : public IPricer {
	public:
		// default Constructor
		GaussianPricer() {};

		// default destructor
		~GaussianPricer() {};

		//constructor
		GaussianPricer(double	q,
					   double	corr,
					   double	R,
					   int		Nb_CDS,
					   double	K1,
					   double	K2);

		double operator()();

		double expected_Loss();

		std::vector<double> expected_LossMC(int N) ;

		double derivativeOfVarTranches(int N, double theta);

		double importance_sampling_tranches(double theta);
		
		double derivativeOfVarCommon(int N, double theta);

		double importance_sampling_common(double theta);

		std::vector<double> expected_LossMCVR(int N, double thetaCommon, double thetaTranches);

		double expected_LossQMC(int N, std::string Qtype);

		double percentage_defaultQMCH(double C, double corr, double R, int Nb_CDS, double K1, double	K2);

		double percentage_defaultQMCK(double C, double corr, double R, int Nb_CDS, double K1, double K2);

	private:
		double				m_q;			// default probability
		double				m_C;
		double				m_corr;			// correlation between CDS
		double				m_R;			// recovery rate
		int					m_Nb_CDS;		// Number of CDS
		double				m_K1;			// lowest default probability of the selected tranche
		double				m_K2;			// highest default probability of the selected tranche
		generator			m_gen;			// uniform generator
		normal_rv	*		m_G;			// normal random variable generator pointer
		Halton2DGauss	*	m_HaltonGauss;		// Halton QMC no. generator
		Kakutani2DGauss	*	m_KakutaniGauss;		// Halton QMC no. generator
		
	
	
	};


}// end namespace MCPROJ


#endif // !_GAUSSIAN_PRICER_H