/*!
Tools used in porject: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 11/03/2017

Version 1.0
*/


#pragma once
#ifndef _MCPROJ_TOOLS_H
#define _MCPROJ_TOOLS_H

namespace MCPROJ {

	static double f(double x, double y, double aprime, double bprime, double rho);

	double NormPDF(double x, double mu, double sigma);

	double NormPDF(double x);

	double NormalCDF(double z);

	double BivariateNormalCDF(double a, double b, double rho);

	double RationalApproximation(double t);

	double normal_CDF_inverse(double p);





}//end namespace MCPROJ

#endif // !_MCPROJ_TOOLS_H