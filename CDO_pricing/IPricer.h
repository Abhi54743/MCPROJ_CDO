/*!
Pricer Interface: CDO Tranche Pricing, Monte Carlo Project

Author: Abhishek MUKHOPADHYAY/ Carlo Pulcini

Date: 11/03/2017

Version 1.0
*/

#include <iostream>

#include<iomanip>
#include "Generators.h"

#pragma once
#ifndef _IPricer_H
#define _IPricer_H

namespace MCPROJ{
	class IPricer {
	public:

		//default constructor
		IPricer() {};

		//virtual destructor
		virtual ~IPricer() {};

		//virtual methods
		
		virtual double operator()()=0;

		virtual double expected_Loss()=0;

		virtual std::vector<double> expected_LossMC(int N) = 0;
		

	};

}//end namepace MCPROJ

#endif // !_IPricer_H
