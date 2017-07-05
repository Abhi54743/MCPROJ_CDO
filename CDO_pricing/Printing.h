#include "Generators.h"

#pragma once
#ifndef _PRINTING_H
#define _PRINTING_H

namespace MCPROJ {

	class Printer 
	{
	
	public:
		void printOutClosedFormula(double res, std::string typeOfVariable, double timeInSeconds);
		void printOutMonteCarlo(std::vector<double> res, std::string typeOfVariable, int N, double timeInSeconds);
	
		std::string seconds2Time(double timeInSeconds);

	};

	
		
		
}
#endif