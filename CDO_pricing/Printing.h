#include "Generators.h"

#pragma once
#ifndef _PRINTING_H
#define _PRINTING_H

namespace MCPROJ {

	class Printer 
	{
	
	public:
		void printOutClosedFormula(double res, std::string typeOfVariable, double q, double corr, double R, int Nb_CDS, double timeInSeconds);
		void printOutMonteCarlo(std::vector<double> res, std::string typeOfVariable, double q, double corr,	double R, int Nb_CDS, int N, double timeInSeconds);
		void printOutQuasiMonteCarlo(double res, std::string typeOfVariable, std::string typeOfPseudo, double q, double corr, double R, int Nb_CDS, int N, double timeInSeconds);

		void printOutStein(std::vector<double> res, std::string typeOfApprox, std::string typeOfVariable, double q, double corr, double R, int Nb_CDS, int N, double timeInSeconds);

		std::string seconds2Time(double timeInSeconds);

	};

	
		
		
}
#endif