#include "Generators.h"
#include "Monte_Carlo.h"
#include "MCPROJ_Tools.h"

#pragma once
#ifndef _NEWTON_H
#define _NEWTON_H

namespace MCPROJ {

	class NewtonSolver {
	public:
		// default Constructor
		NewtonSolver() {};

		// default destructor
		~NewtonSolver() {};

		//constructor
		NewtonSolver(int interMax, double tol);

		double Solve(std::function<double(double)>, double firstGuess);

	private:
		int				m_iterMax;			
		double			m_tol;

	};


}// end namespace MCPROJ


#endif // !_NEWTON_H


