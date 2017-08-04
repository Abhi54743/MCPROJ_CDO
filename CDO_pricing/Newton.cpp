#include "Newton.h"


MCPROJ::NewtonSolver::NewtonSolver(int iterMax, double tol) : m_iterMax(iterMax), m_tol(tol) {};

double MCPROJ::NewtonSolver::Solve(std::function<double(double)> f, double firstGuess)
{
	double eps = 0.001;
	std::function<double(double)> fprime = [f, eps](double theta) {return (f(theta + eps) - f(theta - eps)) / (2*eps); };

	double fx, fx1;
	double x1 = firstGuess;
	double x = x1 + 2 * m_tol;

	int i = 0;
	while (i < m_iterMax && std::fabs(x1 - x) > m_tol)
	{
		x = x1;
		fx = f(x);
		fx1 = fprime(x);
		x1 = x - (fx / fx1);
		i++;
	}

	if (i == m_iterMax)
	{
		throw std::exception("Newton Max Iteration");
	}

	return x1;

};
