#define _USE_MATH_DEFINES

#include "SteinPricer.h"
#include "GaussianPricer.h"
#include "NIGPricer.h"
#include "zero_search.h"
#include "QMCGenerators.h"
#include "Printing.h"
#include "Newton.h"



#include <fstream>
#include <iomanip>
#include <iostream>
#include <boost/algorithm/string.hpp>
using namespace std;

using namespace MCPROJ;

//This is the main

struct initVariables
{
	initVariables(double q, double corr, double R, int Nb_CDS, double K1, double K2, double alpha, double beta, int gridSize)
		: q(q), corr(corr), R(R), Nb_CDS(Nb_CDS), K1(K1), K2(K2), alpha(alpha), beta(beta), gridSize(gridSize) {}

	double q;	//default probability
	double corr;//correlation between common factor and independent factors
	double R;	//recovery rate
	int Nb_CDS;	
	double K1;	//tranche between K1 and K2
	double K2;
	double alpha;//alpha and beta are NIG parameters
	double beta;
	int gridSize;//for NIG
};


int main() {
	/*
	
	*/
	/***************************** STEIN METHOD********************************************/
	/*
	PoissonSteinPricer PoissonStein(std::string("Gaussian"), q, corr, R, Nb_CDS, K1, K2);

	std::vector<double> vect9 = PoissonStein.price(10000);

	std::cout << vect9[0] << "    " << vect9[1] << "    " << vect9[2] << std::endl;

	PoissonSteinPricer PoissonStein2(std::string("NIG"), q, corr, R, Nb_CDS, K1, K2, alpha, beta, gridSize);

	std::vector<double> vect10 = PoissonStein2.price(100);

	std::cout << vect10[0] << "    " << vect10[1] << "    " << vect10[2] << std::endl;

	GaussianSteinPricer GaussStein(std::string("Gaussian"), q, corr, R, Nb_CDS, K1, K2);

	std::vector<double> vect3 = GaussStein.price(10000);

	std::cout << vect3[0] << "    " << vect3[1] << "    " << vect3[2] << std::endl;

	GaussianSteinPricer GaussStein2(std::string("NIG"), q, corr, R, Nb_CDS, K1, K2, alpha, beta, gridSize);

	std::vector<double> vect4 = GaussStein2.price(1000);

	std::cout << vect4[0] << "    " << vect4[1] << "    " << vect4[2] << std::endl;

	*/
	/**********************************************************************************************/

	/**********************************************************************************************/
	/*
	GaussianPricer Gauss(q, corr, R, Nb_CDS, K1, K2);

	std::cout << "Gauss Closed" << std::endl;
	std::cout << Gauss.expected_Loss() << std::endl ;

	std::vector<double> vect = Gauss.expected_LossMC(10000);
	std::cout << "Gauss MonteCarlo" << std::endl;
	std::cout << vect[0] << "    " << vect[1] << "    " << vect[2] << std::endl;

	double rofl = Gauss.expected_LossQMC(100, "Halton");
	std::cout << "Gauss QMC Halton" << std::endl;
	std::cout << rofl << std::endl;

	double rofl1 = Gauss.expected_LossQMC(100, "Kakutani");
	std::cout << "Gauss QMC Kakutani" << std::endl;
	std::cout << rofl1 << std::endl;
	*/

	/*Variance Reduction - Importance Sampling*/

	/*
	int iterMax = 100;
	double tol = 1.e-4;
	int NN = 1.e3;

	NewtonSolver newton(iterMax, tol);

	std::function<double(double)> f1 = [&Gauss, NN](double theta) {return Gauss.derivativeOfVarCommon(NN, theta); };

	double thetaCommon = newton.Solve(f1, 0.0);

	std::function<double(double)> f2 = [&Gauss, NN](double theta) {return Gauss.derivativeOfVarTranches(NN, theta); };

	double thetaTranches = newton.Solve(f2, 0.0);

	std::vector<double> vectVR = Gauss.expected_LossMCVR(10000, thetaCommon, thetaTranches);
	std::cout << "Gauss Variance Reduction" << std::endl;
	std::cout << vectVR[0] << "    " << vectVR[1] << "    " << vectVR[2] << std::endl;


	NIGPricer NIG(q, corr, R, Nb_CDS, K1, K2, alpha, beta, gridSize);

	std::vector<double> vect2 = NIG.expected_LossMC(100);
	std::cout << "NIG MonteCarlo" << std::endl;
	std::cout << vect2[0] << "    " << vect2[1] << "    " << vect2[2] << std::endl;

	double vect_667 = NIG.expected_LossQMC(100, "Halton");
	std::cout << "NIG QMC Halton" << std::endl;
	std::cout << vect_667 << std::endl;

	double vect_668 = NIG.expected_LossQMC(100, "Kakutani");
	std::cout << "NIG QMC Kakutani" << std::endl;
	std::cout << vect_668 << std::endl;
	*/

	/**************************************/




	/**************************************A SORT OF FINAL MAIN**************************************/

	

	Printer Printer;
	//PART 1: CLOSED FORMULA AND MONTECARLO
	//1.A: GAUSSIAN COPULA
	//1.A.I: GUASSIAN COPULA, R = 0

	initVariables init(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	GaussianPricer Gauss(init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);
	
	//close formula
	clock_t begin = clock();

	double res = Gauss.expected_Loss();
	
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "Gauss Closed without Recovery Rate" << std::endl;
	std::cout << res << "\n" << std::endl;

	Printer.printOutClosedFormula(res, "Gaussian", init.q, init.corr, init.R, init.Nb_CDS, elapsed_secs);

	//montecarlo
	std::vector<int> mcIterations;

	int start = 10000;
	int endmc = 20000;
	int step = 10000;

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);
	
	std::cout << "Gauss MonteCarlo without Recovery Rate is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		std::vector<double> res = Gauss.expected_LossMC(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutMonteCarlo(res, "Gaussian", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//1.A.II: GUASSIAN COPULA, R = 0.2
	
	init = initVariables(0.5, 0.3, 0.2, 100, 0.2, 0.7, 1, 0.5, 1000);

	Gauss = GaussianPricer(init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);

	//close formula
	begin = clock();

	res = Gauss.expected_Loss();

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	std::cout << "Gauss Closed with Recovery Rate" << std::endl;
	std::cout << res << "\n" << std::endl;

	Printer.printOutClosedFormula(res, "Gaussian", init.q, init.corr, init.R, init.Nb_CDS, elapsed_secs);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Gauss MonteCarlo with Recovery Rate is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		std::vector<double> res = Gauss.expected_LossMC(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutMonteCarlo(res, "Gaussian", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}
	

	//1.A.III: GUASSIAN COPULA, R = 0 and R!=0, VARIANCE REDUCTION
	//...
	//...aaaaaaaaaaaaaaaaaaaaa
	//...

	//1.B: NIG COPULA
	//1.B.I: NIG COPULA, R = 0

	init = initVariables(0.5, 0.3, 0.0, 100, 0.2, 0.7, 1, 0.5, 1000);

	NIGPricer NIG(init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2, init.alpha, init.beta, init.gridSize);

	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "NIG MonteCarlo without Recovery Rate is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		std::vector<double> res = NIG.expected_LossMC(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutMonteCarlo(res, "NIG", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//1.B.I: NIG COPULA, R = 0.2

	init = initVariables(0.5, 0.3, 0.2, 100, 0.2, 0.7, 1, 0.5, 1000);

	NIG = NIGPricer(init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2, init.alpha, init.beta, init.gridSize);

	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "NIG MonteCarlo with Recovery Rate is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		std::vector<double> res = NIG.expected_LossMC(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutMonteCarlo(res, "NIG", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//1.C: Other parameters effect
	//...
	//...
	//...

	//PART 2 Quasi Monte Carlo

	//To decide if modify in order to generate psuedo-number BEFORE pricing or DURING pricing (now it is during pricing,
	//then QMC is not faster than MC, it is even slower)

	//2.A: QMC for GAUSSIAN COPULA
	//2.A.I: GUASSIAN COPULA, Halton

	init = initVariables(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	Gauss = GaussianPricer(init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Halton Gauss Quasi MonteCarlo is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		double res = Gauss.expected_LossQMC(mcIterations[j], "Halton");

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutQuasiMonteCarlo(res, "Gaussian", "Halton", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//2.A.II: GUASSIAN COPULA, Kakutani
	
	std::cout << "Kakutani Gauss Quasi MonteCarlo is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		double res = Gauss.expected_LossQMC(mcIterations[j], "Kakutani");

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutQuasiMonteCarlo(res, "Gaussian", "Kakutani", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//2.B: QMC for NIG COPULA
	//2.B.I NIG COPULA, Halton

	init = initVariables(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	NIG = NIGPricer(init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2, init.alpha, init.beta, init.gridSize);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Halton NIG Quasi MonteCarlo is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		double res = NIG.expected_LossQMC(mcIterations[j], "Halton");

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutQuasiMonteCarlo(res, "NIG", "Halton", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	std::cout << "Kakutani NIG Quasi MonteCarlo is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		double res = NIG.expected_LossQMC(mcIterations[j], "Kakutani");

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutQuasiMonteCarlo(res, "NIG", "Kakutani", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//PART 3 Stein Method
	
	//3.A: Stein for GAUSSIAN COPULA
	//3.A.I: GUASSIAN COPULA, Gaussian Stein

	init = initVariables(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	GaussianSteinPricer GaussStein(std::string("Gaussian"), init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Stein with Gaussian approx for Gauss copula is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		vector<double> res = GaussStein.price(mcIterations[j]);
		//std::vector<double> res2 = Stein.price(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutStein(res, "GaussApprox", "Gauss", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//3.A.II: GUASSIAN COPULA, Poisson Stein

	init = initVariables(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	PoissonSteinPricer PoissonStein(std::string("Gaussian"), init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Stein with Poisson approx for Gauss copula is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		vector<double> res = PoissonStein.price(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutStein(res, "PoissApprox", "Gauss", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//3.B: Stein for NIG COPULA
	//3.B.I: NIG COPULA, Gaussian Stein

	init = initVariables(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	GaussStein = GaussianSteinPricer(std::string("NIG"), init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Stein with Gaussian approx for NIG copula is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		vector<double> res = GaussStein.price(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutStein(res, "GaussApprox", "NIG", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}

	//3.B.II: NIG COPULA, Poisson Stein

	init = initVariables(0.5, 0.3, 0, 100, 0.2, 0.7, 1, 0.5, 1000);

	PoissonStein = PoissonSteinPricer(std::string("NIG"), init.q, init.corr, init.R, init.Nb_CDS, init.K1, init.K2);

	//montecarlo
	mcIterations.clear();

	for (int i = start; i < endmc + step; i += step)
		mcIterations.push_back(i);

	std::cout << "Stein with Poisson approx for NIG coupla is working, see printed txt" << "\n" << std::endl;
	for (int j = 0; j < mcIterations.size(); j++)
	{
		clock_t begin = clock();

		vector<double> res = PoissonStein.price(mcIterations[j]);

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

		Printer.printOutStein(res, "PoissApprox", "NIG", init.q, init.corr, init.R, init.Nb_CDS, mcIterations[j], elapsed_secs);
	}



	int lol;
	std::cin >> lol;

	return 0;

	

	/*************************************************************************************************************/






}