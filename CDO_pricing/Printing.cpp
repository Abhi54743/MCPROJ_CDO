#include "Printing.h"

void MCPROJ::Printer::printOutClosedFormula(double res, std::string typeOfVariable, double timeInSeconds)
{
	std::ofstream cfOutput;
	//*name of the file
	std::string name = "CloseFormula" + typeOfVariable + ".txt";
	const char* file_name = name.c_str();
	//*
	cfOutput.open(file_name);

	cfOutput << res << "\n";
	cfOutput << seconds2Time(timeInSeconds) << "\n";

	cfOutput.close();
}

void MCPROJ::Printer::printOutMonteCarlo(std::vector<double> res, std::string typeOfVariable, int N, double timeInSeconds)
{
	std::ofstream mcOutput;
	//*name of the file
	std::stringstream NN;
	NN << N;
	std::string name = "MonteCarlo" + typeOfVariable +  NN.str() + "Iterations.txt";
	const char* file_name = name.c_str();
	//*
	mcOutput.open(file_name);

	mcOutput << N << "\n";

	mcOutput << res[0] << "\n";
	mcOutput << res[1] << "\n";
	mcOutput << res[2] << "\n";
	
	mcOutput << seconds2Time(timeInSeconds) << "\n";

	mcOutput.close();
}

std::string MCPROJ::Printer::seconds2Time(double timeInSeconds)
{
	int numberSec = (int)floor(timeInSeconds);
	int numberMin = floor(numberSec / 60);

	int sec = numberSec % 60;
	int min = numberMin % 60;
	int hour = floor(numberMin / 60);

	int cent = floor((timeInSeconds - floor(timeInSeconds)) * 100);

	std::stringstream cc;
	std::stringstream SS;
	std::stringstream MM;
	std::stringstream HH;

	cc << cent;
	SS << sec;
	MM << min;
	HH << hour;

	std::string time = HH.str() + "h" + MM.str() + "m" + SS.str() + "s" + cc.str();
	return time;
}
		
	