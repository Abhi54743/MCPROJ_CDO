#include "Printing.h"

void MCPROJ::Printer::printOutClosedFormula(double res, std::string typeOfVariable, double q,
	double corr,
	double R,
	int Nb_CDS, double timeInSeconds)
{
	std::ofstream cfOutput;
	//*name of the file
	std::stringstream qq;
	qq << q;
	std::stringstream ccorr;
	ccorr << corr;
	std::stringstream RR;
	RR << R;
	std::stringstream NNb_CDS;
	NNb_CDS << Nb_CDS;
	std::string name = "zzCloseFormula" + typeOfVariable + +"q" + qq.str() + "corr" + ccorr.str() + "R" + RR.str() + "NbCDS" + NNb_CDS.str() +  ".txt";
	const char* file_name = name.c_str();
	//*
	cfOutput.open(file_name);

	cfOutput << res << "\n";
	cfOutput << seconds2Time(timeInSeconds) << "\n";

	cfOutput.close();
}

void MCPROJ::Printer::printOutMonteCarlo(std::vector<double> res, 
										 std::string typeOfVariable, 
										 double q,
										 double corr,
										 double R,
										 int Nb_CDS,	
										 int N, double timeInSeconds)
{
	std::ofstream mcOutput;
	//*name of the file
	std::stringstream qq;
	qq << q;
	std::stringstream ccorr;
	ccorr << corr;
	std::stringstream RR;
	RR << R;
	std::stringstream NNb_CDS;
	NNb_CDS << Nb_CDS;
	std::stringstream NN;
	NN << N;
	std::string name = "zzMonteCarlo" + typeOfVariable +  "q" + qq.str() + "corr" + ccorr.str() + "R" + RR.str() + "NbCDS" + NNb_CDS.str() + "mc" + NN.str() + "Iterations.txt";
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

void MCPROJ::Printer::printOutQuasiMonteCarlo(double res, std::string typeOfVariable, std::string typeOfPseudo, double q, double corr, double R, int Nb_CDS, int N, double timeInSeconds)
{
	std::ofstream mcOutput;
	//*name of the file
	std::stringstream qq;
	qq << q;
	std::stringstream ccorr;
	ccorr << corr;
	std::stringstream RR;
	RR << R;
	std::stringstream NNb_CDS;
	NNb_CDS << Nb_CDS;
	std::stringstream NN;
	NN << N;
	std::string name = "zzQuasiMonteCarlo" + typeOfVariable + typeOfPseudo +  "q" + qq.str() + "corr" + ccorr.str() + "R" + RR.str() + "NbCDS" + NNb_CDS.str() + "qmc" + NN.str() + "Iterations.txt";
	const char* file_name = name.c_str();
	//*
	mcOutput.open(file_name);

	mcOutput << N << "\n";

	mcOutput << res << "\n";

	mcOutput << seconds2Time(timeInSeconds) << "\n";

	mcOutput.close();
}

void MCPROJ::Printer::printOutStein(std::vector<double> res, std::string typeOfApprox, std::string typeOfVariable, double q, double corr, double R, int Nb_CDS, int N, double timeInSeconds)
{
	std::ofstream mcOutput;
	//*name of the file
	std::stringstream qq;
	qq << q;
	std::stringstream ccorr;
	ccorr << corr;
	std::stringstream RR;
	RR << R;
	std::stringstream NNb_CDS;
	NNb_CDS << Nb_CDS;
	std::stringstream NN;
	NN << N;
	std::string name = "zzStein" + typeOfApprox + typeOfVariable + "q" + qq.str() + "corr" + ccorr.str() + "R" + RR.str() + "NbCDS" + NNb_CDS.str() + "qmc" + NN.str() + "Iterations.txt";
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

	int microsec = (int)floor((timeInSeconds - numberSec)*1000);
	int sec = numberSec % 60;
	int min = numberMin % 60;
	int hour = floor(numberMin / 60);

	std::stringstream mm;
	std::stringstream SS;
	std::stringstream MM;
	std::stringstream HH;

	mm << microsec;
	SS << sec;
	MM << min;
	HH << hour;

	std::string time = HH.str() + "h" + MM.str() + "m" + SS.str() + "s" + mm.str();
	return time;
}
		
	