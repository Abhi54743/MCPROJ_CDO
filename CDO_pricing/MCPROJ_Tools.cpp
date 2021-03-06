#define _USE_MATH_DEFINES


#include "MCPROJ_Tools.h"

static double MCPROJ::f(double x, double y, double aprime, double bprime, double rho)
{
	double r = aprime * (2 * x - aprime) + bprime * (2 * y - bprime) + 2 * rho * (x - aprime) * (y - bprime);
	return exp(r);
}


double MCPROJ::NormPDF(double x, double mu, double sigma) {
	double xn = (x - mu) / sigma;
	return exp(-0.5 * xn * xn) / (sqrt(2 * M_PI) * sigma);
}

double MCPROJ::NormPDF(double x) {
	return NormPDF(x, 0., 1.);
}

double MCPROJ::NormalCDF(double z)
{
	if (z > 6.0) return 1;
	if (z < -6.0) return 0;

	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;

	double a = abs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-1 * z) * (z / 2.0));
	double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
}


/// <summary>
/// This is the cumulative bivariate normal distribution.
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <param name="rho"></param>
/// <returns></returns>
double MCPROJ::BivariateNormalCDF(double a, double b, double rho) {
	if ((a <= 0.0) && (b <= 0.0) && (rho <= 0.0))
	{
		double aprime = a / sqrt(2.0 * (1.0 - rho * rho));
		double bprime = b / sqrt(2.0 * (1 - rho * rho));
		double A[4] = { 0.3253030, 0.4211071, 0.1334425, 0.006374323 };
		double B[4] = { 0.1337764, 0.6243247, 1.3425378, 2.2626645 };

		double sum = 0.0;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				sum = sum + A[i] * A[j] * f(B[i], B[j], aprime, bprime, rho);
			}
		}
		sum = sum * (sqrt(1.0 - rho * rho) / M_PI);
		return sum;
	}
	else if (a * b * rho <= 0.0)
	{
		try {
			if ((a <= 0.0) && (b >= 0.0) && (rho >= 0.0))
				return NormalCDF(a) - BivariateNormalCDF(a, -1.0 * b, -1.0 * rho);
			else if ((a >= 0.0) && (b <= 0.0) && (rho >= 0.0))
				return NormalCDF(b) - BivariateNormalCDF(-1.0 * a, b, -1.0 * rho);
			else if ((a >= 0.0) && (b >= 0.0) && (rho <= 0.0))
				return NormalCDF(a) + NormalCDF(b) - 1.0 + BivariateNormalCDF(-1 * a, -1 * b, rho);
			else throw;
		}
		catch (double a) { std::cout << "This part of the code should never be reached." << std::endl; };

	}
	else if (a * b * rho >= 0.0)
	{
		double denum = sqrt(a * a - 2 * rho * a * b + b * b);
		double rho1 = ((rho * a - b) * copysign(1, a)) / denum;
		double rho2 = ((rho * b - a) * copysign(1, b)) / denum;
		double delta = (1.0 - copysign(1, a) * copysign(1, b)) / 4.0;
		return BivariateNormalCDF(a, 0.0, rho1) +
			BivariateNormalCDF(b, 0.0, rho2) - delta;
	}
	else
	{
		throw std::exception("This part of the code should never be reached: Bivariate Normal CDF");

	}
}




double MCPROJ::RationalApproximation(double t)
{
	// Abramowitz and Stegun formula 26.2.23.
	// The absolute value of the error should be less than 4.5 e-4.
	double c[] = { 2.515517, 0.802853, 0.010328 };
	double d[] = { 1.432788, 0.189269, 0.001308 };
	return t - ((c[2] * t + c[1])*t + c[0]) /
		(((d[2] * t + d[1])*t + d[0])*t + 1.0);
}

double MCPROJ::normal_CDF_inverse(double p)
{
	if (p <= 0.0 || p >= 1.0)
	{
		std::stringstream os;
		os << "Invalid input argument (" << p
			<< "); must be larger than 0 but less than 1.";
		throw std::invalid_argument(os.str());
	}


	if (p < 0.5)
	{
		// F^-1(p) = - G^-1(p)
		return -RationalApproximation(sqrt(-2.0*log(p)));
	}
	else
	{
		// F^-1(p) = G^-1(1-p)
		return RationalApproximation(sqrt(-2.0*log(1 - p)));
	}
}


void MCPROJ::gauleg(double x1, double x2, std::vector<double> & abs, std::vector<double> & weight) {

	const double EPS = 1.0e-14;			//Relative precision
	size_t n = abs.size(), m = (n+1)/2;

	double xm = 0.5*(x2 + x1);
	double xl = 0.5*(x2 - x1);

	double p1(0.0), p2(0.0), p3(0.0), pp(0.0);
	double z1(0.0), z(0.0);

	for (size_t i = 0; i < m; i++) {
		z = cos(M_PI*(i + 0.75)/(n+0.5));			//starting approximation for root

		do {
			p1 = 1.0;
			p2 = 0.0;

			for (size_t j = 0; j < n; j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0*j+1.0)*z*p2 -j*p3)/ (j+1);
			}

			//p1 is now the required legendre polynomial at z.
			//pp, its derivative, may be found via a standard relationship.

			pp = n*(z*p1 - p2) / (z*z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;					// Newton's method

		} while (std::fabs(z - z1) > EPS);

		abs[i] = xm - xl*z;
		abs[n - 1 - i] = xm + xl*z;
		weight[i] = 2.0*xl / ((1.0 - z*z)*pp*pp);
		weight[n - 1 - i] = weight[i];

	}

}

double MCPROJ::inverseIG(double x, std::vector<double> listofvalues)
{
	if (x < 0 || x > 1) { throw std::exception("wrong input value for inversion"); }
	return listofvalues[floor(x*listofvalues.size())];
}

