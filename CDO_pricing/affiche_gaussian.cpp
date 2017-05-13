#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <chrono>
#include <cmath>
#include <random>

#include "zero_search.h"
#include "ziggurat.h"

using namespace std;
using namespace ziggurat;

struct gaussian_density {
    double operator()(double x) const { return density(x); }
    static double constant() { return _cst; }
    static double density(double x) { 
        return std::exp(-0.5*x*x) / _cst; 
    }
    static double derivative(double x) {
        return -x * density(x);
    }
    static std::pair<double, double> inverse(double y) {
        if (std::fabs(y-max()) < 1e-12) return std::pair<double, double>(1,1);
        double x = std::sqrt(-2 * std::log(_cst * y));
        return std::pair<double, double>(-x, x);
    }
    static double mode() { return 0; }
    static double max() { return 1. / _cst; }
    private:
        static double _cst;
};
// initialisation du champ static _cst
double gaussian_density::_cst = std::sqrt(2*M_PI); 

template <typename TRealUniform = std::uniform_real_distribution<double>>
struct gaussian_queue {
    gaussian_queue(double l = -1, double r = 1)
        : _left(l), _right(l), _F_left(cdf(l)), _F_right(cdf(r)) {} 
    static double cdf(double x) {
        return 0.5 * (1 + std::erf(x / std::sqrt(2)));
    }
    void set_left(double l) { _left = l; _F_left = cdf(l); }
    void set_right(double r) { _right = r; _F_right = cdf(r); }
private:
    double _left, _right, _F_left, _F_right;
    TRealUniform _Ureal;
};

void affiche_gaussian_ziggurat(unsigned n) {
    gaussian_density g;
    gaussian_queue<> q;
    {
        auto Z = mk_ziggurat(g, &gaussian_density::inverse, q, n);
        ostringstream name;
        name << "../data/gaussian_" << n << ".dat";
        ofstream file(name.str());
        file << Z << std::endl;
        file.close();
    }
}
/*
int main() {
    affiche_gaussian_ziggurat(8);
    affiche_gaussian_ziggurat(32);
    affiche_gaussian_ziggurat(256);
    
    return 0;
}*/
