#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <chrono>
#include <cmath>
#include <random>

#include "zero_search.h"
#include "inverse_function.h"
#include "ziggurat.h"
#include <boost/math/special_functions/beta.hpp>

struct beta_density {
    beta_density(double a, double b) 
        : _alpha(a), _beta(b), _cst(boost::math::beta(a,b)) {}
    double operator()(double x) const { 
        return std::pow(x, _alpha-1) * std::pow(1-x, _beta-1) / _cst;
    }
    double derivative(double x) const {
        return ((_alpha-1) * std::pow(x, _alpha-2) * std::pow(1-x, _beta-1)
               - (_beta-1) * std::pow(x, _alpha-1) * std::pow(1-x, _beta-2)) / _cst;
    }
    double mode() const { return (_alpha-1.) / (_alpha+_beta-2.); }
    double max() const { return this->operator()(mode()); }
    private:
        double _alpha, _beta, _cst;
};

using namespace std;
using namespace utils;
using namespace ziggurat;

void affiche_ziggurat(double alpha, double beta, unsigned n) {
    beta_density f(alpha, beta);
    {
        ostringstream name;
        name << "../data/beta_" << alpha << "_" << beta << "_density.dat";
        ofstream file(name.str());
        for (double x = 0; x <= 1; x += 1e-3)
            file << x << "\t" << f(x) << endl;
        file.close();
    }
    auto inv_f = mk_inverse_function(f);
    {
        auto Z = mk_ziggurat(f, inv_f, n);
        ostringstream name;
        name << "../data/beta_" << alpha << "_" << beta << "_" << n << ".dat";
        ofstream file(name.str());
        file << Z << std::endl;
        file.close();
    }
}

int main() {  
    affiche_ziggurat(2, 2, 8);  
    affiche_ziggurat(2, 2, 32);  
    affiche_ziggurat(2, 2, 256);  
    
    affiche_ziggurat(1.5, 2.5, 8);  
    affiche_ziggurat(1.5, 2.5, 32);  
    affiche_ziggurat(1.5, 2.5, 256);  
    
    return 0;
}
