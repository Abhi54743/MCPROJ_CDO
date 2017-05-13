#ifndef INVERSE_FUNCTION_H
#define INVERSE_FUNCTION_H
#include <iostream>

#include "zero_search.h"

namespace useful {

template <typename TUnimodalFunction>
struct inverse_unimodal_function {
    inverse_unimodal_function(
        TUnimodalFunction const & f,
        double eps = 1e3*std::numeric_limits<double>::epsilon()) 
        : _f(f), _eps(eps) {}
    std::pair<double, double> operator()(double y) const; 
    private:
        TUnimodalFunction _f;
        double _eps;
};

template <typename TUnimodalFunction>
std::pair<double, double> inverse_unimodal_function<TUnimodalFunction>::operator()(double y) const {
    if (std::fabs(_f.max()-y) < _eps) 
        return std::pair<double, double>(_f.mode(), _f.mode());
    double x_m, x_p;
    try {
        x_m = zero_search(_f, y, 0, _f.mode(), _eps);
    } catch (std::logic_error & ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        throw ex;
    }
    try {
        x_p = zero_search(_f, y, _f.mode(), 1, _eps);
    } catch (std::logic_error & ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        throw ex;
    }
    return std::pair<double, double>(x_m, x_p);
};

template <typename TUnimodalFunction>
inverse_unimodal_function<TUnimodalFunction> 
mk_inverse_function(TUnimodalFunction const & f) {
    return inverse_unimodal_function<TUnimodalFunction>(f);
};

}
#endif
