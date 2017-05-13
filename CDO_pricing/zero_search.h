#ifndef ZERO_SEARCH_H
#define ZERO_SEARCH_H

namespace useful { 

/** 
* @brief 
* 
* @param f
* @param l
* @param r
* @param eps
* 
* @return 
*/
template <typename TRealFunction>
double dichotomie(TRealFunction const & f, 
                  double y = 0, double l = 0, double r = 1, 
                  double eps = 2*std::numeric_limits<double>::epsilon()) 
{
    #if DEBUG
    std::cout << "#Dichotomie: Search of " << y << " in ]" << l << ";" << r << "[" << std::endl;
    #endif
    double yl = f(l) - y, yr = f(r) - y;
    double x, yx;
    do {
        if (yl * yr > 0) throw std::logic_error("Problem of monotonicity in dichotomie algorithm");
        x = 0.5 * (l + r);
        yx = f(x) - y;
        if (yx * yl < 0) { r = x; yr = yx; }
        else {
            if (yx * yr < 0) { l = x; yl = yx; }
            else return x;
        }
        #if DEBUG 
        std::cout << "#Dichotomie: " << x << "\t" << yx << "\t" << yl << "\t" << yr << std::endl;
        #endif
    } while (std::fabs(yx) > eps);
    return x;
};

/** 
* @brief 
* 
* @param f
* @param df
* @param y
* @param l
* @param r
* @param eps
* 
* @return 
*/
template <typename TRealFunction, typename TDerivativeFunction>
double newton(TRealFunction const & f, 
              TDerivativeFunction const & df,
              double y = 0, double l = 0, double r = 1, 
              double eps = 2.*std::numeric_limits<double>::epsilon())
{
    #if DEBUG
    std::cout << "#Newton: Search of " << y << " in ]" << l << ";" << r << "[" << std::endl;
    #endif
    double x = 0.5 * (l + r);
    unsigned max_iter = 0;
    while (std::fabs(f(x)-y) > eps) {
        double df_x = df(x);
        if (std::fabs(df_x) < 10*eps) { df_x = (df_x < 0) ? -1 : 1; }
        if (std::fabs(df_x) > 1e2) { df_x = (df_x < 0) ? -1e2 : 1e2; }
        x -= (f(x) - y) / df_x;
        x = std::max(std::min(x, r), l);
        #if DEBUG
        std::cout << "#Newton: " << x << "\t" << f(x) << "\t" << df_x << std::endl;
        #endif
        if (++max_iter > 1e3) 
            throw std::logic_error("Non convergence of Newton's algorithm"); 
    };
    return x;
};

/** 
* @brief 
* 
* @param f
* @param y
* @param l
* @param r
* @param eps
* 
* @return 
*/
template <typename TFunctionWithDerivative>
double newton(TFunctionWithDerivative const & f, 
              double y = 0, double l = 0, double r = 1, 
              double eps = 2.*std::numeric_limits<double>::epsilon())
{
    return newton(f, [f](double x) { return f.derivative(x); }, y, l, r, eps);
};


template <typename TFunctionWithDerivative>
double zero_search(TFunctionWithDerivative const & f, 
                   double y = 0, double l = 0, double r = 1, 
                   double eps = 2*std::numeric_limits<double>::epsilon()) 
{
    double x;
    try {
        x = newton(f, y, l, r, eps);
    } catch (std::logic_error & ex) {
        std::cerr << "Exception: " << ex.what() << std::endl;
        try {
            x = dichotomie(f, y, l, r, eps);
        } catch (std::logic_error & ex) {
            std::cerr << "Exception: " << ex.what() << std::endl;
            throw std::logic_error("HUGE PROBLEM");
        }
    }
    return x;
};


} // end namespace useful

#endif // ZERO_SEARCH_H
