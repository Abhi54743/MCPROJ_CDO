#ifndef ZIGGURAT_H
#define ZIGGURAT_H
#include <vector>
#include <random>
#include <algorithm>

#include"zero_search.h"

namespace ziggurat {

template <typename TRealUniform = std::uniform_real_distribution<double>>
class rectangle {
public:
    explicit rectangle(double l = 0, double r = 1, 
                       double b = 0, double t = 1)
        : _l(l), _r(r), _dx(r-l), _b(b), _t(t), _dy(t-b) {}
    double area() const { return _dx * _dy; }
    bool contains_x(double a) const { return (a > _l) && (a < _r); }
    bool contains_y(double a) const { return (a > _b) && (a < _t); }
    double left() const { return _l; }
    double right() const { return _r; }
    double top() const { return _t; }
    double bottom() const { return _b; }
    template <typename TGenerator>
    double generate_x(TGenerator & gen) { return _l + _dx * _Ureal(gen); }
    template <typename TGenerator>
    double generate_y(TGenerator & gen) { return _b + _dy * _Ureal(gen); }
    friend std::ostream & operator<<(std::ostream & o, rectangle<TRealUniform> const & R) {
        return o << "rect from " << R._l << "," << R._b << " to " 
                 << R._r << "," << R._t << " fs empty border fc ls 2";
    };
private:
    double _l, _r, _dx, _b, _t, _dy;
    TRealUniform _Ureal;
};

template <typename TDensity, typename TInverseDensity, 
          typename TIntUniform = std::uniform_int_distribution<int>,
          typename TRealUniform = std::uniform_real_distribution<double>>
class algo_compact {
public:
    algo_compact(TDensity const & f, TInverseDensity const & inv_f, unsigned n);
    double total_area() const;
    algo_compact & build(double a);
    template <typename TGenerator>
    double operator()(TGenerator & gen);
    friend std::ostream & operator<<(std::ostream & o, algo_compact const & az) {
        for (unsigned i = 0; i < az._strates.size(); ++i)
            o << "set object " << i+1 << " " << az._strates[i] << std::endl;
        return o;
    };
protected:
    TDensity _density;
    TInverseDensity _inverse_density;
    TIntUniform _Uint;
    std::vector<rectangle<TRealUniform>> _strates;
};

template <typename TDensity, typename TInverseDensity, 
          typename TIntUniform, typename TRealUniform>
algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform>
    ::algo_compact(TDensity const & f, TInverseDensity const & inv_f, unsigned n) 
        : _density(f), _inverse_density(inv_f), _Uint(0, n-1), _strates(n) {
    double min_area = 1. / (double) n;
    auto support = _inverse_density(0);
    double max_area = _density.max() * (support.second - support.first) / (double) n;
    auto fct = [&](double a) -> double {
        return build(a).total_area() - n*a;
    };
    double area = useful::dichotomie(fct, 0, min_area, max_area);
    build(area);
};

template <typename TDensity, typename TInverseDensity,
          typename TIntUniform, typename TRealUniform>
double algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform>
    ::total_area() const {
    double area = 0;
    for (auto & rect : _strates) {
        area += rect.area();
    }
    return area;
};

template <typename TDensity, typename TInverseDensity,
          typename TIntUniform, typename TRealUniform>
algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform> & 
algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform>
    ::build(double a) {
    double yi = 0, yip1;
    unsigned n = _strates.size();
    for (unsigned i = 0; i < n-1; ++i) {
        auto xi = _inverse_density(yi);
        yip1 = (xi.second > xi.first) ? 
            std::min(yi + a / (xi.second - xi.first), _density.max()) : yi;
        _strates[i] = rectangle<TRealUniform>(xi.first, xi.second, yi, yip1);
        yi = yip1;
    }
    auto xi = _inverse_density(yi);
    _strates[n-1] = rectangle<TRealUniform>(xi.first, xi.second, yi, _density.max());
    return *this;
};

template <typename TDensity, typename TInverseDensity, 
          typename TIntUniform, typename TRealUniform>
template <typename TGenerator>
double algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform>
    ::operator()(TGenerator & gen) {
    while (true) {
        unsigned i = _Uint(gen);
        double x = _strates[i].generate_x(gen);
        if (i+1 < _strates.size() && _strates[i+1].contains_x(x)) return x;
        double y = _strates[i].generate_y(gen);
        if (y < _density(x)) return x;
    }
};

template <typename TDensity, typename TInverseDensity,
          typename TQueueDistribution,  
          typename TIntUniform = std::uniform_int_distribution<int>,
          typename TRealUniform = std::uniform_real_distribution<double>>
class algo_non_compact {
public:
    algo_non_compact(TDensity const & f, TInverseDensity const & inv_f, 
                     TQueueDistribution q, unsigned n);
    double total_area() const;
    double strate0_area(double x1_l, double x1_r, double y1) const {
        return _queue.cdf(x1_l) + (1-_queue.cdf(x1_r)) + (x1_r - x1_l) * y1;
    }
    algo_non_compact & build(double a);
    template <typename TGenerator>
    double operator()(TGenerator & gen);
    friend std::ostream & operator<<(std::ostream & o, algo_non_compact const & az) {
        for (unsigned i = 0; i < az._strates.size(); ++i)
            o << "set object " << i+1 << " " << az._strates[i] << std::endl;
        return o;
    };
protected:
    TDensity _density;
    TInverseDensity _inverse_density;
    TQueueDistribution _queue;
    TIntUniform _Uint;
    TIntUniform _Ureal;
    std::vector<rectangle<TRealUniform>> _strates;
};

template <typename TDensity, typename TInverseDensity, typename TQueueDistribution,  
          typename TIntUniform, typename TRealUniform>
algo_non_compact<TDensity, TInverseDensity, 
                 TQueueDistribution, TIntUniform, TRealUniform>
    ::algo_non_compact(TDensity const & f, TInverseDensity const & inv_f, 
                       TQueueDistribution q, unsigned n)
        : _density(f), _inverse_density(inv_f), 
          _queue(q), _Uint(0, n-1), _strates(n-1) {
    double min_area = 1. / (double) n;
    double y1 = _density.max() / 2.;
    auto x1 = _inverse_density(y1);
    double max_area = ((_density.max()-y1) * (x1.second - x1.first) 
                       + strate0_area(x1.first, x1.second, y1)) / (double) n;
    auto fct = [&](double a) -> double {
        return build(a).total_area() - n*a;
    };
    double area = useful::dichotomie(fct, 0, min_area, max_area);
    build(area);
    _queue.set_left(_strates[0].left());
    _queue.set_right(_strates[0].right());
};

template <typename TDensity, typename TInverseDensity, typename TQueueDistribution,  
          typename TIntUniform, typename TRealUniform>
double algo_non_compact<TDensity, TInverseDensity, 
                        TQueueDistribution, TIntUniform, TRealUniform>
    ::total_area() const {
    double area = 0;
    for (auto & rect : _strates) {
        area += rect.area();
    }
    area += strate0_area(_strates[0].left(), _strates[0].right(), _strates[0].bottom());
    return area;
};

template <typename TDensity, typename TInverseDensity,
          typename TQueueDistribution, typename TIntUniform, typename TRealUniform>
algo_non_compact<TDensity, TInverseDensity, 
                 TQueueDistribution, TIntUniform, TRealUniform> &
algo_non_compact<TDensity, TInverseDensity,
                 TQueueDistribution, TIntUniform, TRealUniform>
    ::build(double a) {
    double yi = 0, yip1;
    auto fct = [&](double y) {
        auto x = _inverse_density(y);
        return strate0_area(x.first, x.second, y);
    };
    unsigned n = _strates.size() + 1;
    yip1 = useful::dichotomie(fct, a, 0, _density.max()) / (double) n;
    yi = yip1;
    for (unsigned i = 0; i < _strates.size()-1; ++i) {
        auto xi = _inverse_density(yi);
        yip1 = (xi.second > xi.first) ? 
            std::min(yi + a / (xi.second - xi.first), _density.max()) : yi;
        _strates[i] = rectangle<TRealUniform>(xi.first, xi.second, yi, yip1);
        yi = yip1;
    }
    auto xi = _inverse_density(yi);
    _strates.back() = rectangle<TRealUniform>(xi.first, xi.second, yi, _density.max());
    return *this;
};

template <typename TDensity, typename TInverseDensity,
          typename TQueueDistribution, typename TIntUniform, typename TRealUniform>
template<typename TGenerator>
double algo_non_compact<TDensity, TInverseDensity,
                        TQueueDistribution, TIntUniform, TRealUniform>
    ::operator()(TGenerator & gen) {
    while (true) {
        unsigned i = _Uint(gen);
        if (i == 0) { return _queue.generate(gen); }
        --i;
        double x = _strates[i].generate_x(gen);
        if (i+1 < _strates.size() && _strates[i+1].contains_x(x)) return x;
        double y = _strates[i].generate_y(gen);
        if (y < _density(x)) return x;
    }
};


template <typename TDensity, typename TInverseDensity, 
          typename TIntUniform = std::uniform_int_distribution<int>,
          typename TRealUniform = std::uniform_real_distribution<double>>
algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform>
mk_ziggurat(TDensity const & f, TInverseDensity const & inv_f, unsigned n) {
    return algo_compact<TDensity, TInverseDensity, TIntUniform, TRealUniform>
           (f, inv_f, n);
};

template <typename TDensity, typename TInverseDensity, 
          typename TQueueDistribution,  
          typename TIntUniform = std::uniform_int_distribution<int>,
          typename TRealUniform = std::uniform_real_distribution<double>>
algo_non_compact<TDensity, TInverseDensity, TQueueDistribution, TIntUniform, TRealUniform>
mk_ziggurat(TDensity const & f, TInverseDensity const & inv_f, 
            TQueueDistribution const & q, unsigned n) {
    return algo_non_compact<TDensity, TInverseDensity, 
                            TQueueDistribution, TIntUniform, TRealUniform>
           (f, inv_f, q, n);
};

} // namespace ziggurat

#endif // ZIGGURAT_H

