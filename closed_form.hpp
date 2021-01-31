#ifndef CLOSED_FORM_HPP
#define CLOSED_FORM_HPP

#include <vector>

namespace dauphine{
    double normal_cdf(double x);
    double vanilla_payoff(double fwd, double strike, bool is_call);
    double bs_time_value(double fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call);
    double bs_time_value(double spot, double strike, double volatility, double maturity, float interest_rate, bool is_call);
    double bs_price(double fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call);
    double bs_price(double spot, double strike, double volatility, double maturity, float interest_rate, bool is_call);

    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call);
    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call);
    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call);
}

#endif
