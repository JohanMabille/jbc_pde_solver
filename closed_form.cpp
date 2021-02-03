#include "closed_form.hpp"
#include <math.h>
#include <algorithm>

namespace dauphine{
    double normal_cdf(double x){
        return 0.5 * (1 + erf(x / sqrt(2.)));
    }

    double vanilla_payoff(double fwd, double strike, bool is_call){
        return std::max(is_call ? fwd - strike : strike - fwd, 0.);
    }

    double bs_time_value(double fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call){
        double discount_factor = exp(-discount_rate*maturity);
        return bs_price(fwd,strike,volatility,maturity,discount_rate,is_call) - discount_factor*vanilla_payoff(fwd,strike,is_call);
    }

    double bs_time_value(double spot, double strike, double volatility, double maturity, float interest_rate, bool is_call){
        return bs_time_value(spot*exp(interest_rate*maturity), strike, volatility, maturity, (double) interest_rate, is_call);
    }

    double bs_price(double fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call){
        // This implmeentation is less stable than the previous one (you can have
        // oscillations on the wings).
        // THe initial implementation was the non discounted function of the forward, it didn't need any fix
        // just call it with the forward instead of the spot, and discount the price.
        // This makes the implementation independent form the rate
        double discount_factor = exp(-discount_rate*maturity);
        if(strike == 0.){
            return discount_factor*vanilla_payoff(fwd,strike,is_call);
        }
        double stddev = volatility * sqrt(maturity);
        if(stddev == 0.){
            return discount_factor*vanilla_payoff(fwd,strike,is_call);
        }
        double tmp = log(fwd / strike) / stddev;
        double d1 = tmp + 0.5 * stddev;
        double d2 = tmp - 0.5 * stddev;
        if(is_call){
            return discount_factor*(fwd * normal_cdf(d1) - strike * normal_cdf(d2));
        }
        return discount_factor*(strike * normal_cdf(-d2) - fwd * normal_cdf(-d1));
    }

    double bs_price(double spot, double strike, double volatility, double maturity, float interest_rate, bool is_call){
        return bs_price(spot*exp(interest_rate*maturity), strike, volatility, maturity, (double) interest_rate, is_call);
    }


    std::vector<double> vanilla_payoff(const std::vector<double>& fwd, double strike, bool is_call){
        std::vector<double> res(fwd.size());
        for(std::size_t i = 0; i < fwd.size(); ++i){
            res[i] = vanilla_payoff(fwd[i], strike, is_call);
        }
        return res;
    }

    std::vector<double> bs_time_value(const std::vector<double>& fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call){
        std::vector<double> res(fwd.size());
        for(std::size_t i = 0; i < fwd.size(); ++i){
            res[i] = bs_time_value(fwd[i], strike, volatility,maturity,discount_rate,is_call);
        }
        return res;
    }

    std::vector<double> bs_price(const std::vector<double>& fwd, double strike, double volatility, double maturity, double discount_rate, bool is_call){
        std::vector<double> res(fwd.size());
        for(std::size_t i = 0; i < fwd.size(); ++i){
            res[i] = bs_price(fwd[i], strike, volatility, maturity, discount_rate, is_call);
        }
        return res;
    }
}

