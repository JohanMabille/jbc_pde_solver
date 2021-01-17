#include "closed_form.hpp"
#include "Solver.h"
#include <iostream>

// Guidelines:
//
// 1] Equation
// Even if the volatility and the rate are constant in the BS model,
// we can prove that, under certain assumptions, we can use the same
// PDE equation with volatility surface and rate curves. We could also
// want to take into account the repo (even if it could theoretically
// be part of the r factor). Therefore, it is more generic to solve
// the following generic equation;
// df/dt = a(x,t)d2f/dx2 + b(x,t)df/dx + c(x, t)f + d(x, t).
// The pricer should ask the coefficients a, b, c and d to an
// abstract class that will be inherited by classes implementing
// different models.
//
// 2] Payoff
// The pricer should be able to price exotic options, such as
// barriers or asian options. To do so, the pricer must be able
// to call an external function between each step. Define an API
// that allows to register a function or an abstract class modeling
// a payoff.
int main(int argc, const char * argv[]){
    double K = 100.;
    double sig = 0.5;
    double matu = 20./365.;
    double r = 0.0;


    for(int i=0;i<-1;i++){
        std::cout <<
        "Forward/Spot price: " << 10.*i << ", call price (fwd): " << dauphine::bs_price(10.*i,K,sig,matu,r,true) <<
        " of which time value (fwd): " << dauphine::bs_time_value(10.*i,K,sig,matu,r,true) <<
        ", call price (spot): " << dauphine::bs_price(10.*i,K,sig,matu,(float) r,true) <<
        " of which time value (spot): " << dauphine::bs_time_value(10.*i,K,sig,matu,(float) r,true) << std::endl;
    }

    Solver s(100.,sig,matu,r,20,50,0.5);
    auto res = s.solve_BS(K, true);
    Solver::print_vector_array(res);

    std::cout << res[26][0] << std::endl;

    return 0;
}
