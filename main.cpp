#include "closed_form.hpp"
#include "Solver.h"
#include "Matrix.h"
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
    double sig = 0.16;
    double matu = 10./365.;
    double r = 0.0;

    for(int i=0;i<-1;i++){
        std::cout <<
        "Forward/Spot price: " << 10.*i << ", call price (fwd): " << dauphine::bs_price(10.*i,K,sig,matu,r,true) <<
        " of which time value (fwd): " << dauphine::bs_time_value(10.*i,K,sig,matu,r,true) <<
        ", call price (spot): " << dauphine::bs_price(10.*i,K,sig,matu,(float) r,true) <<
        " of which time value (spot): " << dauphine::bs_time_value(10.*i,K,sig,matu,(float) r,true) << std::endl;
    }

    /*
    Solver s(100.,sig,matu,r,10,40,0.5);
    auto res = s.solve_BS_theta0(K, true);
    Solver::print_vector_array(res);

    std::cout << res[26][0] << std::endl;
    */

    Matrix<double> iden(6);
    std::cout << iden.to_string() << std::endl;

    std::cout << "elem at 2,2: "<< iden.elem_at(2,2) << std::endl;
    iden.set_elem_at(2,2,8);
    std::cout << "elem at 2,2: "<< iden.elem_at(2,2) << std::endl;

    iden.set_elem_at(2,2,8)->set_elem_at(3,2,-4)->set_elem_at(4,1,3);
    std::cout << iden.to_string() << std::endl;

    Matrix<double> test(4);
    test.set_elem_at(2,3,8.)->set_elem_at(1,2,-4.)->set_elem_at(0,3,5.);

    Matrix<double> test2(4,2);
    test2.set_elem_at(3,1,1.5)->set_elem_at(1,1,-2.35)->set_elem_at(0,0,5.1);
    std::cout << test.to_string() << std::endl;
    std::cout << test2.to_string() << std::endl;
    std::cout << ((test * 2.231) * (test + 1.2)).transpose().to_string() << std::endl;
    std::cout << (test.dot(test2)).to_string() << std::endl;



    return 0;
}
