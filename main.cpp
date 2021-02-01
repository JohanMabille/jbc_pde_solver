#include "Solver.h"
#include <iostream>

double exotic_payout(double S){
    if(S<90.)
        return 10.;
    else if(S>110.)
        return 20.;
    return 0.;
}

int main(int argc, const char * argv[]){
    double K = 100.;
    double sig = 0.2;
    double matu = 0.25;
    double r = 0.02;

    Solver s1(100.,new ConstantVolatility(sig),matu,new ConstantIR(r),100,101,new VanillaCall(K),0.5);
    auto res1 = s1.solve(true);
    s1.displayValues();

    Solver s2(100.,new ConstantVolatility(sig),matu,new ConstantIR(r),20,21,new GeneralPayoff(exotic_payout),0.5);
    auto res2 = s2.solve();
    std::cout << res2.to_string() << std::endl;
    s2.displayValues();

    double temp[]{0.2,0.25,0.3,0.28,0.22,0.18,0.1,0.15,0.5,0.8};
    double temp2[]{0.02,0.025,0.03,0.028,0.022,0.018,0.01,0.015,0.05,0.08};
    Solver s(100.,new GeneralVolatility(temp,10,101),matu,new GeneralIR(temp2,10,101),100,101,new VanillaCall(K),0.5);
    auto res = s.solve(true);
    s.displayValues();

    return 0;
}
