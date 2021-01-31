#include "closed_form.hpp"
#include "solver.hpp"
#include "matrix.hpp"
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


void test_matrix() {
	Matrix<double> M = Matrix<double>(3);

	M.set_elem_at(1, 1, 8.5);
	M.set_elem_at(0, 1, 1.5);
	M.set_elem_at(0, 2, 0.3);
	M.set_elem_at(2, 1, -5);
	M.set_elem_at(2, 2, 8.5);
	M.display();
	
	Matrix<double> N = M.inverse();
	N.display();

	Matrix<double> P = N.dot(M);
	P.display();

	std::cout << "Test matrix is finished" << std::endl;
}

void test_solver() {
	double spot = 100.;
	double volatility = 0.16;
	double maturity = 0.25;
	float interest_rate = (float) 0.04;
	ul spot_points = 16;
	ul time_points = 16;
	double theta = 0.5;
	double strike = 100.;
	bool is_call = true;

	Solver mysolver = Solver(spot, volatility, maturity, interest_rate, time_points, spot_points, theta);
	Matrix<double> M = mysolver.solve(strike, is_call);
}


int main(int argc, const char * argv[]) {
    double K = 100.;
    double sig = 0.16;
    double matu = 10./365.;
    double r = 0.0;

	for (double i = 5.; i <= 15.; i = i + 1.) {
    	std::cout << "Forward/Spot price: " << 10. * i << ", call price (fwd): " << dauphine::bs_price(10. * i, K, sig, matu, r, true)
			  	  << " of which time value (fwd): " << dauphine::bs_time_value(10. * i, K, sig, matu, r, true)
			  	  << ", call price (spot): " << dauphine::bs_price(10. * i, K, sig, matu, (float) r, true)
			  	  << " of which time value (spot): " << dauphine::bs_time_value(10. * i, K, sig, matu,(float) r, true)
			  	  << std::endl;
	}

	test_matrix();

	test_solver();

    return 0;
}
