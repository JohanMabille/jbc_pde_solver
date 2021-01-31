

#include "closed_form.hpp"
#include "Matrix.hpp"
#include "Solver_matrix.hpp"

#include <iostream>


void test_constant()
{
	double spot = 100.;
	double volatility = 0.2;
	double maturity = 0.25;
	float interest_rate = 0.02;
	int points = 16;
	double theta = 0.;
	double strike = 100.;
	bool is_call = true;


	Solver_matrix mysolver = Solver_matrix(spot, volatility, maturity, interest_rate, points, theta);
	mysolver.solve(strike, is_call);
	mysolver.display_end_res();

	double prix = mysolver.price();
	double delta = mysolver.delta();
	double gamma = mysolver.gamma();
	double g_theta = mysolver.theta();
	//double vega = mysolver.vega();

	std::cout << "prix = " << prix << std::endl;
	std::cout << "delta = " << delta << std::endl;
	std::cout << "gamma = " << gamma << std::endl;
	std::cout << "theta = " << g_theta << std::endl;
	//std::cout << "vega = " << vega << std::endl;
}

void test_non_constant()
{
	int points = 21;
	double spot = 100.;
	double maturity = 0.25;
	Matrix<float> interest_rate(points, 1);
	Matrix<double> volatility(points, 1);
	for (int i = 0; i < points; i++)
	{
		interest_rate.set_elem_at(i, 0, 0.02);
		volatility.set_elem_at(i, 0, 0.2);
	}
	double theta = 0.5;
	double strike = 100.;
	bool is_call = true;


	Solver_matrix mysolver = Solver_matrix(spot, volatility, maturity, interest_rate, points, theta);
	mysolver.solve(strike, is_call);
	mysolver.display_end_res();

	double prix = mysolver.price();
	double delta = mysolver.delta();
	double gamma = mysolver.gamma();
	double g_theta = mysolver.theta();
	double vega = mysolver.vega();

	std::cout << "prix = " << prix << std::endl;
	std::cout << "delta = " << delta << std::endl;
	std::cout << "gamma = " << gamma << std::endl;
	std::cout << "theta = " << g_theta << std::endl;
	std::cout << "vega = " << vega << std::endl;
}


int main()
{
	//test_constant();
	test_non_constant();
    return 0;
}

