#include "stdafx.h"
/*

#include "Solver_matrix.h"


Solver_matrix::Solver_matrix(double spot, double volatility, double maturity, float interest_rate, int points, double theta)
{
	Solver(spot, volatility, maturity, interest_rate, points, points, theta);
}


Solver_matrix::~Solver_matrix()
{
}


std::vector<std::vector<double> > Solver::solve_matrix(double strike, bool is_call)
{
	// Compute last column = payoff
	for (int i = 0; i <= m_spot_points; i++)
	{
		m_res[i][m_time_points] = dauphine::vanilla_payoff(exp(m_max_space - i*m_dx), strike, is_call);
		std::cout << "Payoff[" << i << "] = " << m_res[i][m_time_points] << std::endl;
	}
	for (int j = m_time_points - 1; j >= 0; j--)
	{

	}
}

std::vector<double> Solver::compute_vector(double values[])
{

}
*/