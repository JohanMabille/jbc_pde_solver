/*
#pragma once
#include "Solver.h"
#include "Matrix.h"

class Solver_matrix :
	public Solver
{
public:
	Solver_matrix(double spot, double volatility, double maturity, float interest_rate, int points, int spot_points, double theta);
	~Solver_matrix();
	std::vector<double> compute_vector(double values[]);
private:
	Matrix<double> M;
};
*/
