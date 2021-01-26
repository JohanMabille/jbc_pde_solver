
#pragma once
#include "Solver.h"
#include "Matrix.h"

class Solver_matrix : public Solver
{
public:
	explicit Solver_matrix(double spot, double volatility, double maturity, float interest_rate, int points, double theta);
	~Solver_matrix();

	Matrix<double> solving(double strike, bool is_call);
	//Matrix<double> compute_vector(double values[]);

protected:
	Matrix<double>* M;
	double m_coefm[3];
	double m_coefm0[3];
	double m_coefmn[3];
	int m_points;
};
