

#include <iostream>
#include "Matrix.hpp"

class Solver_matrix
{
public:
	explicit Solver_matrix(double spot, double volatility, double maturity, float interest_rate, int points, double theta);
	explicit Solver_matrix(double spot, Matrix<double> volatility, double maturity, Matrix<float> interest_rate, int points, double theta);
	~Solver_matrix();

	void solve(double strike, bool is_call);
	void solve_constant(double strike, bool is_call);
	void solve_non_constant(double strike, bool is_call);
	void display_res();
	void display_end_res();
	
	void compare_with_BS();
	double price();
	double delta();
	double gamma();
	double theta();
	double vega();

protected:
	double m_spot, m_vol, m_matu, m_ir, m_theta, m_strike;
	bool m_is_call;
	double m_min_time, m_max_time, m_dt;
	double m_std_dev;
	double m_min_space, m_max_space, m_dx;

	Matrix<double> M_vol;
	Matrix<float> M_ir;

	Matrix<double>* Res;
	Matrix<double> Spot; 
	Matrix<double> BS_Price;
	Matrix<double> W;

	double m_coefm[3];
	double m_coefm0[3];
	double m_coefmn[3];
	int m_points;
};
