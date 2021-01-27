#ifndef SOLVER_H
#define SOLVER_H

#include "matrix.hpp"

class Solver {
	public:
		explicit Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, double theta);
		~Solver();

		Matrix<double> solve(double strike, bool is_call);

	protected:
		Matrix<double> M;
		Matrix<double> M_INV;
		double m_coefm0, m_coefm1, m_coefm2;
		double m_coefm00, m_coefm01, m_coefm02;
		double m_coefmn0, m_coefmn1, m_coefmn2;
	
		double m_spot, m_vol, m_matu, m_ir, m_theta;
    	double m_min_time, m_max_time, m_dt;
    	double m_std_dev;
    	double m_min_space, m_max_space, m_dx;
    	int m_spot_points, m_time_points;
};

#endif
