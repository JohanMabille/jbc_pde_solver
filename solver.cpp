#include "solver.hpp"
#include "closed_form.hpp"
#include "stdio.h"

using namespace std;

Solver::Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, double theta) :
    m_spot(spot),
	m_vol(volatility),
	m_matu(maturity),
	m_ir(interest_rate),
	m_theta(theta),
	m_spot_points(spot_points),
	m_time_points(time_points) {
	
	m_min_time = 0.;
    m_max_time = m_matu;
    m_dt = (m_max_time - m_min_time) / m_time_points;
    m_std_dev = m_vol * sqrt(m_matu);
    m_min_space = log(m_spot) - 3 * m_std_dev;
    m_max_space = log(m_spot) + 3 * m_std_dev;
    m_dx = (m_max_space - m_min_space) / m_spot_points;

	// Initialize the matrix M and its inverse
	M = Matrix<double>(m_spot_points, m_spot_points);
	Matrix<double> iden = Matrix<double>(m_spot_points);

	// M is full of zeros
	// initialize coefs
	m_coefm0 = -0.5 * pow(m_vol / m_dx, 2) - pow(m_vol, 2) / (4 * m_dx) + 0.5 * m_ir / m_dx;
	m_coefm1 = pow(m_vol / m_dx, 2) + m_ir;
	m_coefm2 = -0.5 * pow(m_vol / m_dx, 2) + pow(m_vol, 2) / (4 * m_dx) - 0.5 * m_ir / m_dx;
	m_coefm00 = -0.5 * pow(m_vol / m_dx, 2) - pow(m_vol, 2) / (2 * m_dx) + m_ir / m_dx + m_ir;
	m_coefm01 = pow(m_vol / m_dx, 2) + 0.5 * pow(m_vol, 2) / m_dx - m_ir / m_dx;
	m_coefm02 = -0.5 * pow(m_vol / m_dx, 2);
	m_coefmn0 = -0.5 * pow(m_vol / m_dx, 2);
	m_coefmn1 = pow(m_vol / m_dx, 2) - 0.5 * pow(m_vol, 2) / m_dx + m_ir / m_dx;
	m_coefmn2 = -0.5 * pow(m_vol / m_dx, 2) + pow(m_vol, 2) / (2 * m_dx) - m_ir / m_dx + m_ir;

	// we start to fill M
	// first row
	M.set_elem_at(0, 0, m_coefm00);
	M.set_elem_at(0, 1, m_coefm01);
	M.set_elem_at(0, 2, m_coefm02);
	// row 1 to last - 1
	for (int i = 1; i < m_spot_points - 1; i++) {
		M.set_elem_at(i, i - 1, m_coefm0);
		M.set_elem_at(i, i, m_coefm1);
		M.set_elem_at(i, i + 1, m_coefm2);
	}
	// last row
	M.set_elem_at(m_spot_points - 1, m_spot_points - 3, m_coefmn0);
	M.set_elem_at(m_spot_points - 1, m_spot_points - 2, m_coefmn1);
	M.set_elem_at(m_spot_points - 1, m_spot_points - 1, m_coefmn2);
	M.display();

	M_INV = (iden + M * m_theta * m_dt).inverse();
	M_INV.display();

	Matrix<double> inverse = M.inverse();
	
	// Unfortunately the inversion doesn't work for some matrices
	inverse.dot(M).display();
}


Matrix<double> Solver::solve(double strike, bool is_call) {
	Matrix<double> res = Matrix<double>(m_spot_points, m_spot_points);
	Matrix<double> iden = Matrix<double>(m_spot_points);
	// Compute last column = payoff
	for (int i = 0; i <= m_spot_points; i++)
		res.set_elem_at(i, m_spot_points - 1, dauphine::vanilla_payoff(exp(m_max_space - i * m_dx), strike, is_call));

	printf("END OF FIRST LOOP\n");

	for (int j = m_time_points - 2; j >= 0; j--)
		res.fill_column(j, (M_INV * (iden + M * (m_theta - 1) * m_dt)).dot(res.column(j + 1)));
	
	printf("END OF SECOND LOOP\n");
	return res;
}

Solver::~Solver() {}
