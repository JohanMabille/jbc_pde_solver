#include "stdafx.h"

#include "Solver_matrix.h"
#include "closed_form.hpp"




Solver_matrix::Solver_matrix(double spot, double volatility, double maturity, float interest_rate, int points, double theta)
	:m_points(points),
	Solver(spot, volatility, maturity, interest_rate, points, points, theta)

{
	// initialize the matrix M 
	M = new Matrix<double>(m_points, m_points);
	// M is full of zeros
	// initialize coefs
	using namespace std;
	m_coefm[0] = -0.5 * pow(m_vol / m_dx, 2) - pow(m_vol, 2) / (4 * m_dx) + 0.5 * m_ir / m_dx;
	m_coefm[1] = pow(m_vol / m_dx, 2) + m_ir;
	m_coefm[2] = -0.5 * pow(m_vol / m_dx, 2) + pow(m_vol, 2) / (4 * m_dx) - 0.5 * m_ir / m_dx;
	m_coefm0[0] = -0.5 * pow(m_vol / m_dx, 2) - pow(m_vol, 2) / (2 * m_dx) + m_ir / m_dx + m_ir;
	m_coefm0[1] = pow(m_vol / m_dx, 2) + 0.5 * pow(m_vol, 2) / m_dx - m_ir / m_dx;
	m_coefm0[2] = -0.5 * pow(m_vol / m_dx, 2);
	m_coefmn[0] = -0.5 * pow(m_vol / m_dx, 2);
	m_coefmn[1] = pow(m_vol / m_dx, 2) - 0.5 * pow(m_vol, 2) / m_dx + m_ir / m_dx;
	m_coefmn[2] = -0.5 * pow(m_vol / m_dx, 2) + pow(m_vol, 2) / (2 * m_dx) - m_ir / m_dx + m_ir;

	// we start to fill M
	// first row
	M->set_elem_at(0, 0, m_coefm0[0]);
	M->set_elem_at(0, 1, m_coefm0[1]);
	M->set_elem_at(0, 2, m_coefm0[2]);
	// row 1 to last - 1
	for (int i = 1; i < m_points - 1; i++)
	{
		M->set_elem_at(i, i-1, m_coefm[0]);
		M->set_elem_at(i, i, m_coefm[1]);
		M->set_elem_at(i, i+1, m_coefm[2]);
	}
	// last row
	M->set_elem_at(m_points - 1, m_points - 3, m_coefmn[0]);
	M->set_elem_at(m_points - 1, m_points - 2, m_coefmn[1]);
	M->set_elem_at(m_points - 1, m_points - 1, m_coefmn[2]);
	M->display();
}


Solver_matrix::~Solver_matrix()
{
}


Matrix<double> Solver_matrix::solving(double strike, bool is_call)
{
	Matrix<double>* Res = new Matrix<double>(m_points, m_points);
	Matrix<double>* Identity = new Matrix<double>(m_points);
	// Compute last column = payoff
	for (int i = 0; i <= m_spot_points; i++)
	{
		Res->set_elem_at(i, m_points-1, dauphine::vanilla_payoff(exp(m_max_space - i*m_dx), strike, is_call));
	}
	//(*Identity + (*M) * (m_theta * m_dt)).display();
	//(((*Identity + (*M) * (m_theta * m_dt)).inverse()) * (*Identity + (*M) * ((m_theta - 1) * m_dt))).dot(Res->column(m_points - 1)).display();
	//Res->fill_column(m_points - 2, (((*Identity + (*M) * (m_theta * m_dt)).inverse()) * (*Identity + (*M) * ((m_theta - 1) * m_dt))).dot(Res->column(m_points - 1)));
	//Res->column(m_points - 2) = (((*Identity + (*M) * (m_theta * m_dt)).inverse()) * (*Identity + (*M) * ((m_theta - 1) * m_dt))).dot(Res->column(m_points - 1));
	
	for (int j = m_time_points - 2; j >= 0; j--)
	{
		Res->fill_column(j, (((*Identity + (*M) * (m_theta * m_dt)).inverse()) * (*Identity + (*M) * ((m_theta - 1) * m_dt))).dot(Res->column(j+1)));
	}
	//Res->column(m_points - 1).display();
	Res->display();
	return *Res;
}

/*
std::vector<double> Solver_matrix::compute_vector(double values[])
{

}
*/