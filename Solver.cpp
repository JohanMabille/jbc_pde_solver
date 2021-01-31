#include "Solver.h"
#include "closed_form.hpp"
#include <math.h>

Solver::Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, double theta):
    m_spot(spot), m_vol(volatility), m_matu(maturity), m_ir(interest_rate), m_theta(theta) , m_spot_points(spot_points), m_time_points(time_points){
    m_min_time = 0.;
    m_max_time = m_matu;
    m_dt = (m_max_time-m_min_time)/double(m_time_points);
    m_std_dev = m_vol * sqrt(m_matu);
    m_min_space = log(m_spot)-3*m_std_dev;
    m_max_space = log(m_spot)+3*m_std_dev;
    m_dx = (m_max_space-m_min_space)/double(m_spot_points);

    for(int i=0;i<=m_spot_points;i++){
        std::vector<double> temp(m_time_points+1,0);
        m_res.push_back(temp);
    }

    double sig2dx2 = m_vol*m_vol/(m_dx*m_dx);
    double sig2dx = m_vol*m_vol/m_dx;
    double thetadt = m_theta*m_dt;
    double invthetadt = (m_theta-1)*m_dt;
    double rdx = m_ir/m_dx;

    double temp1 = 0.5*sig2dx2+0.25*sig2dx-0.5*rdx;
    double temp2 =  -sig2dx2-m_ir;
    double temp3 = 0.5*sig2dx2-0.25*sig2dx+0.5*rdx;
    m_coeffs_theta0[0] = m_dt*temp1;
    m_coeffs_theta0[1] = 1+m_dt*temp2;
    m_coeffs_theta0[2] = m_dt*temp3;
    m_coeffs_theta0_edge[0] = m_dt*(0.5*sig2dx2);
    m_coeffs_theta0_edge[1] = m_dt*(-sig2dx2+0.5*sig2dx+rdx);
    m_coeffs_theta0_edge[2] = 1+m_dt*(0.5*sig2dx2-0.5*sig2dx+rdx-m_ir);

    Matrix<double> m_weights_matrix(m_spot_points+1, m_spot_points+1);
	double coeff_general_0 = -0.5*sig2dx2 - 0.25*sig2dx + 0.5*rdx;
	double coeff_general_1 = sig2dx2 + m_ir;
	double coeff_general_2 = -0.5*sig2dx2 + 0.25*sig2dx2 - 0.5*rdx;
	double coeff_top_0 = -0.5*sig2dx2 - 0.5*sig2dx + rdx + m_ir;
	double coeff_top_1 = sig2dx2 + 0.5*sig2dx - rdx;
	double coeff_top_2 = -0.5*sig2dx2;
	double coeff_bot_0 = -0.5*sig2dx2;
	double coeff_bot_1 = sig2dx2 - 0.5*sig2dx + rdx;
	double coeff_bot_2 = -0.5*sig2dx2 + 0.5*sig2dx - rdx + m_ir;

	// first row
	m_weights_matrix.set_elem_at(0, 0, coeff_bot_0);
	m_weights_matrix.set_elem_at(0, 1, coeff_bot_1);
	m_weights_matrix.set_elem_at(0, 2, coeff_bot_2);
	// row 1 to last - 1
	for (int i=1; i < m_spot_points; i++){
		m_weights_matrix.set_elem_at(i, i-1, coeff_general_0);
		m_weights_matrix.set_elem_at(i, i, coeff_general_1);
		m_weights_matrix.set_elem_at(i, i+1, coeff_general_2);
	}
	// last row
	m_weights_matrix.set_elem_at(m_spot_points, m_spot_points-2, coeff_top_0);
	m_weights_matrix.set_elem_at(m_spot_points, m_spot_points-1, coeff_top_1);
	m_weights_matrix.set_elem_at(m_spot_points, m_spot_points, coeff_top_2);

	Matrix<double> right_mat = (Matrix<double>(m_spot_points+1) + m_weights_matrix * thetadt);
    Matrix<double> right_mat_inv = right_mat.inverse();
    std::cout << right_mat_inv.to_string() << std::endl;
	Matrix<double> left_mat = (Matrix<double>(m_spot_points+1) + m_weights_matrix * invthetadt);
	std::cout << left_mat.to_string() << std::endl;
    m_transition_matrix = right_mat_inv.dot(left_mat);
    std::cout << m_transition_matrix.to_string() << std::endl;
}

double Solver::compute_vertex_theta0(double values[]){
    double tot=0;
    for(int i=0;i<3;i++)
        tot+=m_coeffs_theta0[i] * values[i];
    std::cout << "tot: " << tot << std::endl;
    return tot;
}

double Solver::compute_vertex_theta0_edge(double values[]){
    double tot=0;
    for(int i=0;i<3;i++)
        tot+=m_coeffs_theta0_edge[i] * values[i];
    std::cout << "tot edge: " << tot << std::endl;
    return tot;
}

std::vector<std::vector<double> > Solver::solve_BS_theta0(double strike, bool is_call){
//    std::cout << "min_space, max_space= " << min_space << ", " << max_space << std::endl;

    for(int i=0;i<=m_spot_points;i++)
        m_res[i][m_time_points] = dauphine::vanilla_payoff(exp(m_max_space-i*m_dx),strike,is_call);

    for(int i=m_time_points-1;i>=0;i--){
        double temp[3] {m_res[2][i+1], m_res[1][i+1], m_res[0][i+1]};
        m_res[0][i]  = compute_vertex_theta0_edge(temp);

        for(int j=1;j<m_spot_points;j++){
            double temp2[3] {m_res[j-1][i+1], m_res[j][i+1], m_res[j+1][i+1]};
            m_res[j][i] = compute_vertex_theta0(temp2);
        }
        double temp3[3] {m_res[m_spot_points-2][i+1], m_res[m_spot_points-1][i+1], m_res[m_spot_points][i+1]};
        m_res[m_spot_points][i]  = compute_vertex_theta0_edge(temp3);
    }

    return m_res;
}


Matrix<double> Solver::solve_BS(double strike, bool is_call){
    Matrix<double> res(m_spot_points+1, m_time_points+1);
    for(int i=0;i<=m_spot_points;i++)
        res.set_elem_at(i,m_time_points,dauphine::vanilla_payoff(exp(m_max_space-i*m_dx),strike,is_call));

    for(int j=m_time_points-1;j>=0;j--){
        //std::cout << res.extract_column(j+1).to_string() << std::endl;
        //std::cout << (m_transition_matrix).to_string() << std::endl;
        res.fill_column(j,(m_transition_matrix).dot(res.extract_column(j+1)));
    }

    for(int i=0;i<=m_spot_points;i++){
        std::cout << "S0= " << exp(m_max_space-i*m_dx);
        std::cout << "; BS price: " << dauphine::bs_price(exp(m_max_space-i*m_dx),strike, m_vol, m_matu, m_ir, is_call);
        std::cout << "; Solver price: " << res.elem_at(i,0) << std::endl;
    }

    return res;
}

Solver::~Solver()
{

}
