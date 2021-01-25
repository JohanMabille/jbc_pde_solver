#include "Solver.h"
#include "closed_form.hpp"
#include <math.h>

Solver::Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, double theta):
    m_spot(spot), m_vol(volatility), m_matu(maturity), m_ir(interest_rate), m_time_points(time_points), m_spot_points(spot_points), m_theta(theta) {
    m_min_time = 0.;
    m_max_time = m_matu;
    m_dt = (m_max_time-m_min_time)/m_time_points;
    m_std_dev = m_vol * sqrt(m_matu);
    m_min_space = log(m_spot)-3*m_std_dev;
    m_max_space = log(m_spot)+3*m_std_dev;

    m_dx = (m_max_space-m_min_space)/m_spot_points;
    std::cout << "m_dt: " << m_dt << std::endl;
    std::cout << "m_std_dev: " << m_std_dev << std::endl;
    std::cout << "m_dx: " << m_dx << std::endl;

    for(int i=0;i<=m_spot_points;i++){
        std::vector<double> temp(m_time_points+1,0);
        m_res.push_back(temp);
    }

    double sig2dx2 = m_vol*m_vol/(m_dx*m_dx);
    double sig2dx = m_vol*m_vol/m_dx;
    double thetadt = m_theta*m_dt;
    double invthetadt = (1-m_theta)*m_dt;

    double rdx = m_ir/m_dx;
    double temp1 = -0.5*sig2dx2-0.25*sig2dx+0.5*rdx;
    double temp2 =  m_ir+sig2dx2;
    double temp3 = -0.5*sig2dx2+0.25*sig2dx-0.5*rdx;

    m_coeffs[0][0] = thetadt*temp1;
    m_coeffs[0][1] = invthetadt*temp1;
    m_coeffs[1][0] = 1+thetadt*temp2;
    m_coeffs[1][1] = invthetadt*temp2-1;
    m_coeffs[2][0] = thetadt*temp3;
    m_coeffs[2][1] = invthetadt*temp3;


    temp1 = 0.5*sig2dx2+0.25*sig2dx-0.5*rdx;
    temp2 =  -sig2dx2-m_ir;
    temp3 = 0.5*sig2dx2-0.25*sig2dx+0.5*rdx;

    m_coeffs_theta0[0] = m_dt*temp1;
    m_coeffs_theta0[1] = 1+m_dt*temp2;
    m_coeffs_theta0[2] = m_dt*temp3;

    m_coeffs_theta0_edge[0] = m_dt*(0.5*sig2dx2);
    m_coeffs_theta0_edge[1] = m_dt*(-sig2dx2+0.5*sig2dx+rdx);
    m_coeffs_theta0_edge[2] = 1+m_dt*(0.5*sig2dx2-0.5*sig2dx+rdx-m_ir);

    std::cout << "coeffs: " << std::endl;
    for(int i=0;i<3;i++){
        for(int n=0;n<2;n++){
            std::cout << m_coeffs[i][n] << " , ";
        }
        std::cout << std::endl;
    }

    std::cout << "coeffs theta-0: " << std::endl;
    for(int n=0;n<2;n++){
        std::cout << m_coeffs_theta0[n] << " , ";
    }
    std::cout << std::endl;

    std::cout << "coeffs: " << std::endl;
    for(int i=0;i<3;i++){
        for(int n=0;n<2;n++){
                std::cout << m_coeffs[i][n] << " , ";
        }
        std::cout << std::endl;
    }

}

double Solver::compute_vertex(double values[][2], int di, int dn){
    double tot=0;
    for(int i=0;i<3;i++){
        for(int n=0;n<2;n++){
            if(di!=i || dn!=n){
                std::cout << m_coeffs[i][n] << " * " << values[i][n] << std::endl;
                tot+=m_coeffs[i][n] * values[i][n];
            }
        }
    }
    std::cout << "tot: " << tot << std::endl;
    return -tot/m_coeffs[di][dn];
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

std::vector<std::vector<double> > Solver::solve_BS(double strike, bool is_call){

//    std::cout << "min_space, max_space= " << min_space << ", " << max_space << std::endl;
    for(int i=0;i<=m_spot_points;i++)
        m_res[i][m_time_points] = dauphine::vanilla_payoff(exp(m_max_space-i*m_dx),strike,is_call);

    for(int i=m_time_points-1;i>=0;i--){
        double temp[3] {m_res[2][i+1], m_res[1][i+1], m_res[0][i+1]};
        m_res[0][i]  = compute_vertex_theta0_edge(temp);

        double temp3[3] {m_res[3][i+1], m_res[2][i+1], m_res[1][i+1]};
        m_res[1][i]  = compute_vertex_theta0_edge(temp3);

        for(int j=1;j<m_spot_points;j++){
            double temp2[3][2] {{m_res[j-1][i], m_res[j-1][i+1]}, {m_res[j][i], m_res[j][i+1]}, {0, m_res[j+1][i+1]}};
            m_res[j+1][i] = compute_vertex(temp2,1,0);
        }
    }
    return m_res;
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

Solver::~Solver()
{
    //dtor
}
