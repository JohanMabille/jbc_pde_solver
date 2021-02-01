#include "Solver.h"
#include "closed_form.hpp"
#include <math.h>
#include <stdexcept>
#include <typeinfo>

Solver::Solver(double spot, Volatility* volatility, double maturity, InterestRates* interest_rate, int time_points, int spot_points, Payoff* payoff, double theta):
    m_spot(spot), m_matu(maturity), m_theta(theta) , m_time_points(time_points), m_vol(volatility), m_ir(interest_rate), m_payoff(payoff), m_computed(false){
    if(spot_points%2==0)
        m_spot_points=spot_points;
    else
        m_spot_points=spot_points+1;

    m_min_time = 0.;
    m_max_time = m_matu;
    m_dt = (m_max_time-m_min_time)/double(m_time_points);
    m_std_dev = m_vol->getMeanVol() * sqrt(m_matu);
    m_min_space = log(m_spot)-5*m_std_dev;
    m_max_space = log(m_spot)+5*m_std_dev;
    m_dx = (m_max_space-m_min_space)/double(m_spot_points);
}

Solver::Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, Payoff* payoff, double theta) :
Solver::Solver(spot, new ConstantVolatility(volatility), maturity, new ConstantIR(interest_rate), time_points, spot_points, payoff, theta){

}

Solver::Solver(double spot, Volatility* volatility, double maturity, float interest_rate, int time_points, int spot_points, Payoff* payoff, double theta) :
Solver::Solver(spot, volatility, maturity, new ConstantIR(interest_rate), time_points, spot_points, payoff, theta){

}

Solver::Solver(double spot, double volatility, double maturity, InterestRates* interest_rate, int time_points, int spot_points, Payoff* payoff, double theta) :
Solver::Solver(spot, new ConstantVolatility(volatility), maturity, interest_rate, time_points, spot_points, payoff, theta){

}

Matrix<double> Solver::solve(bool verbose){
    Matrix<double> res(m_spot_points+1, m_time_points+1);
    for(int i=0;i<=m_spot_points;i++)
        res.set_elem_at(i,m_time_points,m_payoff->payout(exp(m_max_space-i*m_dx)));

    for(int j=m_time_points-1;j>=0;j--){
        double vol = m_vol->getVolAt(j);
        double ir = m_ir->getIRAt(j);
        double sig2dx2 = vol*vol/(m_dx*m_dx);
        double sig2dx = vol*vol/m_dx;
        double rdx = ir/m_dx;

        Matrix<double> weights_matrix(m_spot_points+1, m_spot_points+1);

        double coeff_general_0 = -0.5*sig2dx2 - 0.25*sig2dx + 0.5*rdx;
        double coeff_general_1 = sig2dx2 + ir;
        double coeff_general_2 = -0.5*sig2dx2 + 0.25*sig2dx - 0.5*rdx;
        double coeff_top_0 = -0.5*sig2dx2 - 0.5*sig2dx + rdx + ir;
        double coeff_top_1 = sig2dx2 + 0.5*sig2dx - rdx;
        double coeff_top_2 = -0.5*sig2dx2;
        double coeff_bot_0 = -0.5*sig2dx2;
        double coeff_bot_1 = sig2dx2 - 0.5*sig2dx + rdx;
        double coeff_bot_2 = -0.5*sig2dx2 + 0.5*sig2dx - rdx + ir;

        // first row
        weights_matrix.set_elem_at(0, 0, coeff_bot_0);
        weights_matrix.set_elem_at(0, 1, coeff_bot_1);
        weights_matrix.set_elem_at(0, 2, coeff_bot_2);
        // row 1 to last - 1
        for (int i=1; i < m_spot_points; i++){
            weights_matrix.set_elem_at(i, i-1, coeff_general_0);
            weights_matrix.set_elem_at(i, i, coeff_general_1);
            weights_matrix.set_elem_at(i, i+1, coeff_general_2);
        }
        // last row
        weights_matrix.set_elem_at(m_spot_points, m_spot_points-2, coeff_top_0);
        weights_matrix.set_elem_at(m_spot_points, m_spot_points-1, coeff_top_1);
        weights_matrix.set_elem_at(m_spot_points, m_spot_points, coeff_top_2);

        Matrix<double> identity(m_spot_points+1);

        Matrix<double> right_mat = (identity + weights_matrix * m_theta * m_dt);
        Matrix<double> right_mat_inv = right_mat.inverse();
        Matrix<double> left_mat = (identity + weights_matrix * (m_theta-1) * m_dt);
        Matrix<double> transition_matrix = right_mat_inv.dot(left_mat);

        res.fill_column(j,transition_matrix.dot(res.extract_column(j+1)));
    }

    if(verbose){
        double strike(-1);
        bool is_call(false);
        try{
            VanillaCall* tempcall = dynamic_cast<VanillaCall*>(m_payoff);
            strike = tempcall->getStrike();
            is_call = true;
        } catch(const std::bad_cast& e){
            try{
            VanillaPut* tempput = dynamic_cast<VanillaPut*>(m_payoff);
            strike = tempput->getStrike();
            is_call = false;
            } catch(const std::bad_cast& e2){
                  std::cout << "WARNING: payoff is neither vanilla call or put, BS comparison is meaningless" << std::endl;
            }
        }
        for(int i=0;i<=m_spot_points;i++){
            std::cout << "S0= " << exp(m_max_space-i*m_dx);
            std::cout << "; BS price: " << dauphine::bs_price(exp(m_max_space-i*m_dx),strike, m_vol->getMeanVol(), m_matu, m_ir->getMeanIR(), is_call);
            std::cout << "; Solver price: " << res.elem_at(i,0) << std::endl;
        }
        std::cout << std::endl;
    }

    m_res = res;
    m_computed = true;
    return res;
}

double Solver::getPrice(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    return m_res.elem_at((m_spot_points+1)/2,0);
}

double Solver::getDelta(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    int midpoint = (m_spot_points+1)/2;
    return (m_res.elem_at(midpoint + 1, 0) - m_res.elem_at(midpoint - 1, 0))/(exp(m_max_space-(midpoint+1)*m_dx)-exp(m_max_space-(midpoint-1)*m_dx));
}

double Solver::getGamma(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    int midpoint = (m_spot_points+1)/2;
    double midspot = exp(m_max_space-midpoint*m_dx);
    return (m_res.elem_at(midpoint + 1, 0) + m_res.elem_at(midpoint - 1, 0) - 2*m_res.elem_at(midpoint, 0)) /
    ( (exp(m_max_space-(midpoint+1)*m_dx) - midspot) * (midspot - exp(m_max_space-(midpoint-1)*m_dx)) );
}

double Solver::getTheta(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    int midpoint = (m_spot_points+1)/2;
    return (m_res.elem_at(midpoint, 1) - m_res.elem_at(midpoint, 0))/m_dt;
}

double Solver::getVega(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    double volBump = 0.01*m_vol->getMeanVol();
    Solver tempSolver(m_spot,m_vol->bumped(volBump),m_matu,m_ir,m_time_points,m_spot_points,m_payoff,m_theta);
    tempSolver.solve();
    return (tempSolver.getPrice() - getPrice())/volBump;
}

double Solver::getRho(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    double irBump = 0.01*m_ir->getMeanIR();
    Solver tempSolver(m_spot,m_vol,m_matu,m_ir->bumped(irBump),m_time_points,m_spot_points,m_payoff,m_theta);
    tempSolver.solve();
    return (tempSolver.getPrice() - getPrice())/irBump;
}

void Solver::displayValues(){
    if(!m_computed)
        throw std::runtime_error("The solver has not been run yet!");
    std::cout << "Price: " << getPrice() << std::endl;
    std::cout << "Delta: " << getDelta() << std::endl;
    std::cout << "Gamma: " << getGamma() << std::endl;
    std::cout << "Theta: " << getTheta() << std::endl;
    std::cout << "Vega: " << getVega() << std::endl;
    std::cout << "Rho: " << getRho() << std::endl;
    std::cout << std::endl;
}
Solver::~Solver()
{

}
