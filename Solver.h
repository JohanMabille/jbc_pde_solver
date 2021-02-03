#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Volatility.h"
#include "InterestRates.h"
#include "Payoff.h"

class Solver{
    public:
        explicit Solver(double spot, Volatility* volatility, double maturity, InterestRates* interest_rate, int time_points, int spot_points, Payoff* payoff, double theta);
        explicit Solver(double spot, Volatility* volatility, double maturity, float interest_rate, int time_points, int spot_points, Payoff* payoff, double theta);
        explicit Solver(double spot, double volatility, double maturity, InterestRates* interest_rate, int time_points, int spot_points, Payoff* payoff, double theta);
        explicit Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, Payoff* payoff, double theta);

        Matrix<double> solve(bool verbose=false);
        double getPrice();
        double getDelta();
        double getGamma();
        double getTheta();
        double getVega();
        double getRho();

        void displayValues();

        // Why virtual?
        virtual ~Solver();

    protected:
        double m_spot, m_matu, m_theta;
        double m_min_time, m_max_time, m_dt;
        double m_std_dev;
        double m_min_space, m_max_space, m_dx;
        int m_spot_points, m_time_points;

        Volatility* m_vol;
        InterestRates* m_ir;
        Payoff* m_payoff;

        Matrix<double> m_res;

        bool m_computed;

};

#endif // SOLVER_H
