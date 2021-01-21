#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>

class Solver{
    public:
        explicit Solver(double spot, double volatility, double maturity, float interest_rate, int time_points, int spot_points, double theta);

        template<typename T>
        static void print_vector_array(std::vector<std::vector<T> > to_print){
            for(auto line: to_print){
                std::cout << "[";
                std::copy(line.begin(), line.end(), std::ostream_iterator<T>(std::cout, ", "));
                std::cout << "]" << std::endl;
            }
        }
        std::vector<std::vector<double> > solve_BS(double strike, bool is_call);
        std::vector<std::vector<double> > solve_BS_theta0(double strike, bool is_call);
        virtual ~Solver();

    protected:
        double m_spot, m_vol, m_matu, m_ir, m_theta;
        double m_min_time, m_max_time, m_dt;
        double m_std_dev;
        double m_min_space, m_max_space, m_dx;
        int m_spot_points, m_time_points;

        double m_coeffs[3][2];
        double m_coeffs_theta0[3];
        double m_coeffs_theta0_edge[3];

        std::vector<std::vector<double> > m_res;

        double compute_vertex(double values[][2], int di, int dn);
        double compute_vertex_theta0(double values[]);
        double compute_vertex_theta0_edge(double values[]);
};

#endif // SOLVER_H
