// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "closed_form.hpp"
#include "Solver.h"
#include "Matrix.h"
#include "Solver_matrix.h"

#include <iostream>



void solver()
{
	double spot = 100.;
	double volatility = 0.16;
	double maturity = 0.25;
	float interest_rate = 0.01;
	int time_points = 20;
	int spot_points = 20;
	double theta = 0.5;
	double strike = 100.;
	bool iscall = true;

	Solver mysolver = Solver(spot, volatility, maturity, interest_rate, time_points, spot_points, theta);
	std::vector<std::vector<double> > res = mysolver.solve_BS_theta0(strike, iscall);

	for (int i = 0; i <= time_points; i++)
	{
		std::cout << "[";
		for (int j = 0; j <= spot_points; j++)
		{
			std::cout << res[i][j] << ", ";
		}
		std::cout << "]" << std::endl;
	}

	std::cout << "price = " << res[int(time_points / 2)][0] << std::endl;

	double price = dauphine::bs_price(spot, strike, volatility, maturity, interest_rate, iscall);


	std::cout << "price bs = " << price << std::endl;

}

void test()
{
	double spot = 100.;
	double volatility = 0.16;
	double maturity = 0.25;
	float interest_rate = 0.04;
	int points = 16;
	double theta = 0.5;
	double strike = 100.;
	bool is_call = true;


	Solver_matrix mysolver = Solver_matrix(spot, volatility, maturity, interest_rate, points, theta);
	Matrix<double> M = mysolver.solving(strike, is_call);

}

void test_matrice()
{
	Matrix<double> M = Matrix<double>(3);

	//Matrix<double> N = Matrix<double>(2);
	M.set_elem_at(1, 1, 8.5);
	M.set_elem_at(0, 1, 1.5);
	M.set_elem_at(0, 2, 0.3);
	M.set_elem_at(2, 1, -5);
	M.set_elem_at(2, 2, 8.5);
	M.display(); 
	
	
	Matrix<double> N = M.inverse();
	N.display();


	Matrix<double> P = N.dot(M);
	P.display();

	std::cout << "over" << std::endl;
	
}


class mere
{
public:
	mere();
	mere(int j);
	~mere(){}
protected:
	int m_i;
};

class fille : mere
{
public:
	fille();
	~fille(){}
private:
	int w;
};

mere::mere(int j) 
{
	m_i = j; 
}
mere::mere()
{
	m_i = 10;
}
fille::fille()
	:mere(15)
{
	std::cout << m_i << std::endl;
}
void test_class()
{
	fille ma_classe = fille();
}

int main()
{
	test();
	//test_class();
	//test_matrice();
    return 0;
}

