

#include "Solver_matrix.hpp"
#include "closed_form.hpp"

void show_both(Matrix<double> M1, Matrix<double> M2)
{
	int size = M1.get_rows();
	for (int i = 0; i < size; i++)
	{
		std::cout << "[" << M1.elem_at(i, 0) << " ; " << M2.elem_at(i, 0) << "]" << std::endl;
	}
}


Solver_matrix::Solver_matrix(double spot, double volatility, double maturity, float interest_rate, int points, double theta)
	:m_spot(spot),
	m_vol(volatility),
	m_matu(maturity),
	m_ir(interest_rate),
	m_points(points),
	m_theta(theta)

{
	M_vol = Matrix<double>(m_points, 1);
	m_min_time = 0.;
	m_max_time = m_matu;
	m_dt = (m_max_time - m_min_time) / m_points;
	m_std_dev = m_vol * sqrt(m_matu);
	m_min_space = log(m_spot) - 3 * m_std_dev;
	m_max_space = log(m_spot) + 3 * m_std_dev;

	m_dx = (m_max_space - m_min_space) / (m_points - 1);

	// initialize the matrix M 
	Matrix<double>* M = new Matrix<double>(m_points, m_points);
	Res = new Matrix<double>(m_points, m_points);
	Spot = Matrix<double>(m_points, 1);
	BS_Price = Matrix<double>(m_points, 1);

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
	//M->display();

	Matrix<double>* Identity = new Matrix<double>(m_points);
	Matrix<double> WN_1 = (*Identity + (*M) * m_theta * m_dt).inverse();
	Matrix<double> WN1 = *Identity + (*M) * (m_theta - 1) * m_dt;
	WN_1.display();
	WN1.display();
	W = WN_1.dot(WN1);
}


Solver_matrix::Solver_matrix(double spot, Matrix<double> volatility, double maturity, Matrix<float> interest_rate, int points, double theta)
	:m_spot(spot), 
	M_vol(volatility),
	m_matu(maturity),
	M_ir(interest_rate),
	m_points(points),
	m_theta(theta)

{
	m_min_time = 0.;
	m_max_time = m_matu;
	m_dt = (m_max_time - m_min_time) / m_points;
	m_std_dev = M_vol.elem_at(0, 0) * sqrt(m_matu);
	m_min_space = log(m_spot) - 3 * m_std_dev;
	m_max_space = log(m_spot) + 3 * m_std_dev;

	m_dx = (m_max_space - m_min_space) / (m_points - 1);

	Res = new Matrix<double>(m_points, m_points);
	Spot = Matrix<double>(m_points, 1);
	BS_Price = Matrix<double>(m_points, 1);

	// Matrix M can't be initialize since it is now changing at every time step
	
}

Solver_matrix::~Solver_matrix()
{

}

void Solver_matrix::solve(double strike, bool is_call)
{
	m_strike = strike;
	m_is_call = is_call;
	if (M_vol.isnull())
		solve_constant(strike, is_call);
	else
		solve_non_constant(strike, is_call);
}

void Solver_matrix::solve_constant(double strike, bool is_call)
{
	//Matrix<double>* Identity = new Matrix<double>(m_points);
	

	//Matrix<double> WN_1 = (*Identity + (*M) * m_theta * m_dt).inverse();
	//Matrix<double> WN1 = *Identity + (*M) * (m_theta - 1) * m_dt;

	// Compute last column = payoff
	for (int i = 0; i < m_points; i++)
	{
		Res->set_elem_at(i, m_points-1, dauphine::vanilla_payoff(exp(m_max_space - i*m_dx), strike, is_call));
		Spot.set_elem_at(i, 0, exp(m_max_space - i*m_dx));
		BS_Price.set_elem_at(i, 0, dauphine::bs_price(Spot.elem_at(i, 0), strike, m_vol, m_matu, (float) m_ir, true));
	}
	
	// Recursion
	for (int j = m_points - 2; j >= 0; j--)
	{
		//Res->fill_column(j, WN_1.dot(WN1).dot(Res->column(j + 1)));
		Res->fill_column(j, W.dot(Res->column(j + 1)));
	}


	//Res->display_end();
}


void Solver_matrix::solve_non_constant(double strike, bool is_call)
{
	Matrix<double>* Identity = new Matrix<double>(m_points);

	Matrix<double>* M = new Matrix<double>(m_points, m_points);

	//last column is the payoff
	for (int i = 0; i < m_points; i++)
	{
		Res->set_elem_at(i, m_points - 1, dauphine::vanilla_payoff(exp(m_max_space - i*m_dx), strike, is_call));
		Spot.set_elem_at(i, 0, exp(m_max_space - i*m_dx));
	}

	// Recursion
	for (int j = m_points - 2; j >= 0; j--)
	{
		m_ir = M_ir.elem_at(j, 0);
		m_vol = M_vol.elem_at(j, 0);
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
			M->set_elem_at(i, i - 1, m_coefm[0]);
			M->set_elem_at(i, i, m_coefm[1]);
			M->set_elem_at(i, i + 1, m_coefm[2]);
		}
		// last row
		M->set_elem_at(m_points - 1, m_points - 3, m_coefmn[0]);
		M->set_elem_at(m_points - 1, m_points - 2, m_coefmn[1]);
		M->set_elem_at(m_points - 1, m_points - 1, m_coefmn[2]);
		
		// We can then compute the matrix for the recursion
		Matrix<double> WN_1 = (*Identity + (*M) * m_theta * m_dt).inverse();
		Matrix<double> WN1 = *Identity + (*M) * (m_theta - 1) * m_dt;


		Res->fill_column(j, WN_1.dot(WN1).dot(Res->column(j + 1)));
	}

}

void Solver_matrix::display_res()
{
	Res->display();
}

void Solver_matrix::display_end_res()
{
	if (m_points < 15)
		display_res();
	else
	{
		Res->display_end();
	}
}


void Solver_matrix::compare_with_BS()
{
	show_both(Res->column(0), BS_Price);
}

double Solver_matrix::price()
{
	if (m_points % 2 == 1)
	{
		//std::cout << Spot.elem_at((m_points - 1) / 2, 0) << std::endl;
		return Res->elem_at((m_points - 1) / 2, 0);
	}
	else
	{	// On dont have the price for exactly the spot, so we compute the mean
		double pd = Res->elem_at(m_points / 2, 0);
		double pu = Res->elem_at((m_points / 2) - 1, 0);
		//std::cout << "price " << pd << std::endl;
		//std::cout << pu << std::endl;
		return (pd + pu) / 2;
	}
}

double Solver_matrix::delta()
{
	double pu, pd, spotu, spotd;

	if (m_points % 2 == 1)
	{
		spotu = Spot.elem_at(((m_points - 1) / 2) - 1, 0);
		spotd = Spot.elem_at(((m_points - 1) / 2) + 1, 0);
		pu = Res->elem_at(((m_points - 1) / 2) - 1, 0);
		pd = Res->elem_at(((m_points - 1) / 2) + 1, 0);
	}
	else
	{	
		spotu = Spot.elem_at((m_points / 2) - 1, 0);
		spotd = Spot.elem_at(m_points / 2, 0);
		pu = Res->elem_at((m_points / 2) - 1, 0);
		pd = Res->elem_at(m_points / 2, 0);
	}
	//std::cout << pu << ", " << pd << std::endl;
	//std::cout << spotu << ", " << spotd << std::endl;
	return (pu - pd) / (spotu - spotd);
}

double Solver_matrix::gamma()
{
	double spotu = Spot.elem_at(((m_points - 1) / 2) - 1, 0);
	double spotd = Spot.elem_at(((m_points - 1) / 2) + 1, 0);
	double dx = (spotu - spotd) / 2;

	double pu, p, pd;

	if (m_points % 2 == 1)
	{
		pu = Res->elem_at(((m_points - 1) / 2) - 1, 0);
		p = Res->elem_at(((m_points - 1) / 2), 0);
		pd = Res->elem_at(((m_points - 1) / 2) + 1, 0);
		
		
	}
	else
	{
		pu = (Res->elem_at((m_points / 2) - 2, 0) + Res->elem_at((m_points / 2) - 1, 0)) / 2;
		p = (Res->elem_at(m_points / 2, 0) + Res->elem_at((m_points / 2) - 1, 0)) / 2;
		pd = (Res->elem_at(m_points / 2, 0) + Res->elem_at((m_points / 2) + 1, 0)) / 2;

	}
	//std::cout << "p " << pu << ", " << p << ", " << pd << std::endl;
	return (pu - 2 * p + pd) / (dx * dx);
}

double Solver_matrix::theta()
{
	double pt, pt1;
	
	if (m_points % 2 == 1)
	{
		pt = Res->elem_at((m_points - 1) / 2, 0);
		pt1 = Res->elem_at((m_points - 1) / 2, 1);
	}
	else
	{	
		pt = (Res->elem_at(m_points / 2, 0) + Res->elem_at((m_points / 2) - 1, 0)) / 2;
		pt1 = (Res->elem_at(m_points / 2, 1) + Res->elem_at((m_points / 2) - 1, 1)) / 2;
	}
	return (pt1 - pt) / m_dt;
}

double Solver_matrix::vega()
{
	double start_price = price();
	double shift = 0.01;
	double end_price;
	if (M_vol.isnull())
	{
		Solver_matrix solver_shifted(m_spot, m_vol + shift, m_matu, m_ir, m_points, m_theta);
		solver_shifted.solve(m_strike, m_is_call);
		end_price = solver_shifted.price();
	}
	else
	{
		Solver_matrix solver_shifted(m_spot, M_vol + shift, m_matu, M_ir, m_points, m_theta);
		solver_shifted.solve(m_strike, m_is_call);
		end_price = solver_shifted.price();
	}
	return (end_price - start_price) / shift;
}