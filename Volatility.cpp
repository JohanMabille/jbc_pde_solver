#include "Volatility.h"

Volatility::Volatility(){

}

Volatility::~Volatility(){

}

ConstantVolatility::ConstantVolatility(double sigma) : m_vol(sigma){

}

double ConstantVolatility::getVolAt(int j){
    return m_vol;
}

double ConstantVolatility::getMeanVol(){
    return m_vol;
}

ConstantVolatility* ConstantVolatility::bumped(double bump){
    return (new ConstantVolatility(m_vol+bump));
}

ConstantVolatility::~ConstantVolatility(){

}

GeneralVolatility::GeneralVolatility(double vols[], int points_available, int points_grid) : m_size(points_available), m_points_grid(points_grid){
    m_vols = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_vols[i] = vols[i];
}

GeneralVolatility::GeneralVolatility(double vols[], int points_available) : m_size(points_available), m_points_grid(points_available){
    m_vols = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_vols[i] = vols[i];
}

GeneralVolatility::GeneralVolatility(std::vector<double> vols, int points_grid) : m_size(vols.size()), m_points_grid(points_grid){
    m_vols = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_vols[i] = vols[i];
}

GeneralVolatility::GeneralVolatility(std::vector<double> vols) : m_size(vols.size()), m_points_grid(vols.size()){
    m_vols = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_vols[i] = vols[i];
}

double GeneralVolatility::getVolAt(int j){
    int closest_point = int(double(m_size)*j/m_points_grid+0.499);
    return m_vols[closest_point];
}

double GeneralVolatility::getMeanVol(){
    double tot=0.;
    for(int i=0;i<m_points_grid;i++)
        tot += getVolAt(i);
    return tot/m_points_grid;
}

GeneralVolatility* GeneralVolatility::bumped(double bump){
    double *tempvol = new double[m_size];
    for(int i=0;i<m_size;i++)
        tempvol[i] = m_vols[i] + bump;

    return (new GeneralVolatility(tempvol,m_size,m_points_grid));
}

GeneralVolatility::~GeneralVolatility(){
    delete[] m_vols;
}
