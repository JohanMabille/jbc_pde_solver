#include "InterestRates.h"

InterestRates::InterestRates(){

}

InterestRates::~InterestRates(){

}

ConstantIR::ConstantIR(double ir) : m_ir(ir){

}

double ConstantIR::getIRAt(int j){
    return m_ir;
}

double ConstantIR::getMeanIR(){
    return m_ir;
}

ConstantIR* ConstantIR::bumped(double bump){
    return (new ConstantIR(m_ir+bump));
}

ConstantIR::~ConstantIR(){

}

GeneralIR::GeneralIR(double irs[], int points_available, int points_grid) : m_size(points_available), m_points_grid(points_grid){
    m_irs = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_irs[i] = irs[i];
}

GeneralIR::GeneralIR(double irs[], int points_available) : m_size(points_available), m_points_grid(points_available){
    m_irs = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_irs[i] = irs[i];
}

GeneralIR::GeneralIR(std::vector<double> irs, int points_grid) : m_size(irs.size()), m_points_grid(points_grid){
    m_irs = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_irs[i] = irs[i];
}

GeneralIR::GeneralIR(std::vector<double> irs) : m_size(irs.size()), m_points_grid(irs.size()){
    m_irs = new double[m_size];
    for(int i=0;i<m_size;i++)
        m_irs[i] = irs[i];
}

double GeneralIR::getIRAt(int j){
    int closest_point = int(double(m_size)*j/m_points_grid+0.499);
    return m_irs[closest_point];
}

double GeneralIR::getMeanIR(){
    double tot=0.;
    for(int i=0;i<m_points_grid;i++)
        tot += getIRAt(i);
    return tot/m_points_grid;
}

GeneralIR* GeneralIR::bumped(double bump){
    double *tempir = new double[m_size];
    for(int i=0;i<m_size;i++)
        tempir[i] = m_irs[i] + bump;

    return (new GeneralIR(tempir,m_size,m_points_grid));
}

GeneralIR::~GeneralIR(){
    delete[] m_irs;
}
