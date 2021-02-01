#include "Payoff.h"
#include <math.h>
#include <algorithm>

Payoff::Payoff(){

}

Payoff::~Payoff(){

}

VanillaCall::VanillaCall(double strike) : m_strike(strike){
}

double VanillaCall::payout(double spot) const{
    return std::max(spot - m_strike , 0.);
}

double VanillaCall::getStrike() const{
    return m_strike;
}

VanillaCall::~VanillaCall(){

}

VanillaPut::VanillaPut(double strike) : m_strike(strike){

}

double VanillaPut::payout(double spot) const{
    return std::max(m_strike - spot, 0.);
}

double VanillaPut::getStrike() const{
    return m_strike;
}

VanillaPut::~VanillaPut(){

}

GeneralPayoff::GeneralPayoff(std::function<double(double)> fpayout): m_payout(fpayout){

}

double GeneralPayoff::payout(double spot) const{
    return m_payout(spot);
}

GeneralPayoff::~GeneralPayoff(){

}
