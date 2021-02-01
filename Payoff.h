#ifndef PAYOFF_H
#define PAYOFF_H

#include <functional>

class Payoff
{
    public:
        Payoff();
        virtual double payout(double spot) const = 0;
        virtual ~Payoff();

    protected:

};

class VanillaCall : public Payoff
{
    public:
        VanillaCall(double strike);
        double payout(double spot) const;
        double getStrike() const;
        virtual ~VanillaCall();

    protected:
        double m_strike;
};

class VanillaPut : public Payoff
{
    public:
        VanillaPut(double strike);
        double payout(double spot) const;
        double getStrike() const;
        virtual ~VanillaPut();

    protected:
        double m_strike;
};

class GeneralPayoff : public Payoff
{
    public:
        GeneralPayoff(std::function<double(double)> fpayout);
        double payout(double spot) const;
        virtual ~GeneralPayoff();

    protected:
        std::function<double(double)> m_payout;
};

#endif // PAYOFF_H
