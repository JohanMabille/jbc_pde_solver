#ifndef PAYOFF_H
#define PAYOFF_H

#include <functional>

class Payoff
{
    public:
        // entity semantics: the constructor can be protected
        Payoff();
        virtual double payout(double spot) const = 0;
        virtual ~Payoff();

        // Entity semantics:
        Payoff(const Payoff&) = delete;
        Payoff& operator=(const Payoff&) = delete;
        Payoff(Payoff&&) = delete;
        Payoff& operator=(Payoff&&) = delete;

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
        // Consider passing by cont ref to avoid copy
        GeneralPayoff(std::function<double(double)> fpayout);
        double payout(double spot) const;
        virtual ~GeneralPayoff();

    protected:
        std::function<double(double)> m_payout;
};

#endif // PAYOFF_H
