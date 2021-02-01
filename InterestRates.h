#ifndef INTERESTRATES_H
#define INTERESTRATES_H
#include <vector>

class InterestRates{
    public:
        InterestRates();
        virtual double getIRAt(int j) = 0;
        virtual double getMeanIR() = 0;
        virtual InterestRates* bumped(double bump) = 0;
        virtual ~InterestRates();

    protected:

};

class ConstantIR: public InterestRates{
    public:
        explicit ConstantIR(double ir);
        double getIRAt(int j);
        double getMeanIR();
        ConstantIR* bumped(double bump);
        virtual ~ConstantIR();

    protected:
        double m_ir;
};

class GeneralIR: public InterestRates{
    public:
        explicit GeneralIR(double irs[], int points_available);
        explicit GeneralIR(double irs[], int points_available, int points_grid);
        explicit GeneralIR(std::vector<double> irs);
        explicit GeneralIR(std::vector<double> irs, int points_grid);
        double getIRAt(int j);
        double getMeanIR();
        GeneralIR* bumped(double bump);
        virtual ~GeneralIR();

    protected:
        double* m_irs;
        int m_size;
        int m_points_grid;
};

#endif // INTERESTRATES_H
