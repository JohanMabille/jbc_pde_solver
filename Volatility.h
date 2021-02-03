#ifndef VOLATILITY_H
#define VOLATILITY_H
#include <vector>

class Volatility{
    public:
        // Entity semantics: declare copy and move constructors
        // and assign operators as deleted functions
        // Trivial constructor and destructor can be declared
        // with = default so that their implementation can be
        // omitter in the source file.
        Volatility();
        virtual double getVolAt(int j) = 0;
        virtual double getMeanVol() = 0;
        virtual Volatility* bumped(double bump) = 0;
        virtual ~Volatility();

    protected:
};


class ConstantVolatility: public Volatility{
    public:
        explicit ConstantVolatility(double sigma);
        double getVolAt(int j);
        double getMeanVol();
        ConstantVolatility* bumped(double bump);
        virtual ~ConstantVolatility();

    protected:
        double m_vol;
};

class GeneralVolatility: public Volatility{
    public:
        // Consider passing the arrays / vector by constant references
        // to avoid copies
        explicit GeneralVolatility(double vols[], int points_available);
        explicit GeneralVolatility(double vols[], int points_available, int points_grid);
        explicit GeneralVolatility(std::vector<double> vols);
        explicit GeneralVolatility(std::vector<double> vols, int points_grid);
        double getVolAt(int j);
        double getMeanVol();
        GeneralVolatility* bumped(double bump);
        virtual ~GeneralVolatility();

    protected:
        // why not using std::vector?
        double* m_vols;
        int m_size;
        int m_points_grid;
};

#endif // VOLATILITY_H
