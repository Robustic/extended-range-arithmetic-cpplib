#ifndef FUKUSHIMA_H_
#define FUKUSHIMA_H_

#include <vector>

namespace floatingExp2Integer
{
    class Fukushima {
        private:
            double scnfcnd;
            std::int64_t exp;
            inline void xnorm();
            // Float64ExtendedExp(double sicnificand, std::int64_t exponent);
            // inline void checkRuleForScale();
        public:
            Fukushima();
            Fukushima(double dbl);
            Fukushima(double dbl, std::int64_t exponent);
            Fukushima(const std::vector<floatingExp2Integer::Fukushima>& vector);
            void doubleToFukushima(double dbl);
            // void log2ToFloat64ExtendedExp(double log2);
            // double float64ExtendedExpToLog2() const;
            double sicnificand() { return scnfcnd; };
            std::int64_t exponent() { return exp; };
            double asDouble() const;

            Fukushima& operator+=(Fukushima z);
            // Float64ExtendedExp& operator*=(Float64ExtendedExp z);
    };

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
