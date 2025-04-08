#ifndef FUKUSHIMA_H_
#define FUKUSHIMA_H_

#include <vector>

namespace floatingExp2Integer
{
    class Fukushima {
        private:
            inline void xnorm();

            static constexpr int EXP_MULTIPLIER = 960;
            static constexpr double BIG = 0x1p960;
            static constexpr double BIGI = 0x1p-960;
            static constexpr double BIGS = 0x1p480;
            static constexpr double BIGSI = 0x1p-480;
        public:
            double scnfcnd;
            int64_t exp;
            Fukushima();
            Fukushima(double dbl);
            Fukushima(double dbl, int64_t exponent);
            void sum(const std::vector<floatingExp2Integer::Fukushima>& vector);
            void multiply(const std::vector<floatingExp2Integer::Fukushima>& vector);
            void doubleToFukushima(double dbl);
            // void log2ToFukushima(double log2);
            double fukushimaToLog2() const;
            double sicnificand() { return scnfcnd; };
            int64_t exponent() { return exp; };
            double asDouble() const;

            Fukushima& operator+=(Fukushima z);
            Fukushima& operator*=(Fukushima z);
    };

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
