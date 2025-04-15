#ifndef FUKUSHIMA_H_
#define FUKUSHIMA_H_

#include <vector>

namespace floatingExp2Integer
{
    class Fukushima {
        private:
            double scnfcnd;
            int64_t exp;

            Fukushima(double dbl, int64_t exponent);
            inline void xnorm();

            static constexpr int64_t EXP_MULTIPLIER = 960;
            static constexpr double BIG = 0x1p960;
            static constexpr double BIGI = 0x1p-960;
            static constexpr double BIGS = 0x1p480;
            static constexpr double BIGSI = 0x1p-480;
        public:
            Fukushima();
            Fukushima(double dbl);

            void log2_to(const double from);
            static void log2s_to(const std::vector<double>& from, std::vector<floatingExp2Integer::Fukushima>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            void sum(const std::vector<floatingExp2Integer::Fukushima>& vector);
            void multiply(const std::vector<floatingExp2Integer::Fukushima>& vector);

            double as_double() const;
            double as_log2() const;

            double sicnificand() const { return scnfcnd; };
            int64_t exponent() const { return exp; };

            Fukushima& operator+=(Fukushima z);
            Fukushima& operator*=(Fukushima z);
    };

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
