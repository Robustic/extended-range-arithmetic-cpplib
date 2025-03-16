#ifndef FLOAT64_EXTENDEDEXP_H_
#define FLOAT64_EXTENDEDEXP_H_

#include <vector>

namespace floatingExp2Integer
{
    class Float64ExtendedExp {
        private:
            double encoded;
            Float64ExtendedExp(double sicnificand, std::int64_t exponent);
            inline void checkRuleForScale();
            inline void scale(double sicnificand, std::int64_t exponent);
        public:
            Float64ExtendedExp();
            Float64ExtendedExp(double dbl);
            Float64ExtendedExp(std::vector<floatingExp2Integer::Float64ExtendedExp>& vector);
            void doubleToFloat64ExtendedExp(double dbl);
            void log2ToFloat64ExtendedExp(double log2);
            double float64ExtendedExpToLog2() const;
            double sicnificand();
            std::int64_t exponent();
            double asDouble() const;

            Float64ExtendedExp& operator+=(Float64ExtendedExp z);
            Float64ExtendedExp& operator*=(Float64ExtendedExp z);
    };

    Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
