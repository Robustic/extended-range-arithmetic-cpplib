#ifndef FLOAT64_EXTENDEDEXP_H_
#define FLOAT64_EXTENDEDEXP_H_

#include <vector>

namespace floatingExp2Integer
{
    class Float64ExtendedExp {
        private:
            double scnfcnd;
            std::int64_t exp;
            Float64ExtendedExp(double sicnificand, std::int64_t exponent);
            // inline void checkRuleForScale();
            inline void encode_double(double sicnificand, std::int64_t exponent);
            inline void decode_double(double& sicnificand, std::int64_t& exponent);
            inline void encode(double sgnfcnd_bits, std::int64_t exponent);
            inline void decode(std::uint64_t& sicnificand, std::int64_t& exponent, double enco);
            void printBinary(std::string message, uint64_t value) const;
            void print(std::string message, int64_t value) const;
            void print(std::string message, double value) const;
        public:
            double encoded;
            Float64ExtendedExp();
            Float64ExtendedExp(double dbl);
            double sumFloat64ExtendedExp(unsigned int n, double* vector_in);
            double decodeFloat64ExtendedExp(double data);
            void doubleToFloat64ExtendedExp(unsigned int n, const double* vector_in, double* vector_out);
            void doubleToFloat64ExtendedExp(double dbl);
            void encodedToFloat64ExtendedExp(double dbl);
            // void log2ToFloat64ExtendedExp(double log2);
            // double float64ExtendedExpToLog2() const;
            double sicnificand();
            std::int64_t exponent();
            double asDouble();

            Float64ExtendedExp& operator+=(Float64ExtendedExp z);
            Float64ExtendedExp& operator+=(double encoded);
            // Float64ExtendedExp& operator*=(Float64ExtendedExp z);
    };

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
