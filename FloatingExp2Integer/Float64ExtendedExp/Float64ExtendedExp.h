#ifndef FLOAT64_EXTENDEDEXP_H_
#define FLOAT64_EXTENDEDEXP_H_

#include <vector>

namespace floatingExp2Integer
{
    class Float64ExtendedExp {
        private:
            double encoded;
            Float64ExtendedExp(double sicnificand, std::int64_t exponent);
            // inline void checkRuleForScale();
            inline void encode_double(double sicnificand, std::int64_t exponent);
            inline void decode_double(double& sicnificand, std::int64_t& exponent);
            inline void encode(double sgnfcnd_bits, std::int64_t exponent);
            inline void decode(std::uint64_t& sicnificand, std::int64_t& exponent);
            void printBinary(std::string message, uint64_t value) const;
            void print(std::string message, int64_t value) const;
            void print(std::string message, double value) const;
            inline void convert_double_to_uint64(__m128d encoded, uint64_t* exponent, uint64_t* significand);
            inline void convert_double_to_uint64(__m512d encoded_values, __m512i exponent, __m512i sicnificand);
        public:
            Float64ExtendedExp();
            Float64ExtendedExp(double dbl);
            double sumFloat64ExtendedExp(unsigned int n, double* vector_in);
            double decodeFloat64ExtendedExp(double data);
            void doubleToFloat64ExtendedExp(unsigned int n, const double* vector_in, double* vector_out);
            void doubleToFloat64ExtendedExp(double dbl);
            // void log2ToFloat64ExtendedExp(double log2);
            // double float64ExtendedExpToLog2() const;
            double sicnificand();
            std::int64_t exponent();
            double asDouble();

            Float64ExtendedExp& operator+=(Float64ExtendedExp z);
            // Float64ExtendedExp& operator*=(Float64ExtendedExp z);
    };

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
