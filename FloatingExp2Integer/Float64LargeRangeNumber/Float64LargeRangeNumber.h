#ifndef FLOAT64_LARGERANGENUMBER_H_
#define FLOAT64_LARGERANGENUMBER_H_

#include <vector>
#include <iostream>

namespace floatingExp2Integer
{
    class Float64LargeRangeNumber {
        private:
            inline static void decode_fast_int(double& encoded_value, int64_t& exponent, double& sicnificand_d) {
                uint64_t dbl_bits = std::bit_cast<uint64_t>(encoded_value);
                int32_t cutter = ((dbl_bits >> 52) & 0x7FF) - 1023;
                int64_t full_mantissa = (dbl_bits & 0x000FFFFFFFFFFFFFull) | 0x0010000000000000ull;
                uint64_t sicnificand_int64;

                if (cutter >= 0) {
                    exponent = ((int64_t)full_mantissa >> (52 - cutter));
                    sicnificand_int64 = (dbl_bits << cutter) & 0x000FFFFFFFFFFFFFull;
                }
                else {
                    exponent = 0;
                    sicnificand_int64 = full_mantissa >> (-1 * cutter);
                }

                exponent = exponent * (encoded_value < 0.0 ? -1LL : 1LL);
                exponent = (encoded_value < 0.0 && sicnificand_int64 > 0ULL) ? exponent - 1 : exponent;

                sicnificand_int64 = encoded_value < 0.0 ? (0x0020000000000000ull - sicnificand_int64) : sicnificand_int64;
                sicnificand_int64 = sicnificand_int64 | 0x3FF0000000000000ull;
                sicnificand_d = std::bit_cast<double>(sicnificand_int64);
            }
        public:
            double encoded;
            Float64LargeRangeNumber() { encoded = 1.0; }
            Float64LargeRangeNumber(double enc) { encoded = enc; }
            static double double_to(double dbl);
            static void doubles_to(const std::vector<double>& in, std::vector<double>& out);
            static double as_double(double dbl);
            double as_double();

            void log2_to(double dbl);
            static void log2s_to(const std::vector<double>& in, std::vector<floatingExp2Integer::Float64LargeRangeNumber>& out) {
                for (size_t i = 0; i < out.size(); i++) {
                    out[i].log2_to(in[i]);
                }
            }
            static double as_log2(double dbl);

            inline static double sum(double lrn1, double lrn2) {
                //if (lrn1 - lrn2 > 0x1p53) {
                //    return lrn1;
                //}
                //else if (lrn2 - lrn1 > 0x1p53) {
                //    return lrn2;
                //}

                int64_t exponent_1;
                double sicnificand_1;

                int64_t exponent_2;
                double sicnificand_2;

                decode_fast_int(lrn1, exponent_1, sicnificand_1);
                decode_fast_int(lrn2, exponent_2, sicnificand_2);

                int64_t exp_diff = (int64_t)(exponent_1 - exponent_2);

                int64_t sicnificand_int64_1 = std::bit_cast<int64_t>(sicnificand_1);
                int64_t sicnificand_int64_2 = std::bit_cast<int64_t>(sicnificand_2);

                if (exp_diff > 0) {
                    if (exp_diff > 63) {
                        return lrn1;
                    }
                    sicnificand_int64_2 -= exp_diff << 52;
                }
                else if (exp_diff < 0) {
                    if (exp_diff < -63) {
                        return lrn2;
                    }
                    sicnificand_int64_1 -= (-exp_diff) << 52;
                    exponent_1 = exponent_2;
                }

                sicnificand_1 = std::bit_cast<double>(sicnificand_int64_1);
                sicnificand_2 = std::bit_cast<double>(sicnificand_int64_2);

                sicnificand_1 += sicnificand_2;

                if (sicnificand_1 >= 2.0) {
                    sicnificand_1 *= 0.5;
                    return (double)exponent_1 + sicnificand_1;
                }
                else {
                    return (double)exponent_1 + sicnificand_1 - 1.0;
                }
            }
            static double multiply(double lrn1, double lrn2);

            static double sum(std::vector<double>& lrns);
            static double multiply(std::vector<double>& lrns);
    };
}

#endif
