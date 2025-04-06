#include <cstdint>
#include <cmath>
#include <bit>
#include <vector>
#include <cstring>
#include <immintrin.h>
#include <iostream>
#include <bitset>
#include <cstdint>

#include "Float64LargeRangeNumber.h"

namespace floatingExp2Integer
{
    double Float64LargeRangeNumber::double_to_largeRangeNumber(double dbl) {
        uint64_t dbl_bits = std::bit_cast<std::uint64_t>(dbl);
        int64_t exponent = (((dbl_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
        dbl_bits &= 0x800FFFFFFFFFFFFFull;
        dbl_bits |= 0x3FF0000000000000ull;

        double sicnificand = std::bit_cast<double>(dbl_bits);
        return (double)exponent + (sicnificand - 1.0);
    }

    void Float64LargeRangeNumber::doubles_to_largeRangeNumbers(const std::vector<double>& in, std::vector<double>& out) {
        for (size_t i = 0; i < in.size(); i++) {
            out[i] = double_to_largeRangeNumber(in[i]);
        }
    }

    double Float64LargeRangeNumber::largeRangeNumber_to_double(double large_range_number) {
        double floored_exponent = std::floor(large_range_number);
        std::int64_t exponent = (std::int64_t)floored_exponent;
        double sicnificand = large_range_number - floored_exponent + 1.0;
        return sicnificand * std::pow(2.0, (int)exponent);
    }

    double Float64LargeRangeNumber::log2_to_largeRangeNumber(double log2_dbl) {
        double floored_exponent = std::floor(log2_dbl);
        double sicnificand = log2_dbl - floored_exponent + 1.0;
        return std::exp2(sicnificand) + floored_exponent;
    }

    void Float64LargeRangeNumber::log2s_to_largeRangeNumbers(const std::vector<double>& in, std::vector<double>& out) {
        for (size_t i = 0; i < in.size(); i++) {
            out[i] = log2_to_largeRangeNumber(in[i]);
        }
    }

    double Float64LargeRangeNumber::largeRangeNumber_to_log2(double large_range_number) {
        double floored_exponent = std::floor(large_range_number);
        double sicnificand = large_range_number - floored_exponent + 1.0;
        return std::log2(sicnificand) + floored_exponent;
    }

    inline void decode_fast(double& encoded_value, std::int64_t& exponent, double& sicnificand_d) {
        uint64_t dbl_bits = std::bit_cast<uint64_t>(encoded_value);
        std::int32_t cutter = ((dbl_bits >> 52) & 0x7FF) - 1023;
        std::int64_t full_mantissa = (dbl_bits & 0x000FFFFFFFFFFFFFull) | 0x0010000000000000ull;
        std::uint64_t sicnificand_int64;

        if (cutter >= 0) {
            exponent = ((std::int64_t)full_mantissa >> (52 - cutter));
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

    double Float64LargeRangeNumber::sum_largeRangeNumbers(double lrn1, double lrn2) {
        std::int64_t exponent_1;
        double sicnificand_1;

        std::int64_t exponent_2;
        double sicnificand_2;

        decode_fast(lrn1, exponent_1, sicnificand_1);
        decode_fast(lrn2, exponent_2, sicnificand_2);
        
        std::int64_t exp_diff = (std::int64_t)(exponent_1 - exponent_2);

        std::int64_t sicnificand_int64_1 = std::bit_cast<std::int64_t>(sicnificand_1);
        std::int64_t sicnificand_int64_2 = std::bit_cast<std::int64_t>(sicnificand_2);

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

    double Float64LargeRangeNumber::sum_largeRangeNumbers(std::vector<double>& large_range_numbers) {
        if (large_range_numbers.size() < 16) {
            double sum_small = large_range_numbers[0] + (double)large_range_numbers[0];
            for (size_t k = 1; k < large_range_numbers.size(); k++) {
                sum_small = sum_largeRangeNumbers(sum_small, large_range_numbers[k]);
            }
            return sum_small;
        }

        __m512d encoded_values = _mm512_loadu_pd(&large_range_numbers[0]);

        __m512d integer_decimals = _mm512_floor_pd(encoded_values);
        __m512i exponent_sum = _mm512_cvtpd_epi64(integer_decimals);
        __m512d one_d = _mm512_set1_pd(1.0);
        __m512d two_d = _mm512_set1_pd(2.0);
        __m512i one_i = _mm512_set1_epi64(1);
        __m512d half_d = _mm512_set1_pd(0.5);
        __m512d zero_d = _mm512_setzero_pd();
        __m512d sicnificand_sum = encoded_values - integer_decimals + one_d;

        size_t i = 8;
        for (i = 8; i + 7 < large_range_numbers.size(); i += 8) {
            encoded_values = _mm512_loadu_pd(&large_range_numbers[i]);

            integer_decimals = _mm512_floor_pd(encoded_values);
            __m512i exponent = _mm512_cvttpd_epi64(integer_decimals);
            __m512d sicnificand = encoded_values - integer_decimals + one_d;

            __m512i exp_diff = exponent_sum - exponent;

            __m512i zero_i = _mm512_setzero_si512();

            __m512i sicnificand_bits = _mm512_castpd_si512(sicnificand);
            sicnificand_bits = exp_diff < zero_i ? sicnificand_bits : sicnificand_bits - (exp_diff << 52);
            sicnificand = _mm512_castsi512_pd(sicnificand_bits);

            __m512i m1_i = _mm512_set1_epi64(-1);
            __m512i sicnificand_sum_bits = _mm512_castpd_si512(sicnificand_sum);
            sicnificand_sum_bits = exp_diff >= zero_i ? sicnificand_sum_bits : sicnificand_sum_bits - ((m1_i * exp_diff) << 52);
            sicnificand_sum = _mm512_castsi512_pd(sicnificand_sum_bits);

            __m512i p63_i = _mm512_set1_epi64(63);
            __m512i m63_i = _mm512_set1_epi64(-63);
            sicnificand_sum = exp_diff < m63_i ? zero_d : sicnificand_sum;
            sicnificand = exp_diff > p63_i ? zero_d : sicnificand;

            exponent_sum = exp_diff < zero_i ? exponent : exponent_sum;

            sicnificand_sum = sicnificand_sum + sicnificand;

            exponent_sum = sicnificand_sum >= two_d ? exponent_sum + one_i : exponent_sum;
            sicnificand_sum = sicnificand_sum >= two_d ? sicnificand_sum * half_d : sicnificand_sum;
        }

        double sum = sicnificand_sum[0] + (double)exponent_sum[0] - 1;

        for (int k = 1; k < 8; k++) {
            sum = sum_largeRangeNumbers(sum, sicnificand_sum[k] + (double)exponent_sum[k] - 1);
        }

        for (i = i; i < large_range_numbers.size(); i++) {
            sum = sum_largeRangeNumbers(sum, large_range_numbers[i]);
        }

        return sum;
    }

    void printBinary(std::string message, uint64_t value) {
        if (message.length() > 0) {
            std::cout << message << std::endl;
        }
        for (int i = 56; i >= 0; i -= 8) {
            std::bitset<8> block((value >> i) & 0xFF); // Extract 8 bits
            std::cout << block << " ";                 // Print block with a space
        }
        std::cout << std::endl;
    }

    void print(std::string message, int64_t value) {
        std::cout << message << value << std::endl;
    }

    void print(std::string message, double value) {
        std::cout << message << value << std::endl;
    }
}

