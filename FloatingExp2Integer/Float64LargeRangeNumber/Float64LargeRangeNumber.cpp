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
    double Float64LargeRangeNumber::double_to(double dbl) {
        uint64_t dbl_bits = std::bit_cast<uint64_t>(dbl);
        int64_t exponent = (((dbl_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
        dbl_bits &= 0x800FFFFFFFFFFFFFull;
        dbl_bits |= 0x3FF0000000000000ull;

        double sicnificand = std::bit_cast<double>(dbl_bits);
        return (double)exponent + (sicnificand - 1.0);
    }

    void Float64LargeRangeNumber::doubles_to(const std::vector<double>& in, std::vector<double>& out) {
        for (size_t i = 0; i < in.size(); i++) {
            out[i] = double_to(in[i]);
        }
    }

    double Float64LargeRangeNumber::as_double(double large_range_number) {
        double floored_exponent = std::floor(large_range_number);
        int64_t exponent = (int64_t)floored_exponent;
        double sicnificand = large_range_number - floored_exponent + 1.0;
        return sicnificand * std::pow(2.0, (int)exponent);
    }

    double Float64LargeRangeNumber::as_double() {
        double floored_exponent = std::floor(encoded);
        int64_t exponent = (int64_t)floored_exponent;
        double sicnificand = encoded - floored_exponent + 1.0;
        return sicnificand * std::pow(2.0, (int)exponent);
    }

    void Float64LargeRangeNumber::log2_to(double log2_dbl) {
        double floored_exponent = std::floor(log2_dbl);
        double sicnificand = log2_dbl - floored_exponent;
        encoded = std::exp2(sicnificand) + floored_exponent - 1;
    }

    double Float64LargeRangeNumber::as_log2(double large_range_number) {
        double floored_exponent = std::floor(large_range_number);
        double sicnificand = large_range_number - floored_exponent + 1.0;
        return std::log2(sicnificand) + floored_exponent;
    }

    

    

    //inline void decode_fast_double(double& encoded_value, double& exponent, double& sicnificand_d) {
    //    uint64_t dbl_bits = std::bit_cast<uint64_t>(encoded_value);
    //    int32_t cutter = ((dbl_bits >> 52) & 0x7FF) - 1023;
    //    uint64_t exponent_int64;

    //    if (cutter >= 0) {
    //        exponent_int64 = dbl_bits & (~(0x000FFFFFFFFFFFFFull >> cutter));
    //    }
    //    else {
    //        exponent_int64 = 0ULL;
    //    }

    //    exponent = std::bit_cast<double>(exponent_int64);
    //    sicnificand_d = encoded_value - exponent;

    //    if (sicnificand_d < 0.0) {
    //        exponent -= 1.0;
    //        sicnificand_d += 2.0;
    //    }
    //    else {
    //        sicnificand_d += 1.0;
    //    }
    //}

    double Float64LargeRangeNumber::multiply(double lrn1, double lrn2) {
        double exponent_1 = std::floor(lrn1);
        double sicnificand_1 = lrn1 - exponent_1 + 1;

        double exponent_2 = std::floor(lrn2);
        double sicnificand_2 = lrn2 - exponent_2 + 1;

        //double exponent_1;
        //double sicnificand_1;

        //double exponent_2;
        //double sicnificand_2;

        //decode_fast_double(lrn1, exponent_1, sicnificand_1);
        //decode_fast_double(lrn2, exponent_2, sicnificand_2);

        sicnificand_1 *= sicnificand_2;
        exponent_1 += exponent_2;

        if (sicnificand_1 >= 2.0) {
            //int64_t sicnificand_int64_1 = std::bit_cast<int64_t>(sicnificand_1);
            //sicnificand_int64_1 -= 1ULL << 52;
            //sicnificand_1 = std::bit_cast<double>(sicnificand_int64_1);
            sicnificand_1 *= 0.5;
            return (double)exponent_1 + sicnificand_1;
        }
        else {
            return (double)exponent_1 + sicnificand_1 - 1.0;
        }
    }

    double Float64LargeRangeNumber::sum(std::vector<double>& large_range_numbers) {
        if (large_range_numbers.size() < 16) {
            double sum_small = large_range_numbers[0];
            for (size_t k = 1; k < large_range_numbers.size(); k++) {
                sum_small = sum(sum_small, large_range_numbers[k]);
            }
            return sum_small;
        }

        __m512d zero_d = _mm512_setzero_pd();
        __m512d half_d = _mm512_set1_pd(0.5);
        __m512d one_d = _mm512_set1_pd(1.0);
        __m512d two_d = _mm512_set1_pd(2.0);
        __m512i m1_i = _mm512_set1_epi64(-1);
        __m512i zero_i = _mm512_setzero_si512();
        __m512i one_i = _mm512_set1_epi64(1);
        __m512i p63_i = _mm512_set1_epi64(63);
        __m512i m63_i = _mm512_set1_epi64(-63);

        constexpr size_t n_parallel = 5;

        __m512d encoded_values[n_parallel];
        __m512d integer_decimals[n_parallel];
        __m512i exponent_sum[n_parallel];
        __m512d sicnificand_sum[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[8 * m]);
            integer_decimals[m] = _mm512_floor_pd(encoded_values[m]);
            exponent_sum[m] = _mm512_cvtpd_epi64(integer_decimals[m]);
            sicnificand_sum[m] = encoded_values[m] - integer_decimals[m] + one_d;
        }

        __m512i exponent[n_parallel];
        __m512d sicnificand[n_parallel];
        __m512i exp_diff[n_parallel];
        __m512i sicnificand_bits[n_parallel];
        __m512i sicnificand_sum_bits[n_parallel];

        size_t i = 8;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < large_range_numbers.size(); i += 8 * n_parallel) {
            for (size_t m = 0; m < n_parallel; m++) {
                encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[i + 8 * m]);

                integer_decimals[m] = _mm512_floor_pd(encoded_values[m]);
                exponent[m] = _mm512_cvttpd_epi64(integer_decimals[m]);
                sicnificand[m] = encoded_values[m] - integer_decimals[m] + one_d;

                exp_diff[m] = exponent_sum[m] - exponent[m];

                sicnificand_bits[m] = _mm512_castpd_si512(sicnificand[m]);
                sicnificand_bits[m] = exp_diff[m] < zero_i ? sicnificand_bits[m] : sicnificand_bits[m] - (exp_diff[m] << 52);
                sicnificand[m] = _mm512_castsi512_pd(sicnificand_bits[m]);

                sicnificand_sum_bits[m] = _mm512_castpd_si512(sicnificand_sum[m]);
                sicnificand_sum_bits[m] = exp_diff[m] >= zero_i ? sicnificand_sum_bits[m] : sicnificand_sum_bits[m] - ((m1_i[m] * exp_diff[m]) << 52);
                sicnificand_sum[m] = _mm512_castsi512_pd(sicnificand_sum_bits[m]);

                sicnificand_sum[m] = exp_diff[m] < m63_i ? zero_d : sicnificand_sum[m];
                sicnificand[m] = exp_diff[m] > p63_i ? zero_d : sicnificand[m];

                exponent_sum[m] = exp_diff[m] < zero_i ? exponent[m] : exponent_sum[m];

                sicnificand_sum[m] = sicnificand_sum[m] + sicnificand[m];

                exponent_sum[m] = sicnificand_sum[m] >= two_d ? exponent_sum[m] + one_i : exponent_sum[m];
                sicnificand_sum[m] = sicnificand_sum[m] >= two_d ? sicnificand_sum[m] * half_d : sicnificand_sum[m];
            }
        }

        double total_sum = -0x1p52;

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                total_sum = sum(total_sum, sicnificand_sum[m][k] + (double)exponent_sum[m][k] - 1);
            }
        }

        for (i = i; i < large_range_numbers.size(); i++) {
            total_sum = sum(total_sum, large_range_numbers[i]);
        }

        return total_sum;
    }

    double Float64LargeRangeNumber::multiply(std::vector<double>& large_range_numbers) {
        if (large_range_numbers.size() < 16) {
            double sum_small = large_range_numbers[0];
            for (size_t k = 1; k < large_range_numbers.size(); k++) {
                sum_small = multiply(sum_small, large_range_numbers[k]);
            }
            return sum_small;
        }

        __m512d half_d = _mm512_set1_pd(0.5);
        __m512d one_d = _mm512_set1_pd(1.0);
        __m512d two_d = _mm512_set1_pd(2.0);

        constexpr size_t n_parallel = 5;

        __m512d encoded_values[n_parallel];
        __m512d exponent_floored_sum[n_parallel];
        __m512d sicnificand_sum[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[8 * m]);
            exponent_floored_sum[m] = _mm512_floor_pd(encoded_values[m]);
            sicnificand_sum[m] = encoded_values[m] - exponent_floored_sum[m] + one_d;
        }

        __m512d exponent_floored[n_parallel];
        __m512d sicnificand[n_parallel];

        size_t i = 8;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < large_range_numbers.size(); i += 8 * n_parallel) {
            for (size_t m = 0; m < n_parallel; m++) {
                encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[i + 8 * m]);

                exponent_floored[m] = _mm512_floor_pd(encoded_values[m]);
                sicnificand[m] = encoded_values[m] - exponent_floored[m] + one_d;

                sicnificand_sum[m] = sicnificand_sum[m] * sicnificand[m];
                exponent_floored_sum[m] = exponent_floored_sum[m] + exponent_floored[m];

                exponent_floored_sum[m] = sicnificand_sum[m] >= two_d ? exponent_floored_sum[m] + one_d : exponent_floored_sum[m];
                sicnificand_sum[m] = sicnificand_sum[m] >= two_d ? sicnificand_sum[m] * half_d : sicnificand_sum[m];
            }
        }

        // exponent + (sicnificand - 1.0);
        double sum = 0.0 + (1.0 - 1.0);

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                sum = multiply(sum, sicnificand_sum[m][k] + exponent_floored_sum[m][k] - 1);
            }
        }

        for (i = i; i < large_range_numbers.size(); i++) {
            sum = multiply(sum, large_range_numbers[i]);
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

