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

    double Float64LargeRangeNumber::log2_to(double log2_dbl) {
        double floored_exponent = std::floor(log2_dbl);
        double sicnificand = log2_dbl - floored_exponent;
        return std::exp2(sicnificand) + floored_exponent - 1;
    }

    double Float64LargeRangeNumber::as_log2(double large_range_number) {
        double floored_exponent = std::floor(large_range_number);
        double sicnificand = large_range_number - floored_exponent + 1.0;
        return std::log2(sicnificand) + floored_exponent;
    }

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

    double Float64LargeRangeNumber::sum(double lrn1, double lrn2) {
        if (lrn1 - lrn2 > 54.0) {
            return lrn1;
        }
        else if (lrn1 - lrn2 < -54.0) {
            return lrn2;
        }

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
            sicnificand_int64_2 -= exp_diff << 52;
        }
        else if (exp_diff < 0) {
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

    Float64LargeRangeNumber& Float64LargeRangeNumber::operator+=(double encoded_2) {
        double diff = (double)exp - encoded_2;
        if (diff > 54.0) {
            return *this;
        }
        else if (diff < -54.0) {
            int64_t exponent_temp;
            double sicnificand_temp;
            decode_fast_int(encoded_2, exponent_temp, sicnificand_temp);
            exp = exponent_temp;
            scnfcnd = sicnificand_temp;
            encoded = encoded_2;
            return *this;
        }

        int64_t exponent_2;
        double sicnificand_2;

        decode_fast_int(encoded_2, exponent_2, sicnificand_2);

        int64_t exp_diff = exp - exponent_2;

        int64_t sicnificand_int64_1 = std::bit_cast<int64_t>(scnfcnd);
        int64_t sicnificand_int64_2 = std::bit_cast<int64_t>(sicnificand_2);

        if (exp_diff > 0) {
            sicnificand_int64_2 -= exp_diff << 52;
        }
        else if (exp_diff < 0) {
            sicnificand_int64_1 -= (-exp_diff) << 52;
            exp = exponent_2;
        }

        scnfcnd = std::bit_cast<double>(sicnificand_int64_1);
        sicnificand_2 = std::bit_cast<double>(sicnificand_int64_2);

        scnfcnd += sicnificand_2;

        if (scnfcnd >= 2.0) {
            scnfcnd *= 0.5;
            encoded = (double)exp + scnfcnd;
            exp++;
        }
        else {
            encoded = (double)exp + scnfcnd - 1.0;
        }
        return *this;
    }

    double Float64LargeRangeNumber::sum(const std::vector<double>& large_range_numbers) {
        constexpr size_t n_parallel = 5;

        if (large_range_numbers.size() < 2 * 8 * n_parallel) {
            double sum_small = large_range_numbers[0];
            for (size_t k = 1; k < large_range_numbers.size(); k++) {
                sum_small = sum(sum_small, large_range_numbers[k]);
            }
            return sum_small;
        }

        __m512d zero_d = _mm512_setzero_pd();
        __m512d one_d = _mm512_set1_pd(1.0);
        __m512i m1_i = _mm512_set1_epi64(-1);
        __m512i zero_i = _mm512_setzero_si512();
        __m512i plimit_i = _mm512_set1_epi64(70);
        __m512i mlimit_i = _mm512_set1_epi64(-70);

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

        size_t counter = 0;
        size_t i = 8 * n_parallel;
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

                sicnificand_sum[m] = exp_diff[m] < mlimit_i ? zero_d : sicnificand_sum[m];
                sicnificand[m] = exp_diff[m] > plimit_i ? zero_d : sicnificand[m];

                exponent_sum[m] = exp_diff[m] < zero_i ? exponent[m] : exponent_sum[m];

                sicnificand_sum[m] = sicnificand_sum[m] + sicnificand[m];
            }

            counter++;
            if (counter > 500) {
                counter = 0;
                for (size_t m = 0; m < n_parallel; m++) {
                    __m512i sgnfcnd_bits = std::bit_cast<__m512i>(sicnificand_sum[m]);
                    exponent_sum[m] = exponent_sum[m] + ((sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1023LL;
                    sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
                    sgnfcnd_bits |= 0x3FF0000000000000ull;
                    sicnificand_sum[m] = std::bit_cast<__m512d>(sgnfcnd_bits);
                }
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            __m512i sgnfcnd_bits = std::bit_cast<__m512i>(sicnificand_sum[m]);
            exponent_sum[m] = exponent_sum[m] + ((sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1023LL;
            sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
            sgnfcnd_bits |= 0x3FF0000000000000ull;
            sicnificand_sum[m] = std::bit_cast<__m512d>(sgnfcnd_bits);
        }

        double total_sum = sicnificand_sum[0][0] + (double)exponent_sum[0][0] - 1;

        for (size_t k = 1; k < 8; k++) {
            total_sum = sum(total_sum, sicnificand_sum[0][k] + (double)exponent_sum[0][k] - 1);
        }

        for (size_t m = 1; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                total_sum = sum(total_sum, sicnificand_sum[m][k] + (double)exponent_sum[m][k] - 1);
            }
        }

        for (i = i; i < large_range_numbers.size(); i++) {
            total_sum = sum(total_sum, large_range_numbers[i]);
        }

        return total_sum;
    }

    Float64LargeRangeNumber& Float64LargeRangeNumber::operator*=(double encoded_2) {
        double exponent_2 = std::floor(encoded_2);
        double sicnificand_2 = encoded_2 - exponent_2 + 1;

        scnfcnd *= sicnificand_2;
        exp += exponent_2;

        if (scnfcnd >= 2.0) {
            scnfcnd *= 0.5;
            encoded = (double)exp + scnfcnd;
            exp++;
        }
        else {
            encoded = (double)exp + scnfcnd - 1.0;
        }

        return *this;
    }

    double Float64LargeRangeNumber::multiply(double lrn1, double lrn2) {
        double exponent_1 = std::floor(lrn1);
        double sicnificand_1 = lrn1 - exponent_1 + 1;

        double exponent_2 = std::floor(lrn2);
        double sicnificand_2 = lrn2 - exponent_2 + 1;

        sicnificand_1 *= sicnificand_2;
        exponent_1 += exponent_2;

        if (sicnificand_1 >= 2.0) {
            sicnificand_1 *= 0.5;
            return (double)exponent_1 + sicnificand_1;
        }
        else {
            return (double)exponent_1 + sicnificand_1 - 1.0;
        }
    }

    double Float64LargeRangeNumber::multiply(const std::vector<double>& large_range_numbers) {
        constexpr size_t n_parallel = 5;

        if (large_range_numbers.size() < 2 * 8 * n_parallel) {
            double sum_small = large_range_numbers[0];
            for (size_t k = 1; k < large_range_numbers.size(); k++) {
                sum_small = multiply(sum_small, large_range_numbers[k]);
            }
            return sum_small;
        }

        __m512d half_d = _mm512_set1_pd(0.5);
        __m512d one_d = _mm512_set1_pd(1.0);
        __m512d two_d = _mm512_set1_pd(2.0);

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

        size_t counter = 0;
        size_t i = 8 * n_parallel;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < large_range_numbers.size(); i += 8 * n_parallel) {
            counter++;
            for (size_t m = 0; m < n_parallel; m++) {
                encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[i + 8 * m]);

                exponent_floored[m] = _mm512_floor_pd(encoded_values[m]);
                sicnificand[m] = encoded_values[m] - exponent_floored[m] + one_d;

                sicnificand_sum[m] = sicnificand_sum[m] * sicnificand[m];
                exponent_floored_sum[m] = exponent_floored_sum[m] + exponent_floored[m];

                exponent_floored_sum[m] = sicnificand_sum[m] >= two_d ? exponent_floored_sum[m] + one_d : exponent_floored_sum[m];
                sicnificand_sum[m] = sicnificand_sum[m] >= two_d ? sicnificand_sum[m] * half_d : sicnificand_sum[m];
            }

            if (counter > 490) {
                counter = 0;
                for (size_t m = 0; m < n_parallel; m++) {
                    __m512i sgnfcnd_bits = std::bit_cast<__m512i>(sicnificand_sum[m]);
                    exponent_floored_sum[m] = exponent_floored_sum[m] + _mm512_castsi512_pd(((sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
                    sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
                    sgnfcnd_bits |= 0x3FF0000000000000ull;
                    sicnificand_sum[m] = std::bit_cast<__m512d>(sgnfcnd_bits);
                }
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            __m512i sgnfcnd_bits = std::bit_cast<__m512i>(sicnificand_sum[m]);
            exponent_floored_sum[m] = exponent_floored_sum[m] + _mm512_castsi512_pd(((sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
            sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
            sgnfcnd_bits |= 0x3FF0000000000000ull;
            sicnificand_sum[m] = std::bit_cast<__m512d>(sgnfcnd_bits);
        }

        // exponent + (sicnificand - 1.0);
        double sum = 0.0 + (1.0 - 1.0);

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                sum = multiply(sum, sicnificand_sum[m][k] + exponent_floored_sum[m][k] - 1.0);
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

