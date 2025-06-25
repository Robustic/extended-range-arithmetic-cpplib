#include <cstdint>
#include <cmath>
#include <bit>
#include <vector>
#include <cstring>
#include <immintrin.h>
#include <iostream>
#include <bitset>
#include <cstdint>

#include "WideRangeNumber64.h"

namespace extended_range_arithmetic
{
    double WideRangeNumber64::double_to(double dbl) {
        uint64_t dbl_bits = std::bit_cast<uint64_t>(dbl);
        int64_t auxiliary_index = (((dbl_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
        dbl_bits &= 0x800FFFFFFFFFFFFFull;
        dbl_bits |= 0x3FF0000000000000ull;

        double principal_part = std::bit_cast<double>(dbl_bits);
        return (double)auxiliary_index + (principal_part - 1.0);
    }

    void WideRangeNumber64::doubles_to(const std::vector<double>& in, std::vector<double>& out) {
        for (size_t i = 0; i < in.size(); i++) {
            out[i] = double_to(in[i]);
        }
    }

    double WideRangeNumber64::as_double(double large_range_number) {
        double floored_auxiliary_index = std::floor(large_range_number);
        int64_t auxiliary_index = (int64_t)floored_auxiliary_index;
        double principal_part = large_range_number - floored_auxiliary_index + 1.0;
        return principal_part * std::pow(2.0, (int)auxiliary_index);
    }

    double WideRangeNumber64::log2_to(double log2_dbl) {
        double floored_auxiliary_index = std::floor(log2_dbl);
        double principal_part = log2_dbl - floored_auxiliary_index;
        return std::exp2(principal_part) + floored_auxiliary_index - 1;
    }

    double WideRangeNumber64::as_log2(double large_range_number) {
        double floored_auxiliary_index = std::floor(large_range_number);
        double principal_part = large_range_number - floored_auxiliary_index + 1.0;
        return std::log2(principal_part) + floored_auxiliary_index;
    }

    inline static void decode_fast_int(double& encoded_value, int64_t& auxiliary_index, double& principal_part_d) {
        uint64_t dbl_bits = std::bit_cast<uint64_t>(encoded_value);
        int32_t cutter = ((dbl_bits >> 52) & 0x7FF) - 1023;
        int64_t full_mantissa = (dbl_bits & 0x000FFFFFFFFFFFFFull) | 0x0010000000000000ull;
        uint64_t principal_part_int64;

        if (cutter >= 0) {
            auxiliary_index = ((int64_t)full_mantissa >> (52 - cutter));
            principal_part_int64 = (dbl_bits << cutter) & 0x000FFFFFFFFFFFFFull;
        }
        else {
            auxiliary_index = 0;
            principal_part_int64 = full_mantissa >> (-1 * cutter);
        }

        auxiliary_index = auxiliary_index * (encoded_value < 0.0 ? -1LL : 1LL);
        auxiliary_index = (encoded_value < 0.0 && principal_part_int64 > 0ULL) ? auxiliary_index - 1 : auxiliary_index;

        principal_part_int64 = encoded_value < 0.0 ? (0x0020000000000000ull - principal_part_int64) : principal_part_int64;
        principal_part_int64 = principal_part_int64 | 0x3FF0000000000000ull;
        principal_part_d = std::bit_cast<double>(principal_part_int64);
    }

    double WideRangeNumber64::sum(double lrn1, double lrn2) {
        if (lrn1 - lrn2 > 54.0) {
            return lrn1;
        }
        else if (lrn1 - lrn2 < -54.0) {
            return lrn2;
        }

        int64_t auxiliary_index_1;
        double principal_part_1;

        int64_t auxiliary_index_2;
        double principal_part_2;

        decode_fast_int(lrn1, auxiliary_index_1, principal_part_1);
        decode_fast_int(lrn2, auxiliary_index_2, principal_part_2);

        int64_t aux_diff = (int64_t)(auxiliary_index_1 - auxiliary_index_2);

        int64_t principal_part_int64_1 = std::bit_cast<int64_t>(principal_part_1);
        int64_t principal_part_int64_2 = std::bit_cast<int64_t>(principal_part_2);

        if (aux_diff > 0) {
            principal_part_int64_2 -= aux_diff << 52;
        }
        else if (aux_diff < 0) {
            principal_part_int64_1 -= (-aux_diff) << 52;
            auxiliary_index_1 = auxiliary_index_2;
        }

        principal_part_1 = std::bit_cast<double>(principal_part_int64_1);
        principal_part_2 = std::bit_cast<double>(principal_part_int64_2);

        principal_part_1 += principal_part_2;

        if (principal_part_1 >= 2.0) {
            principal_part_1 *= 0.5;
            return (double)auxiliary_index_1 + principal_part_1;
        }
        else {
            return (double)auxiliary_index_1 + principal_part_1 - 1.0;
        }
    }

    WideRangeNumber64& WideRangeNumber64::operator+=(double encoded_2) {
        double diff = (double)aux - encoded_2;
        if (diff > 54.0) {
            return *this;
        }
        else if (diff < -54.0) {
            int64_t auxiliary_index_temp;
            double principal_part_temp;
            decode_fast_int(encoded_2, auxiliary_index_temp, principal_part_temp);
            aux = auxiliary_index_temp;
            principal = principal_part_temp;
            encoded = encoded_2;
            return *this;
        }

        int64_t auxiliary_index_2;
        double principal_part_2;

        decode_fast_int(encoded_2, auxiliary_index_2, principal_part_2);

        int64_t aux_diff = aux - auxiliary_index_2;

        int64_t principal_part_int64_1 = std::bit_cast<int64_t>(principal);
        int64_t principal_part_int64_2 = std::bit_cast<int64_t>(principal_part_2);

        if (aux_diff > 0) {
            principal_part_int64_2 -= aux_diff << 52;
        }
        else if (aux_diff < 0) {
            principal_part_int64_1 -= (-aux_diff) << 52;
            aux = auxiliary_index_2;
        }

        principal = std::bit_cast<double>(principal_part_int64_1);
        principal_part_2 = std::bit_cast<double>(principal_part_int64_2);

        principal += principal_part_2;

        if (principal >= 2.0) {
            principal *= 0.5;
            encoded = (double)aux + principal;
            aux++;
        }
        else {
            encoded = (double)aux + principal - 1.0;
        }
        return *this;
    }

    double WideRangeNumber64::sum(const std::vector<double>& large_range_numbers) {
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
        __m512i auxiliary_index_sum[n_parallel];
        __m512d principal_part_sum[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[8 * m]);
            integer_decimals[m] = _mm512_floor_pd(encoded_values[m]);
            auxiliary_index_sum[m] = _mm512_cvtpd_epi64(integer_decimals[m]);
            principal_part_sum[m] = encoded_values[m] - integer_decimals[m] + one_d;
        }

        __m512i auxiliary_index[n_parallel];
        __m512d principal_part[n_parallel];
        __m512i aux_diff[n_parallel];
        __m512i principal_part_bits[n_parallel];
        __m512i principal_part_sum_bits[n_parallel];

        size_t counter = 0;
        size_t i = 8 * n_parallel;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < large_range_numbers.size(); i += 8 * n_parallel) {
            for (size_t m = 0; m < n_parallel; m++) {
                encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[i + 8 * m]);

                integer_decimals[m] = _mm512_floor_pd(encoded_values[m]);
                auxiliary_index[m] = _mm512_cvttpd_epi64(integer_decimals[m]);
                principal_part[m] = encoded_values[m] - integer_decimals[m] + one_d;

                aux_diff[m] = auxiliary_index_sum[m] - auxiliary_index[m];

                principal_part_bits[m] = _mm512_castpd_si512(principal_part[m]);
                principal_part_bits[m] = aux_diff[m] < zero_i ? principal_part_bits[m] : principal_part_bits[m] - (aux_diff[m] << 52);
                principal_part[m] = _mm512_castsi512_pd(principal_part_bits[m]);

                principal_part_sum_bits[m] = _mm512_castpd_si512(principal_part_sum[m]);
                principal_part_sum_bits[m] = aux_diff[m] >= zero_i ? principal_part_sum_bits[m] : principal_part_sum_bits[m] - ((m1_i[m] * aux_diff[m]) << 52);
                principal_part_sum[m] = _mm512_castsi512_pd(principal_part_sum_bits[m]);

                principal_part_sum[m] = aux_diff[m] < mlimit_i ? zero_d : principal_part_sum[m];
                principal_part[m] = aux_diff[m] > plimit_i ? zero_d : principal_part[m];

                auxiliary_index_sum[m] = aux_diff[m] < zero_i ? auxiliary_index[m] : auxiliary_index_sum[m];

                principal_part_sum[m] = principal_part_sum[m] + principal_part[m];
            }

            counter++;
            if (counter > 500) {
                counter = 0;
                for (size_t m = 0; m < n_parallel; m++) {
                    __m512i principal_bits = std::bit_cast<__m512i>(principal_part_sum[m]);
                    auxiliary_index_sum[m] = auxiliary_index_sum[m] + ((principal_bits & 0x7FF0000000000000ull) >> 52) - 1023LL;
                    principal_bits &= 0x800FFFFFFFFFFFFFull;
                    principal_bits |= 0x3FF0000000000000ull;
                    principal_part_sum[m] = std::bit_cast<__m512d>(principal_bits);
                }
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            __m512i principal_bits = std::bit_cast<__m512i>(principal_part_sum[m]);
            auxiliary_index_sum[m] = auxiliary_index_sum[m] + ((principal_bits & 0x7FF0000000000000ull) >> 52) - 1023LL;
            principal_bits &= 0x800FFFFFFFFFFFFFull;
            principal_bits |= 0x3FF0000000000000ull;
            principal_part_sum[m] = std::bit_cast<__m512d>(principal_bits);
        }

        double total_sum = principal_part_sum[0][0] + (double)auxiliary_index_sum[0][0] - 1;

        for (size_t k = 1; k < 8; k++) {
            total_sum = sum(total_sum, principal_part_sum[0][k] + (double)auxiliary_index_sum[0][k] - 1);
        }

        for (size_t m = 1; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                total_sum = sum(total_sum, principal_part_sum[m][k] + (double)auxiliary_index_sum[m][k] - 1);
            }
        }

        for (i = i; i < large_range_numbers.size(); i++) {
            total_sum = sum(total_sum, large_range_numbers[i]);
        }

        return total_sum;
    }

    WideRangeNumber64& WideRangeNumber64::operator*=(double encoded_2) {
        double auxiliary_index_2 = std::floor(encoded_2);
        double principal_part_2 = encoded_2 - auxiliary_index_2 + 1;

        principal *= principal_part_2;
        aux += auxiliary_index_2;

        if (principal >= 2.0) {
            principal *= 0.5;
            encoded = (double)aux + principal;
            aux++;
        }
        else {
            encoded = (double)aux + principal - 1.0;
        }

        return *this;
    }

    double WideRangeNumber64::multiply(double lrn1, double lrn2) {
        double auxiliary_index_1 = std::floor(lrn1);
        double principal_part_1 = lrn1 - auxiliary_index_1 + 1;

        double auxiliary_index_2 = std::floor(lrn2);
        double principal_part_2 = lrn2 - auxiliary_index_2 + 1;

        principal_part_1 *= principal_part_2;
        auxiliary_index_1 += auxiliary_index_2;

        if (principal_part_1 >= 2.0) {
            principal_part_1 *= 0.5;
            return (double)auxiliary_index_1 + principal_part_1;
        }
        else {
            return (double)auxiliary_index_1 + principal_part_1 - 1.0;
        }
    }

    double WideRangeNumber64::multiply(const std::vector<double>& large_range_numbers) {
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
        __m512d auxiliary_index_floored_sum[n_parallel];
        __m512d principal_part_sum[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[8 * m]);
            auxiliary_index_floored_sum[m] = _mm512_floor_pd(encoded_values[m]);
            principal_part_sum[m] = encoded_values[m] - auxiliary_index_floored_sum[m] + one_d;
        }

        __m512d auxiliary_index_floored[n_parallel];
        __m512d principal_part[n_parallel];

        size_t counter = 0;
        size_t i = 8 * n_parallel;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < large_range_numbers.size(); i += 8 * n_parallel) {
            counter++;
            for (size_t m = 0; m < n_parallel; m++) {
                encoded_values[m] = _mm512_loadu_pd(&large_range_numbers[i + 8 * m]);

                auxiliary_index_floored[m] = _mm512_floor_pd(encoded_values[m]);
                principal_part[m] = encoded_values[m] - auxiliary_index_floored[m] + one_d;

                principal_part_sum[m] = principal_part_sum[m] * principal_part[m];
                auxiliary_index_floored_sum[m] = auxiliary_index_floored_sum[m] + auxiliary_index_floored[m];

                auxiliary_index_floored_sum[m] = principal_part_sum[m] >= two_d ? auxiliary_index_floored_sum[m] + one_d : auxiliary_index_floored_sum[m];
                principal_part_sum[m] = principal_part_sum[m] >= two_d ? principal_part_sum[m] * half_d : principal_part_sum[m];
            }

            if (counter > 490) {
                counter = 0;
                for (size_t m = 0; m < n_parallel; m++) {
                    __m512i principal_bits = std::bit_cast<__m512i>(principal_part_sum[m]);
                    auxiliary_index_floored_sum[m] = auxiliary_index_floored_sum[m] + _mm512_castsi512_pd(((principal_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
                    principal_bits &= 0x800FFFFFFFFFFFFFull;
                    principal_bits |= 0x3FF0000000000000ull;
                    principal_part_sum[m] = std::bit_cast<__m512d>(principal_bits);
                }
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            __m512i principal_bits = std::bit_cast<__m512i>(principal_part_sum[m]);
            auxiliary_index_floored_sum[m] = auxiliary_index_floored_sum[m] + _mm512_castsi512_pd(((principal_bits & 0x7FF0000000000000ull) >> 52) - 1023LL);
            principal_bits &= 0x800FFFFFFFFFFFFFull;
            principal_bits |= 0x3FF0000000000000ull;
            principal_part_sum[m] = std::bit_cast<__m512d>(principal_bits);
        }

        // auxiliary_index + (principal_part - 1.0);
        double sum = 0.0 + (1.0 - 1.0);

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                sum = multiply(sum, principal_part_sum[m][k] + auxiliary_index_floored_sum[m][k] - 1.0);
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

