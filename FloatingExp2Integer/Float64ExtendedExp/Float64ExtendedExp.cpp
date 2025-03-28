#include <cstdint>
#include <cmath>
#include <bit>
#include <vector>
#include <cstring>
#include <immintrin.h>
#include <iostream>
#include <bitset>
#include <cstdint>

#include "Float64ExtendedExp.h"

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
typedef uint64_t uint64_4_t __attribute__((vector_size(4 * sizeof(uint64_t))));
typedef int64_t int64_4_t __attribute__((vector_size(4 * sizeof(int64_t))));

constexpr int64_4_t int64_4_1023 { 1023LL, 1023LL, 1023LL, 1023LL };
constexpr double4_t double4_0 { 0.0, 0.0, 0.0, 0.0 };

namespace floatingExp2Integer
{
    Float64ExtendedExp::Float64ExtendedExp() {
        encode_double(1.0, 0LL);
    }

    Float64ExtendedExp::Float64ExtendedExp(double dbl) {
        encode_double(dbl, 0LL);
    }

    Float64ExtendedExp::Float64ExtendedExp(double sicnificand, std::int64_t exponent) {
        encode_double(sicnificand, exponent);
    }

    double Float64ExtendedExp::decodeFloat64ExtendedExp(double data) {
        std::int64_t exponent = (std::int64_t)std::floor(data);
        double sicnificand = data - exponent + 1.0;
        return sicnificand * std::pow(2.0, (int)exponent);
    }

    double Float64ExtendedExp::sumFloat64ExtendedExp(unsigned int n, double* vector_in) {
        __m512d encoded_values = _mm512_loadu_pd(&vector_in[0]);

        __m512d integer_decimals = _mm512_floor_pd(encoded_values);
        __m512i exponent_sum = _mm512_cvtpd_epi64(integer_decimals);
        __m512d one_d = _mm512_set1_pd(1.0);
        __m512d two_d = _mm512_set1_pd(2.0);
        __m512i one_i = _mm512_set1_epi64(1);
        __m512d half_d = _mm512_set1_pd(0.5);
        __m512d zero_d = _mm512_setzero_pd();
        __m512d sicnificand_sum = encoded_values - integer_decimals + one_d;

        //convert_double_to_uint64(encoded_values, exponent, sicnificand);

        for (unsigned int i = 8; i + 7 < n; i += 8) {
            //unsigned int k = sicnificand_sum[0] > 1000000000.0 ? i - 1 : i;
            encoded_values = _mm512_loadu_pd(&vector_in[i]);

            integer_decimals = _mm512_floor_pd(encoded_values);
            __m512i exponent = _mm512_cvttpd_epi64(integer_decimals);
            __m512d sicnificand = encoded_values - integer_decimals + one_d;

            //__m512i exponent = _mm512_cvttpd_epi64(encoded_values); // function truncatei: round towards zero
            //__m512d sicnificand = encoded_values - _mm512_cvtepi64_pd(exponent);
            //exponent = sicnificand < zero_d ? exponent - one_i : exponent;
            //sicnificand = sicnificand < zero_d ? sicnificand + two_d : sicnificand + one_d;

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

            //print("double 0 ", vector_in[0]);
            //print("exponent ", exponent_sum[0]);
            //print("sicnificand ", sicnificand_sum[0]);
            //print("double 8 ", vector_in[0]);
            //print("exponent ", exponent[0]);
            //print("sicnificand ", sicnificand[0]);

            sicnificand_sum = sicnificand_sum + sicnificand;

            exponent_sum = sicnificand_sum >= two_d ? exponent_sum + one_i : exponent_sum;
            sicnificand_sum = sicnificand_sum >= two_d ? sicnificand_sum * half_d : sicnificand_sum;
        }

        //for (int i = 0; i < 8; i++) {
        //    print("double ", vector_in[i]);
        //    print("exponent ", exponent_sum[i]);
        //    print("sicnificand ", sicnificand_sum[i]);
        //}

        floatingExp2Integer::Float64ExtendedExp collector(sicnificand_sum[0], exponent_sum[0]);

        for (int i = 1; i < 8; i++) {
            floatingExp2Integer::Float64ExtendedExp item(sicnificand_sum[i], exponent_sum[i]);
            collector += item;
        }

        //print("exponent ", collector.exponent());
        //print("sicnificand ", collector.sicnificand());
        //print("collector.asDouble ", collector.asDouble());

        return collector.encoded;
    }

    void Float64ExtendedExp::doubleToFloat64ExtendedExp(double dbl) {
        encode_double(dbl, 0LL);
    }

    void Float64ExtendedExp::encodedToFloat64ExtendedExp(double dbl) {
        encoded = dbl;
        double floored = std::floor(dbl);
        exp = (std::int64_t)floored;
        scnfcnd = dbl - floored + 1.0;
    }

    void Float64ExtendedExp::doubleToFloat64ExtendedExp(unsigned int n, const double* vector_in, double* vector_out) {
        for (unsigned int i = 0; i < n; i++) {
            uint64_t dbl_bits = std::bit_cast<std::uint64_t>(vector_in[i]);
            int64_t exponent = ((((dbl_bits & 0x7FF0000000000000ull)) >> 52) - 1023LL);
            dbl_bits &= 0x800FFFFFFFFFFFFFull;
            dbl_bits |= 0x3FF0000000000000ull;
            double sicnificand = std::bit_cast<double>(dbl_bits);
            vector_out[i] = (double)exponent + (sicnificand - 1.0);
        }
    }

    // void Float64ExtendedExp::log2ToFloat64ExtendedExp(double log2) {
    //     std::int64_t exponent = (std::int64_t)log2;
    //     double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
    //     this->scale(sicnificand, exponent);
    // }

    // double Float64ExtendedExp::float64ExtendedExpToLog2() const {
    //     std::int64_t exponent = (std::int64_t)encoded;
    //     double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
    //     return std::log2(sicnificand) + exponent;
    // }

    double Float64ExtendedExp::sicnificand() {
        std::int64_t exponent;
        double sicnificand;
        decode_double(sicnificand, exponent); 
        double sicnificand_return = sicnificand;
        return sicnificand_return; 
    }

    std::int64_t Float64ExtendedExp::exponent() {
        std::int64_t exponent;
        double sicnificand;
        decode_double(sicnificand, exponent);
        return exponent;
    }

    double Float64ExtendedExp::asDouble() {
        std::int64_t exponent;
        double sicnificand;
        decode_double(sicnificand, exponent);
        return sicnificand * std::pow(2.0, (int)exponent);
    }

    Float64ExtendedExp& Float64ExtendedExp::operator+=(double z_encoded) {
        //std::uint64_t sicnificand;
        //std::int64_t exponent;
        //decode(sicnificand, exponent);

        //std::uint64_t sicnificand_z;
        //std::int64_t exponent_z;
        //decode(sicnificand_z, exponent_z, z_encoded);

        //****************

        //std::int64_t exponent_z = (std::int64_t)std::floor(z_encoded);
        //double sicnificandZ_d = z_encoded - exponent_z + 1.0;

        //std::int64_t exponent_z = (std::int64_t)z_encoded;
        //double sicnificandZ_d = z_encoded - exponent_z;

        //if (sicnificandZ_d < 0.0) {
        //    exponent_z -= 1LL;
        //    sicnificandZ_d += 2.0;
        //}
        //else
        //{
        //    sicnificandZ_d += 1.0;
        //}

        std::int64_t exponent_z = (std::int64_t)z_encoded - 1;
        double sicnificandZ_d = z_encoded - exponent_z;

        if (sicnificandZ_d >= 1.0) {
            exponent_z += 1LL;
            sicnificandZ_d -= 1.0;
        }

        //****************

         //__m128d encoded_2 = { encoded, z.encoded };
         //uint64_t exponent_2[2];
         //uint64_t sicnificand_2[2];

         //convert_double_to_uint64(encoded_2, exponent_2, sicnificand_2);
         //uint64_t exponent = exponent_2[0];
         //uint64_t exponent_z = exponent_2[1];
         //uint64_t sicnificand = sicnificand_2[0];
         //uint64_t sicnificand_z = sicnificand_2[1];

        //****************

        double sicnificand_d = scnfcnd;
        std::int64_t exponent = exp;

        //if (encoded + z_encoded >= 1.79769e+100) {
        //    sicnificand_d = 1000000;
        //    exponent = 10000000;
        //    z_encoded = 23453253;
        //}

        //double sicnificand_d;
        //std::int64_t exponent;
        //double floored_1 = std::floor(encoded);
        //exponent = (std::int64_t)floored_1;
        //sicnificand_d = encoded - floored_1 + 1.0;

        //double sicnificandZ_d;
        //std::int64_t exponent_z;
        //double floored_2 = std::floor(z_encoded);
        //exponent_z = (std::int64_t)floored_2;
        //sicnificandZ_d = z_encoded - floored_2 + 1.0;

        std::int64_t sicnificand = std::bit_cast<std::int64_t>(sicnificand_d);
        std::int64_t sicnificand_z = std::bit_cast<std::int64_t>(sicnificandZ_d);

        //****************

        std::int64_t exp_diff = (std::int64_t)(exponent - exponent_z);

        if (exp_diff > 0) {
            if (exp_diff > 63) {
                return *this;
            }
            sicnificand_z -= exp_diff << 52;
        }
        else if (exp_diff < 0) {
            if (exp_diff < -63) {
                encoded = z_encoded;
                double sicnificandZ_d = std::bit_cast<double>(sicnificand_z);
                scnfcnd = sicnificandZ_d;
                exp = exponent_z;
                return *this;
            }
            sicnificand += exp_diff << 52;
            exponent = exponent_z;
        }

        sicnificand_d = std::bit_cast<double>(sicnificand);
        sicnificandZ_d = std::bit_cast<double>(sicnificand_z);

        //sicnificand_d = std::bit_cast<double>(sicnificand);
        //sicnificandZ_d = std::bit_cast<double>(sicnificand_z);

        sicnificand_d += sicnificandZ_d;

        if (sicnificand_d >= 2.0) {
            sicnificand_d *= 0.5;
            exponent++;
        }

        encode(sicnificand_d, exponent);
        scnfcnd = sicnificand_d;
        exp = exponent;
        return *this;
    }

    Float64ExtendedExp& Float64ExtendedExp::operator+=(Float64ExtendedExp z) {
        std::uint64_t sicnificand;
        std::int64_t exponent;
        decode(sicnificand, exponent, 0.0);

        std::uint64_t sicnificand_z;
        std::int64_t exponent_z;
        z.decode(sicnificand_z, exponent_z, 0.0);

        //****************

         //__m128d encoded_2 = { encoded, z.encoded };
         //uint64_t exponent_2[2];
         //uint64_t sicnificand_2[2];

         //convert_double_to_uint64(encoded_2, exponent_2, sicnificand_2);
         //uint64_t exponent = exponent_2[0];
         //uint64_t exponent_z = exponent_2[1];
         //uint64_t sicnificand = sicnificand_2[0];
         //uint64_t sicnificand_z = sicnificand_2[1];

        //****************

        //double sicnificand_d;
        //std::int64_t exponent;
        //exponent = (std::int64_t)std::floor(encoded);
        //sicnificand_d = encoded - exponent + 1.0;

        //double sicnificandZ_d;
        //std::int64_t exponent_z;
        //exponent_z = (std::int64_t)std::floor(z.encoded);
        //sicnificandZ_d = z.encoded - exponent_z + 1.0;

        //std::int64_t sicnificand = std::bit_cast<std::int64_t>(sicnificand_d);
        //std::int64_t sicnificand_z = std::bit_cast<std::int64_t>(sicnificandZ_d);

        //****************

        std::int64_t exp_diff = (std::int64_t)(exponent - exponent_z);

        if (exp_diff > 0) {
            if (exp_diff > 63) {
                return *this;
            }
            sicnificand_z -= exp_diff << 52;
        }
        else if (exp_diff < 0){
            if (exp_diff < -63) {
                encoded = z.encoded;
                return *this;
            }
            sicnificand += exp_diff << 52;
            exponent = exponent_z;
        }
        
        double sicnificand_d = std::bit_cast<double>(sicnificand);
        double sicnificandZ_d = std::bit_cast<double>(sicnificand_z);
        
        //sicnificand_d = std::bit_cast<double>(sicnificand);
        //sicnificandZ_d = std::bit_cast<double>(sicnificand_z);

        sicnificand_d += sicnificandZ_d;

        if (sicnificand_d >= 2.0) {
            sicnificand_d *= 0.5;
            exponent++;
        }

        encode(sicnificand_d, exponent);
        return *this;
    }

    // Float64ExtendedExp& Float64ExtendedExp::operator*=(Float64ExtendedExp z) {
    //     std::int64_t exponent = (std::int64_t)encoded;
    //     std::int64_t exponentZ = (std::int64_t)z.encoded;
    //     double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
    //     double sicnificandZ = exponentZ < 0 ? -1*(z.encoded - exponentZ) : z.encoded - exponentZ;
    //     exponent += exponentZ;
    //     sicnificand *= sicnificandZ;

    //     if (sicnificand < 0.5) {
    //         sicnificand *= 2;
    //         exponent--;
    //     }

    //     encoded = exponent < 0 ? (double)exponent - sicnificand : (double)exponent + sicnificand;

    //     return *this;
    // }

    void Float64ExtendedExp::printBinary(std::string message, uint64_t value) const {
        if (message.length() > 0) {
            std::cout << message << std::endl;
        }
        for (int i = 56; i >= 0; i -= 8) {
            std::bitset<8> block((value >> i) & 0xFF); // Extract 8 bits
            std::cout << block << " ";                  // Print block with a space
        }
        std::cout << std::endl;
    }

    void Float64ExtendedExp::print(std::string message, int64_t value) const {        
        std::cout << message << value << std::endl;
    }

    void Float64ExtendedExp::print(std::string message, double value) const {        
        std::cout << message << value << std::endl;
    }

    inline void Float64ExtendedExp::encode(double sicnificand, std::int64_t exponent) {
        encoded = (double)exponent + (sicnificand - 1.0);
    }

    inline void Float64ExtendedExp::decode(std::uint64_t& sicnificand, std::int64_t& exponent, double enco) {
        //exponent = (std::int64_t)std::floor(encoded);
        //sicnificand = encoded - exponent + 1.0;

        double encoded_value = enco == 0.0 ? encoded : enco;

        uint64_t dbl_bits = std::bit_cast<uint64_t>(encoded_value);
        std::int32_t cutter = ((dbl_bits >> 52) & 0x7FF) - 1023;
        std::int64_t full_mantissa = (dbl_bits & 0x000FFFFFFFFFFFFFull) | 0x0010000000000000ull;

        if (cutter >= 0) {
            exponent = ((std::int64_t)full_mantissa >> (52 - cutter));
            sicnificand = (dbl_bits << cutter) & 0x000FFFFFFFFFFFFFull;
        }
        else {
            exponent = 0;
            sicnificand = full_mantissa >> (-1 * cutter);
        }

        exponent = exponent * (encoded_value < 0.0 ? -1LL : 1LL);
        exponent = (encoded_value < 0.0 && sicnificand > 0ULL) ? exponent - 1 : exponent;

        sicnificand = encoded_value < 0.0 ? (0x0020000000000000ull - sicnificand) : sicnificand;
        sicnificand = sicnificand | 0x3FF0000000000000ull;
    }

    inline void Float64ExtendedExp::encode_double(double dbl, std::int64_t exponent) {
        uint64_t* dbl_bits = reinterpret_cast<std::uint64_t*>(&dbl);
        exponent += ((((*dbl_bits & 0x7FF0000000000000ull)) >> 52) - 1023LL);
        *dbl_bits &= 0x800FFFFFFFFFFFFFull;
        *dbl_bits |= 0x3FF0000000000000ull;

        encode(dbl, exponent);
    }

    inline void Float64ExtendedExp::decode_double(double& sicnificand, std::int64_t& exponent) {
        std::uint64_t scnfcnd;
        std::int64_t exp;
        decode(scnfcnd, exp, 0.0);

        sicnificand = std::bit_cast<double>(scnfcnd);
        exponent = exp;

        //exponent = (std::int64_t)std::floor(encoded);
        //sicnificand = encoded - exponent + 1.0;
    }

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b) { return a+=b; }
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b) { return a*=b; }
}

