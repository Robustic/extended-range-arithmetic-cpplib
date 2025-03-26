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
        __m512d sicnificand_sum = encoded_values - integer_decimals + one_d;

        //convert_double_to_uint64(encoded_values, exponent, sicnificand);

        for (int i = 8; i + 7 < n; i += 8) {
            encoded_values = _mm512_loadu_pd(&vector_in[i]);

            integer_decimals = _mm512_floor_pd(encoded_values);
            __m512i exponent = _mm512_cvtpd_epi64(integer_decimals);
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
            __m512d zero_d = _mm512_setzero_pd();
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

            __m512i one_i = _mm512_set1_epi64(1);
            __m512d half_d = _mm512_set1_pd(0.5);
            __m512d two_d = _mm512_set1_pd(2.0);
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

    Float64ExtendedExp& Float64ExtendedExp::operator+=(Float64ExtendedExp z) {
        std::uint64_t sicnificand;
        std::int64_t exponent;
        decode(sicnificand, exponent);

        std::uint64_t sicnificand_z;
        std::int64_t exponent_z;
        z.decode(sicnificand_z, exponent_z);

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

    inline void Float64ExtendedExp::convert_double_to_uint64(__m512d encoded_values, __m512i exponent, __m512i sicnificand) {
        //__m512i dbl_bits = _mm512_castpd_si512(encoded_values);  // Bit-cast double to integer representation

        //__m512i exp_mask = _mm512_set1_epi64(0x7FF0000000000000ull);  // Extract exponent mask
        //__m512i mant_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFull);  // Extract mantissa mask
        //__m512i bias = _mm512_set1_epi64(1023);

        //__m512i exp_bits = _mm512_and_si512(dbl_bits, exp_mask);
        //exp_bits = _mm512_srli_epi64(exp_bits, 52);  // Shift exponent into position
        //__m512i cutter = _mm512_sub_epi64(exp_bits, bias);  // Subtract bias to get true exponent

        //// Extract sign bit to determine negative numbers
        //__m512i sign_mask = _mm512_set1_epi64(0x8000000000000000ull);
        //__m512i sign = _mm512_and_si512(dbl_bits, sign_mask);

        //__m512i mantissa = _mm512_and_si512(dbl_bits, mant_mask);
        //__m512i implicit_bit = _mm512_set1_epi64(0x0010000000000000ull);
        //__m512i full_mantissa = _mm512_or_si512(mantissa, implicit_bit);  // Restore implicit bit

        //// Compute exponent and sicnificand based on cutter value
        //__m512i zero = _mm512_setzero_si512();
        //__m512i i52 = _mm512_set1_epi64(52);
        //__m512i im1 = _mm512_set1_epi64(-1);
        //__m512i ip1 = _mm512_set1_epi64(1);
        //__m512i sign_multiplier = _mm512_mask_blend_epi64(_mm512_cmpeq_epu64_mask(sign, sign_mask), im1, ip1);

        //__m512i shifted_mantissa_left = _mm512_sllv_epi64(full_mantissa, cutter);
        //__m512i shifted_mantissa_right = _mm512_srlv_epi64(full_mantissa, zero - cutter);
        //__m512i shifted_mantissa = cutter < zero ? shifted_mantissa_right : shifted_mantissa_left;
        //shifted_mantissa = _mm512_and_si512(shifted_mantissa, mant_mask);

        //__m512i shifted_exponent = _mm512_srlv_epi64(full_mantissa, i52 - cutter);
        //shifted_exponent = _mm512_max_epi64(shifted_exponent, zero); // !!! change to _mm_max_epi64
        //shifted_exponent = shifted_exponent * sign_multiplier;

        //__m512i condition = _mm512_and_si512(_mm512_cmpgt_epi64_mask(shifted_mantissa, _mm512_setzero_si512()), _mm512_cmpeq_epu64_mask(sign, sign_mask));
        //shifted_exponent = _mm_sub_epi64(shifted_exponent, _mm_and_si128(condition, ip1));

        //__m512i negative_offset = _mm_set1_epi64x(0x0020000000000000ull);

        //shifted_mantissa = _mm512_mask_blend_epi64(_mm512_cmpeq_epu64_mask(sign, sign_mask), negative_offset - shifted_mantissa, shifted_mantissa);
        //__m512i zero_exp_mask = _mm512_set1_epi64(0x3FF0000000000000ull);
        //shifted_mantissa = _mm512_or_si512(shifted_mantissa, zero_exp_mask);

        //// Store results
        //_mm512_storeu_si512(sicnificand, shifted_mantissa);
        //_mm512_storeu_si512(exponent, shifted_exponent);
    }

    inline void Float64ExtendedExp::convert_double_to_uint64(__m128d encoded_values, uint64_t* exponent, uint64_t* sicnificand) {
        __m128i dbl_bits = _mm_castpd_si128(encoded_values);  // Bit-cast double to integer representation
    
        __m128i exp_mask = _mm_set1_epi64x(0x7FF0000000000000ull);  // Extract exponent mask
        __m128i mant_mask = _mm_set1_epi64x(0x000FFFFFFFFFFFFFull);  // Extract mantissa mask
        __m128i bias = _mm_set1_epi64x(1023);
    
        __m128i exp_bits = _mm_and_si128(dbl_bits, exp_mask);
        exp_bits = _mm_srli_epi64(exp_bits, 52);  // Shift exponent into position
        __m128i cutter = _mm_sub_epi64(exp_bits, bias);  // Subtract bias to get true exponent
    
        // Extract sign bit to determine negative numbers
        __m128i sign_mask = _mm_set1_epi64x(0x8000000000000000ull);
        __m128i sign = _mm_and_si128(dbl_bits, sign_mask);
    
        __m128i mantissa = _mm_and_si128(dbl_bits, mant_mask);
        __m128i implicit_bit = _mm_set1_epi64x(0x0010000000000000ull);  
        __m128i full_mantissa = _mm_or_si128(mantissa, implicit_bit);  // Restore implicit bit
    
        // Compute exponent and sicnificand based on cutter value
        __m128i zero = _mm_setzero_si128();
        __m128i i52 = _mm_set1_epi64x(52);
        __m128i im1 = _mm_set1_epi64x(-1);
        __m128i ip1 = _mm_set1_epi64x(1);
        __m128i sign_multiplier = _mm_blendv_epi8(ip1, im1, _mm_cmpeq_epi64(sign, sign_mask));

        __m128i shifted_mantissa_left = _mm_sllv_epi64(full_mantissa, cutter);
        __m128i shifted_mantissa_right = _mm_srlv_epi64(full_mantissa, zero - cutter);
        __m128i shifted_mantissa = cutter < zero ? shifted_mantissa_right : shifted_mantissa_left;
        shifted_mantissa = _mm_and_si128(shifted_mantissa, mant_mask);

        __m128i shifted_exponent = _mm_srlv_epi64(full_mantissa, i52 - cutter);
        shifted_exponent = _mm_max_epi64(shifted_exponent, zero); // !!! change to _mm_max_epi64
        shifted_exponent = shifted_exponent * sign_multiplier;
    
        __m128i condition = _mm_and_si128(_mm_cmpgt_epi64(shifted_mantissa, _mm_setzero_si128()), _mm_cmpeq_epi64(sign, sign_mask));
        shifted_exponent = _mm_sub_epi64(shifted_exponent, _mm_and_si128(condition, ip1));

        __m128i negative_offset = _mm_set1_epi64x(0x0020000000000000ull);

        shifted_mantissa = _mm_blendv_epi8(shifted_mantissa, negative_offset - shifted_mantissa, _mm_cmpeq_epi64(sign, sign_mask));
        __m128i zero_exp_mask = _mm_set1_epi64x(0x3FF0000000000000ull);
        shifted_mantissa = _mm_or_si128(shifted_mantissa, zero_exp_mask);
    
        // Store results
        _mm_storeu_si128(reinterpret_cast<__m128i*>(sicnificand), shifted_mantissa);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(exponent), shifted_exponent);
    }
    

    inline void Float64ExtendedExp::decode(std::uint64_t& sicnificand, std::int64_t& exponent) {
        //exponent = (std::int64_t)std::floor(encoded);
        //sicnificand = encoded - exponent + 1.0;

        uint64_t dbl_bits = std::bit_cast<uint64_t>(encoded);
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

        exponent = exponent * (encoded < 0.0 ? -1LL : 1LL);
        exponent = (encoded < 0.0 && sicnificand > 0ULL) ? exponent - 1 : exponent;

        sicnificand = encoded < 0.0 ? (0x0020000000000000ull - sicnificand) : sicnificand;
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
        decode(scnfcnd, exp);

        sicnificand = std::bit_cast<double>(scnfcnd);
        exponent = exp;

        //exponent = (std::int64_t)std::floor(encoded);
        //sicnificand = encoded - exponent + 1.0;
    }

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b) { return a+=b; }
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b) { return a*=b; }
}

