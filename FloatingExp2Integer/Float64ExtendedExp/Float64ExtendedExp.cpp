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

    // Float64ExtendedExp::Float64ExtendedExp(double sicnificand, std::int64_t exponent) {
    //     this->scale(sicnificand, exponent);
    // }

    Float64ExtendedExp::Float64ExtendedExp(const std::vector<floatingExp2Integer::Float64ExtendedExp>& vector) {
    //     // if (vector.size() < 4 * 2) {
    //     //     Float64ExtendedExp a(vector[0].encoded);
    //     //     for (unsigned int i = 0; i < vector.size(); i++) {
    //     //         a += vector[i];
    //     //     }    
    //     //     encoded = a.encoded;
    //     // }

    //     // double4_t sa1;
    //     // int64_4_t ea1;
    //     // double4_t sa2;
    //     // int64_4_t ea2;

    //     // double4_t sb1;
    //     // int64_4_t eb1;
    //     // double4_t sb2;
    //     // int64_4_t eb2;
    
    //     // for (int k = 0; k < 4; k++) {
    //     //     sa1[k] = vector[k].encoded;
    //     //     sa2[k] = vector[k + 4].encoded;
    //     //     ea1[k] = vector[k].exp;
    //     //     ea2[k] = vector[k + 4].exp;
    //     // }

    //     // const auto* vec_ptr = vector.data();
    
    //     // int counter = 0;
    //     // unsigned int i;
    //     // for (i = 8; i + 7 < vector.size(); i += 8) {
    
    //     //     for (int k = 0; k < 4; k++) {
    //     //         int ik = i + k;
    //     //         int ik4 = ik + 4;
    //     //         sb1[k] = vector[ik].encoded;
    //     //         sb2[k] = vector[ik4].encoded;
    //     //         eb1[k] = vector[ik].exp;
    //     //         eb2[k] = vector[ik4].exp;
    //     //     }
    
    //     //     int64_4_t exp_diff1 = ea1 - eb1;
    //     //     int64_4_t exp_diff2 = ea2 - eb2;

    //     //     uint64_4_t& ib1 = reinterpret_cast<uint64_4_t&>(sb1); 
    //     //     uint64_4_t& ib2 = reinterpret_cast<uint64_4_t&>(sb2);
    //     //     ib1 = ((((ib1 & 0x7FF0000000000000ull) >> 52) - exp_diff1) << 52) | (ib1 & 0x800FFFFFFFFFFFFFull);
    //     //     ib2 = ((((ib2 & 0x7FF0000000000000ull) >> 52) - exp_diff2) << 52) | (ib2 & 0x800FFFFFFFFFFFFFull);
    
    //     //     sa1 = exp_diff1 <= -64LL ? double4_0 : sa1;
    //     //     sb1 = exp_diff1 >= 64LL ? double4_0 : sb1;
    //     //     sa2 = exp_diff2 <= -64LL ? double4_0 : sa2;
    //     //     sb2 = exp_diff2 >= 64LL ? double4_0 : sb2;
    
    //     //     sa1 += sb1;
    //     //     sa2 += sb2;

    //     //     counter++;
    //     //     if (counter > 30) {
    //     //         counter = 0;
    //     //         uint64_4_t& ia1 = reinterpret_cast<uint64_4_t&>(sa1); 
    //     //         uint64_4_t& ia2 = reinterpret_cast<uint64_4_t&>(sa2);
    //     //         ea1 = (ea1 - int64_4_1023) + (int64_4_t)((ia1 & 0x7FF0000000000000ull) >> 52);
    //     //         ea2 = (ea2 - int64_4_1023) + (int64_4_t)((ia2 & 0x7FF0000000000000ull) >> 52);
    //     //         ia1 &= 0x800FFFFFFFFFFFFFull;
    //     //         ia1 |= 0x3FF0000000000000ull;
    //     //         ia2 &= 0x800FFFFFFFFFFFFFull;
    //     //         ia2 |= 0x3FF0000000000000ull;
    //     //     }
    //     // }        
        
    //     // Float64ExtendedExp a1(sa1[0], ea1[0]);
    //     // Float64ExtendedExp a2(sa2[0], ea2[0]);
    //     // for (int k = 1; k < 4; k++) {
    //     //     Float64ExtendedExp z1(sa1[k], ea1[k]);
    //     //     Float64ExtendedExp z2(sa2[k], ea2[k]);
    //     //     a1 += z1;
    //     //     a2 += z2;
    //     // }
    //     // a1 += a2;

    //     // for (i = i; i < vector.size(); i++) {
    //     //     a1 += vector[i];
    //     // }

    //     // encoded = a1.encoded;
    }

    void Float64ExtendedExp::doubleToFloat64ExtendedExp(double dbl) {
        encode_double(dbl, 0LL);
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

        // __m128d encoded_2 = { encoded, z.encoded };
        // uint64_t exponent_2[2];
        // uint64_t sicnificand_2[2];

        // convert_double_to_uint64(encoded_2, exponent_2, sicnificand_2);
        // uint64_t exponent = exponent_2[0];
        // uint64_t exponent_z = exponent_2[1];
        // uint64_t sicnificand = sicnificand_2[0];
        // uint64_t sicnificand_z = sicnificand_2[1];

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
        shifted_exponent = _mm_max_epi32(shifted_exponent, zero); // !!! change to _mm_max_epi64
        shifted_exponent = shifted_exponent * sign_multiplier;
    
        __m128i condition = _mm_and_si128(_mm_cmpgt_epi64(shifted_mantissa, _mm_setzero_si128()), _mm_cmpeq_epi64(sign, sign_mask));
        shifted_exponent = _mm_sub_epi64(shifted_exponent, _mm_and_si128(condition, ip1));

        __m128i negative_offset = _mm_set1_epi64x(0x0020000000000000ull);

        shifted_mantissa = _mm_blendv_epi8(negative_offset - shifted_mantissa, shifted_mantissa, _mm_cmpeq_epi64(sign, sign_mask));
        __m128i zero_exp_mask = _mm_set1_epi64x(0x3FF0000000000000ull);
        shifted_mantissa = _mm_or_si128(shifted_mantissa, zero_exp_mask);
    
        // Store results
        _mm_storeu_si128(reinterpret_cast<__m128i*>(exponent), shifted_exponent);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(sicnificand), shifted_mantissa);
    }
    

    inline void Float64ExtendedExp::decode(std::uint64_t& sicnificand, std::int64_t& exponent) {
        // exponent = (std::int64_t)std::floor(encoded);
        // sicnificand = encoded - exponent + 1.0;

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
    }

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b) { return a+=b; }
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b) { return a*=b; }
}

