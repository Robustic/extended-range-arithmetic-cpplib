#include <cstdint>
#include <cmath>
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
        this->encode_double(1.0, 0LL);
    }

    Float64ExtendedExp::Float64ExtendedExp(double dbl) {
        this->encode_double(dbl, 0LL);
    }

    // Float64ExtendedExp::Float64ExtendedExp(double sicnificand, std::int64_t exponent) {
    //     this->scale(sicnificand, exponent);
    // }

    // Float64ExtendedExp::Float64ExtendedExp(std::vector<floatingExp2Integer::Float64ExtendedExp>& vector) {
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
    // }

    void Float64ExtendedExp::doubleToFloat64ExtendedExp(double dbl) {
        this->encode_double(dbl, 0LL);
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
        this->decode_double(sicnificand, exponent);        
        return sicnificand; 
    }

    std::int64_t Float64ExtendedExp::exponent() {
        std::int64_t exponent;
        double sicnificand;
        this->decode_double(sicnificand, exponent);
        return exponent;
    }

    double Float64ExtendedExp::asDouble() const {
        std::int64_t exponent;
        double sicnificand;
        this->decode_double(sicnificand, exponent);
        return sicnificand * std::pow(2.0, (int)exponent);
    }

    Float64ExtendedExp& Float64ExtendedExp::operator+=(Float64ExtendedExp z) {
        std::uint64_t sicnificand;
        std::int64_t exponent;
        decode(sicnificand, exponent);

        std::uint64_t sicnificandZ;
        std::int64_t exponentZ;
        z.decode(sicnificandZ, exponentZ);

        std::int64_t exp_diff = (std::int64_t)(exponent - exponentZ);

        if (exp_diff > 0) {
            if (exp_diff > 63) {
                return *this;
            }
            sicnificandZ = sicnificandZ >> exp_diff;
        }
        else if (exp_diff < 0){
            if (exp_diff < -63) {
                encoded = z.encoded;
                return *this;
            }
            sicnificand = sicnificand >> (-exp_diff);
            exponent = exponentZ;
        }

        sicnificand += sicnificandZ;
        if (sicnificand >= 0x8000000000000000) {
            sicnificand = sicnificand >> 1;
            exponent++;
        }

        encode(sicnificand, exponent);
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

    inline void Float64ExtendedExp::encode(std::uint64_t sgnfcnd_bits, std::int64_t exponent) {
        uint64_t sign_of_exponent = 0x0ULL;
        //print("exponent: ", exponent);
        if (exponent < 0) {
            exponent = -1 * exponent;
            sign_of_exponent = 0x8000000000000000ULL;
        }
        //print("exponent: ", exponent);

        std::uint32_t length_of_exponent = 70 - __builtin_clzll(exponent);

        uint64_t return_bits = 0x0ULL;
        return_bits |= exponent << 6;
        //printBinary("", return_bits);

        return_bits |= sign_of_exponent;
        //printBinary("", return_bits);
        
        //print("length_of_exponent: ", length_of_exponent);
        return_bits |= length_of_exponent;
        //printBinary("", return_bits);
        
        sgnfcnd_bits &= 0xFFFFFFFFFFFFFFC0ull << length_of_exponent;
        //printBinary("", sgnfcnd_bits);

        return_bits |= sgnfcnd_bits;
        //printBinary("return_bits: ", return_bits);

        encoded = return_bits;
    }

    inline void Float64ExtendedExp::decode(std::uint64_t& significand_bits, std::int64_t& exponent) const {
        //printBinary("decode", encoded);
        // __m256i encoded_4 = _mm256_set1_epi64x(encoded);
        // __m256i mask_4 = _mm256_set_epi64x(
        //     0ULL,
        //     0x7FFFFFFFFFFFFFFFull,
        //     0x3FULL,
        //     0x8000000000000000ULL
        // );
        // __m256i result = _mm256_and_si256(encoded_4, mask_4);

        uint64_t sign_of_exponent = encoded & 0x8000000000000000ULL;
        uint32_t separator_location = encoded & 0x3FULL;
        uint64_t temp = encoded & 0x7FFFFFFFFFFFFFFFull;
        uint64_t cut_mirror = 0xFFFFFFFFFFFFFFC0ull << separator_location;
        //print("separator_location: ", separator_location);
        exponent = (encoded & (~cut_mirror)) >> 6;
        if (sign_of_exponent != 0) {
            exponent = -1 * exponent;
        }
        //print("exponent: ", exponent);
        significand_bits = temp & cut_mirror;
        //printBinary("decoded_significand", significand_bits);
    }

    inline void Float64ExtendedExp::encode_double(double sicnificand, std::int64_t exponent) {
        //std::cout << "Double: " << sicnificand << std::endl;
        uint64_t sgnfcnd_bits = *reinterpret_cast<std::uint64_t*>(&sicnificand);
        exponent += ((((sgnfcnd_bits & 0x7FF0000000000000ull)) >> 52) - 1023LL);
        //print("Double: exponent: ", exponent);
        //printBinary("original_significand", sgnfcnd_bits);
        sgnfcnd_bits &= 0x000FFFFFFFFFFFFFull; // significand sign is forced to zero
        sgnfcnd_bits = sgnfcnd_bits << 10;
        sgnfcnd_bits |= 0x4000000000000000ull;
        //printBinary("sgnfcnd_bits_to_encode:", sgnfcnd_bits);

        encode(sgnfcnd_bits, exponent);
        //printBinary("encoded_bits:", encoded);
    }

    inline void Float64ExtendedExp::decode_double(double& sicnificand, std::int64_t& exponent) const {
        uint64_t sgnfcnd_bits;
        //printBinary("bits to Decode: ", encoded);
        decode(sgnfcnd_bits, exponent);
        sgnfcnd_bits &= 0x3FFFFFFFFFFFFFFFull;
        sgnfcnd_bits = sgnfcnd_bits >> 10;
        sgnfcnd_bits |= 0x3FF0000000000000ull;
        uint64_t print2 = sgnfcnd_bits;
        //printBinary("sgnfcnd_bits_before_return: ", print2);
        sicnificand = *reinterpret_cast<double*>(&sgnfcnd_bits);
        //std::cout << "Double (decode_double): " << sicnificand << std::endl;
    }

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b) { return a+=b; }
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b) { return a*=b; }
}

