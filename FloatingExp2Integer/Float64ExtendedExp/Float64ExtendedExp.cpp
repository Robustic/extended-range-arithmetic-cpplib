#include <cstdint>
#include <cmath>
#include <vector>
#include <cstring>
#include "Float64ExtendedExp.h"

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
typedef uint64_t uint64_4_t __attribute__((vector_size(4 * sizeof(uint64_t))));
typedef int64_t int64_4_t __attribute__((vector_size(4 * sizeof(int64_t))));

constexpr int64_4_t int64_4_1023 { 1023LL, 1023LL, 1023LL, 1023LL };
constexpr double4_t double4_0 { 0.0, 0.0, 0.0, 0.0 };

namespace floatingExp2Integer
{
    Float64ExtendedExp::Float64ExtendedExp() {
        this->scale(0.5, 0LL);
    }

    Float64ExtendedExp::Float64ExtendedExp(double dbl) {
        this->scale(dbl, 0LL);
    }

    Float64ExtendedExp::Float64ExtendedExp(double sicnificand, std::int64_t exponent) {
        this->scale(sicnificand, exponent);
    }

    Float64ExtendedExp::Float64ExtendedExp(std::vector<floatingExp2Integer::Float64ExtendedExp>& vector) {
        // if (vector.size() < 4 * 2) {
        //     Float64ExtendedExp a(vector[0].encoded);
        //     for (unsigned int i = 0; i < vector.size(); i++) {
        //         a += vector[i];
        //     }    
        //     encoded = a.encoded;
        // }

        // double4_t sa1;
        // int64_4_t ea1;
        // double4_t sa2;
        // int64_4_t ea2;

        // double4_t sb1;
        // int64_4_t eb1;
        // double4_t sb2;
        // int64_4_t eb2;
    
        // for (int k = 0; k < 4; k++) {
        //     sa1[k] = vector[k].encoded;
        //     sa2[k] = vector[k + 4].encoded;
        //     ea1[k] = vector[k].exp;
        //     ea2[k] = vector[k + 4].exp;
        // }

        // const auto* vec_ptr = vector.data();
    
        // int counter = 0;
        // unsigned int i;
        // for (i = 8; i + 7 < vector.size(); i += 8) {
    
        //     for (int k = 0; k < 4; k++) {
        //         int ik = i + k;
        //         int ik4 = ik + 4;
        //         sb1[k] = vector[ik].encoded;
        //         sb2[k] = vector[ik4].encoded;
        //         eb1[k] = vector[ik].exp;
        //         eb2[k] = vector[ik4].exp;
        //     }
    
        //     int64_4_t exp_diff1 = ea1 - eb1;
        //     int64_4_t exp_diff2 = ea2 - eb2;

        //     uint64_4_t& ib1 = reinterpret_cast<uint64_4_t&>(sb1); 
        //     uint64_4_t& ib2 = reinterpret_cast<uint64_4_t&>(sb2);
        //     ib1 = ((((ib1 & 0x7FF0000000000000ull) >> 52) - exp_diff1) << 52) | (ib1 & 0x800FFFFFFFFFFFFFull);
        //     ib2 = ((((ib2 & 0x7FF0000000000000ull) >> 52) - exp_diff2) << 52) | (ib2 & 0x800FFFFFFFFFFFFFull);
    
        //     sa1 = exp_diff1 <= -64LL ? double4_0 : sa1;
        //     sb1 = exp_diff1 >= 64LL ? double4_0 : sb1;
        //     sa2 = exp_diff2 <= -64LL ? double4_0 : sa2;
        //     sb2 = exp_diff2 >= 64LL ? double4_0 : sb2;
    
        //     sa1 += sb1;
        //     sa2 += sb2;

        //     counter++;
        //     if (counter > 30) {
        //         counter = 0;
        //         uint64_4_t& ia1 = reinterpret_cast<uint64_4_t&>(sa1); 
        //         uint64_4_t& ia2 = reinterpret_cast<uint64_4_t&>(sa2);
        //         ea1 = (ea1 - int64_4_1023) + (int64_4_t)((ia1 & 0x7FF0000000000000ull) >> 52);
        //         ea2 = (ea2 - int64_4_1023) + (int64_4_t)((ia2 & 0x7FF0000000000000ull) >> 52);
        //         ia1 &= 0x800FFFFFFFFFFFFFull;
        //         ia1 |= 0x3FF0000000000000ull;
        //         ia2 &= 0x800FFFFFFFFFFFFFull;
        //         ia2 |= 0x3FF0000000000000ull;
        //     }
        // }        
        
        // Float64ExtendedExp a1(sa1[0], ea1[0]);
        // Float64ExtendedExp a2(sa2[0], ea2[0]);
        // for (int k = 1; k < 4; k++) {
        //     Float64ExtendedExp z1(sa1[k], ea1[k]);
        //     Float64ExtendedExp z2(sa2[k], ea2[k]);
        //     a1 += z1;
        //     a2 += z2;
        // }
        // a1 += a2;

        // for (i = i; i < vector.size(); i++) {
        //     a1 += vector[i];
        // }

        // encoded = a1.encoded;
    }

    void Float64ExtendedExp::doubleToFloat64ExtendedExp(double dbl) {
        this->scale(dbl, 0LL);
    }

    void Float64ExtendedExp::log2ToFloat64ExtendedExp(double log2) {
        std::int64_t exponent = (std::int64_t)log2;
        double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
        this->scale(sicnificand, exponent);
    }

    double Float64ExtendedExp::float64ExtendedExpToLog2() const {
        std::int64_t exponent = (std::int64_t)encoded;
        double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
        return std::log2(sicnificand) + exponent;
    }

    double Float64ExtendedExp::sicnificand() {
        std::int64_t exponent = (std::int64_t)encoded;
        double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
        return sicnificand; 
    }

    std::int64_t Float64ExtendedExp::exponent() {
        std::int64_t exponent = (std::int64_t)encoded;
        return exponent;
    }

    double Float64ExtendedExp::asDouble() const {
        std::int64_t exponent = (std::int64_t)encoded;
        double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
        return sicnificand * std::pow(2.0, exponent);
    }

    Float64ExtendedExp& Float64ExtendedExp::operator+=(Float64ExtendedExp z) {
        std::int64_t exponent = (std::int64_t)encoded;
        double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;

        std::int64_t exponentZ = (std::int64_t)z.encoded;
        double sicnificandZ = exponentZ < 0 ? -1*(z.encoded - exponentZ) : z.encoded - exponentZ;

        std::int64_t exp_diff = (std::int64_t)(exponent - exponentZ);

        std::uint64_t* sgnfcnd_bits = reinterpret_cast<std::uint64_t*>(&sicnificand);
        std::uint64_t* sgnfcnd_bits_z = reinterpret_cast<std::uint64_t*>(&sicnificandZ);

        if (exp_diff > 0) {
            if (exp_diff > 52) {
                return *this;
            }
            *sgnfcnd_bits_z -= exp_diff << 52;
        }
        else if (exp_diff < 0){
            if (exp_diff < -52) {
                encoded = z.encoded;
                return *this;
            }
            *sgnfcnd_bits -= (-exp_diff) << 52;
            exponent = exponentZ;
        }

        sicnificand += sicnificandZ;
        if (sicnificand >= 1.0) {
            sicnificand /= 2.0;
            exponent++;
        }

        encoded = exponent < 0 ? (double)exponent - sicnificand : (double)exponent + sicnificand;
        return *this;
    }

    Float64ExtendedExp& Float64ExtendedExp::operator*=(Float64ExtendedExp z) {
        std::int64_t exponent = (std::int64_t)encoded;
        std::int64_t exponentZ = (std::int64_t)z.encoded;
        double sicnificand = exponent < 0 ? -1*(encoded - exponent) : encoded - exponent;
        double sicnificandZ = exponentZ < 0 ? -1*(z.encoded - exponentZ) : z.encoded - exponentZ;
        exponent += exponentZ;
        sicnificand *= sicnificandZ;

        if (sicnificand < 0.5) {
            sicnificand *= 2;
            exponent--;
        }

        encoded = exponent < 0 ? (double)exponent - sicnificand : (double)exponent + sicnificand;

        return *this;
    }

    inline void Float64ExtendedExp::scale(double sicnificand, std::int64_t exponent) {
        uint64_t* sgnfcnd_bits = reinterpret_cast<uint64_t*>(&sicnificand);
        exponent += ((*sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1022;
        *sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
        *sgnfcnd_bits |= 0x3FE0000000000000ull;
        sicnificand = *reinterpret_cast<double*>(sgnfcnd_bits);
        encoded = exponent < 0 ? (double)exponent - sicnificand : (double)exponent + sicnificand;
    }

    Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b) { return a+=b; }
    Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b) { return a*=b; }
}

