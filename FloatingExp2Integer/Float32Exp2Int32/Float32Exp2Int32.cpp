#include <cstdint>
#include <cmath>
#include <cstring>
// #include <limits>
#include <iostream>
#include "Float32Exp2Int32.h"

namespace floatingExp2Integer
{
    Float32Exp2Int32::Float32Exp2Int32() {
        scnfcnd = 1.0;
        exp = 0;
        this->scale();
    }

    Float32Exp2Int32::Float32Exp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->scale();
    }

    Float32Exp2Int32::Float32Exp2Int32(float sicnificand, std::int32_t exponent) {
        scnfcnd = sicnificand;
        exp = exponent;
        this->scale();
    }

    void Float32Exp2Int32::log2ToFloat32Exp2Int32(float logarithm2) {
        exp = (std::int32_t)logarithm2;
        scnfcnd = std::exp2(logarithm2 - exp);
        this->scale();
    }

    float Float32Exp2Int32::float32Exp2Int32ToLog2() {
        return std::log2(scnfcnd) + exp ;
    }


    float Float32Exp2Int32::sicnificand() { 
        this->scale();
        return scnfcnd; 
    }

    std::int32_t Float32Exp2Int32::exponent() { 
        this->scale();
        return exp;
    }

    float Float32Exp2Int32::asFloat() const
    { 
        return scnfcnd * std::pow(2, exp);
    }

    Float32Exp2Int32& Float32Exp2Int32::operator+=(Float32Exp2Int32 z) {
        if (z.scnfcnd == 0) {
            return *this;
        }
        if (scnfcnd == 0) {
            scnfcnd = z.scnfcnd;
            exp = z.exp;
            return *this;
        }

        std::int32_t exp_diff = exp - z.exp;
        std::uint32_t* sgnfcndBitsZ = reinterpret_cast<std::uint32_t*>(&z.scnfcnd);

        if (exp_diff >= 0) {
            if (exp_diff > 31) {
                return *this;
            }
            *sgnfcndBitsZ -= exp_diff << 23;
        }
        else {
            if (exp_diff < -31) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            *sgnfcndBitsZ += exp_diff << 23;
        }
        
        scnfcnd += z.scnfcnd;
        this->checkLimitForScale();    
        return *this;
    }

    Float32Exp2Int32& Float32Exp2Int32::operator*=(Float32Exp2Int32 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkLimitForScale();
        return *this;
    }

    Float32Exp2Int32& Float32Exp2Int32::operator*=(float& flt) {
        Float32Exp2Int32 z(flt);
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkLimitForScale();
        return *this;
    }

    inline void Float32Exp2Int32::checkLimitForScale() {
        if (16 <= scnfcnd || (-0.0625 <= scnfcnd && scnfcnd <= 0.0625) || -16 >= scnfcnd) {
            this->scale();
        }
    }

    inline void Float32Exp2Int32::scale() {
        if (std::numeric_limits<float>::min() > scnfcnd && scnfcnd > -std::numeric_limits<float>::min()) {
            scnfcnd = 0.0;
            exp = 0;
            return;
        }

        int32_t sgnfcndBits = *reinterpret_cast<int32_t*>(&scnfcnd);
        exp += ((sgnfcndBits & 0x7F800000) >> 23) - 127;
        sgnfcndBits &= 0x807FFFFF;
        sgnfcndBits |= 0x3F800000;
        scnfcnd = *reinterpret_cast<float*>(&sgnfcndBits);
    }

    Float32Exp2Int32 operator+(Float32Exp2Int32 a, const Float32Exp2Int32 b) { return a+=b; }
    Float32Exp2Int32 operator*(Float32Exp2Int32 a, const Float32Exp2Int32 b) { return a*=b; }
}

