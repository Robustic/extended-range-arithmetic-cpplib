#include <cstdint>
#include <cmath>
#include <cstring>
// #include <limits>
#include <iostream>
#include <iomanip>
#include "Int32PosExp2Int32.h"

namespace floatingExp2Integer
{
    Int32PosExp2Int32::Int32PosExp2Int32() {
        scnfcnd = 1.0;
        exp = 0;
        this->scale();
    }

    Int32PosExp2Int32::Int32PosExp2Int32(float flt) {
        std::uint32_t* fltBits = reinterpret_cast<std::uint32_t*>(&flt);
        scnfcnd = *fltBits & 0x007FFFFF;
        scnfcnd |= 0x00800000;
        exp = ((*fltBits & 0x7F800000) >> 23) - 127 - 23;
        this->scale();
    }

    Int32PosExp2Int32::Int32PosExp2Int32(std::uint32_t sicnificand, std::int32_t exponent) {
        scnfcnd = sicnificand;
        exp = exponent;
        this->scale();
    }

    void Int32PosExp2Int32::log2ToInt32Exp2Int32(float logarithm2) {
        exp = (std::int32_t)logarithm2;
        float flt = std::exp2(logarithm2 - exp);
        std::uint32_t* fltBits = reinterpret_cast<std::uint32_t*>(&flt);
        scnfcnd = *fltBits & 0x807FFFFF;
        scnfcnd |= 0x00800000;
        exp += ((*fltBits & 0x7F800000) >> 23) - 127 - 23;
        this->scale();
    }

    float Int32PosExp2Int32::int32Exp2Int32ToLog2() {
        return std::log2(scnfcnd) + exp ;
    }


    std::int32_t Int32PosExp2Int32::sicnificand() { 
        this->scale();
        return scnfcnd; 
    }

    std::int32_t Int32PosExp2Int32::exponent() { 
        this->scale();
        return exp;
    }

    float Int32PosExp2Int32::asFloat() const
    { 
        return (float)scnfcnd * std::pow(2, exp);
    }

    Int32PosExp2Int32& Int32PosExp2Int32::operator+=(Int32PosExp2Int32 z) {
        std::int32_t exp_diff = exp - z.exp;

        if (exp_diff >= 0) {
            if (exp_diff > 31) {
                return *this;
            }
            scnfcnd += z.scnfcnd >> exp_diff;
        }
        else {
            if (exp_diff < -31) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            exp -= exp_diff;
            scnfcnd = scnfcnd >> (-exp_diff);
            scnfcnd += z.scnfcnd;
        }
        
        
        this->checkLimitForScale();    
        return *this;
    }

    Int32PosExp2Int32& Int32PosExp2Int32::operator*=(Int32PosExp2Int32 z) {
        std::int32_t offset = 16 - __builtin_clz(scnfcnd);
        std::int32_t offsetZ = 16 - __builtin_clz(z.scnfcnd);
        scnfcnd = (scnfcnd >> offset) * (z.scnfcnd >> offsetZ);
        exp += z.exp + offset + offsetZ;
        this->checkLimitForScale();
        return *this;
    }

    inline void Int32PosExp2Int32::checkLimitForScale() {
        if (2147483648 <= scnfcnd) {
            this->scale();
        }
    }

    inline void Int32PosExp2Int32::scale() {
        std::int32_t offset = 8 - __builtin_clz(scnfcnd);
        exp += offset;
        scnfcnd = scnfcnd >> offset;
    }

    Int32PosExp2Int32 operator+(Int32PosExp2Int32 a, const Int32PosExp2Int32 b) { return a+=b; }
    Int32PosExp2Int32 operator*(Int32PosExp2Int32 a, const Int32PosExp2Int32 b) { return a*=b; }
}
