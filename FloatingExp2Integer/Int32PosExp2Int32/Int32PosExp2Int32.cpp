#include <cstdint>
#include <cmath>
#include "Int32PosExp2Int32.h"

namespace floatingExp2Integer
{
    Int32PosExp2Int32::Int32PosExp2Int32() {
        this->fromFloat(1.0f);
        this->scale();
    }

    Int32PosExp2Int32::Int32PosExp2Int32(float flt) {
        this->fromFloat(flt);
        this->scale();
    }

    void Int32PosExp2Int32::floatToInt32PosExp2Int32(float flt) {
        this->fromFloat(flt);
        this->scale();
    }

    void Int32PosExp2Int32::log2ToInt32PosExp2Int32(float log2) {
        std::int32_t exponent = (std::int32_t)log2;
        float flt = std::exp2f(log2 - exponent);
        this->fromFloat(flt);
        exp += exponent;
        this->scale();
    }

    float Int32PosExp2Int32::int32PosExp2Int32ToLog2() const {
        return std::log2(scnfcnd) + exp;
    }

    std::uint64_t Int32PosExp2Int32::sicnificand() { 
        this->scale();
        return scnfcnd; 
    }

    std::int32_t Int32PosExp2Int32::exponent() { 
        this->scale();
        return exp;
    }

    float Int32PosExp2Int32::asFloat() const { 
        return (float)scnfcnd * std::pow(2.0f, exp);
    }

    Int32PosExp2Int32& Int32PosExp2Int32::operator+=(Int32PosExp2Int32 z) {
        std::int64_t exp_diff = exp - z.exp;

        if (exp_diff >= 0) {
            if (exp_diff > 63) {
                return *this;
            }
            scnfcnd += z.scnfcnd >> exp_diff;
        }
        else {
            if (exp_diff < -63) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            exp -= exp_diff;
            scnfcnd = scnfcnd >> (-exp_diff);
            scnfcnd += z.scnfcnd;
        }        
        
        this->checkRuleForScale();    
        return *this;
    }

    Int32PosExp2Int32& Int32PosExp2Int32::operator*=(Int32PosExp2Int32 z) {
        std::uint32_t offset = 32 - __builtin_clzll(scnfcnd);
        std::uint32_t offset_z = 32 - __builtin_clzll(z.scnfcnd);
        scnfcnd = (scnfcnd >> offset) * (z.scnfcnd >> offset_z);
        exp += z.exp + offset + offset_z;
        this->checkRuleForScale();
        return *this;
    }

    inline void Int32PosExp2Int32::fromFloat(float flt) {
        std::uint32_t* sgnfcnd_bits = reinterpret_cast<std::uint32_t*>(&flt);
        scnfcnd = ((std::uint64_t)(*sgnfcnd_bits & 0x007FFFFFu)) << 29;
        scnfcnd |= 0x0010000000000000ull;
        exp = ((*sgnfcnd_bits & 0x7F800000) >> 23) - (127 + 52);
    }

    inline void Int32PosExp2Int32::checkRuleForScale() {
        if (0x1000000000000000ull <= scnfcnd) {
            this->scale();
        }
    }

    inline void Int32PosExp2Int32::scale() {
        std::uint32_t offset = 11 - __builtin_clzll(scnfcnd);
        exp += offset;
        scnfcnd = scnfcnd >> offset;
    }

    Int32PosExp2Int32 operator+(Int32PosExp2Int32 a, const Int32PosExp2Int32 b) { return a+=b; }
    Int32PosExp2Int32 operator*(Int32PosExp2Int32 a, const Int32PosExp2Int32 b) { return a*=b; }
}
