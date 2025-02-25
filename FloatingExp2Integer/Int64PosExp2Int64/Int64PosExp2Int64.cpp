#include <cstdint>
#include <cmath>
#include "Int64PosExp2Int64.h"

namespace floatingExp2Integer
{
    Int64PosExp2Int64::Int64PosExp2Int64() {
        this->fromDouble(1.0);
        this->scale();
    }

    Int64PosExp2Int64::Int64PosExp2Int64(double dbl) {
        this->fromDouble(dbl);
        this->scale();
    }

    void Int64PosExp2Int64::doubleToInt64PosExp2Int64(double dbl) {
        this->fromDouble(dbl);
        this->scale();
    }

    void Int64PosExp2Int64::log2ToInt64PosExp2Int64(double log2) {
        std::int64_t exponent = (std::int64_t)log2;
        double dbl = std::exp2(log2 - exponent);
        this->fromDouble(dbl);
        exp += exponent;
        this->scale();
    }

    double Int64PosExp2Int64::int64PosExp2Int64ToLog2() const {
        return std::log2(scnfcnd) + exp;
    }

    std::uint64_t Int64PosExp2Int64::sicnificand() { 
        this->scale();
        return scnfcnd; 
    }

    std::int64_t Int64PosExp2Int64::exponent() { 
        this->scale();
        return exp;
    }

    double Int64PosExp2Int64::asDouble() const { 
        return (double)scnfcnd * std::pow(2.0, exp);
    }

    Int64PosExp2Int64& Int64PosExp2Int64::operator+=(Int64PosExp2Int64 z) {
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

    Int64PosExp2Int64& Int64PosExp2Int64::operator*=(Int64PosExp2Int64 z) {
        std::uint32_t offset = 32 - __builtin_clzll(scnfcnd);
        std::uint32_t offset_z = 32 - __builtin_clzll(z.scnfcnd);
        scnfcnd = (scnfcnd >> offset) * (z.scnfcnd >> offset_z);
        exp += z.exp + offset + offset_z;
        this->checkRuleForScale();
        return *this;
    }

    inline void Int64PosExp2Int64::fromDouble(double dbl) {
        std::uint64_t* sgnfcnd_bits = reinterpret_cast<std::uint64_t*>(&dbl);
        scnfcnd = *sgnfcnd_bits & 0x000FFFFFFFFFFFFFull;
        scnfcnd |= 0x0010000000000000ull;
        exp = ((*sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - (1023 + 52);
    }

    inline void Int64PosExp2Int64::checkRuleForScale() {
        if (0x1000000000000000ull <= scnfcnd) {
            this->scale();
        }
    }

    inline void Int64PosExp2Int64::scale() {
        std::uint32_t offset = 11 - __builtin_clzll(scnfcnd);
        exp += offset;
        scnfcnd = scnfcnd >> offset;
    }

    Int64PosExp2Int64 operator+(Int64PosExp2Int64 a, const Int64PosExp2Int64 b) { return a+=b; }
    Int64PosExp2Int64 operator*(Int64PosExp2Int64 a, const Int64PosExp2Int64 b) { return a*=b; }
}
