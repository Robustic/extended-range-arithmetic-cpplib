#include <cstdint>
#include <cmath>
#include "Float32PosExp2Int32.h"

namespace floatingExp2Integer
{
    Float32PosExp2Int32::Float32PosExp2Int32() {
        scnfcnd = 1.0;
        exp = 0;
        this->scale();
    }

    Float32PosExp2Int32::Float32PosExp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->scale();
    }

    void Float32PosExp2Int32::floatToFloat32PosExp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->scale();
    }

    void Float32PosExp2Int32::log2ToFloat32PosExp2Int32(float logarithm2) {
        exp = (std::int32_t)logarithm2;
        scnfcnd = std::exp2(logarithm2 - exp);
        this->scale();
    }

    float Float32PosExp2Int32::float32PosExp2Int32ToLog2() {
        return std::log2(scnfcnd) + exp ;
    }


    float Float32PosExp2Int32::sicnificand() { 
        this->scale();
        return scnfcnd; 
    }

    std::int32_t Float32PosExp2Int32::exponent() { 
        this->scale();
        return exp;
    }

    float Float32PosExp2Int32::asFloat() const
    { 
        return scnfcnd * std::pow(2, exp);
    }

    Float32PosExp2Int32& Float32PosExp2Int32::operator+=(Float32PosExp2Int32 z) {
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
            *sgnfcndBitsZ += (-exp_diff) << 23;
        }
        
        scnfcnd += z.scnfcnd;
        this->checkLimitForScale();    
        return *this;
    }

    Float32PosExp2Int32& Float32PosExp2Int32::operator*=(Float32PosExp2Int32 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkLimitForScale();
        return *this;
    }

    Float32PosExp2Int32& Float32PosExp2Int32::operator*=(float& flt) {
        Float32PosExp2Int32 z(flt);
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkLimitForScale();
        return *this;
    }

    inline void Float32PosExp2Int32::checkLimitForScale() {
        if (16 <= scnfcnd || scnfcnd <= 0.0625) {
            this->scale();
        }
    }

    inline void Float32PosExp2Int32::scale() {
        int32_t* sgnfcndBits = reinterpret_cast<int32_t*>(&scnfcnd);
        exp += ((*sgnfcndBits & 0x7F800000) >> 23) - 127;
        *sgnfcndBits &= 0x807FFFFF;
        *sgnfcndBits |= 0x3F800000;
        scnfcnd = *reinterpret_cast<float*>(sgnfcndBits);
    }

    Float32PosExp2Int32 operator+(Float32PosExp2Int32 a, const Float32PosExp2Int32 b) { return a+=b; }
    Float32PosExp2Int32 operator*(Float32PosExp2Int32 a, const Float32PosExp2Int32 b) { return a*=b; }
}

