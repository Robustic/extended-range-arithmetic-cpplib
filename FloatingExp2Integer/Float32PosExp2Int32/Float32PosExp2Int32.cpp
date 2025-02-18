#include <cstdint>
#include <cmath>
#include "Float32PosExp2Int32.h"

namespace floatingExp2Integer
{
    Float32PosExp2Int32::Float32PosExp2Int32() {
        scnfcnd = 1.0f;
        exp = 0;
        this->checkRuleForScale();
    }

    Float32PosExp2Int32::Float32PosExp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->checkRuleForScale();
    }

    void Float32PosExp2Int32::floatToFloat32PosExp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->checkRuleForScale();
    }

    void Float32PosExp2Int32::log2ToFloat32PosExp2Int32(float log2) {
        exp = (std::int32_t)log2;
        scnfcnd = std::exp2f(log2 - exp);
        this->checkRuleForScale();
    }

    float Float32PosExp2Int32::float32PosExp2Int32ToLog2() const {
        return std::log2f(scnfcnd) + exp;
    }

    float Float32PosExp2Int32::sicnificand() {
        this->scale();
        return scnfcnd; 
    }

    std::int32_t Float32PosExp2Int32::exponent() {
        this->scale();
        return exp;
    }

    float Float32PosExp2Int32::asFloat() const {
        return scnfcnd * std::pow(2.0f, exp);
    }

    Float32PosExp2Int32& Float32PosExp2Int32::operator+=(Float32PosExp2Int32 z) {
        std::int32_t exp_diff = exp - z.exp;
        std::uint32_t* sgnfcnd_bits_z = reinterpret_cast<std::uint32_t*>(&z.scnfcnd);

        if (exp_diff >= 0) {
            if (exp_diff > 31) {
                return *this;
            }
            *sgnfcnd_bits_z -= exp_diff << 23;
        }
        else {
            if (exp_diff < -31) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            *sgnfcnd_bits_z += (-exp_diff) << 23;
        }

        scnfcnd += z.scnfcnd;
        this->checkRuleForScale();    
        return *this;
    }

    Float32PosExp2Int32& Float32PosExp2Int32::operator*=(Float32PosExp2Int32 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    Float32PosExp2Int32& Float32PosExp2Int32::operator*=(float& flt) {
        Float32PosExp2Int32 z(flt);
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    inline void Float32PosExp2Int32::checkRuleForScale() {
        if (16 <= scnfcnd || scnfcnd <= 0.0625) {
            this->scale();
        }
    }

    inline void Float32PosExp2Int32::scale() {
        uint32_t* sgnfcnd_bits = reinterpret_cast<uint32_t*>(&scnfcnd);
        exp += ((*sgnfcnd_bits & 0x7F800000) >> 23) - 127;
        *sgnfcnd_bits &= 0x807FFFFF;
        *sgnfcnd_bits |= 0x3F800000;
        scnfcnd = *reinterpret_cast<float*>(sgnfcnd_bits);
    }

    Float32PosExp2Int32 operator+(Float32PosExp2Int32 a, const Float32PosExp2Int32 b) { return a+=b; }
    Float32PosExp2Int32 operator*(Float32PosExp2Int32 a, const Float32PosExp2Int32 b) { return a*=b; }
}

