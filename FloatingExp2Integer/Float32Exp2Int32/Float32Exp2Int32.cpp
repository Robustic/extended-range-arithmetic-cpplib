#include <cstdint>
#include <cmath>
#include "Float32Exp2Int32.h"

namespace floatingExp2Integer
{
    Float32Exp2Int32::Float32Exp2Int32() {
        scnfcnd = 1.0f;
        exp = 0;
        this->scaleIfNotZero();
    }

    Float32Exp2Int32::Float32Exp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->scaleIfNotZero();
    }

    void Float32Exp2Int32::floatToFloat32Exp2Int32(float flt) {
        scnfcnd = flt;
        exp = 0;
        this->scaleIfNotZero();
    }

    void Float32Exp2Int32::log2ToFloat32Exp2Int32(float log2) {
        exp = (std::int32_t)log2;
        scnfcnd = std::exp2f(log2 - exp);
        this->scaleIfNotZero();
    }

    float Float32Exp2Int32::float32Exp2Int32ToLog2() const {
        return std::log2f(scnfcnd) + exp;
    }

    float Float32Exp2Int32::sicnificand() {
        this->scaleIfNotZero();
        return scnfcnd; 
    }

    std::int32_t Float32Exp2Int32::exponent() {
        this->scaleIfNotZero();
        return exp;
    }

    float Float32Exp2Int32::asFloat() const {
        if (scnfcnd == 0.0f) {
            return 0.0f;
        }
        return scnfcnd * std::pow(2.0f, exp);
    }

    Float32Exp2Int32& Float32Exp2Int32::operator+=(Float32Exp2Int32 z) {
        std::int32_t exp_diff = exp - z.exp;
        std::uint32_t* sgnfcnd_bits_z = reinterpret_cast<std::uint32_t*>(&z.scnfcnd);

        if (exp_diff >= 0) {
            if (exp_diff > 35) {
                return *this;
            }
            *sgnfcnd_bits_z -= exp_diff << 23;
        }
        else {
            if (exp_diff < -35) {
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

    Float32Exp2Int32& Float32Exp2Int32::operator*=(Float32Exp2Int32 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    Float32Exp2Int32& Float32Exp2Int32::operator*=(float& flt) {
        Float32Exp2Int32 z(flt);
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    inline void Float32Exp2Int32::checkRuleForScale() {
        if (0x1p-6f < scnfcnd) {
            if (0x1p6f <= scnfcnd) {
                this->scale();
            }
        }
        else if (scnfcnd < -0x1p-6f) {
            if (scnfcnd <= -0x1p6f) {
                this->scale();
            }
        }
        else {
            this->scaleIfNotZero();
        }
    }

    inline void Float32Exp2Int32::scaleIfNotZero() {
        if (scnfcnd <= -0x1p-126f || 0x1p-126f <= scnfcnd) {
            this->scale();
        }
        else {
            scnfcnd = 0.0f;
            std::int32_t exp_for_zero = 0x80000000;
            exp = exp_for_zero;
        }
    }

    inline void Float32Exp2Int32::scale() {
        uint32_t* sgnfcnd_bits = reinterpret_cast<uint32_t*>(&scnfcnd);
        exp += ((*sgnfcnd_bits & 0x7F800000u) >> 23) - 127;
        *sgnfcnd_bits &= 0x807FFFFFu;
        *sgnfcnd_bits |= 0x3F800000u;
        scnfcnd = *reinterpret_cast<float*>(sgnfcnd_bits);
    }

    Float32Exp2Int32 operator+(Float32Exp2Int32 a, const Float32Exp2Int32 b) { return a+=b; }
    Float32Exp2Int32 operator*(Float32Exp2Int32 a, const Float32Exp2Int32 b) { return a*=b; }
}

