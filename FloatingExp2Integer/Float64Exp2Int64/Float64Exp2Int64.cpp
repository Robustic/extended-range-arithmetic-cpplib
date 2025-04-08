#include <cstdint>
#include <cmath>
#include "Float64Exp2Int64.h"

namespace floatingExp2Integer
{
    Float64Exp2Int64::Float64Exp2Int64() {
        scnfcnd = 1.0;
        exp = 0;
        this->scaleIfNotZero();
    }

    Float64Exp2Int64::Float64Exp2Int64(double dbl) {
        scnfcnd = dbl;
        exp = 0;
        this->scaleIfNotZero();
    }

    void Float64Exp2Int64::doubleToFloat64Exp2Int64(double dbl) {
        scnfcnd = dbl;
        exp = 0;
        this->scaleIfNotZero();
    }

    void Float64Exp2Int64::log2ToFloat64Exp2Int64(double log2) {
        exp = (int64_t)log2;
        scnfcnd = std::exp2(log2 - exp);
        this->scaleIfNotZero();
    }

    double Float64Exp2Int64::float64Exp2Int64ToLog2() const {
        return std::log2(scnfcnd) + exp;
    }

    double Float64Exp2Int64::sicnificand() {
        this->scaleIfNotZero();
        return scnfcnd; 
    }

    int64_t Float64Exp2Int64::exponent() {
        this->scaleIfNotZero();
        return exp;
    }

    double Float64Exp2Int64::asDouble() const {
        if (scnfcnd == 0.0) {
            return 0.0;
        }
        return scnfcnd * std::pow(2.0, exp);
    }

    Float64Exp2Int64& Float64Exp2Int64::operator+=(Float64Exp2Int64 z) {
        int64_t exp_diff = exp - z.exp;
        uint64_t* sgnfcnd_bits_z = reinterpret_cast<uint64_t*>(&z.scnfcnd);

        if (exp_diff >= 0) {
            if (exp_diff > 64) {
                return *this;
            }
            *sgnfcnd_bits_z -= exp_diff << 52;
            scnfcnd += z.scnfcnd;
        }
        else {
            if (exp_diff < -64) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            *sgnfcnd_bits_z += (-exp_diff) << 52;
            scnfcnd += z.scnfcnd;
        }

        this->checkRuleForScale();    
        return *this;
    }

    Float64Exp2Int64& Float64Exp2Int64::operator*=(Float64Exp2Int64 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    Float64Exp2Int64& Float64Exp2Int64::operator*=(double& dbl) {
        Float64Exp2Int64 z(dbl);
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    inline void Float64Exp2Int64::checkRuleForScale() {
        if (0x1p-6 < scnfcnd) {
            if (0x1p6 <= scnfcnd) {
                this->scale();
            }
        }
        else if (scnfcnd < -0x1p-6) {
            if (scnfcnd <= -0x1p6) {
                this->scale();
            }
        }
        else {
            this->scaleIfNotZero();
        }
    }

    inline void Float64Exp2Int64::scaleIfNotZero() {
        if (scnfcnd <= -0x1p-1022 || 0x1p-1022 <= scnfcnd) {
            this->scale();
        }
        else {
            scnfcnd = 0.0;
            int64_t exp_for_zero = 0x8000000000000000ll;
            exp = exp_for_zero;
        }
    }

    inline void Float64Exp2Int64::scale() {
        uint64_t* sgnfcnd_bits = reinterpret_cast<uint64_t*>(&scnfcnd);
        exp += ((*sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1023;
        *sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
        *sgnfcnd_bits |= 0x3FF0000000000000ull;
        scnfcnd = *reinterpret_cast<double*>(sgnfcnd_bits);
    }

    Float64Exp2Int64 operator+(Float64Exp2Int64 a, const Float64Exp2Int64 b) { return a+=b; }
    Float64Exp2Int64 operator*(Float64Exp2Int64 a, const Float64Exp2Int64 b) { return a*=b; }
}

