#include <iostream>
#include <cstring>
#include <cmath>
#include <limits>
#include "Float64PosExp2Int64.h"

namespace floatingExp2Integer
{
    Float64PosExp2Int64::Float64PosExp2Int64(double dbl) {
        mant = dbl;
        exp = 0;
        this->scale();
    }

    Float64PosExp2Int64::Float64PosExp2Int64(double dbl, std::int64_t ex) {
        mant = dbl;
        exp = ex;
        this->scale();
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator+=(Float64PosExp2Int64 z) {
        if (z.mant == 0) {
            return *this;
        }
        if (mant == 0) {
            mant = z.mant;
            exp = z.exp;
            return *this;
        }

        std::int64_t exp_diff = exp - z.exp;
        if (exp_diff > 53) {
            return *this;
        }
        if (exp_diff < -53) {
            mant = z.mant;
            exp = z.exp;
            return *this;
        }

        std::int64_t parts[1];
        std::memcpy(&parts, &z.mant, sizeof(std::int64_t) * 1);
        std::int64_t zMantOldExp = parts[0] & 0x7FF0000000000000;
        parts[0] &= 0x800FFFFFFFFFFFFF;
        parts[0] |= zMantOldExp + ((z.exp - this->exp) << 52);
        std::memcpy(&z.mant, &parts, sizeof(double));
        mant += z.mant;
        this->scale();
        return *this;
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator*=(Float64PosExp2Int64 z) {
        mant *= z.mant;
        exp += z.exp;
        this->scale();
        return *this;
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator*=(double& z) {
        mant *= z;
        this->scale();
        return *this;
    }

    inline void Float64PosExp2Int64::scale() {
        if (std::numeric_limits<double>::min() > mant && mant > -std::numeric_limits<double>::min()) {
            mant = 0.0;
            exp = 0;
            return;
        }
        std::int64_t parts[1];
        std::memcpy(&parts, &mant, sizeof(std::int64_t) * 1);
        exp += ((parts[0] & 0x7FF0000000000000) >> 52) - 1023;
        parts[0] &= 0x800FFFFFFFFFFFFF;
        parts[0] |= 0x3FF0000000000000;
        std::memcpy(&mant, &parts, sizeof(double));
    }

    Float64PosExp2Int64 operator+(Float64PosExp2Int64 a, const Float64PosExp2Int64 b) { return a+=b; }
    Float64PosExp2Int64 operator*(Float64PosExp2Int64 a, const Float64PosExp2Int64 b) { return a*=b; }

    double Float64PosExp2Int64::asDouble()
    { 
        return mant * std::pow(2, exp);
    }
}
