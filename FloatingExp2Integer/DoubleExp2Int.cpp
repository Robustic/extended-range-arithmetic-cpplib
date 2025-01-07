#include <iostream>
#include <cstring>
#include <cmath>
#include "DoubleExp2Int.h"

namespace floatingExp2Integer
{
    DoubleExp2Int::DoubleExp2Int(double dbl) {
        if (dbl == 0) {
            mant = 0.0;
            exp = 0;
            return;
        }

        std::int32_t parts[2];
        std::memcpy(&parts, &dbl, sizeof(std::int32_t) * 2);
        exp = ((parts[1] & 0x7FF00000) >> 20) - 1023;
        parts[1] &= 0x800FFFFF;
        parts[1] |= 0x3FF00000;
        std::memcpy(&mant, &parts, sizeof(double));
    }

    DoubleExp2Int& DoubleExp2Int::operator+=(DoubleExp2Int z) {
        std::int32_t parts[2];
        std::memcpy(&parts, &z.mant, sizeof(std::int32_t) * 2);
        std::int32_t zMantOldExp = parts[1] & 0x7FF00000;
        parts[1] &= 0x800FFFFF;
        parts[1] |= zMantOldExp + ((z.exp - this->exp) << 20);
        std::memcpy(&z.mant, &parts, sizeof(double));
        mant += z.mant;
        this->scale();
        return *this;
    }

    void DoubleExp2Int::scale() {
        if (mant == 0) {
            mant = 0.0;
            exp = 0;
            return;
        }

        std::int32_t parts[2];
        std::memcpy(&parts, &mant, sizeof(std::int32_t) * 2);
        exp += ((parts[1] & 0x7FF00000) >> 20) - 1023;
        parts[1] &= 0x800FFFFF;
        parts[1] |= 0x3FF00000;
        std::memcpy(&mant, &parts, sizeof(double));
    }

    DoubleExp2Int operator+(DoubleExp2Int a, const DoubleExp2Int b) { return a+=b; }

    double DoubleExp2Int::asDouble()
    { 
        return mant * std::pow(2, exp);
    }
}
