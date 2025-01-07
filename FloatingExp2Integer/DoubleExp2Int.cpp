#include <iostream>
#include <cstring>
#include <cmath>
#include "DoubleExp2Int.h"

namespace floatingExp2Integer
{
    DoubleExp2Int::DoubleExp2Int(double number) {
        if (number == 0) {
            mant = 0.0;
            exp = 0;
            return;
        }

        std::int32_t value[2];
        std::memcpy(&value, &number, sizeof(std::int32_t) * 2);
        exp = ((value[1] & 0x7FF00000) >> 20) - 1023;
        value[1] &= 0x800FFFFF;
        value[1] |= 0x3FF00000;
        std::memcpy(&mant, &value, sizeof(double));
    }

    double DoubleExp2Int::asDouble()
    { 
        return mant * std::pow(2, exp);
    }
}
