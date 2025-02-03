#ifndef DOUBLE_EXP2_INT_H_
#define DOUBLE_EXP2_INT_H_

namespace floatingExp2Integer
{
    class DoubleExp2Int {
        private:
            double mant;
            std::int32_t exp;
        public:
            DoubleExp2Int(double number);
            double mantissa() const { return mant; }
            int exponent() const { return exp; }
            double asDouble();
            inline void scale();

            DoubleExp2Int& operator+=(DoubleExp2Int z);
            DoubleExp2Int& operator*=(DoubleExp2Int z);
            DoubleExp2Int& operator*=(double& z);
    };

    DoubleExp2Int operator+(DoubleExp2Int a, DoubleExp2Int b);
    DoubleExp2Int operator*(DoubleExp2Int a, DoubleExp2Int b);
}

#endif
