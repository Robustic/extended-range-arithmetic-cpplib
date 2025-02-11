#ifndef DOUBLE_EXP2_INT_H_
#define DOUBLE_EXP2_INT_H_

namespace floatingExp2Integer
{
    class Float64PosExp2Int64 {
        private:
            double mant;
            std::int64_t exp;
        public:
            Float64PosExp2Int64();
            Float64PosExp2Int64(double number);
            Float64PosExp2Int64(double dbl, std::int64_t ex);
            void Log2ToFloat64Exp2Int64(double logarithm2);
            double Float64Exp2Int64ToLog2();
            double mantissa() const { return mant; }
            int exponent() const { return exp; }
            double asDouble();
            inline void scale();

            Float64PosExp2Int64& operator+=(Float64PosExp2Int64 z);
            Float64PosExp2Int64& operator*=(Float64PosExp2Int64 z);
            Float64PosExp2Int64& operator*=(double& z);
    };

    Float64PosExp2Int64 operator+(Float64PosExp2Int64 a, Float64PosExp2Int64 b);
    Float64PosExp2Int64 operator*(Float64PosExp2Int64 a, Float64PosExp2Int64 b);
}

#endif
