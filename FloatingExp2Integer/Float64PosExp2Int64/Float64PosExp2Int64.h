#ifndef FLOAT64POS_EXP2INT64_H_
#define FLOAT64POS_EXP2INT64_H_

namespace floatingExp2Integer
{
    class Float64PosExp2Int64 {
        private:
            double scnfcnd;
            std::int64_t exp;
            inline void checkRuleForScale();
            inline void scale();
        public:
            Float64PosExp2Int64();
            Float64PosExp2Int64(double dbl);
            void doubleToFloat64PosExp2Int64(double dbl);
            void log2ToFloat64PosExp2Int64(double log2);
            double float64PosExp2Int64ToLog2() const;
            double sicnificand();
            std::int64_t exponent();
            double asDouble() const;

            Float64PosExp2Int64& operator+=(Float64PosExp2Int64 z);
            Float64PosExp2Int64& operator*=(Float64PosExp2Int64 z);
            Float64PosExp2Int64& operator*=(double& dbl);
    };

    Float64PosExp2Int64 operator+(Float64PosExp2Int64 a, const Float64PosExp2Int64 b);
    Float64PosExp2Int64 operator*(Float64PosExp2Int64 a, const Float64PosExp2Int64 b);
}

#endif
