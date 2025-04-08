#ifndef FLOAT64POS_EXP2INT64_H_
#define FLOAT64POS_EXP2INT64_H_

#include <vector>

namespace floatingExp2Integer
{
    class Float64PosExp2Int64 {
        private:
            Float64PosExp2Int64(double sicnificand, std::int64_t exponent);
            inline void checkRuleForScale();
            inline void scale();
        public:
            double scnfcnd;
            std::int64_t exp;
            Float64PosExp2Int64();
            Float64PosExp2Int64(double dbl);
            void sum(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector);
            void multiply(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector);
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
