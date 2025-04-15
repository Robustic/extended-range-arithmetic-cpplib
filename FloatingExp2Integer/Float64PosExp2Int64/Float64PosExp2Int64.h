#ifndef FLOAT64POS_EXP2INT64_H_
#define FLOAT64POS_EXP2INT64_H_

#include <vector>

namespace floatingExp2Integer
{
    class Float64PosExp2Int64 {
        private:
            double scnfcnd;
            int64_t exp;

            Float64PosExp2Int64(double sicnificand, int64_t exponent);
            inline void checkRuleForScale();
            inline void scale();
        public:
            Float64PosExp2Int64();
            Float64PosExp2Int64(double dbl);

            void double_to(double dbl);
            void log2_to(double log2);
            static void log2s_to(const std::vector<double>& from, std::vector<floatingExp2Integer::Float64PosExp2Int64>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            void sum(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector);
            void multiply(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector);

            double as_double() const;
            double as_log2() const;

            double sicnificand() const { return scnfcnd; };
            int64_t exponent() const { return exp; };

            Float64PosExp2Int64& operator+=(Float64PosExp2Int64 z);
            Float64PosExp2Int64& operator*=(Float64PosExp2Int64 z);
            Float64PosExp2Int64& operator*=(double& dbl);
    };

    Float64PosExp2Int64 operator+(Float64PosExp2Int64 a, const Float64PosExp2Int64 b);
    Float64PosExp2Int64 operator*(Float64PosExp2Int64 a, const Float64PosExp2Int64 b);
}

#endif
