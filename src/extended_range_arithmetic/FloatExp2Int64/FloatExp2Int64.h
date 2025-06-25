#ifndef FLOAT_EXP2_INT_64_H_
#define FLOAT_EXP2_INT_64_H_

#include <vector>

namespace extended_range_arithmetic
{
    class FloatExp2Int64 {
        private:
            double principal;
            int64_t aux;

            FloatExp2Int64(double principal_part, int64_t auxiliary_index);
            inline void checkRuleForScale();
            inline void scale();
        public:
            FloatExp2Int64();
            FloatExp2Int64(double dbl);

            void double_to(double dbl);
            void log2_to(double log2);
            static void log2s_to(const std::vector<double>& from, std::vector<extended_range_arithmetic::FloatExp2Int64>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            void sum(const std::vector<extended_range_arithmetic::FloatExp2Int64>& vector);
            void multiply(const std::vector<extended_range_arithmetic::FloatExp2Int64>& vector);

            double as_double() const;
            double as_log2() const;

            double principal_part() const { return principal; };
            int64_t auxiliary_index() const { return aux; };

            FloatExp2Int64& operator+=(FloatExp2Int64 z);
            FloatExp2Int64& operator*=(FloatExp2Int64 z);
            FloatExp2Int64& operator*=(double& dbl);
    };

    FloatExp2Int64 operator+(FloatExp2Int64 a, const FloatExp2Int64 b);
    FloatExp2Int64 operator*(FloatExp2Int64 a, const FloatExp2Int64 b);
}

#endif
