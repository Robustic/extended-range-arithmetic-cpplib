#ifndef INT_EXP2_INT_64_H_
#define INT_EXP2_INT_64_H_

#include <vector>

namespace extended_range_arithmetic
{
    class IntExp2Int64 {
        private:
            uint64_t principal;
            int64_t aux;

            IntExp2Int64(uint64_t principal_part, int64_t auxiliary_index);
            inline void fromDouble(double dbl);
            inline void checkRuleForScale();
            inline void scale();
        public:
            IntExp2Int64();
            IntExp2Int64(double dbl);

            void double_to(double dbl);
            void log2_to(double log2);
            static void log2s_to(const std::vector<double>& from, std::vector<extended_range_arithmetic::IntExp2Int64>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            void sum(const std::vector<extended_range_arithmetic::IntExp2Int64>& vector);
            void multiply(const std::vector<extended_range_arithmetic::IntExp2Int64>& vector);

            double as_log2() const;
            double as_double() const;

            uint64_t principal_part() const { return principal; };
            int64_t auxiliary_index() const { return aux; };

            IntExp2Int64& operator+=(IntExp2Int64 z);
            IntExp2Int64& operator*=(IntExp2Int64 z);
    };

    IntExp2Int64 operator+(IntExp2Int64 a, const IntExp2Int64 b);
    IntExp2Int64 operator*(IntExp2Int64 a, const IntExp2Int64 b);
}

#endif
