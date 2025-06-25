#ifndef XNUMBER_64_H_
#define XNUMBER_64_H_

#include <vector>

namespace extended_range_arithmetic
{
    class Xnumber64 {
        private:
            double principal;
            int64_t aux;

            Xnumber64(double principal_part, int64_t auxiliary_index);
            inline void xnorm();

            static constexpr int64_t EXP_MULTIPLIER = 960;
            static constexpr double BIG = 0x1p960;
            static constexpr double BIGI = 0x1p-960;
            static constexpr double BIGS = 0x1p480;
            static constexpr double BIGSI = 0x1p-480;
        public:
            Xnumber64();
            Xnumber64(double dbl);

            void log2_to(const double from);
            static void log2s_to(const std::vector<double>& from, std::vector<extended_range_arithmetic::Xnumber64>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            void sum(const std::vector<extended_range_arithmetic::Xnumber64>& vector);
            void multiply(const std::vector<extended_range_arithmetic::Xnumber64>& vector);

            double as_double() const;
            double as_log2() const;

            double principal_part() const { return principal; };
            int64_t auxiliary_index() const { return aux; };

            Xnumber64& operator+=(Xnumber64 z);
            Xnumber64& operator*=(Xnumber64 z);
    };

    // Float64ExtendedExp operator+(Float64ExtendedExp a, const Float64ExtendedExp b);
    // Float64ExtendedExp operator*(Float64ExtendedExp a, const Float64ExtendedExp b);
}

#endif
