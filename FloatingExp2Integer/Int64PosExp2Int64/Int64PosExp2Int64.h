#ifndef INT64POS_EXP2INT64_H_
#define INT64POS_EXP2INT64_H_

namespace floatingExp2Integer
{
    class Int64PosExp2Int64 {
        private:
            Int64PosExp2Int64(double significand, uint64_t exponent);
            inline void fromDouble(double dbl);
            inline void checkRuleForScale();
            inline void scale();
        public:
            uint64_t scnfcnd;
            int64_t exp;

            Int64PosExp2Int64();
            Int64PosExp2Int64(double dbl);

            void double_to(double dbl);
            void log2_to(double log2);
            static void log2s_to(const std::vector<double>& from, std::vector<floatingExp2Integer::Int64PosExp2Int64>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            void sum(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& vector);
            void multiply(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& vector);

            double as_log2() const;
            double as_double() const;

            uint64_t sicnificand() const { return scnfcnd; };
            int64_t exponent() const { return exp; };

            Int64PosExp2Int64& operator+=(Int64PosExp2Int64 z);
            Int64PosExp2Int64& operator*=(Int64PosExp2Int64 z);
    };

    Int64PosExp2Int64 operator+(Int64PosExp2Int64 a, const Int64PosExp2Int64 b);
    Int64PosExp2Int64 operator*(Int64PosExp2Int64 a, const Int64PosExp2Int64 b);
}

#endif
