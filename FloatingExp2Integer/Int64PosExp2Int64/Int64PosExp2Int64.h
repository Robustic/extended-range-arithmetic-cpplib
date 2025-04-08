#ifndef INT64POS_EXP2INT64_H_
#define INT64POS_EXP2INT64_H_

namespace floatingExp2Integer
{
    class Int64PosExp2Int64 {
        private:
            inline void fromDouble(double dbl);
            inline void checkRuleForScale();
            inline void scale();
        public:
            std::uint64_t scnfcnd;
            std::int64_t exp;
            Int64PosExp2Int64();
            Int64PosExp2Int64(double dbl);
            Int64PosExp2Int64(double significand, std::uint64_t exponent);
            void doubleToInt64PosExp2Int64(double dbl);
            void log2ToInt64PosExp2Int64(double log2);
            double int64PosExp2Int64ToLog2() const;
            std::uint64_t sicnificand();
            std::int64_t exponent();
            double asDouble() const;

            void sum(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& vector);
            void multiply(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& vector);

            Int64PosExp2Int64& operator+=(Int64PosExp2Int64 z);
            Int64PosExp2Int64& operator*=(Int64PosExp2Int64 z);
    };

    Int64PosExp2Int64 operator+(Int64PosExp2Int64 a, const Int64PosExp2Int64 b);
    Int64PosExp2Int64 operator*(Int64PosExp2Int64 a, const Int64PosExp2Int64 b);
}

#endif
