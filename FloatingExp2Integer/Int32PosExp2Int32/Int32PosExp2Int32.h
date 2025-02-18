#ifndef INT32POS_EXP2INT32_H_
#define INT32POS_EXP2INT32_H_

namespace floatingExp2Integer
{
    class Int32PosExp2Int32 {
        private:
            std::uint32_t scnfcnd;
            std::int32_t exp;
            inline void checkLimitForScale();
            inline void scale();
        public:
            Int32PosExp2Int32();
            Int32PosExp2Int32(float flt);
            void floatToInt32PosExp2Int32(float flt);
            void log2ToInt32Exp2Int32(float logarithm2);
            float int32Exp2Int32ToLog2();
            std::int32_t sicnificand();
            std::int32_t exponent();
            float asFloat() const;

            Int32PosExp2Int32& operator+=(Int32PosExp2Int32 z);
            Int32PosExp2Int32& operator*=(Int32PosExp2Int32 z);
    };

    Int32PosExp2Int32 operator+(Int32PosExp2Int32 a, Int32PosExp2Int32 b);
    Int32PosExp2Int32 operator*(Int32PosExp2Int32 a, Int32PosExp2Int32 b);
}

#endif
