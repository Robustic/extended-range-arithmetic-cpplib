#ifndef INT32POS_EXP2INT32_H_
#define INT32POS_EXP2INT32_H_

namespace floatingExp2Integer
{
    class Int32PosExp2Int32 {
        private:
            uint32_t scnfcnd;
            int32_t exp;
            inline void fromFloat(float flt);
            inline void checkRuleForScale();
            inline void scale();
        public:
            Int32PosExp2Int32();
            Int32PosExp2Int32(float flt);
            void floatToInt32PosExp2Int32(float flt);
            void log2ToInt32PosExp2Int32(float log2);
            float int32PosExp2Int32ToLog2() const;
            uint32_t sicnificand();
            int32_t exponent();
            float asFloat() const;

            Int32PosExp2Int32& operator+=(Int32PosExp2Int32 z);
            Int32PosExp2Int32& operator*=(Int32PosExp2Int32 z);
    };

    Int32PosExp2Int32 operator+(Int32PosExp2Int32 a, const Int32PosExp2Int32 b);
    Int32PosExp2Int32 operator*(Int32PosExp2Int32 a, const Int32PosExp2Int32 b);
}

#endif
