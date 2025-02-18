#ifndef FLOAT32POS_EXP2INT32_H_
#define FLOAT32POS_EXP2INT32_H_

namespace floatingExp2Integer
{
    class Float32PosExp2Int32 {
        private:
            float scnfcnd;
            std::int32_t exp;
            inline void checkRuleForScale();
            inline void scale();
        public:
            Float32PosExp2Int32();
            Float32PosExp2Int32(float flt);
            void floatToFloat32PosExp2Int32(float flt);
            void log2ToFloat32PosExp2Int32(float log2);
            float float32PosExp2Int32ToLog2() const;
            float sicnificand();
            std::int32_t exponent();
            float asFloat() const;

            Float32PosExp2Int32& operator+=(Float32PosExp2Int32 z);
            Float32PosExp2Int32& operator*=(Float32PosExp2Int32 z);
            Float32PosExp2Int32& operator*=(float& flt);
    };

    Float32PosExp2Int32 operator+(Float32PosExp2Int32 a, const Float32PosExp2Int32 b);
    Float32PosExp2Int32 operator*(Float32PosExp2Int32 a, const Float32PosExp2Int32 b);
}

#endif

