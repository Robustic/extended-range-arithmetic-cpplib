#ifndef FLOAT32_EXP2INT32_H_
#define FLOAT32_EXP2INT32_H_

namespace floatingExp2Integer
{
    class Float32Exp2Int32 {
        private:
            float scnfcnd;
            std::int32_t exp;
            inline void checkLimitForScale();
            inline void scale();
        public:
            Float32Exp2Int32();
            Float32Exp2Int32(float flt);
            Float32Exp2Int32(float sicnificand, std::int32_t exponent);
            void log2ToFloat32Exp2Int32(float logarithm2);
            float float32Exp2Int32ToLog2();
            float sicnificand();
            std::int32_t exponent();
            float asFloat() const;

            Float32Exp2Int32& operator+=(Float32Exp2Int32 z);
            Float32Exp2Int32& operator*=(Float32Exp2Int32 z);
            Float32Exp2Int32& operator*=(float& flt);
    };

    Float32Exp2Int32 operator+(Float32Exp2Int32 a, Float32Exp2Int32 b);
    Float32Exp2Int32 operator*(Float32Exp2Int32 a, Float32Exp2Int32 b);
}

#endif
