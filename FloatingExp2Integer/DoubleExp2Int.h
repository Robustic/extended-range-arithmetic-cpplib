#ifndef DOUBLE_EXP2_INT_H_
#define DOUBLE_EXP2_INT_H_

namespace floatingExp2Integer
{
    class DoubleExp2Int {
        private:
            double mant;
            std::int32_t exp;        
        public:
            DoubleExp2Int(double number);
            double mantissa() const { return mant; }
            int exponent() const { return exp; }
            double asDouble();
    };
}

#endif
