#ifndef DOUBLE_EXP2_INT_H_
#define DOUBLE_EXP2_INT_H_

namespace floatingExp2Integer
{
    class DoubleExp2Int {
        private:
            double mant;
            int exp;        
        public:
            DoubleExp2Int(double number) :mant{number}, exp{0} {}
            double mantissa() const { return mant; }
            int exponent() const { return exp; }
    };
}

#endif
