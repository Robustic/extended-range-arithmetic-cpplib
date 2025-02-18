#ifndef DBL2_H_
#define DBL2_H_

namespace floatingExp2Integer
{
    class Dbl2 {
        private:
            double dbl1;
            double dbl2;
        public:
            Dbl2() { dbl1 = 0; dbl2 = 0; }
            Dbl2(double d) { dbl1 = d; dbl2 = -d; }
            double asDouble() { return 2 * dbl1 + dbl2; }

            Dbl2& operator+=(Dbl2 d) {
                dbl1 += d.dbl1;
                dbl2 += d.dbl2;
                return *this;
            }
    };
}

#endif
