#ifndef DBL_H_
#define DBL_H_

namespace floatingExp2Integer
{
    class Dbl {
        private:
            double dbl;
        public:
            Dbl() { dbl = 0; }
            Dbl(double d) { dbl = d; }
            double asDouble() { return dbl; }

            Dbl& operator+=(Dbl d) { 
                dbl += d.dbl;
                return *this;
            }
    };
}

#endif
