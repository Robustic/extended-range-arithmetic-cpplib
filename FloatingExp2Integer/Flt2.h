#ifndef FLT2_H_
#define FLT2_H_

namespace floatingExp2Integer
{
    class Flt2 {
        private:
            float flt1;
            float flt2;
        public:
            Flt2() { flt1 = 0; flt2 = 0; }
            Flt2(float f) { flt1 = f; flt2 = -f; }
            float asFloat() { return 2 * flt1 + flt2; }

            Flt2& operator+=(Flt2 f) { 
                flt1 += f.flt1;
                flt2 += f.flt2;
                return *this;
            }
    };
}

#endif
