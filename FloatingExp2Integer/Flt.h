#ifndef FLT_H_
#define FLT_H_

namespace floatingExp2Integer
{
    class Flt {
        private:
            float flt;
        public:
            Flt() { flt = 0; }
            Flt(float f) { flt = f; }
            float asFloat() { return flt; }

            Flt& operator+=(Flt f) { 
                flt += f.flt;
                return *this;
            }
    };
}

#endif
