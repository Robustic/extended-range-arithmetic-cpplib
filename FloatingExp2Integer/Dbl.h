#ifndef DBL_H_
#define DBL_H_

#include <vector>

namespace floatingExp2Integer
{
    class Dbl {
        private:
            double dbl;
        public:
            Dbl() { dbl = 0; }
            Dbl(double d) { dbl = d; }
            Dbl(const std::vector<floatingExp2Integer::Dbl>& dblValues) { 
                floatingExp2Integer::Dbl dblSum1 = 0.0;
                floatingExp2Integer::Dbl dblSum2 = 0.0;
                floatingExp2Integer::Dbl dblSum3 = 0.0;
                floatingExp2Integer::Dbl dblSum4 = 0.0;
                unsigned int i;
                for (i = 0; i + 3 < dblValues.size(); i += 4) {
                    dblSum1.dbl += dblValues[i].dbl;
                    dblSum2.dbl += dblValues[i+1].dbl;
                    dblSum3.dbl += dblValues[i+2].dbl;
                    dblSum4.dbl += dblValues[i+3].dbl;
                }
                for (i = i; i < dblValues.size(); i++) {
                    dblSum1.dbl += dblValues[i].dbl;
                }
                dbl = dblSum1.dbl + dblSum2.dbl + dblSum3.dbl + dblSum4.dbl;
            }
            double asDouble() { return dbl; }

            Dbl& operator+=(Dbl d) { 
                dbl += d.dbl;
                return *this;
            }
    };
}

#endif
