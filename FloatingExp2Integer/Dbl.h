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
            double asDouble() { return dbl; }

            Dbl& operator+=(Dbl d) { 
                dbl += d.dbl;
                return *this;
            }

            Dbl& operator*=(Dbl d) {
                dbl *= d.dbl;
                return *this;
            }

            void sumDbl(const std::vector<floatingExp2Integer::Dbl>& dblValues) {
                floatingExp2Integer::Dbl dblSum1 = 0.0;
                floatingExp2Integer::Dbl dblSum2 = 0.0;
                floatingExp2Integer::Dbl dblSum3 = 0.0;
                floatingExp2Integer::Dbl dblSum4 = 0.0;
                floatingExp2Integer::Dbl dblSum5 = 0.0;
                floatingExp2Integer::Dbl dblSum6 = 0.0;
                floatingExp2Integer::Dbl dblSum7 = 0.0;
                floatingExp2Integer::Dbl dblSum8 = 0.0;
                unsigned int i;
                for (i = 0; i + 7 < dblValues.size(); i += 8) {
                    dblSum1.dbl += dblValues[i].dbl;
                    dblSum2.dbl += dblValues[i + 1].dbl;
                    dblSum3.dbl += dblValues[i + 2].dbl;
                    dblSum4.dbl += dblValues[i + 3].dbl;
                    dblSum5.dbl += dblValues[i + 4].dbl;
                    dblSum6.dbl += dblValues[i + 5].dbl;
                    dblSum7.dbl += dblValues[i + 6].dbl;
                    dblSum8.dbl += dblValues[i + 7].dbl;
                }
                for (i = i; i < dblValues.size(); i++) {
                    dblSum1.dbl += dblValues[i].dbl;
                }
                dbl = dblSum1.dbl + dblSum2.dbl + dblSum3.dbl + dblSum4.dbl + dblSum5.dbl + dblSum6.dbl + dblSum7.dbl + dblSum8.dbl;
            }

            void multiplyDbl(const std::vector<floatingExp2Integer::Dbl>& dblValues) {
                floatingExp2Integer::Dbl dblSum1 = 0.0;
                floatingExp2Integer::Dbl dblSum2 = 0.0;
                floatingExp2Integer::Dbl dblSum3 = 0.0;
                floatingExp2Integer::Dbl dblSum4 = 0.0;
                floatingExp2Integer::Dbl dblSum5 = 0.0;
                floatingExp2Integer::Dbl dblSum6 = 0.0;
                floatingExp2Integer::Dbl dblSum7 = 0.0;
                floatingExp2Integer::Dbl dblSum8 = 0.0;
                unsigned int i;
                for (i = 0; i + 7 < dblValues.size(); i += 8) {
                    dblSum1.dbl *= dblValues[i].dbl;
                    dblSum2.dbl *= dblValues[i + 1].dbl;
                    dblSum3.dbl *= dblValues[i + 2].dbl;
                    dblSum4.dbl *= dblValues[i + 3].dbl;
                    dblSum5.dbl *= dblValues[i + 4].dbl;
                    dblSum6.dbl *= dblValues[i + 5].dbl;
                    dblSum7.dbl *= dblValues[i + 6].dbl;
                    dblSum8.dbl *= dblValues[i + 7].dbl;
                }
                for (i = i; i < dblValues.size(); i++) {
                    dblSum1.dbl *= dblValues[i].dbl;
                }
                dbl = dblSum1.dbl * dblSum2.dbl * dblSum3.dbl * dblSum4.dbl * dblSum5.dbl * dblSum6.dbl * dblSum7.dbl * dblSum8.dbl;
            }
    };
}

#endif
