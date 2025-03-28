#ifndef DBL2_H_
#define DBL2_H_

#include <vector>

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

            Dbl2& operator*=(Dbl2 d) {
                dbl1 *= d.dbl1;
                dbl2 *= d.dbl2;
                return *this;
            }

            void sumDbl2(const std::vector<floatingExp2Integer::Dbl2>& dbl2Values) {
                floatingExp2Integer::Dbl2 dblSum1 = 0.0;
                floatingExp2Integer::Dbl2 dblSum2 = 0.0;
                floatingExp2Integer::Dbl2 dblSum3 = 0.0;
                floatingExp2Integer::Dbl2 dblSum4 = 0.0;
                floatingExp2Integer::Dbl2 dblSum5 = 0.0;
                floatingExp2Integer::Dbl2 dblSum6 = 0.0;
                floatingExp2Integer::Dbl2 dblSum7 = 0.0;
                floatingExp2Integer::Dbl2 dblSum8 = 0.0;
                unsigned int i;
                for (i = 0; i + 7 < dbl2Values.size(); i += 8) {
                    dblSum1.dbl1 += dbl2Values[i].dbl1;
                    dblSum2.dbl1 += dbl2Values[i + 1].dbl1;
                    dblSum3.dbl1 += dbl2Values[i + 2].dbl1;
                    dblSum4.dbl1 += dbl2Values[i + 3].dbl1;
                    dblSum5.dbl1 += dbl2Values[i + 4].dbl1;
                    dblSum6.dbl1 += dbl2Values[i + 5].dbl1;
                    dblSum7.dbl1 += dbl2Values[i + 6].dbl1;
                    dblSum8.dbl1 += dbl2Values[i + 7].dbl1;

                    dblSum1.dbl2 += dbl2Values[i].dbl2;
                    dblSum2.dbl2 += dbl2Values[i + 1].dbl2;
                    dblSum3.dbl2 += dbl2Values[i + 2].dbl2;
                    dblSum4.dbl2 += dbl2Values[i + 3].dbl2;
                    dblSum5.dbl2 += dbl2Values[i + 4].dbl2;
                    dblSum6.dbl2 += dbl2Values[i + 5].dbl2;
                    dblSum7.dbl2 += dbl2Values[i + 6].dbl2;
                    dblSum8.dbl2 += dbl2Values[i + 7].dbl2;
                }
                for (i = i; i < dbl2Values.size(); i++) {
                    dblSum1.dbl1 += dbl2Values[i].dbl1;
                    dblSum1.dbl2 += dbl2Values[i].dbl2;
                }
                dbl1 = dblSum1.dbl1 + dblSum2.dbl1 + dblSum3.dbl1 + dblSum4.dbl1 + dblSum5.dbl1 + dblSum6.dbl1 + dblSum7.dbl1 + dblSum8.dbl1;
                dbl2 = dblSum1.dbl2 + dblSum2.dbl2 + dblSum3.dbl2 + dblSum4.dbl2 + dblSum5.dbl2 + dblSum6.dbl2 + dblSum7.dbl2 + dblSum8.dbl2;
            }

            void multiplyDbl2(const std::vector<floatingExp2Integer::Dbl2>& dbl2Values) {
                floatingExp2Integer::Dbl2 dblSum1 = 0.0;
                floatingExp2Integer::Dbl2 dblSum2 = 0.0;
                floatingExp2Integer::Dbl2 dblSum3 = 0.0;
                floatingExp2Integer::Dbl2 dblSum4 = 0.0;
                floatingExp2Integer::Dbl2 dblSum5 = 0.0;
                floatingExp2Integer::Dbl2 dblSum6 = 0.0;
                floatingExp2Integer::Dbl2 dblSum7 = 0.0;
                floatingExp2Integer::Dbl2 dblSum8 = 0.0;
                unsigned int i;
                for (i = 0; i + 7 < dbl2Values.size(); i += 8) {
                    dblSum1.dbl1 *= dbl2Values[i].dbl1;
                    dblSum2.dbl1 *= dbl2Values[i + 1].dbl1;
                    dblSum3.dbl1 *= dbl2Values[i + 2].dbl1;
                    dblSum4.dbl1 *= dbl2Values[i + 3].dbl1;
                    dblSum5.dbl1 *= dbl2Values[i + 4].dbl1;
                    dblSum6.dbl1 *= dbl2Values[i + 5].dbl1;
                    dblSum7.dbl1 *= dbl2Values[i + 6].dbl1;
                    dblSum8.dbl1 *= dbl2Values[i + 7].dbl1;

                    dblSum1.dbl2 *= dbl2Values[i].dbl2;
                    dblSum2.dbl2 *= dbl2Values[i + 1].dbl2;
                    dblSum3.dbl2 *= dbl2Values[i + 2].dbl2;
                    dblSum4.dbl2 *= dbl2Values[i + 3].dbl2;
                    dblSum5.dbl2 *= dbl2Values[i + 4].dbl2;
                    dblSum6.dbl2 *= dbl2Values[i + 5].dbl2;
                    dblSum7.dbl2 *= dbl2Values[i + 6].dbl2;
                    dblSum8.dbl2 *= dbl2Values[i + 7].dbl2;
                }
                for (i = i; i < dbl2Values.size(); i++) {
                    dblSum1.dbl1 *= dbl2Values[i].dbl1;
                    dblSum1.dbl2 *= dbl2Values[i].dbl2;
                }
                dbl1 = dblSum1.dbl1 * dblSum2.dbl1 * dblSum3.dbl1 * dblSum4.dbl1 * dblSum5.dbl1 * dblSum6.dbl1 * dblSum7.dbl1 * dblSum8.dbl1;
                dbl2 = dblSum1.dbl2 * dblSum2.dbl2 * dblSum3.dbl2 * dblSum4.dbl2 * dblSum5.dbl2 * dblSum6.dbl2 * dblSum7.dbl2 * dblSum8.dbl2;
            }
    };
}

#endif
