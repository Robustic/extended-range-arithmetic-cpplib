#ifndef DBL3_H_
#define DBL3_H_

#include <vector>

namespace floatingExp2Integer
{
    class Dbl3 {
        private:
            double dbl1;
            double dbl2;
            double dbl3;
        public:
            Dbl3() { dbl1 = 0; dbl2 = 0; dbl3 = 0; }
            Dbl3(double d) { dbl1 = d; dbl2 = -d; dbl3 = d; }
            Dbl3(std::vector<floatingExp2Integer::Dbl3>& dbl3Values) { 
                floatingExp2Integer::Dbl3 dblSum1 = 0.0;
                floatingExp2Integer::Dbl3 dblSum2 = 0.0;
                floatingExp2Integer::Dbl3 dblSum3 = 0.0;
                unsigned int i;
                for (i = 0; i + 2 < dbl3Values.size(); i += 3) {
                    dblSum1.dbl1 += dbl3Values[i].dbl1;
                    dblSum2.dbl1 += dbl3Values[i+1].dbl1;
                    dblSum3.dbl1 += dbl3Values[i+2].dbl1;

                    dblSum1.dbl2 += dbl3Values[i].dbl2;
                    dblSum2.dbl2 += dbl3Values[i+1].dbl2;
                    dblSum3.dbl2 += dbl3Values[i+2].dbl2;

                    dblSum1.dbl3 += dbl3Values[i].dbl3;
                    dblSum2.dbl3 += dbl3Values[i+1].dbl3;
                    dblSum3.dbl3 += dbl3Values[i+2].dbl3;
                }
                for (i = i; i < dbl3Values.size(); i++) {
                    dblSum1.dbl1 += dbl3Values[i].dbl1;
                    dblSum1.dbl2 += dbl3Values[i].dbl2;
                    dblSum1.dbl3 += dbl3Values[i].dbl3;
                }
                dbl1 = dblSum1.dbl1 + dblSum2.dbl1 + dblSum3.dbl1;
                dbl2 = dblSum1.dbl2 + dblSum2.dbl2 + dblSum3.dbl2;
                dbl3 = dblSum1.dbl3 + dblSum2.dbl3 + dblSum3.dbl3;
            }
            double asDouble() { return dbl1 + dbl2 + dbl3; }

            Dbl3& operator+=(Dbl3 d) {
                dbl1 += d.dbl1;
                dbl2 += d.dbl2;
                dbl3 += d.dbl3;
                return *this;
            }
    };
}

#endif
