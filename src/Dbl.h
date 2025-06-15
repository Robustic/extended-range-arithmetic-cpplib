#ifndef DBL_H_
#define DBL_H_

#include <vector>

namespace dbl
{
    class Dbl {
        private:
        public:
            double dbl;
            Dbl() { dbl = 1.0; }
            Dbl(double d) { dbl = d; }

            void log2_to(const double from) { dbl = std::exp2(from); }
            static void log2s_to(const std::vector<double>& from, std::vector<dbl::Dbl>& to) {
                for (size_t i = 0; i < to.size(); i++) {
                    to[i].log2_to(from[i]);
                }
            }

            double as_double() { return dbl; }
            double as_log2() { return std::log2(dbl); }

            Dbl& operator+=(Dbl d) { 
                dbl += d.dbl;
                return *this;
            }

            Dbl& operator*=(Dbl d) {
                dbl *= d.dbl;
                return *this;
            }

            void sum(const std::vector<dbl::Dbl>& dblValues) {
                dbl::Dbl dblSum1(0.0);
                dbl::Dbl dblSum2(0.0);
                dbl::Dbl dblSum3(0.0);
                dbl::Dbl dblSum4(0.0);
                dbl::Dbl dblSum5(0.0);
                dbl::Dbl dblSum6(0.0);
                dbl::Dbl dblSum7(0.0);
                dbl::Dbl dblSum8(0.0);
                size_t i = 0;
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

            void multiply(const std::vector<dbl::Dbl>& dblValues) {
                dbl::Dbl dblSum1(1.0);
                dbl::Dbl dblSum2(1.0);
                dbl::Dbl dblSum3(1.0);
                dbl::Dbl dblSum4(1.0);
                dbl::Dbl dblSum5(1.0);
                dbl::Dbl dblSum6(1.0);
                dbl::Dbl dblSum7(1.0);
                dbl::Dbl dblSum8(1.0);
                size_t i = 0;
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
