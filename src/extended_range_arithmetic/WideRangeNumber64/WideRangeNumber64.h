#ifndef WIDERANGENUMBER_64_H_
#define WIDERANGENUMBER_64_H_

#include <vector>
#include <cmath>
#include <iostream>

namespace extended_range_arithmetic
{
    class WideRangeNumber64 {
        private:
            double principal;
            int64_t aux;
            double encoded;
        public:
            WideRangeNumber64() {
                aux = 0;
                principal = 1;
                encoded = 0;
            }

            void set_encoded(double encoded_value) {
                double floored = std::floor(encoded_value);
                aux = (int64_t)floored;
                principal = encoded_value - floored + 1;
                encoded = encoded_value;
            }

            static double double_to(double dbl);
            static void doubles_to(const std::vector<double>& in, std::vector<double>& out);
            static double as_double(double dbl);

            static double log2_to(double dbl);
            static void log2s_to(const std::vector<double>& in, std::vector<double>& out) {
                for (size_t i = 0; i < out.size(); i++) {
                    out[i] = log2_to(in[i]);
                }
            }
            static double as_log2(double dbl);
            double as_log2() {
                return WideRangeNumber64::as_log2(encoded);
            }

            static double sum(double lrn1, double lrn2);
            static double multiply(double lrn1, double lrn2);

            static double sum(const std::vector<double>& lrns);
            static double multiply(const std::vector<double>& lrns);

            WideRangeNumber64& operator+=(double encoded_2);
            WideRangeNumber64& operator*=(double encoded_2);
    };
}

#endif
