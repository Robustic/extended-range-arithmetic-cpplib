#ifndef FLOAT64_LARGERANGENUMBER_H_
#define FLOAT64_LARGERANGENUMBER_H_

#include <vector>
#include <cmath>
#include <iostream>

namespace floatingExp2Integer
{
    class WideRangeNumber64 {
        private:
            double scnfcnd;
            int64_t exp;
            double encoded;
        public:
            WideRangeNumber64() {
                exp = 0;
                scnfcnd = 1;
                encoded = 0;
            }

            void set_encoded(double encoded_value) {
                double floored = std::floor(encoded_value);
                exp = (int64_t)floored;
                scnfcnd = encoded_value - floored + 1;
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
