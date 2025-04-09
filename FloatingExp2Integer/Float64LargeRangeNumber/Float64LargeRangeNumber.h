#ifndef FLOAT64_LARGERANGENUMBER_H_
#define FLOAT64_LARGERANGENUMBER_H_

#include <vector>
#include <iostream>

namespace floatingExp2Integer
{
    class Float64LargeRangeNumber {
        private:
        public:
            double encoded;
            Float64LargeRangeNumber() { encoded = 1.0; }
            Float64LargeRangeNumber(double enc) { encoded = enc; }
            static double double_to(double dbl);
            static void doubles_to(const std::vector<double>& in, std::vector<double>& out);
            static double as_double(double dbl);
            double as_double();

            void log2_to(double dbl);
            static void log2s_to(const std::vector<double>& in, std::vector<floatingExp2Integer::Float64LargeRangeNumber>& out) {
                for (size_t i = 0; i < out.size(); i++) {
                    out[i].log2_to(in[i]);
                }
            }
            static double as_log2(double dbl);

            static double sum(double lrn1, double lrn2);
            static double multiply(double lrn1, double lrn2);

            static double sum(std::vector<double>& lrns);
            static double multiply(std::vector<double>& lrns);
    };
}

#endif
