#ifndef FLOAT64_LARGERANGENUMBER_H_
#define FLOAT64_LARGERANGENUMBER_H_

#include <vector>

namespace floatingExp2Integer
{
    class Float64LargeRangeNumber {
        public:
            static double double_to_largeRangeNumber(double dbl);
            static void doubles_to_largeRangeNumbers(const std::vector<double>& in, std::vector<double>& out);
            static double largeRangeNumber_to_double(double dbl);

            static double log2_to_largeRangeNumber(double dbl);
            static void log2s_to_largeRangeNumbers(const std::vector<double>& in, std::vector<double>& out);
            static double largeRangeNumber_to_log2(double dbl);

            static double sum_largeRangeNumbers(double lrn1, double lrn2);
            static double multiply_largeRangeNumbers(double lrn1, double lrn2);

            static double sum_largeRangeNumbers(std::vector<double>& lrns);
            static double multiply_largeRangeNumbers(std::vector<double>& lrns);
    };
}

#endif
