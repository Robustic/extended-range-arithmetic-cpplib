#include <iostream>
#include <limits>
#include <functional>
#include <vector>
#include <immintrin.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>
#include "Timer.h"
#include "Dbl.h"
#include "Dbl2.h"
#include "Log2Scale.h"
#include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
#include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
#include "./Float64LargeRangeNumber/Float64LargeRangeNumber.h"
#include "./Fukushima/Fukushima.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec, double min_log2, double max_log2) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(min_log2, max_log2);
    for (unsigned int i = 0; i < vec.size(); i++) {
        double a_random_doable = unif(rng);
        vec[i] = a_random_doable;
    }
}

constexpr size_t n[] = { 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000 };
constexpr size_t n_rounds[] = { 100000, 30000, 10000, 3000, 1000, 300, 100, 30, 10, 3, 1 };
//constexpr size_t n[] = { 1 };
//constexpr size_t n_rounds[] = { 1 };
constexpr size_t n_count = sizeof(n) / sizeof(n[0]);

struct ResultCollection {
    std::string header;
    std::array<double, n_count> results;
    std::array<int64_t, n_count> times;
};

template<typename T>
void calculate_avg(const std::vector<double>& values, std::vector<ResultCollection>& resultCollections, std::function<int64_t(const std::vector<T>&, double&, std::string&)> function) {

    ResultCollection rc;

    for (size_t i = 0; i < n_count; i++) {
        size_t n_current = n[i];
        size_t n_rounds_current = n_rounds[i];

        std::vector<T> values_converted(n_current);
        T::log2s_to(values, values_converted);

        int64_t time_sum = 0;

        double result;
        std::string header;
        for (size_t i = 0; i < n_rounds_current; i++) {
            time_sum += function(values_converted, result, header);
        }

        rc.header = header;
        rc.results[i] = result;
        rc.times[i] = time_sum / n_rounds_current;
    }
    resultCollections.push_back(rc);
}

// ******************************************************************************************************



void DoubleToInt64PosExp2Int64Values(const std::vector<double>& doubleValues,
    std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values) {
    for (unsigned int i = 0; i < int64PosExp2Int64Values.size(); i++) {
        int64PosExp2Int64Values[i].doubleToInt64PosExp2Int64(doubleValues[i]);
    }
}

void DoubleToFloat64PosExp2Int64Values(const std::vector<double>& doubleValues,
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    for (unsigned int i = 0; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Values[i].doubleToFloat64PosExp2Int64(doubleValues[i]);
    }
}

/////////////////////////////////

/// double sum and multiply

int64_t calculate_array_sum_dbl(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0; 
    if (16 <= values_converted.size()) {
        __m512d vsum = _mm512_loadu_pd(&values_converted[0]);
        for (i = 8; i + 7 < values.size(); i += 8) {
            __m512d vb = _mm512_loadu_pd(&values_converted[i]);
            vsum = vsum + vb;
        }

        for (unsigned int k = 0; k < 8; k++) {
            sum += vsum[k];
        }
    }

    for (i = i; i < values_converted.size(); i++) {
        sum += values_converted[i];
    }
    timer.stop();
    result = std::log2(sum);
    return timer.time();
}

int64_t calculate_array_multiply_dbl(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double multiplied = 1.0;
    if (16 <= values_converted.size()) {
        __m512d vmul = _mm512_set1_pd(1.0);
        for (i = 0; i + 7 < values_converted.size(); i += 8) {
            __m512d vb = _mm512_loadu_pd(&values_converted[i]);
            vmul = vmul * vb;
        }

        for (unsigned int k = 0; k < 8; k++) {
            multiplied *= vmul[k];
        }
    }

    for (i = i; i < values_converted.size(); i++) {
        multiplied *= values_converted[i];
    }
    timer.stop();
    result = std::log2(multiplied);
    return timer.time();
}

/// Log2 sum trick and (sequential and array)

int64_t calculate_sequential_sum_log2(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double sum = values_converted[0];
    for (unsigned int i = 1; i < values_converted.size(); i++) {
        double value = values_converted[i];
        if (sum > value) {
            sum = sum + std::log2(1 + std::exp2(value - sum));
        }
        else {
            sum = value + std::log2(std::exp2(sum - value) + 1);
        }
    }
    timer.stop();
    result = sum;
    return timer.time();
}

int64_t calculate_array_sum_log2(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0;
    double max = *std::max_element(values_converted.begin(), values_converted.end());
    for (i = 0; i < values_converted.size(); i++) {
        sum += std::exp2(values_converted[i] - max);
    }
    timer.stop();
    result = max + std::log2(sum);
    return timer.time();
}

int64_t calculate_array_multiply_log2(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0;
    for (i = 0; i < values_converted.size(); i++) {
        sum += values_converted[i];
    }
    timer.stop();
    result = sum;
    return timer.time();
}

/// Dbl1 and Dbl2 sum



int64_t calculate_sequential_sum_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum = 0.0;
    for (size_t i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t calculate_sequential_sum_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum = 0.0;
    for (size_t i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t calculate_array_sum_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum = 0.0;
    sum.sum(values);
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t  calculate_array_sum_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum = 0.0;
    sum.sumDbl2(values);
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

/// Dbl1 and Dbl2 multiply

int64_t  calculate_sequential_multiply_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();
    result = multiplied.as_log2();
    return timer.time();
}

int64_t  calculate_sequential_multiply_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();
    result = multiplied.as_log2();
    return timer.time();
}

int64_t  calculate_array_multiply_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied = 1.0;
    multiplied.multiply(values);
    timer.stop();
    result = multiplied.as_log2();
    return timer.time();
}

int64_t  calculate_array_multiply_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied = 1.0;
    multiplied.multiplyDbl2(values);
    timer.stop();
    result = multiplied.as_log2();
    return timer.time();
}

///

int64_t  calculate_sequential_sum_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima sum = 0.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        sum += values[i];

        if (sum.scnfcnd > 0.1) {
            i++;
        }
    }
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t  calculate_array_sum_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.sum(values);
    timer.stop();
    result = collector.as_log2();
    return timer.time();
}

int64_t  calculate_sequential_multiply_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima res = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        res *= values[i];

        if (res.scnfcnd > 0x1.999p479) {
            i++;
        }
    }
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  calculate_array_multiply_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.multiply(values);
    timer.stop();
    result = collector.as_log2();
    return timer.time();
}

int64_t  calculate_array_sum_Float64LargeRangeNumber(std::vector<double> values, double& result) {
    floatingExp2Integer::Timer timer;
    double sum = floatingExp2Integer::Float64LargeRangeNumber::sum_largeRangeNumbers(values);
    result = floatingExp2Integer::Float64LargeRangeNumber::largeRangeNumber_to_double(sum);
    timer.stop();
    return timer.time();
}

int64_t  calculate_avg_array_sum_Float64LargeRangeNumber(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64LargeRangeNumber_array_sum:";

    std::vector<double> values_converted(values.size());
    floatingExp2Integer::Float64LargeRangeNumber::doubles_to_largeRangeNumbers(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Float64LargeRangeNumber(values_converted, result);
    }

    return time_sum / n_rounds;
}

int64_t  calculate_sequential_sum_Float64LargeRangeNumber(std::vector<double> values, double& result) {
    floatingExp2Integer::Timer timer;
    double sum = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        double value = values[i];
        sum = floatingExp2Integer::Float64LargeRangeNumber::sum_largeRangeNumbers(sum, value);
        if (sum > 1.0) {
            i++;
        }
    }
    result = floatingExp2Integer::Float64LargeRangeNumber::largeRangeNumber_to_double(sum);
    timer.stop();
    return timer.time();
}

int64_t  calculate_avg_sequential_sum_Float64LargeRangeNumber(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64LargeRangeNumber_sequential_sum:";

    std::vector<double> values_converted(values.size());
    floatingExp2Integer::Float64LargeRangeNumber::doubles_to_largeRangeNumbers(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Float64LargeRangeNumber(values_converted, result);
    }

    return time_sum / n_rounds;
}

int64_t  calculate_array_multiply_Float64LargeRangeNumber(std::vector<double> values, double& result) {
    floatingExp2Integer::Timer timer;
    double res = floatingExp2Integer::Float64LargeRangeNumber::multiply_largeRangeNumbers(values);
    result = res;
    timer.stop();
    return timer.time();
}

int64_t  calculate_avg_array_multiply_Float64LargeRangeNumber(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64LargeRangeNumber_array_multiply:";

    std::vector<double> values_converted(values.size());
    floatingExp2Integer::Float64LargeRangeNumber::doubles_to_largeRangeNumbers(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_Float64LargeRangeNumber(values_converted, result);
    }

    return time_sum / n_rounds;
}

int64_t  calculate_sequential_multiply_Float64LargeRangeNumber(std::vector<double> values, double& result) {
    floatingExp2Integer::Timer timer;
    double res = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        double value = values[i];
        res = floatingExp2Integer::Float64LargeRangeNumber::multiply_largeRangeNumbers(res, value);
        if (res > 1.0) {
            i++;
        }
    }
    result = floatingExp2Integer::Float64LargeRangeNumber::largeRangeNumber_to_log2(res);
    timer.stop();
    return timer.time();
}

int64_t  calculate_avg_sequential_multiply_Float64LargeRangeNumber(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64LargeRangeNumber_sequential_multiply:";

    std::vector<double> values_converted(values.size());
    floatingExp2Integer::Float64LargeRangeNumber::doubles_to_largeRangeNumbers(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_multiply_Float64LargeRangeNumber(values_converted, result);
    }

    return time_sum / n_rounds;
}

int64_t  calculate_sequential_sum_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res = values[0];
    for (unsigned int i = 1; i < values.size(); i++) {
        res += values[i];

        if (res.scnfcnd < 1000ull) {
            i++;
        }
    }
    timer.stop();
    result = res.asDouble();
    return timer.time();
}

int64_t  calculate_avg_sequential_sum_Int64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Int64PosExp2Int64_sequential_sum:";

    std::vector<floatingExp2Integer::Int64PosExp2Int64> values_converted(values.size());
    DoubleToInt64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Int64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int64_t  calculate_array_sum_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res;
    res.sum(values);
    timer.stop();
    result = res.asDouble();
    return timer.time();
}

int64_t  calculate_avg_array_sum_Int64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Int64PosExp2Int64_array_sum:";

    std::vector<floatingExp2Integer::Int64PosExp2Int64> values_converted(values.size());
    DoubleToInt64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Int64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int64_t  calculate_sequential_multiply_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        res *= values[i];

        if (res.exp > -3ll) {
            i++;
        }
    }
    timer.stop();
    result = res.int64PosExp2Int64ToLog2();
    return timer.time();
}

int64_t  calculate_avg_sequential_multiply_Int64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Int64PosExp2Int64_sequential_multiply:";

    std::vector<floatingExp2Integer::Int64PosExp2Int64> values_converted(values.size());
    DoubleToInt64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_multiply_Int64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int64_t  calculate_array_multiply_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res;
    res.multiply(values);
    timer.stop();
    result = res.int64PosExp2Int64ToLog2();
    return timer.time();
}

int64_t  calculate_avg_array_multiply_Int64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Int64PosExp2Int64_array_multiply:";

    std::vector<floatingExp2Integer::Int64PosExp2Int64> values_converted(values.size());
    DoubleToInt64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_Int64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

///////////////////////////////////////////////

int64_t  calculate_sequential_sum_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res = values[0];
    for (unsigned int i = 1; i < values.size(); i++) {
        res += values[i];

        if (res.scnfcnd < 0.5) {
            i++;
        }
    }
    timer.stop();
    result = res.asDouble();
    return timer.time();
}

int64_t  calculate_avg_sequential_sum_Float64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64PosExp2Int64_sequential_sum:";

    std::vector<floatingExp2Integer::Float64PosExp2Int64> values_converted(values.size());
    DoubleToFloat64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Float64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int64_t  calculate_array_sum_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res;
    res.sum(values);
    timer.stop();
    result = res.asDouble();
    return timer.time();
}

int64_t  calculate_avg_array_sum_Float64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64PosExp2Int64_array_sum:";

    std::vector<floatingExp2Integer::Float64PosExp2Int64> values_converted(values.size());
    DoubleToFloat64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Float64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int64_t  calculate_sequential_multiply_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        res *= values[i];

        if (res.exp > -3ll) {
            i++;
        }
    }
    timer.stop();
    result = res.float64PosExp2Int64ToLog2();
    return timer.time();
}

int64_t  calculate_avg_sequential_multiply_Float64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64PosExp2Int64_sequential_multiply:";

    std::vector<floatingExp2Integer::Float64PosExp2Int64> values_converted(values.size());
    DoubleToFloat64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_multiply_Float64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int64_t  calculate_array_multiply_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res;
    res.multiply(values);
    timer.stop();
    result = res.float64PosExp2Int64ToLog2();
    return timer.time();
}

int64_t  calculate_avg_array_multiply_Float64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64PosExp2Int64_array_multiply:";

    std::vector<floatingExp2Integer::Float64PosExp2Int64> values_converted(values.size());
    DoubleToFloat64PosExp2Int64Values(values, values_converted);

    int64_t  time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_Float64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

int main() {
    constexpr double min_log2 = -33;
    constexpr double max_log2 = -30;



    //std::vector<std::function<int64_t (std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result, std::function<int64_t(void)>)>> functions;
    //std::vector < std::function<int64_t(void)> > functions2;
    //functions.push_back(calculate_avg<floatingExp2Integer::Dbl>);
    //functions2.push_back(calculate_sequential_sum_Dbl1);
    //functions.push_back(calculate_avg_sequential_sum_Dbl2);
    //functions.push_back(calculate_avg_sequential_sum_log2);
    //functions.push_back(calculate_avg_sequential_sum_Fukushima);
    //functions.push_back(calculate_avg_sequential_sum_Float64LargeRangeNumber);
    //functions.push_back(calculate_avg_sequential_sum_Int64PosExp2Int64);
    //functions.push_back(calculate_avg_sequential_sum_Float64PosExp2Int64);
    //functions.push_back(calculate_avg_array_sum_dbl);
    //functions.push_back(calculate_avg_array_sum_Dbl1);
    //functions.push_back(calculate_avg_array_sum_Dbl2);
    //functions.push_back(calculate_avg_array_sum_log2);
    //functions.push_back(calculate_avg_array_sum_Fukushima);
    //functions.push_back(calculate_avg_array_sum_Float64LargeRangeNumber);
    //functions.push_back(calculate_avg_array_sum_Int64PosExp2Int64);
    //functions.push_back(calculate_avg_array_sum_Float64PosExp2Int64);
    //functions.push_back(calculate_avg_sequential_multiply_Dbl1);
    //functions.push_back(calculate_avg_sequential_multiply_Dbl2);
    //functions.push_back(calculate_avg_sequential_multiply_Fukushima);
    //functions.push_back(calculate_avg_sequential_multiply_Float64LargeRangeNumber);
    //functions.push_back(calculate_avg_sequential_multiply_Int64PosExp2Int64);
    //functions.push_back(calculate_avg_sequential_multiply_Float64PosExp2Int64);
    //functions.push_back(calculate_avg_array_multiply_dbl);
    //functions.push_back(calculate_avg_array_multiply_Dbl1);
    //functions.push_back(calculate_avg_array_multiply_Dbl2);
    //functions.push_back(calculate_avg_array_multiply_log2);
    //functions.push_back(calculate_avg_array_multiply_Fukushima);
    //functions.push_back(calculate_avg_array_multiply_Float64LargeRangeNumber);
    //functions.push_back(calculate_avg_array_multiply_Int64PosExp2Int64);
    //functions.push_back(calculate_avg_array_multiply_Float64PosExp2Int64);



    std::vector<ResultCollection> resultCollections;

    constexpr unsigned int n_max = *std::max_element(n, n + n_count);
    std::vector<double> double_values(n_max);
    InitializeRandomNumbers(double_values, min_log2, max_log2);


    calculate_avg<floatingExp2Integer::Log2Scale>(double_values, resultCollections, calculate_sequential_sum_log2);
    calculate_avg<floatingExp2Integer::Dbl>(double_values, resultCollections, calculate_sequential_sum_Dbl1);
    calculate_avg<floatingExp2Integer::Dbl2>(double_values, resultCollections, calculate_sequential_sum_Dbl2);
    calculate_avg<floatingExp2Integer::Fukushima>(double_values, resultCollections, calculate_sequential_sum_Fukushima);

    calculate_avg<floatingExp2Integer::Dbl>(double_values, resultCollections, calculate_array_sum_dbl);
    calculate_avg<floatingExp2Integer::Dbl>(double_values, resultCollections, calculate_array_sum_Dbl1);
    calculate_avg<floatingExp2Integer::Dbl2>(double_values, resultCollections, calculate_array_sum_Dbl2);
    calculate_avg<floatingExp2Integer::Log2Scale>(double_values, resultCollections, calculate_array_sum_log2);
    calculate_avg<floatingExp2Integer::Fukushima>(double_values, resultCollections, calculate_array_sum_Fukushima);

    calculate_avg<floatingExp2Integer::Dbl>(double_values, resultCollections, calculate_sequential_multiply_Dbl1);
    calculate_avg<floatingExp2Integer::Dbl2>(double_values, resultCollections, calculate_sequential_multiply_Dbl2);
    calculate_avg<floatingExp2Integer::Fukushima>(double_values, resultCollections, calculate_sequential_multiply_Fukushima);

    calculate_avg<floatingExp2Integer::Dbl>(double_values, resultCollections, calculate_array_multiply_dbl);
    calculate_avg<floatingExp2Integer::Log2Scale>(double_values, resultCollections, calculate_array_multiply_log2);
    calculate_avg<floatingExp2Integer::Dbl>(double_values, resultCollections, calculate_array_multiply_Dbl1);
    calculate_avg<floatingExp2Integer::Dbl2>(double_values, resultCollections, calculate_array_multiply_Dbl2);
    calculate_avg<floatingExp2Integer::Fukushima>(double_values, resultCollections, calculate_array_multiply_Fukushima);


    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    std::cout << std::endl << "Result" << std::endl << "n ";
    for (size_t i = 0; i < n_count; i++) {
        std::cout << n[i] << " ";
    }

    std::cout << std::endl;
    for (size_t f = 0; f < resultCollections.size(); f++) {
        ResultCollection& current = resultCollections[f];
        std::cout << current.header << " ";
        for (size_t i = 0; i < n_count; i++) {
            std::cout << current.results[i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "Time" << std::endl << "n ";
    for (size_t i = 0; i < n_count; i++) {
        std::cout << n[i] << " ";
    }

    std::cout << std::endl;
    for (size_t f = 0; f < resultCollections.size(); f++) {
        ResultCollection& current = resultCollections[f];
        std::cout << current.header << " ";
        for (size_t i = 0; i < n_count; i++) {
            std::cout << current.times[i] << " ";
        }
        std::cout << std::endl;
    }

    //std::cout << "Press any key to exit... ";
    //std::cin.get();
    //std::cout << std::endl;

    resultCollections.clear();

    return 0;
}
