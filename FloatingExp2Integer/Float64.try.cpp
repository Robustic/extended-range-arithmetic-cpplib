#include <iostream>
#include <fstream>
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

// *******  PARAMETERS  *******

// CHANGE THESE VALUES TO CALCULATE WITH DIFFERENT VALUE RANGE AND DIFFERENT COUNT OF NUMBERS

constexpr double min_log2 = -10000;
constexpr double max_log2 = -10;

//constexpr size_t n[] = { 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000 };
//constexpr size_t n_rounds[] = { 100000, 30000, 10000, 3000, 1000, 300, 100, 30, 10, 3, 1 };
constexpr size_t n[] = { 10000000 };
constexpr size_t n_rounds[] = { 1 };

// *******  COMMON FUNCTIONS  *******

constexpr size_t n_count = sizeof(n) / sizeof(n[0]);

void save_vector_as_binary(const std::vector<double>& vec, const std::string& filename) {
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    for (const double& value : vec) {
        outFile.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }

    outFile.close();
}

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec, double min_log2, double max_log2) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(min_log2, max_log2);
    for (size_t i = 0; i < vec.size(); i++) {
        double a_random_double = unif(rng);
        vec[i] = a_random_double;
    }

    // save_vector_as_binary(vec, "data.bin");
}

// *******  TEMPLATE  *******

struct ResultCollection {
    std::string header;
    std::array<double, n_count> results;
    std::array<int64_t, n_count> times;
};

template<typename T>
void calc_perf(const std::vector<double>& values, std::vector<ResultCollection>& resultCollections,
        std::function<int64_t(const std::vector<T>&, double&, std::string&)> function) {

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

// *******  CASES  *******

// *******  double  *******

int64_t sum_sequential_dbl(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double sum = 0.0;
    for (size_t i = 0; i < values_converted.size(); i++) {
        sum += values_converted[i];
        if (sum > 0x1p999) {
            i = values_converted.size();
        }
    }
    timer.stop();

    result = std::log2(sum);
    return timer.time();
}

int64_t sum_parallel_dbl(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    size_t i = 0;
    double sum = 0.0;
    if (16 <= values_converted.size()) {
        __m512d vsum = _mm512_loadu_pd(&values_converted[0]);
        for (i = 8; i + 7 < values_converted.size(); i += 8) {
            __m512d vb = _mm512_loadu_pd(&values_converted[i]);
            vsum = vsum + vb;
        }

        for (size_t k = 0; k < 8; k++) {
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

int64_t multiply_sequential_dbl(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double multiplied = 1.0;
    for (size_t i = 0; i < values_converted.size(); i++) {
        multiplied *= values_converted[i];
        if (multiplied > 0x1p999) {
            i = values_converted.size();
        }
    }
    timer.stop();

    result = std::log2(multiplied);
    return timer.time();
}

int64_t multiply_parallel_dbl(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    size_t i = 0;
    double multiplied = 1.0;
    if (16 <= values_converted.size()) {
        __m512d vmul = _mm512_set1_pd(1.0);
        for (i = 0; i + 7 < values_converted.size(); i += 8) {
            __m512d vb = _mm512_loadu_pd(&values_converted[i]);
            vmul = vmul * vb;
        }

        for (size_t k = 0; k < 8; k++) {
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

// *******  log2scale  *******

int64_t sum_sequential_log2scale(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
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

int64_t sum_parallel_log2scale(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
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

int64_t multiply_sequential_log2scale(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0;
    for (i = 0; i < values_converted.size(); i++) {
        sum += values_converted[i];
        if (sum > 1e100) {
            i++;
        }
    }
    timer.stop();
    result = sum;
    return timer.time();
}

int64_t multiply_parallel_log2scale(const std::vector<floatingExp2Integer::Log2Scale>& values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    size_t i = 0;
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
    result = sum;
    return timer.time();
}

// *******  Dbl1  *******

int64_t sum_sequential_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
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

int64_t sum_parallel_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum = 0.0;
    sum.sum(values);
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
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

int64_t  multiply_parallel_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied = 1.0;
    multiplied.multiply(values);
    timer.stop();
    result = multiplied.as_log2();
    return timer.time();
}

// *******  Dbl2  *******

int64_t sum_sequential_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
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

int64_t  sum_parallel_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum = 0.0;
    sum.sumDbl2(values);
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
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

int64_t  multiply_parallel_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied = 1.0;
    multiplied.multiplyDbl2(values);
    timer.stop();
    result = multiplied.as_log2();
    return timer.time();
}

// *******  Fukushima  *******

int64_t  sum_sequential_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima sum = 0.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        sum += values[i];

        if (sum.scnfcnd > 0x1p970) {
            i++;
        }
    }
    timer.stop();
    result = sum.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.sum(values);
    timer.stop();
    result = collector.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
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

int64_t  multiply_parallel_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.multiply(values);
    timer.stop();
    result = collector.as_log2();
    return timer.time();
}

// *******  Float64LargeRangeNumber  *******

int64_t  sum_sequential_Float64LargeRangeNumber(std::vector<floatingExp2Integer::Float64LargeRangeNumber> values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double sum = values_converted[0];
    for (size_t i = 1; i < values_converted.size(); i++) {
        double value = values_converted[i];
        sum = floatingExp2Integer::Float64LargeRangeNumber::sum(sum, value);
        if (sum > 1.0) {
            i++;
        }
    }
    result = floatingExp2Integer::Float64LargeRangeNumber::as_log2(sum);
    timer.stop();
    return timer.time();
}

int64_t  sum_parallel_Float64LargeRangeNumber(std::vector<floatingExp2Integer::Float64LargeRangeNumber> values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double sum = floatingExp2Integer::Float64LargeRangeNumber::sum(values_converted);
    result = floatingExp2Integer::Float64LargeRangeNumber::as_log2(sum);
    timer.stop();
    return timer.time();
}

int64_t  multiply_sequential_Float64LargeRangeNumber(std::vector<floatingExp2Integer::Float64LargeRangeNumber> values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double res = values_converted[0];
    for (size_t i = 1; i < values_converted.size(); i++) {
        double value = values_converted[i];
        res = floatingExp2Integer::Float64LargeRangeNumber::multiply(res, value);
        if (res > 1.0) {
            i++;
        }
    }
    result = floatingExp2Integer::Float64LargeRangeNumber::as_log2(res);
    timer.stop();
    return timer.time();
}

int64_t  multiply_parallel_Float64LargeRangeNumber(std::vector<floatingExp2Integer::Float64LargeRangeNumber> values, double& result, std::string& header) {
    header = __func__;
    auto dblArray = std::bit_cast<double*>(values.data());
    std::vector<double> values_converted(dblArray, dblArray + values.size());

    floatingExp2Integer::Timer timer;
    double res = floatingExp2Integer::Float64LargeRangeNumber::multiply(values_converted);
    result = floatingExp2Integer::Float64LargeRangeNumber::as_log2(res);
    timer.stop();
    return timer.time();
}

// *******  Int64PosExp2Int64  *******

int64_t  sum_sequential_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res = values[0];
    for (unsigned int i = 1; i < values.size(); i++) {
        res += values[i];

        if (res.scnfcnd < 1000ull) {
            i++;
        }
    }
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res;
    res.sum(values);
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        res *= values[i];

        if (res.exp > -3ll) {
            i++;
        }
    }
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res;
    res.multiply(values);
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

// *******  Float64PosExp2Int64  *******

int64_t  sum_sequential_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res = values[0];
    for (unsigned int i = 1; i < values.size(); i++) {
        res += values[i];

        if (res.scnfcnd < 0.5) {
            i++;
        }
    }
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res;
    res.sum(values);
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        res *= values[i];

        if (res.exp > -3ll) {
            i++;
        }
    }
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result, std::string& header) {
    header = __func__;
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res;
    res.multiply(values);
    timer.stop();
    result = res.as_log2();
    return timer.time();
}

// *******  MAIN  *******

int main() {
    std::cout << "START" << std::endl;

    std::vector<ResultCollection> resultCollections;

    constexpr unsigned int n_max = *std::max_element(n, n + n_count);
    std::vector<double> double_values(n_max);
    InitializeRandomNumbers(double_values, min_log2, max_log2);

    std::cout << "RANDOM NUMBERS DEFINED" << std::endl;


    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, sum_sequential_dbl);
    calc_perf<floatingExp2Integer::Log2Scale>(double_values, resultCollections, sum_sequential_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, sum_sequential_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, sum_sequential_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, sum_sequential_Fukushima);
    calc_perf<floatingExp2Integer::Float64LargeRangeNumber>(double_values, resultCollections, sum_sequential_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, sum_sequential_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, sum_sequential_Float64PosExp2Int64);

    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, sum_parallel_dbl);
    calc_perf<floatingExp2Integer::Log2Scale>(double_values, resultCollections, sum_parallel_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, sum_parallel_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, sum_parallel_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, sum_parallel_Fukushima);
    calc_perf<floatingExp2Integer::Float64LargeRangeNumber>(double_values, resultCollections, sum_parallel_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, sum_parallel_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, sum_parallel_Float64PosExp2Int64);

    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, multiply_sequential_dbl);
    calc_perf<floatingExp2Integer::Log2Scale>(double_values, resultCollections, multiply_sequential_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, multiply_sequential_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, multiply_sequential_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, multiply_sequential_Fukushima);
    calc_perf<floatingExp2Integer::Float64LargeRangeNumber>(double_values, resultCollections, multiply_sequential_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, multiply_sequential_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, multiply_sequential_Float64PosExp2Int64);

    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, multiply_parallel_dbl);
    calc_perf<floatingExp2Integer::Log2Scale>(double_values, resultCollections, multiply_parallel_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, multiply_parallel_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, multiply_parallel_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, multiply_parallel_Fukushima);
    calc_perf<floatingExp2Integer::Float64LargeRangeNumber>(double_values, resultCollections, multiply_parallel_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, multiply_parallel_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, multiply_parallel_Float64PosExp2Int64);


    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    std::cout << std::endl << "Results" << std::endl << "n ";
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

    std::cout << std::endl << "Calculation time" << std::endl << "n ";
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

    std::cout << "END" << std::endl;

    return 0;
}
