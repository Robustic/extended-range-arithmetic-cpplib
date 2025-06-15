#include <cmath>
#include <iostream>
#include "IntExp2Int64.h"
#include <gtest/gtest.h>

namespace {
    double round_to(double value, double precision = 1.0)
    {
        return std::round(value / precision) * precision;
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        extended_range_arithmetic::IntExp2Int64 num {1.0};

        EXPECT_EQ(num.sicnificand(), std::exp2(52));
        EXPECT_EQ(num.exponent(), -52);
        EXPECT_EQ(num.as_double(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        extended_range_arithmetic::IntExp2Int64 num {2.0};

        EXPECT_EQ(num.sicnificand(), std::exp2(52));
        EXPECT_EQ(num.exponent(), -51);
        EXPECT_EQ(num.as_double(), 2.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
        extended_range_arithmetic::IntExp2Int64 num1 {0.5};
        extended_range_arithmetic::IntExp2Int64 num2 {2.0};

        EXPECT_EQ(num1.sicnificand(), std::exp2(52));
        EXPECT_EQ(num1.exponent(), -53);
        EXPECT_EQ(num1.as_double(), 0.5);

        EXPECT_EQ(num2.sicnificand(), std::exp2(52));
        EXPECT_EQ(num2.exponent(), -51);
        EXPECT_EQ(num2.as_double(), 2.0);

        extended_range_arithmetic::IntExp2Int64 result {num1 + num2};

        EXPECT_EQ(result.sicnificand(), std::exp2(52) + std::exp2(50));
        EXPECT_EQ(result.exponent(), -51);
        EXPECT_EQ(result.as_double(), 2.5);
    }

    TEST(DoubleExp2Int, OperatorMultiplyWorksWith_0_5and2) {
        extended_range_arithmetic::IntExp2Int64 num1 {0.5};
        extended_range_arithmetic::IntExp2Int64 num2 {2.0};

        EXPECT_EQ(num1.sicnificand(), std::exp2(52));
        EXPECT_EQ(num1.exponent(), -53);
        EXPECT_EQ(num1.as_double(), 0.5);

        EXPECT_EQ(num2.sicnificand(), std::exp2(52));
        EXPECT_EQ(num2.exponent(), -51);
        EXPECT_EQ(num2.as_double(), 2.0);

        extended_range_arithmetic::IntExp2Int64 result {num1 * num2};

        EXPECT_EQ(result.sicnificand(), std::exp2(52));
        EXPECT_EQ(result.exponent(), -52);
        EXPECT_EQ(result.as_double(), 1.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
        extended_range_arithmetic::IntExp2Int64 num1 {6.345634};
        extended_range_arithmetic::IntExp2Int64 num2 {27.2728};

        EXPECT_EQ(num1.sicnificand(), 1.5864085 * std::exp2(52));
        EXPECT_EQ(num1.exponent(), -50);
        EXPECT_EQ(num1.as_double(), 6.345634);

        EXPECT_EQ(num2.sicnificand(), 1.70455 * std::exp2(52));
        EXPECT_EQ(num2.exponent(), -48);
        EXPECT_EQ(num2.as_double(), 27.2728);

        extended_range_arithmetic::IntExp2Int64 result {num1 + num2};

        EXPECT_EQ(round_to(result.sicnificand(), 10), round_to(1.0505760625 * std::exp2(52), 10));
        EXPECT_EQ(result.exponent(), -47);
        EXPECT_EQ(round_to(result.as_double(), 0.00001), round_to(33.618434, 0.00001));
    }
}
