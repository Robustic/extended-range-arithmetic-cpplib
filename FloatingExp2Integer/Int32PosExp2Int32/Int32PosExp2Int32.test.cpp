#include <cmath>
#include <iostream>
#include "Int32PosExp2Int32.h"
#include <gtest/gtest.h>

namespace {
    double round_to(double value, double precision = 1.0)
    {
        return std::round(value / precision) * precision;
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        floatingExp2Integer::Int32PosExp2Int32 num {1.0};

        EXPECT_EQ(num.sicnificand(), std::exp2(23));
        EXPECT_EQ(num.exponent(), -23);
        EXPECT_EQ(num.asFloat(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        floatingExp2Integer::Int32PosExp2Int32 num {2.0};

        EXPECT_EQ(num.sicnificand(), std::exp2(23));
        EXPECT_EQ(num.exponent(), -22);
        EXPECT_EQ(num.asFloat(), 2.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
        floatingExp2Integer::Int32PosExp2Int32 num1 {0.5};
        floatingExp2Integer::Int32PosExp2Int32 num2 {2.0};

        EXPECT_EQ(num1.sicnificand(), std::exp2(23));
        EXPECT_EQ(num1.exponent(), -24);
        EXPECT_EQ(num1.asFloat(), 0.5);

        EXPECT_EQ(num2.sicnificand(), std::exp2(23));
        EXPECT_EQ(num2.exponent(), -22);
        EXPECT_EQ(num2.asFloat(), 2.0);

        floatingExp2Integer::Int32PosExp2Int32 result {num1 + num2};

        EXPECT_EQ(result.sicnificand(), std::exp2(23) + std::exp2(21));
        EXPECT_EQ(result.exponent(), -22);
        EXPECT_EQ(result.asFloat(), 2.5);
    }

    TEST(DoubleExp2Int, OperatorMultiplyWorksWith_0_5and2) {
        floatingExp2Integer::Int32PosExp2Int32 num1 {0.5};
        floatingExp2Integer::Int32PosExp2Int32 num2 {2.0};

        EXPECT_EQ(num1.sicnificand(), std::exp2(23));
        EXPECT_EQ(num1.exponent(), -24);
        EXPECT_EQ(num1.asFloat(), 0.5);

        EXPECT_EQ(num2.sicnificand(), std::exp2(23));
        EXPECT_EQ(num2.exponent(), -22);
        EXPECT_EQ(num2.asFloat(), 2.0);

        floatingExp2Integer::Int32PosExp2Int32 result {num1 * num2};

        EXPECT_EQ(result.sicnificand(), std::exp2(23));
        EXPECT_EQ(result.exponent(), -23);
        EXPECT_EQ(result.asFloat(), 1.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
        floatingExp2Integer::Int32PosExp2Int32 num1 {6.345634};
        floatingExp2Integer::Int32PosExp2Int32 num2 {27.2728};

        EXPECT_EQ(num1.sicnificand(), (float)1.5864085 * std::exp2(23));
        EXPECT_EQ(num1.exponent(), -21);
        EXPECT_EQ(num1.asFloat(), (float)6.345634);

        EXPECT_EQ(num2.sicnificand(), (float)1.70455 * std::exp2(23));
        EXPECT_EQ(num2.exponent(), -19);
        EXPECT_EQ(num2.asFloat(), (float)27.2728);

        floatingExp2Integer::Int32PosExp2Int32 result {num1 + num2};

        EXPECT_EQ(round_to(result.sicnificand(), 10), round_to((float)1.0505760625 * std::exp2(23), 10));
        EXPECT_EQ(result.exponent(), -18);
        EXPECT_EQ(round_to(result.asFloat(), 0.00001), round_to((float)33.618434, 0.00001));
    }
}
