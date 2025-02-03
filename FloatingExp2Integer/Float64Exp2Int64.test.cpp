#include <iostream>
#include "Float64Exp2Int64.h"
#include <gtest/gtest.h>

namespace {
    TEST(DoubleExp2Int, ConstructorWorksWith_0) {
        floatingExp2Integer::Float64Exp2Int64 num {0.0};

        EXPECT_EQ(num.mantissa(), 0.0);
        EXPECT_EQ(num.exponent(), 0);
        EXPECT_EQ(num.asDouble(), 0.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        floatingExp2Integer::Float64Exp2Int64 num {1.0};

        EXPECT_EQ(num.mantissa(), 1.0);
        EXPECT_EQ(num.exponent(), 0);
        EXPECT_EQ(num.asDouble(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        floatingExp2Integer::Float64Exp2Int64 num {2.0};

        EXPECT_EQ(num.mantissa(), 1.0);
        EXPECT_EQ(num.exponent(), 1);
        EXPECT_EQ(num.asDouble(), 2.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
        floatingExp2Integer::Float64Exp2Int64 num1 {0.5};
        floatingExp2Integer::Float64Exp2Int64 num2 {2.0};

        EXPECT_EQ(num1.mantissa(), 1.0);
        EXPECT_EQ(num1.exponent(), -1);
        EXPECT_EQ(num1.asDouble(), 0.5);

        EXPECT_EQ(num2.mantissa(), 1.0);
        EXPECT_EQ(num2.exponent(), 1);
        EXPECT_EQ(num2.asDouble(), 2.0);

        floatingExp2Integer::Float64Exp2Int64 result {num1 + num2};

        EXPECT_EQ(result.mantissa(), 1.25);
        EXPECT_EQ(result.exponent(), 1);
        EXPECT_EQ(result.asDouble(), 2.5);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
        floatingExp2Integer::Float64Exp2Int64 num1 {6.345634};
        floatingExp2Integer::Float64Exp2Int64 num2 {27.2728};

        EXPECT_EQ(num1.mantissa(), 1.5864085);
        EXPECT_EQ(num1.exponent(), 2);
        EXPECT_EQ(num1.asDouble(), 6.345634);

        EXPECT_EQ(num2.mantissa(), 1.70455);
        EXPECT_EQ(num2.exponent(), 4);
        EXPECT_EQ(num2.asDouble(), 27.2728);

        floatingExp2Integer::Float64Exp2Int64 result {num1 + num2};

        EXPECT_EQ(result.mantissa(), 1.0505760625);
        EXPECT_EQ(result.exponent(), 5);
        EXPECT_EQ(result.asDouble(), 33.618434);
    }
}
