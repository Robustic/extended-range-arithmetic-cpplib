#include <iostream>
#include "DoubleExp2Int.h"
#include <gtest/gtest.h>

namespace {
    TEST(DoubleExp2Int, ConstructorWorksWith_0) {
        floatingExp2Integer::DoubleExp2Int num {0.0};

        EXPECT_EQ(num.mantissa(), 0.0);
        EXPECT_EQ(num.exponent(), 0);
        EXPECT_EQ(num.asDouble(), 0.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        floatingExp2Integer::DoubleExp2Int num {1.0};

        EXPECT_EQ(num.mantissa(), 1.0);
        EXPECT_EQ(num.exponent(), 0);
        EXPECT_EQ(num.asDouble(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        floatingExp2Integer::DoubleExp2Int num {2.0};

        EXPECT_EQ(num.mantissa(), 1.0);
        EXPECT_EQ(num.exponent(), 1);
        EXPECT_EQ(num.asDouble(), 2.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
        floatingExp2Integer::DoubleExp2Int num1 {0.5};
        floatingExp2Integer::DoubleExp2Int num2 {2.0};

        EXPECT_EQ(num1.mantissa(), 1.0);
        EXPECT_EQ(num1.exponent(), -1);
        EXPECT_EQ(num1.asDouble(), 0.5);

        EXPECT_EQ(num2.mantissa(), 1.0);
        EXPECT_EQ(num2.exponent(), 1);
        EXPECT_EQ(num2.asDouble(), 2.0);

        floatingExp2Integer::DoubleExp2Int result {num1 + num2};

        EXPECT_EQ(result.mantissa(), 1.25);
        EXPECT_EQ(result.exponent(), 1);
        EXPECT_EQ(result.asDouble(), 2.5);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
        floatingExp2Integer::DoubleExp2Int num1 {6.345634};
        floatingExp2Integer::DoubleExp2Int num2 {27.2728};

        // EXPECT_EQ(num1.mantissa(), 1.0);
        // EXPECT_EQ(num1.exponent(), -1);
        EXPECT_EQ(num1.asDouble(), 6.345634);

        // EXPECT_EQ(num2.mantissa(), 1.0);
        // EXPECT_EQ(num2.exponent(), 1);
        EXPECT_EQ(num2.asDouble(), 27.2728);

        floatingExp2Integer::DoubleExp2Int result {num1 + num2};

        // EXPECT_EQ(result.mantissa(), 1.25);
        // EXPECT_EQ(result.exponent(), 1);
        EXPECT_EQ(result.asDouble(), 33.618434);
    }
}
