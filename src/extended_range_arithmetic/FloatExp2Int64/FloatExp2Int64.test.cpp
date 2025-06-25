#include <iostream>
#include "FloatExp2Int64.h"
#include <gtest/gtest.h>

namespace {
    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        extended_range_arithmetic::FloatExp2Int64 num {1.0};

        EXPECT_EQ(num.principal_part(), 1.0);
        EXPECT_EQ(num.auxiliary_index(), 0);
        EXPECT_EQ(num.as_double(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        extended_range_arithmetic::FloatExp2Int64 num {2.0};

        EXPECT_EQ(num.principal_part(), 1.0);
        EXPECT_EQ(num.auxiliary_index(), 1);
        EXPECT_EQ(num.as_double(), 2.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
        extended_range_arithmetic::FloatExp2Int64 num1 {0.5};
        extended_range_arithmetic::FloatExp2Int64 num2 {2.0};

        EXPECT_EQ(num1.principal_part(), 1.0);
        EXPECT_EQ(num1.auxiliary_index(), -1);
        EXPECT_EQ(num1.as_double(), 0.5);

        EXPECT_EQ(num2.principal_part(), 1.0);
        EXPECT_EQ(num2.auxiliary_index(), 1);
        EXPECT_EQ(num2.as_double(), 2.0);

        extended_range_arithmetic::FloatExp2Int64 result {num1 + num2};

        EXPECT_EQ(result.principal_part(), 1.25);
        EXPECT_EQ(result.auxiliary_index(), 1);
        EXPECT_EQ(result.as_double(), 2.5);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
        extended_range_arithmetic::FloatExp2Int64 num1 {6.345634};
        extended_range_arithmetic::FloatExp2Int64 num2 {27.2728};

        EXPECT_EQ(num1.principal_part(), 1.5864085);
        EXPECT_EQ(num1.auxiliary_index(), 2);
        EXPECT_EQ(num1.as_double(), 6.345634);

        EXPECT_EQ(num2.principal_part(), 1.70455);
        EXPECT_EQ(num2.auxiliary_index(), 4);
        EXPECT_EQ(num2.as_double(), 27.2728);

        extended_range_arithmetic::FloatExp2Int64 result {num1 + num2};

        EXPECT_EQ(result.principal_part(), 1.0505760625);
        EXPECT_EQ(result.auxiliary_index(), 5);
        EXPECT_EQ(result.as_double(), 33.618434);
    }
}

