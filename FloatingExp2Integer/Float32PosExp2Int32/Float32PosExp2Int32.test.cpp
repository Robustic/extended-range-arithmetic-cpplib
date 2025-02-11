#include <iostream>
#include "Float32PosExp2Int32.h"
#include <gtest/gtest.h>

namespace {
    // TEST(DoubleExp2Int, ConstructorWorksWith_0) {
    //     floatingExp2Integer::Float32PosExp2Int32 num {0.0};

    //     EXPECT_EQ(num.sicnificand(), 0.0);
    //     EXPECT_EQ(num.exponent(), 0);
    //     EXPECT_EQ(num.asFloat(), 0.0);
    // }

    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        floatingExp2Integer::Float32PosExp2Int32 num {1.0};

        EXPECT_EQ(num.sicnificand(), 1.0);
        EXPECT_EQ(num.exponent(), 0);
        EXPECT_EQ(num.asFloat(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        floatingExp2Integer::Float32PosExp2Int32 num {2.0};

        EXPECT_EQ(num.sicnificand(), 1.0);
        EXPECT_EQ(num.exponent(), 1);
        EXPECT_EQ(num.asFloat(), 2.0);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
        floatingExp2Integer::Float32PosExp2Int32 num1 {0.5};
        floatingExp2Integer::Float32PosExp2Int32 num2 {2.0};

        EXPECT_EQ(num1.sicnificand(), 1.0);
        EXPECT_EQ(num1.exponent(), -1);
        EXPECT_EQ(num1.asFloat(), 0.5);

        EXPECT_EQ(num2.sicnificand(), 1.0);
        EXPECT_EQ(num2.exponent(), 1);
        EXPECT_EQ(num2.asFloat(), 2.0);

        floatingExp2Integer::Float32PosExp2Int32 result {num1 + num2};

        EXPECT_EQ(result.sicnificand(), 1.25);
        EXPECT_EQ(result.exponent(), 1);
        EXPECT_EQ(result.asFloat(), 2.5);
    }

    TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
        floatingExp2Integer::Float32PosExp2Int32 num1 {6.345634};
        floatingExp2Integer::Float32PosExp2Int32 num2 {27.2728};

        EXPECT_EQ(num1.sicnificand(), (float)1.5864085);
        EXPECT_EQ(num1.exponent(), 2);
        EXPECT_EQ(num1.asFloat(), (float)6.345634);

        EXPECT_EQ(num2.sicnificand(), (float)1.70455);
        EXPECT_EQ(num2.exponent(), 4);
        EXPECT_EQ(num2.asFloat(), (float)27.2728);

        floatingExp2Integer::Float32PosExp2Int32 result {num1 + num2};

        EXPECT_EQ(result.sicnificand(), (float)1.0505760625);
        EXPECT_EQ(result.exponent(), 5);
        EXPECT_EQ(result.asFloat(), (float)33.618434);
    }
}
