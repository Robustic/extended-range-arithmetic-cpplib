#include <iostream>
#include "DoubleExp2Int.h"
#include <gtest/gtest.h>

namespace {
    TEST(DoubleExp2Int, ConstructorWorksWith_0) {
        floatingExp2Integer::DoubleExp2Int test0 {0.0};

        std::cout << test0.mantissa() << std::endl;

        EXPECT_EQ(test0.mantissa(), 0.0);
        EXPECT_EQ(test0.exponent(), 0);
        EXPECT_EQ(test0.asDouble(), 0.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_1) {
        floatingExp2Integer::DoubleExp2Int test1 {1.0};

        std::cout << test1.mantissa() << std::endl;

        EXPECT_EQ(test1.mantissa(), 1.0);
        EXPECT_EQ(test1.exponent(), 0);
        EXPECT_EQ(test1.asDouble(), 1.0);
    }

    TEST(DoubleExp2Int, ConstructorWorksWith_2) {
        floatingExp2Integer::DoubleExp2Int test2 {2.0};

        std::cout << test2.mantissa() << std::endl;

        EXPECT_EQ(test2.mantissa(), 1.0);
        EXPECT_EQ(test2.exponent(), 1);
        EXPECT_EQ(test2.asDouble(), 2.0);
    }
}
