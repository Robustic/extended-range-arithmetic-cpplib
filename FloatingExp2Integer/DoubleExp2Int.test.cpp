#include "DoubleExp2Int.h"
#include <gtest/gtest.h>

namespace {
    TEST(DoubleExp2Int, ConstructionWorks) {
        floatingExp2Integer::DoubleExp2Int d2i {3.14};

        EXPECT_EQ(d2i.mantissa(), 3.14);
        EXPECT_EQ(d2i.exponent(), 0);
    }
}
