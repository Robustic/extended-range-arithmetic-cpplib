// #include <iostream>
// #include <cmath>
// #include "Float64ExtendedExp.h"
// #include <gtest/gtest.h>

// namespace {
//     double roundToDecimal(double value, int decimalPlaces) {
//         double factor = std::pow(10, decimalPlaces);
//         return std::round(value * factor) / factor;
//     }

//     TEST(DoubleExp2Int, ConstructorWorksWith_1) {
//         floatingExp2Integer::Float64ExtendedExp num {1.0};

//         EXPECT_EQ(num.sicnificand(), 0.5);
//         EXPECT_EQ(num.exponent(), 1);
//         EXPECT_EQ(num.asDouble(), 1.0);
//     }

//     TEST(DoubleExp2Int, ConstructorWorksWith_2) {
//         floatingExp2Integer::Float64ExtendedExp num {2.0};

//         EXPECT_EQ(num.sicnificand(), 0.5);
//         EXPECT_EQ(num.exponent(), 2);
//         EXPECT_EQ(num.asDouble(), 2.0);
//     }

//     TEST(DoubleExp2Int, OperatorPlusWorksWith_0_5and2) {
//         floatingExp2Integer::Float64ExtendedExp num1 {0.5};
//         floatingExp2Integer::Float64ExtendedExp num2 {2.0};

//         EXPECT_EQ(num1.sicnificand(), 0.5);
//         EXPECT_EQ(num1.exponent(), 0);
//         EXPECT_EQ(num1.asDouble(), 0.5);

//         EXPECT_EQ(num2.sicnificand(), 0.5);
//         EXPECT_EQ(num2.exponent(), 2);
//         EXPECT_EQ(num2.asDouble(), 2.0);

//         floatingExp2Integer::Float64ExtendedExp result {num1 + num2};

//         EXPECT_EQ(result.sicnificand(), 0.625);
//         EXPECT_EQ(result.exponent(), 2);
//         EXPECT_EQ(result.asDouble(), 2.5);
//     }

//     TEST(DoubleExp2Int, OperatorPlusWorksWith_6_345634and27_2728) {
//         floatingExp2Integer::Float64ExtendedExp num1 {6.345634};
//         floatingExp2Integer::Float64ExtendedExp num2 {27.2728};

//         EXPECT_EQ(roundToDecimal(num1.sicnificand(), 8), 0.79320425);
//         EXPECT_EQ(num1.exponent(), 3);
//         EXPECT_EQ(num1.asDouble(), 6.345634);

//         EXPECT_EQ(roundToDecimal(num2.sicnificand(), 8), 0.852275);
//         EXPECT_EQ(num2.exponent(), 5);
//         EXPECT_EQ(roundToDecimal(num2.asDouble(), 4), 27.2728);

//         floatingExp2Integer::Float64ExtendedExp result {num1 + num2};

//         EXPECT_EQ(roundToDecimal(result.sicnificand(), 11), 0.52528803125);
//         EXPECT_EQ(result.exponent(), 6);
//         EXPECT_EQ(roundToDecimal(result.asDouble(), 8), 33.618434);
//     }
// }

