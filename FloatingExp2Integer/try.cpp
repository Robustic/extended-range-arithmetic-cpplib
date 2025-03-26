#include <iostream>
#include <vector>
#include <immintrin.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>
// #include "exp_table.h"
#include "Timer.h"
// #include "Dbl.h"
// #include "Dbl2.h"
// #include "Dbl3.h"
// #include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
// #include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
// #include "./Float64Exp2Int64/Float64Exp2Int64.h"
#include "./Float64ExtendedExp/Float64ExtendedExp.h"
#include "./Fukushima/Fukushima.h"

void print_Float64ExtendedExp(double number) {
    std::cout << "Test: " << number << std::endl;
    floatingExp2Integer::Float64ExtendedExp nmbr(number);
    // std::cout << "nmbr.sicnificand(): " << std::to_string(nmbr.sicnificand()) << std::endl;
    // std::cout << "nmbr.exponent(): " << std::to_string(nmbr.exponent()) << std::endl;
    // std::cout << "nmbr.asDouble(): " << std::to_string(nmbr.asDouble()) << std::endl;
    std::cout << std::endl;
}

void print_Fukushima(double number) {
    std::cout << "Test: " << number << std::endl;
    floatingExp2Integer::Fukushima nmbr(number);
    std::cout << "nmbr.sicnificand(): " << nmbr.sicnificand() << std::endl;
    std::cout << "nmbr.exponent(): " << nmbr.exponent() << std::endl;
    std::cout << "nmbr.asDouble(): " << nmbr.asDouble() << std::endl;
    std::cout << std::endl;
}

int main() {
    const unsigned int n = 16;
    double input[n] = { 0.25, 0.375, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 3.0, 2.0, 1.5, 1.0, 0.75, 0.5, 0.375, 0.25 };
    double encoded[n];
    floatingExp2Integer::Float64ExtendedExp extendedExp(0.2625);
    extendedExp.doubleToFloat64ExtendedExp(n, input, encoded);

    std::cout << "input[0] " << input[0] << std::endl;
    std::cout << "encoded[0] " << encoded[0] << std::endl;

    double sum = extendedExp.sumFloat64ExtendedExp(n, encoded);
    std::cout << "sum " << extendedExp.decodeFloat64ExtendedExp(sum) << std::endl;


    //print_Float64ExtendedExp(0.25);
    //print_Float64ExtendedExp(0.375);
    //print_Float64ExtendedExp(0.5);
    //print_Float64ExtendedExp(0.75);
    //print_Float64ExtendedExp(1.0);
    //print_Float64ExtendedExp(1.5);
    //print_Float64ExtendedExp(2.0);
    //print_Float64ExtendedExp(3.0);
    //print_Float64ExtendedExp(4.0);
    //print_Float64ExtendedExp(4.5);
    // print_Float64ExtendedExp(2.0/3.0);
    // print_Float64ExtendedExp(4503599627370494);
    // print_Float64ExtendedExp(4503599627370495);
    // print_Float64ExtendedExp(4503599627370496);

    //floatingExp2Integer::Float64ExtendedExp nmbr1(4.5);
    //floatingExp2Integer::Float64ExtendedExp nmbr2(2.0);
    //nmbr1 += nmbr2;
    //std::cout << "nmbr.sicnificand(): " << nmbr1.sicnificand() << std::endl;
    //std::cout << "nmbr.exponent(): " << nmbr1.exponent() << std::endl;
    //std::cout << "nmbr.asDouble(): " << nmbr1.asDouble() << std::endl;



    // floatingExp2Integer::Float64ExtendedExp nmbr3(0.2625);
    // floatingExp2Integer::Float64ExtendedExp nmbr4(0.75);
    // nmbr3 += nmbr4;
    // std::cout << "nmbr.sicnificand(): " << nmbr3.sicnificand() << std::endl;
    // std::cout << "nmbr.exponent(): " << nmbr3.exponent() << std::endl;
    // std::cout << "nmbr.asDouble(): " << nmbr3.asDouble() << std::endl;

    // print_Fukushima(0.25);
    // print_Fukushima(1.0);
    // print_Fukushima(2.0/3.0);
    // print_Fukushima(4503599627370494);
    // print_Fukushima(4503599627370495);
    // print_Fukushima(4503599627370496);

    //floatingExp2Integer::Fukushima nmbr3(0.08044567);
    //floatingExp2Integer::Fukushima nmbr4(0.3456388);
    //nmbr3 += nmbr4;
    //std::cout << "nmbr.sicnificand(): " << nmbr3.sicnificand() << std::endl;
    //std::cout << "nmbr.exponent(): " << nmbr3.exponent() << std::endl;
    //std::cout << "nmbr.asDouble(): " << nmbr3.asDouble() << std::endl;
}
 


