#include <iostream>
#include "DoubleExp2Int.h"

int main() {
    union {
        unsigned long long num;
        double fp;
    } pun;
    
    floatingExp2Integer::DoubleExp2Int test0 {0};
    std::cout << std::dec << test0.mantissa() << " " << test0.exponent() << std::endl;
    std::cout << std::dec << test0.asDouble() << std::endl;

    floatingExp2Integer::DoubleExp2Int test1 {1};
    std::cout << std::dec << test1.mantissa() << " " << test1.exponent() << std::endl;
    std::cout << std::dec << test1.asDouble() << std::endl;
    
    floatingExp2Integer::DoubleExp2Int test2 {2};
    std::cout << std::dec << test2.mantissa() << " " << test2.exponent() << std::endl;
    std::cout << std::dec << test2.asDouble() << std::endl;
}
