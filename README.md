# extended-range-arithmetic-cpplib

## License

[extended-range-arithmetic-cpplib](https://github.com/Robustic/extended-range-arithmetic-cpplib) © 1999 by [Juha Malinen](https://github.com/Robustic) is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)<img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" style="max-width: 1.5em; max-height: 1.5em; margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" style="max-width: 1.5em; max-height: 1.5em; margin-left: .2em;">


## Automatic tests

Run the following commands in the cloned repository root:

```
cmake --preset debug
// some output
cmake --preset release
// some output
cmake --build build_debug
// some output
cmake --build build_release
// some output
./build_debug/FloatExp2Int64.test
Running main() from /extended-range-arithmetic-cpplib/build_debug/_deps/googletest-src/googletest/src/gtest_main.cc
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 4 tests from DoubleExp2Int
[ RUN      ] DoubleExp2Int.ConstructorWorksWith_1
[       OK ] DoubleExp2Int.ConstructorWorksWith_1 (0 ms)
[ RUN      ] DoubleExp2Int.ConstructorWorksWith_2
[       OK ] DoubleExp2Int.ConstructorWorksWith_2 (0 ms)
[ RUN      ] DoubleExp2Int.OperatorPlusWorksWith_0_5and2
[       OK ] DoubleExp2Int.OperatorPlusWorksWith_0_5and2 (0 ms)
[ RUN      ] DoubleExp2Int.OperatorPlusWorksWith_6_345634and27_2728
[       OK ] DoubleExp2Int.OperatorPlusWorksWith_6_345634and27_2728 (0 ms)
[----------] 4 tests from DoubleExp2Int (16 ms total)

[----------] Global test environment tear-down
[==========] 4 tests from 1 test suite ran. (25 ms total)
[  PASSED  ] 4 tests.
```

In Windows 11 you can give high priority with:

```
start "build_release/PerformanceTest-ArraySize.exe" /high
 ```

