echo "Running myProgram.exe with 0"
./build_release/Float64.try.exe -1000 0

echo "Running myProgram.exe with -1 000"
./build_release/Float64.try.exe -2000 -1000

echo "Running myProgram.exe with -2 000"
./build_release/Float64.try.exe -3000 -2000

echo "Running myProgram.exe with -5 000"
./build_release/Float64.try.exe -6000 -5000

echo "Running myProgram.exe with -10 000"
./build_release/Float64.try.exe -11000 -10000

echo "Running myProgram.exe with -20 000"
./build_release/Float64.try.exe -21000 -20000

echo "Running myProgram.exe with -50 000"
./build_release/Float64.try.exe -51000 -50000

echo "Running myProgram.exe with -100 000"
./build_release/Float64.try.exe -101000 -100000

echo "Running myProgram.exe with -200 000"
./build_release/Float64.try.exe -201000 -200000

echo "Running myProgram.exe with -500 000"
./build_release/Float64.try.exe -501000 -500000

echo "Running myProgram.exe with -1 000 000"
./build_release/Float64.try.exe -1001000 -1000000

echo All tasks completed.
