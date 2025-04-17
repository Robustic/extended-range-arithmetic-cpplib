import time
import struct
import math
import torch

# Read the binary file
def read_binary_file(filename):
    with open(filename, 'rb') as file:
        data = file.read()
        num_doubles = len(data) // struct.calcsize('d')
        doubles = struct.unpack('d' * num_doubles, data)
        return list(doubles)
    
def ln_to_log2(ln_value):
    return ln_value / math.log(2)

# Calculate log-sum-exp using torch
def calculate_log_sum_exp(doubles, repeat_count):
    tensor = torch.tensor(doubles, dtype=torch.float64)
    
    start_time = time.perf_counter()
    for _ in range(repeat_count):
        log_sum_exp = torch.logsumexp(tensor, dim=0)
    end_time = time.perf_counter()
    
    execution_time_nanoseconds = (end_time - start_time) * 1_000_000_000 / repeat_count;

    print("Count:", tensor.numel())
    print("Log-Sum-Exp:", ln_to_log2(log_sum_exp.item()))
    print(f"Execution time: {execution_time_nanoseconds} nanoseconds")

    return log_sum_exp

def log2_to_ln(log2_value):
    return log2_value * math.log(2)

# Read doubles from binary file
doubles_log2 = read_binary_file('log2range_-100_to_0.bin')
doubles = [log2_to_ln(value) for value in doubles_log2]

print("torch.get_num_threads():", torch.get_num_threads())

for _ in range(1):
    result = calculate_log_sum_exp(doubles[:1000], 10000)
    result = calculate_log_sum_exp(doubles[:3000], 10000)
    result = calculate_log_sum_exp(doubles[:10000], 10000)
    result = calculate_log_sum_exp(doubles[:30000], 3000)
    result = calculate_log_sum_exp(doubles[:100000], 1000)
    result = calculate_log_sum_exp(doubles[:300000], 300)
    result = calculate_log_sum_exp(doubles[:1000000], 100)
    result = calculate_log_sum_exp(doubles[:3000000], 30)
    result = calculate_log_sum_exp(doubles[:10000000], 10)
    result = calculate_log_sum_exp(doubles[:30000000], 3)
    result = calculate_log_sum_exp(doubles, 1)

# C:\msys64\home\FloatingExp2Integer\Python>python logsumexp.py
# torch.get_num_threads(): 8
# Count: 1000
# Log-Sum-Exp: -21.3121011101092
# Execution time: 1013099.9726243317 nanoseconds
# Count: 3000
# Log-Sum-Exp: -19.708301267289112
# Execution time: 1056200.0097706914 nanoseconds
# Count: 10000
# Log-Sum-Exp: -17.964701494069118
# Execution time: 270899.96729046106 nanoseconds
# Count: 30000
# Log-Sum-Exp: -16.380524377550557
# Execution time: 238600.01238062978 nanoseconds
# Count: 100000
# Log-Sum-Exp: -14.63813971690725
# Execution time: 336299.99961704016 nanoseconds
# Count: 300000
# Log-Sum-Exp: -13.055326874980866
# Execution time: 750200.0080421567 nanoseconds
# Count: 1000000
# Log-Sum-Exp: -11.317468364916614
# Execution time: 1617599.9771803617 nanoseconds
# Count: 3000000
# Log-Sum-Exp: -9.732198940653394
# Execution time: 6065300.025511533 nanoseconds
# Count: 10000000
# Log-Sum-Exp: -7.995498936018165
# Execution time: 21973999.973852187 nanoseconds
# Count: 30000000
# Log-Sum-Exp: -6.4105928230173745
# Execution time: 61859800.01464486 nanoseconds
# Count: 100000000
# Log-Sum-Exp: -4.673451135115633
# Execution time: 190469200.02438128 nanoseconds

# C:\msys64\home\FloatingExp2Integer\Python>
