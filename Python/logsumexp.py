import time
import struct
import torch

# Read the binary file
def read_binary_file(filename):
    with open(filename, 'rb') as file:
        data = file.read()
        num_doubles = len(data) // struct.calcsize('d')
        doubles = struct.unpack('d' * num_doubles, data)
        return list(doubles)

# Calculate log-sum-exp using torch
def calculate_log_sum_exp(doubles):
    tensor = torch.tensor(doubles, dtype=torch.float64)
    
    start_time = time.perf_counter()
    log_sum_exp = torch.logsumexp(tensor, dim=0)
    end_time = time.perf_counter()
    
    execution_time_nanoseconds = (end_time - start_time) * 1_000_000_000

    print("Count:", tensor.numel())
    print("Log-Sum-Exp:", log_sum_exp.item())
    print(f"Execution time: {execution_time_nanoseconds} nanoseconds")

    return log_sum_exp

# Read doubles from binary file
doubles = read_binary_file('data.bin')

result = calculate_log_sum_exp(doubles[:1000])
result = calculate_log_sum_exp(doubles[:3000])
result = calculate_log_sum_exp(doubles[:10000])
result = calculate_log_sum_exp(doubles[:30000])
result = calculate_log_sum_exp(doubles[:100000])
result = calculate_log_sum_exp(doubles[:300000])
result = calculate_log_sum_exp(doubles[:1000000])
result = calculate_log_sum_exp(doubles[:3000000])
result = calculate_log_sum_exp(doubles[:10000000])
result = calculate_log_sum_exp(doubles[:30000000])
result = calculate_log_sum_exp(doubles)

# C:\msys64\home\FloatingExp2Integer\Python>python logsumexp.py
# Count: 1000
# Log-Sum-Exp: -24.277303072587276
# Execution time: 1116200.0009790063 nanoseconds
# Count: 3000
# Log-Sum-Exp: -23.155393202897304
# Execution time: 1283499.994315207 nanoseconds
# Count: 10000
# Log-Sum-Exp: -21.943974860976112
# Execution time: 133500.02700462937 nanoseconds
# Count: 30000
# Log-Sum-Exp: -20.84549363731455
# Execution time: 217799.9704144895 nanoseconds
# Count: 100000
# Log-Sum-Exp: -19.636003929565195
# Execution time: 480799.9939657748 nanoseconds
# Count: 300000
# Log-Sum-Exp: -18.53908126562819
# Execution time: 703800.0039756298 nanoseconds
# Count: 1000000
# Log-Sum-Exp: -17.334269991359015
# Execution time: 1836600.0149399042 nanoseconds
# Count: 3000000
# Log-Sum-Exp: -16.235441172484787
# Execution time: 6711399.997584522 nanoseconds
# Count: 10000000
# Log-Sum-Exp: -15.031710085573613
# Execution time: 21436500.013805926 nanoseconds
# Count: 30000000
# Log-Sum-Exp: -13.9331482950458
# Execution time: 59958000.02710894 nanoseconds
# Count: 100000000
# Log-Sum-Exp: -12.729030034082715
# Execution time: 172586700.00592247 nanoseconds

# C:\msys64\home\FloatingExp2Integer\Python>
