# %%
import time
from concurrent.futures import ProcessPoolExecutor

num_workers = 4  # Adjust based on your system's cores and resources

# Define a CPU-intensive task
def sum_of_squares(n):
    return sum(i * i for i in range(n))

# Single-process execution
def single_process_demo(task_inputs):
    results = []
    for n in task_inputs:
        results.append(sum_of_squares(n))
    return results

# Multi-process execution using ProcessPoolExecutor
def multi_process_demo(task_inputs, num_workers):
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(sum_of_squares, task_inputs))
    return results

# Main function to compare performance
def main():
    # Task inputs: list of large numbers
    task_inputs = [10**7, 10**7, 10**7, 10**7]

    # Determine the number of workers (logical cores)
    import os
    num_workers = os.cpu_count()
    num_workers = 1
    print(f"Number of workers: {num_workers}")

    # Single-process execution
    print("Running single-process execution...")
    start_time = time.time()
    single_process_demo(task_inputs)
    single_duration = time.time() - start_time
    print(f"Single-process duration: {single_duration:.2f} seconds")

    # Multi-process execution
    print("Running multi-process execution...")
    start_time = time.time()
    multi_process_demo(task_inputs, num_workers)
    multi_duration = time.time() - start_time
    print(f"Multi-process duration: {multi_duration:.2f} seconds")

    # Compare results
    print("\nPerformance comparison:")
    print(f"Single-process time: {single_duration:.2f} seconds")
    print(f"Multi-process time: {multi_duration:.2f} seconds")
    print(f"Speedup: {single_duration / multi_duration:.2f}x")

if __name__ == "__main__":
    main()
# %%

from multiprocessing import Pool
import time

# Define a CPU-intensive task
def sum_of_squares(n):
    return sum(i * i for i in range(n))

def main():
    task_inputs = [10**7, 10**7, 10**7, 10**7]
    
    # Single-process execution
    print("Running single-process execution...")
    start_time = time.time()
    results = [sum_of_squares(n) for n in task_inputs]
    single_duration = time.time() - start_time
    print(f"Single-process duration: {single_duration:.2f} seconds")

    # Multi-process execution
    print("Running multi-process execution...")
    start_time = time.time()
    with Pool(processes=4) as pool:  # Use 4 worker processes
        results = pool.map(sum_of_squares, task_inputs)
    multi_duration = time.time() - start_time
    print(f"Multi-process duration: {multi_duration:.2f} seconds")

    print(f"Speedup: {single_duration / multi_duration:.2f}x")

if __name__ == "__main__":
    main()
# %%
from joblib import Parallel, delayed
import time

# Define a CPU-intensive task
def sum_of_squares(n):
    return sum(i * i for i in range(n))

def main():
    task_inputs = [10**7, 10**7, 10**7, 10**7]
    
    # Single-process execution
    print("Running single-process execution...")
    start_time = time.time()
    results = [sum_of_squares(n) for n in task_inputs]
    single_duration = time.time() - start_time
    print(f"Single-process duration: {single_duration:.2f} seconds")

    # Multi-process execution
    print("Running multi-process execution...")
    start_time = time.time()
    results = Parallel(n_jobs=10)(delayed(sum_of_squares)(n) for n in task_inputs)
    multi_duration = time.time() - start_time
    print(f"Multi-process duration: {multi_duration:.2f} seconds")

    print(f"Speedup: {single_duration / multi_duration:.2f}x")

if __name__ == "__main__":
    main()
# %%
# %%

import ipyparallel as ipp
import time

# Create a client to connect to the ipyparallel cluster
client = ipp.Client()
dview = client[:]

# Define a function to compute the sum of squares for a range of numbers
def compute_sum_of_squares(n):
    return sum(i**2 for i in range(n))

# Test parameters
n_values = [10**6] * 20  # List of 20 tasks, each computing sum of squares for 1 million numbers

# Serial computation
start_serial = time.time()
serial_results = [compute_sum_of_squares(n) for n in n_values]
end_serial = time.time()
print(f"Serial computation time: {end_serial - start_serial:.2f} seconds")

# Parallel computation
start_parallel = time.time()
parallel_results = dview.map_sync(compute_sum_of_squares, n_values)
end_parallel = time.time()
print(f"Parallel computation time: {end_parallel - start_parallel:.2f} seconds")

# Verify the results are the same
assert serial_results == parallel_results, "Results don't match!"

# Print speedup
speedup = (end_serial - start_serial) / (end_parallel - start_parallel)
print(f"Speedup: {speedup:.2f}x")
# %%
