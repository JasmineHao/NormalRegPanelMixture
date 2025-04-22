# %%
from MixtureTestFunctions import *
import time
from numba import set_num_threads, get_num_threads

set_num_threads(12)

# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
alphaset = [np.array([1/3,1/3,1/3])]
muset = [np.array([-0.5, 0 , 1.5]), np.array([-1.5, 0 , 1.5])]
sigmaset = [np.array([1.,1.,1.])]
beta = np.array([[1.0],[1.0],[1.0]])
gamma = np.array([0.0])
alpha_bounds = [0.01,0.05,0.1]
Nset = [200, 400]
# Test panel mixture
N = 200
T = 3
M = 3
p = 0
q = 0
nrep = 100
BB = 199

alpha = alphaset[0]
mu = muset[0]
sigma = sigmaset[0]
n_grid=3
r_test=2
# Generate data
weights_equal = np.full(N, 1 / N)


 # Determine the total number of parameter combinations
total_combinations = len(Nset) * len(muset) * len(alpha_bounds)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 7))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0

# Loop over parameters
for N in Nset:
    for mu in muset:  
        Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]
            # Nonparametric test 
        data_nopar = [data[0] for data in Data]
        for alpha_bound in alpha_bounds:           
            start_time = time.time()
            # Nonparametric test
            result_rk_each = NonParTestParallel(data_nopar, N, T, 2, p, q, nrep, n_grid, BB, r_test)
            result_lr_each = LRTestParallel(Data, N, T, 2, p, q, nrep, BB = BB, alpha_bound=alpha_bound)
                        
            # Compute the mean results across replications and store them
            simulation_result_matrix[count, :2] = result_rk_each.mean(axis=0)
            simulation_result_matrix[count, 2:7] = result_lr_each.mean(axis=0)[:5]
            
             # Print execution time for this parameter combination
            print(f"Execution time for N={N}, mu={mu}: {time.time() - start_time:.2f} seconds")
            print(simulation_result_matrix)
            parameter_combinations.append((N, mu, alpha_bound))
            # Increment the counter
            count += 1
# %%
            
# Print the final result matrix
print("Simulation Result Matrix:")
print(simulation_result_matrix)


result_freq_table = pd.DataFrame(100 * simulation_result_matrix, index=parameter_combinations, columns = ['rk_mean','rk_max', "lr 10%", "lr 5%", "lr 1%", "aic", "bic" ])

# Save to CSV
result_freq_table.to_csv("test_power_table2.csv")

# %%
