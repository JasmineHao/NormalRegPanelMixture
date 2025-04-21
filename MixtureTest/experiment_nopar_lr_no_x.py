from MixtureTestFunctions import *
import time

# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
alphaset = [np.array([0.5, 0.5]), np.array([0.2, 0.8])]
muset = [np.array([-1.0, 1.0]), np.array([-0.5, 0.5])]
sigmaset = [np.array([0.8, 1.2])]
beta = np.array([[1.0],[1.0]])
gamma = np.array([0.0])

# Test panel mixture
N = 200
T = 3
M = 2
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
total_combinations = len(alphaset) * len(muset)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 2))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0
# Loop over parameters
for alpha in alphaset:
    for mu in muset:       
        start_time = time.time()
        result_rk_each = np.zeros((nrep,2))
        Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]
        data_nopar = [data[0] for data in Data]
        # Nonparametric test
        result_rk_each = NonParTestParallel(data_nopar, N, T, M, p, q, nrep, n_grid, BB, r_test)
        
        # Compute the mean results across replications and store them
        simulation_result_matrix[count, :] = result_rk_each.mean(axis=0)
        
         # Print execution time for this parameter combination
        print(f"Execution time for alpha={alpha}, mu={mu}: {time.time() - start_time:.2f} seconds")

        # Increment the counter
        count += 1

# Print the final result matrix
print("Simulation Result Matrix:")
print(simulation_result_matrix)

np.savetxt("TestStatNonparNoCovariates.txt", 100 * simulation_result_matrix, delimiter=',', fmt='%.1f')
