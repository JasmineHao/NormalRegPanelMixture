# %%

# %%
from MixtureTestFunctions import *
import time

# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
alphaset = [np.array([0.5, 0.5]), np.array([0.2, 0.8])]
muset = [np.array([-1.0, 1.0]), np.array([-0.5, 0.5])]
sigmaset = [np.array([0.8, 1.2])]
betaset = [np.array([[0.0],[0.0]]), np.array([[1.0],[1.0]])]
gamma = np.zeros(0)

# Test panel mixture
N = 200
T = 3
M = 2
p = 0
q = 1
nrep = 100
BB = 199

n_grid=3
r_test=2
# Generate data
weights_equal = np.full(N, 1 / N)


# %%
test = True
if test:
    alpha = alphaset[0]
    mu = muset[0]
    sigma = sigmaset[0]
    beta = betaset[1]
    gamma = np.array([2.0])
    p = 1
    
    data = generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) 
    y = data[0]
    x = data[1]
    z = data[2]
    out_h0 = regpanelmixPMLE(y,x,z, p, q, M)
    print(out_h0)

# %%

 # Determine the total number of parameter combinations
total_combinations = len(alphaset) * len(muset) * len(betaset)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 3))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0
# Loop over parameters
for alpha in alphaset:
    for mu in muset: 
        for beta in betaset: 
            start_time = time.time()
            result_rk_each = np.zeros((nrep,2))
            Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]
            
            result_rk_each = NonParTestParallelCovariate(Data, N, T, M, p, q, nrep, n_grid=n_grid, BB=BB, r_test=r_test, n_bins=2)
            result_lr_each = LRTestParallel(Data, N, T, M, p, q, nrep, BB = BB)
            
            # Compute the mean results across replications and store them
            simulation_result_matrix[count, :] = np.append(result_rk_each.mean(axis=0), result_lr_each.mean())
            print(simulation_result_matrix[count, :])
            # Print execution time for this parameter combination
            print(f"Execution time for alpha={alpha}, mu={mu}, beta={beta}: {time.time() - start_time:.2f} seconds")

            parameter_combinations.append(f"alpha={alpha}, mu={mu}, beta={beta}")
            # Increment the counter
            count += 1

# %%
# Print the final result matrix
import pandas as pd
simulation_result_matrix_pd = pd.DataFrame(simulation_result_matrix, index=parameter_combinations, columns=['rk_mean','rk_max','lrt'])
(100 * simulation_result_matrix_pd).to_csv('NonparTestWithCovariates.csv')
# np.savetxt("NonparTestWithCovariates.txt", 100 * simulation_result_matrix, delimiter=',', fmt='%.1f')

# %%
