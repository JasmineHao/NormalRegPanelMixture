# %%
from MixtureTestFunctions import *
import time

# Set print options
# Define parameters

# alpha = np.array([0.5, 0.5])  # Category probabilities (for M categories)
alphaset = [np.array([0.5, 0.5]), np.array([0.2, 0.8])]
rhoset = [np.array([0.5, 0.7]), np.array([0.5, 0.9])]
# rho = np.array([0.5, 0.7])  # Mean values for each subcategory (M x K)
mu = np.array([-1.0, 1.0])  # Mean values for each subcategory (M x K)
sigma = np.array([0.8,1.2])  # Standard deviation for each category (length M)


# beta = np.array([[0.0], [0.0]])  # Coefficients for q covariates (M x q)

N = 200  # Number of individuals
T = 5  # Number of time periods
M = 2  # Number of categories
p = 0  # Number of covariates in z
q = 1  # Number of covariates in x
nrep = 10
BB = 19

beta = np.zeros((M,q))
gamma = np.zeros(p) # Coefficients for p covariates (length p)
beta = np.array([[1.0], [-1.0]])  # Coefficients for q covariates (M x q)


# %%
# Call the function

 # Determine the total number of parameter combinations
total_combinations = len(alphaset) * len(rhoset)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 2))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0
# Loop over parameters
for alpha in alphaset:
    for rho in rhoset:       
        # alpha = alphaset[0]
        # rho = rhoset[0]
        start_time = time.time()
        # Initial period distribution
        mu_0 = mu / (1 - rho)
        beta_0 = np.zeros(beta.shape)
        # beta 
        sigma_0_sq = sigma**2 / (1 - rho**2)
        for mm in prange(M):
            for j in prange(q):
                beta_0[mm, j] = beta[mm, j] / (1- rho[mm])
                sigma_0_sq[mm] += beta[mm, j]**2 / (1- rho[mm]**2)
        sigma_0 = np.sqrt(sigma_0_sq)
        gamma_0 = gamma
        
        result_rk_each = np.zeros((nrep,2))
        Data = [generate_data_ar1(alpha, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0, N, T, M, p, q) for _ in prange(nrep)]

        # Nonparametric test
        result_rk_each = LRTestAR1Parallel(Data, N, T, M, p, q, nrep, BB = BB)
        
        # Compute the mean results across replications and store them
        simulation_result_matrix[count, :] = result_rk_each.mean(axis=0)
        
         # Print execution time for this parameter combination
        print(f"Execution time for alpha={alpha}, mu={mu}: {time.time() - start_time:.2f} seconds")

        # Increment the counter
        count += 1

# Print the final result matrix
print("Simulation Result Matrix:")
print(simulation_result_matrix)

np.savetxt("TestStatAR1.txt", 100 * simulation_result_matrix, delimiter=',', fmt='%.1f')




# %%
