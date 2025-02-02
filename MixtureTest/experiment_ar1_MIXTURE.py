# %%
from MixtureTestFunctions import *
import time

# Set print options
# Define parameters

alphaset = [np.array([0.5, 0.5]), np.array([0.2, 0.8])]
rhoset = [np.array([0.5, 0.7]), np.array([0.5, 0.9])]
tau = np.array([[0.5, 0.5], [0.4, 0.6]])  # Subcategory probabilities for each M 
muset = [np.array([[-4.0, -2.0], [2.0, 4.0]]) ]
sigma = np.array([0.5,0.5])  # Standard deviation for each category (length M)

# beta = np.array([[0.0], [0.0]])  # Coefficients for q covariates (M x q)

N = 200  # Number of individuals
T = 5  # Number of time periods
M = 2  # Number of categories
K = 2  # Number of sub-categories 
p = 0  # Number of covariates in z
q = 0  # Number of covariates in x
nrep = 10
BB = 19

beta = np.zeros((M,q))
# beta = np.array([[1.0], [-1.0]])  # Coefficients for q covariates (M x q)

gamma = np.zeros(p) # Coefficients for p covariates (length p)
# beta = [[1, 0.5], [0.2, 0.3]]  # Coefficients for q covariates (M x q)



# Call the function

 # Determine the total number of parameter combinations
total_combinations = len(alphaset) * len(rhoset)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 2))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0

mu = muset[0]
alpha = alphaset[0]
rho = rhoset[0]
# %%
# Loop over parameters
for alpha in alphaset:
    for rho in rhoset:       
        start_time = time.time()
        # Initial period distribution
        
        #mubeta
        mu_0 = np.zeros(mu.shape) 
        beta_0 = np.zeros(beta.shape)
        for mm in range(M):
            mu_0[mm] = mu[mm] / (1 - rho[mm])
            beta_0[mm] = beta[mm] / (1 - rho[mm])
        
        # sigma
        sigma_0_sq = sigma**2 / (1 - rho**2)
        for mm in range(M):
            for j in range(q):
                sigma_0_sq[mm] += beta[mm, j]**2 / (1- rho[mm]**2)
        sigma_0 = np.sqrt(sigma_0_sq)
        gamma_0 = gamma
        
        result_rk_each = np.zeros((nrep,2))
        Data = [generate_data_ar1_mixture(alpha, tau, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0,  N, T, M, K, p, q) for _ in prange(nrep)]

        # Test
        # ----
        # ii = 0
        # data = Data[ii]
        # y = data[0]
        # x = data[1]
        # z = data[2]
        # LRTestAR1Mixture(y, x, z, p, q, M, K, N, T, bootstrap = True, BB= 199)

        # Nonparametric test
        result_rk_each = LRTestAR1MixtureParallel(Data, N, T, M, p, q, nrep, BB = BB)
        
        # Compute the mean results across replications and store them
        simulation_result_matrix[count, :] = result_rk_each.mean(axis=0)
        
         # Print execution time for this parameter combination
        print(f"Execution time for alpha={alpha}, mu={mu}: {time.time() - start_time:.2f} seconds")

        # Increment the counter
        count += 1

# Print the final result matrix
print("Simulation Result Matrix:")
print(simulation_result_matrix)

np.savetxt("TestStatAR1Mixture.txt", 100 * simulation_result_matrix, delimiter=',', fmt='%.1f')




# %%
