# %%
from MixtureTestFunctions import *
import time

# Set print options
np.set_printoptions(
    precision=3,
    threshold=None,
    edgeitems=None,
    linewidth=100,
    suppress=True,
    nanstr=None,
    infstr=None,
    formatter=None,
    sign=None,
    floatmode=None,
    legacy=None
)

# Define parameters
alpha = np.array([0.5, 0.5])  # Category probabilities (for M categories)
alphaset = [np.array([0.5, 0.5]), np.array([0.2, 0.8])]
mu = np.array([[-2.0, -1.0], [1.0, 2.0]])  # Mean values for each subcategory (M x K)
muset = [np.array([[-2.0, -1.0], [1.0, 2.0]]) ]
sigmaset = [np.array([0.5,0.5]), np.array([1.0,1.0])]
# sigma = np.array([1.0,1.0])  # Standard deviation for each category (length M)
# sigma = np.array([0.5,0.5])  # Standard deviation for each category (length M)
tau = np.array([[0.5, 0.5], [0.4, 0.6]])  # Subcategory probabilities for each M (M x K)
beta = None  # Coefficients for q covariates (M x q)
# beta = [[1, 0.5], [0.2, 0.3]]  # Coefficients for q covariates (M x q)
gam = None  # Coefficients for p covariates (length p)

N = 200  # Number of individuals
T = 5  # Number of time periods
M = 2  # Number of categories
K = 2  # Number of subcategories per category
p = 0  # Number of covariates in z
q = 0  # Number of covariates in x
nrep = 100
BB = 199

 # Determine the total number of parameter combinations
total_combinations = len(alphaset) * len(muset) * len(sigmaset)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 2))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0
# Loop over parameters
for alpha in alphaset:
    for mu in muset:
        for sigma in sigmaset:
            start_time = time.time()
            # Call the function
            Data = [generate_data_mixture(alpha, mu, sigma, tau, N, T, M, K, p, q) for _ in prange(nrep)]
            result_lr_each = LRTestMixtureParallel(Data, N, T, M, K, p, q, nrep, BB = BB)
            
            # Compute the mean results across replications and store them
            simulation_result_matrix[count, :] = result_lr_each.mean(axis=0)
            
            # Print execution time for this parameter combination
            print(f"Execution time for alpha={alpha}, mu={mu}: {time.time() - start_time:.2f} seconds")

            # Increment the counter
            count += 1

# Print the final result matrix
print("Simulation Result Matrix:")
print(simulation_result_matrix)

np.savetxt("TestStatMixture.txt", 100 * simulation_result_matrix, delimiter=',', fmt='%.1f')

        
# %%
data = generate_data_mixture(alpha, mu, sigma, tau, N, T, M, K, p, q) 
y = data[0]
x = data[1]
z = data[2]
estim_result = regpanelmixmixturePMLE(y,x,z, p, q, M, K)

print(estim_result["alpha_hat"])
print(estim_result["tau_hat"])
print(estim_result["sigma_hat"])
print(estim_result["mubeta_hat"])

# %%
