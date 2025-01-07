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
mu = np.array([[-4.0, -2.0], [2.0, 4.0]])  # Mean values for each subcategory (M x K)
sigma = np.array([1.0,1.0])  # Standard deviation for each category (length M)
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

# Call the function
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
