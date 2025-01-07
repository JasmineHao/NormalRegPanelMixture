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
rho = np.array([0.5, 0.7])  # Mean values for each subcategory (M x K)
mu = np.array([-1.0, 1.0])  # Mean values for each subcategory (M x K)
sigma = np.array([0.8,1.2])  # Standard deviation for each category (length M)

beta = np.array([[0.0], [0.0]])  # Coefficients for q covariates (M x q)
# beta = [[1, 0.5], [0.2, 0.3]]  # Coefficients for q covariates (M x q)
gamma = np.array([0.0]) # Coefficients for p covariates (length p)

N = 200  # Number of individuals
T = 5  # Number of time periods
M = 2  # Number of categories
p = 0  # Number of covariates in z
q = 0  # Number of covariates in x

# Call the function
data = generate_data_ar1(alpha, rho, mu, sigma, gamma, beta, N, T, M, p, q)
y = data[0]
x = data[1]
z = data[2]
# Test the estimation