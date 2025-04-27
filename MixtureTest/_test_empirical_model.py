# %%
# Import necessary libraries and functions
from MixtureTestFunctions import *
import time
from scipy.optimize import minimize
from numdifftools import Hessian
import numpy as np
import pandas as pd

# %%
# Set print options and define parameters
# ----------------------------------------

# Mixture model parameters
alpha = np.array([0.5, 0.5])  # Category probabilities (for M categories)
rho = np.array([0.5, 0.7])    # Mean values for each subcategory (M x K)
mu = np.array([-1.0, 1.0])    # Mean values for each subcategory (M x K)
sigma = np.array([0.8, 1.2])  # Standard deviation for each category (length M)

# Regression coefficients
beta = np.array([[0.0, 1.0], [1.0, 0.0]])  # Coefficients for q covariates (M x q)
gamma = np.array([0.])                     # Coefficients for p covariates (length p)

# Model dimensions
N = 200  # Number of individuals
T = 5    # Number of time periods
M = 2    # Number of categories

# %%
# Test case: Stationary normal mixture models
# -------------------------------------------

p = 1    # Number of covariates in z
q = 2    # Number of covariates in x

# Simulation results storage
simulation_results = []

# Combine true parameters into a single array for comparison
true_parameters = np.concatenate([alpha[:-1], mu, beta.flatten(), sigma, gamma.flatten()])

parameter_names = (
    [f"alpha_{i+1}" for i in range(len(alpha) - 1)] +
    [f"mu_{i+1}" for i in range(len(mu))] +
    [f"beta_{i+1}{j+1}" for i in range(beta.shape[0]) for j in range(beta.shape[1])] +
    [f"sigma_{i+1}" for i in range(len(sigma))] +
    [f"gamma_{i+1}" for i in range(len(gamma))]
)

# Run multiple simulations
for iteration in range(100):
    print(f"Running simulation {iteration + 1}...")  # Display progress

    # Generate simulated data
    data = generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q)
    y, x, z = data

    # Fit the model and extract parameters
    model_output = regpanelmixPMLE(y, x, z, p, q, M)
    params_dict, params_array = get_params_stationary_normal(model_output)
    
    score = score_stationary_normal(data, params_dict)
    hessian = score.T @ score     
    hessian_inv = np.linalg.inv(hessian)

    standard_errors = np.sqrt(np.diag(hessian_inv))

    # Check if parameters are within the confidence interval
    within_confidence_interval = (
        np.abs(params_array - true_parameters) < 1.96 * standard_errors
    )
    simulation_results.append(within_confidence_interval)

# Convert results to a matrix and compute the mean
results_matrix = np.array(simulation_results)
confidence_interval_coverage = np.mean(results_matrix, axis=0)
print("Confidence Interval Coverage:", confidence_interval_coverage)

# %%
# Test case: Stationary mixture models with covariates
# ----------------------------------------------------

M = 2
K = 2
p = 1
q = 2
N = 200  # Number of individuals
T = 5    # Number of time periods

# Mixture model parameters
alpha = np.array([0.5, 0.5])  # Category probabilities (for M categories)
mu = np.array([[-2.0, -1.0], [1.0, 2.0]])  # Mean values for each subcategory 
tau = np.array([[0.5, 0.5], [0.4, 0.6]])  # Subcategory probabilities for each 
sigma = np.array([0.5, 0.5])  # Standard deviation for each category (length M)
gamma = np.array([0.5])
beta = np.array([[1, 0.5], [0.2, 0.3]])  # Coefficients for q covariates (M x q)

# Simulation results storage
simulation_results = []

# Combine true parameters into a single array for comparison
true_parameters = np.concatenate([alpha[:-1], tau[:, :-1].flatten(), mu.flatten(), beta.T.flatten(), sigma, gamma.flatten()])

# Run multiple simulations
for iteration in range(100):
    start_time = time.time()
    
    data = generate_data_mixture(alpha, mu, beta, sigma, tau, gamma, N, T, M, K, p, q)
    y, x, z = data

    # Fit the stationary mixture model
    out_h0 = regpanelmixmixturePMLE(y, x, z, p, q, M, K)
    params_dict, params_array = get_params_stationary_mixture(out_h0)

    # Compute the score and Hessian for stationary model
    score, hessian = score_stationary_mixture(data, params_dict) 
    
    # Compute the sandwich formula for standard errors
    bread = np.linalg.pinv(np.mean(hessian, axis=0))  # Inverse of the Hessian (bread matrix)
    meat = score.T @ score  # Outer product of the score (meat matrix)
    sandwich = bread @ meat @ bread  # Sandwich formula
    stationary_standard_errors = np.sqrt(np.abs(np.diag(sandwich))) / N
    
    stationary_within_confidence_interval = (
        np.abs(params_array - true_parameters) < 1.96 * stationary_standard_errors
    )
    simulation_results.append(stationary_within_confidence_interval)

    end_time = time.time()
    print(f"Iteration {iteration + 1} duration: {end_time - start_time} seconds")

# Convert results to a matrix and compute the mean
results_matrix = np.array(simulation_results)
confidence_interval_coverage = np.mean(results_matrix, axis=0)

# Print results for stationary mixture model
print("Stationary Model Parameters:", true_parameters)
print("Stationary Model Standard Errors:", stationary_standard_errors)
print("Stationary Model Confidence Interval Coverage:", confidence_interval_coverage)

# %%
# Test case: Stationary mixture models with no covariates
# -------------------------------------------------------

M = 2
K = 2
p = 0
q = 0
N = 200  # Number of individuals
T = 5    # Number of time periods

# Mixture model parameters
alpha = np.array([0.5, 0.5])  # Category probabilities (for M categories)
mu = np.array([[-2.0, -1.0], [1.0, 2.0]])  # Mean values for each subcategory 
tau = np.array([[0.5, 0.5], [0.4, 0.6]])  # Subcategory probabilities for each 
sigma = np.array([0.5, 0.5])  # Standard deviation for each category (length M)

beta = np.zeros((M, q))  # No covariates
gamma = np.zeros(p)  # No covariates

# Simulation results storage
simulation_results = []

# Combine true parameters into a single array for comparison
true_parameters = np.concatenate([alpha[:-1], tau[:, :-1].flatten(), mu.flatten(), beta.T.flatten(), sigma, gamma.flatten()])

# Run multiple simulations
for iteration in range(100):
    start_time = time.time()
    
    data = generate_data_mixture(alpha, mu, beta, sigma, tau, gamma, N, T, M, K, p, q)
    y, x, z = data

    # Fit the stationary mixture model
    out_h0 = regpanelmixmixturePMLE(y, x, z, p, q, M, K)
    params_dict, params_array = get_params_stationary_mixture(out_h0)

    # Compute the score and Hessian for stationary model
    score, hessian = score_stationary_mixture(data, params_dict) 
    
    # Compute the sandwich formula for standard errors
    bread = np.linalg.pinv(np.mean(hessian, axis=0))  # Inverse of the Hessian (bread matrix)
    meat = score.T @ score  # Outer product of the score (meat matrix)
    sandwich = bread @ meat @ bread  # Sandwich formula
    stationary_standard_errors = np.sqrt(np.abs(np.diag(sandwich))) / N
    
    stationary_within_confidence_interval = (
        np.abs(params_array - true_parameters) < 1.96 * stationary_standard_errors
    )
    simulation_results.append(stationary_within_confidence_interval)

    end_time = time.time()
    print(f"Iteration {iteration + 1} duration: {end_time - start_time} seconds")

# Convert results to a matrix and compute the mean
results_matrix = np.array(simulation_results)
confidence_interval_coverage = np.mean(results_matrix, axis=0)

# Print results for stationary mixture model
print("Stationary Model Parameters:", true_parameters)
print("Stationary Model Standard Errors:", stationary_standard_errors)
print("Stationary Model Confidence Interval Coverage:", confidence_interval_coverage)

# %%
# Test case: AR(1) mixture error model with covariates
# ----------------------------------------------------

# Parameters for AR(1) mixture model
alpha = np.array([0.5, 0.5])  # Category probabilities
rho = np.array([0.5, 0.7])    # AR(1) coefficients
mu = np.array([[-4.0, -2.0], [2.0, 4.0]])
tau = np.array([[0.5, 0.5], [0.4, 0.6]])  # Subcategory probabilities for each M 
sigma = np.array([0.1, 0.1])  # Standard deviations
beta = np.array([[0.0, 1.0], [1.0, 0.0]])  # Regression coefficients
gamma = np.array([0.])        # Coefficients for covariates
N = 200  # Number of individuals
T = 5    # Number of time periods
M = 2    # Number of categories
K = 2    # Number of sub-categories 

p = 1    # Number of covariates in z
q = 2    # Number of covariates in x

mu_0 = mu
beta_0 = beta        
sigma_0 = np.array([0.1, 0.1])
gamma_0 = gamma

# Simulation results storage
simulation_results = []

# Combine true parameters into a single array for comparison
true_parameters = np.concatenate([
    alpha[:-1],
    tau[:, :-1].flatten(),
    rho,
    mu.flatten(),
    beta.T.flatten(),
    sigma,
    gamma.flatten(),
    mu_0.flatten(),
    beta_0.flatten(),
    sigma_0,
    gamma_0
])

parameter_names = (
    [f"alpha_{i+1}" for i in range(len(alpha) - 1)] +
    [f"tau_{m+1}{k+1}" for m in range(tau.shape[0]) for k in range(tau.shape[1] - 1)] +
    [f"rho_{m+1}" for m in range(len(rho))] +
    [f"mu_{m+1}{k+1}" for m in range(mu.shape[0]) for k in range(mu.shape[1])] +
    [f"beta_{m+1}{j+1}" for j in range(beta.shape[1]) for m in range(beta.shape[0])] +
    [f"sigma_{m+1}" for m in range(len(sigma))] +
    [f"gamma_{i+1}" for i in range(len(gamma))] +
    [f"mu_0_{m+1}{k+1}" for m in range(mu_0.shape[0]) for k in range(mu_0.shape[1])] +
    [f"beta_0_{m+1}{j+1}" for j in range(beta_0.shape[1]) for m in range(beta_0.shape[0])] +
    [f"sigma_0_{m+1}" for m in range(len(sigma_0))] +
    [f"gamma_0_{i+1}" for i in range(len(gamma_0))]
)

# Run multiple simulations
for iteration in range(100):
    start_time = time.time()
    
    data = generate_data_ar1_mixture(alpha, tau, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0, N, T, M, K, p, q)
    y, x, z = data
    out_h0 = regpanelmixAR1mixturePMLE(y, x, z, p, q, M, K)
    
    params_dict, params_array = get_params_ar1_mixture(out_h0)
    score, hessian = score_ar1_mixture(data, params_dict)
    
    # Compute the sandwich formula for standard errors
    bread = np.linalg.pinv(np.mean(hessian, axis=0))  # Inverse of the Hessian (bread matrix)
    meat = score.T @ score  # Outer product of the score (meat matrix)
    sandwich = bread @ meat @ bread  # Sandwich formula
    ar1_standard_errors = np.sqrt(np.abs(np.diag(sandwich))) / N
    
    ar1_within_confidence_interval = (
        np.abs(params_array - true_parameters) < 1.96 * ar1_standard_errors
    )
    simulation_results.append(ar1_within_confidence_interval)

    end_time = time.time()
    print(f"Iteration {iteration + 1} duration: {end_time - start_time} seconds")

# Convert results to a matrix and compute the mean
results_matrix = np.array(simulation_results)
confidence_interval_coverage = np.mean(results_matrix, axis=0)

# Print results for AR(1) mixture model
# Create a DataFrame to tabulate results
results_df = pd.DataFrame({
    "Parameter Name": parameter_names,
    "True Parameter": true_parameters,
    "Standard Error": ar1_standard_errors,
    "Confidence Interval Coverage": confidence_interval_coverage
})

# Print the table
print(results_df)
print("AR(1) Mixture Model Confidence Interval Coverage:", confidence_interval_coverage)

# %%
# Test case: AR(1) normal error model
# -----------------------------------

rho = np.array([0.5, 0.7])    # AR(1) coefficients
mu = np.array([-1.0, 1.0])    # Mean values for each subcategory (M x K)
sigma = np.array([0.8, 1.2])  # Standard deviation for each category (length M)

# Derived parameters
mu_0 = mu / (1 - rho)
beta_0 = np.zeros(beta.shape)
sigma_0_sq = sigma**2 / (1 - rho**2)

# Adjust beta_0 and sigma_0_sq based on beta and rho
for mm in range(M):
    for j in range(q):
        beta_0[mm, j] = beta[mm, j] / (1 - rho[mm])
        sigma_0_sq[mm] += beta[mm, j]**2 / (1 - rho[mm]**2)

sigma_0 = np.sqrt(sigma_0_sq)
gamma_0 = gamma

# Generate data
data = generate_data_ar1(alpha, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0, N, T, M, p, q)
y, x, z = data

# Fit the AR(1) normal error model
out_h0 = regpanelmixAR1PMLE(y, x, z, p, q, M)

# Directly print the fields from the result dictionary
print("alpha_hat:", out_h0['alpha_hat'], "True alpha:", alpha)
print("rho_hat:", out_h0['rho_hat'], "True rho:", rho)
print("mu_hat:", out_h0['mu_hat'], "True mu:", mu)
print("sigma_hat:", out_h0['sigma_hat'], "True sigma:", sigma)
print("mubeta_hat:", out_h0['mubeta_hat'], "True beta:", beta)
print("gamma_hat:", out_h0['gamma_hat'], "True gamma:", gamma)
print("sigma_0_hat:", out_h0['sigma_0_hat'], "True sigma_0:", sigma_0)
print("mubeta_0_hat:", out_h0['mubeta_0_hat'], "True beta_0:", beta_0)
print("gamma_0_hat:", out_h0['gamma_0_hat'], "True gamma_0:", gamma_0)
# %%
