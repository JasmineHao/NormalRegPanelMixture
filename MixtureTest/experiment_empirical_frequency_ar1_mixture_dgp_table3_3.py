# Example DataFrame
# %%
from MixtureTestFunctions import *
import time
import pyreadr
import pandas as pd
import numpy as np
import math
import time
from numba import njit, prange
from numba import set_num_threads, get_num_threads

import os
# os.environ["OPENBLAS_NUM_THREADS"] = "64"

# Your code here
# Set the number of threads you want Numba to use
# set_num_threads(64)

# Verify the number of threads
print(f"Numba is using {get_num_threads()} threads.")


result = pyreadr.read_r('ChileanClean.rds')  # Load the RDS file
df = result[None]  # Extract the dataframe

# Call the function
processed_data = process_chilean_data(df, each_code=321, T=3, p=0)

# Access the results
y = processed_data['processed_y']
x = processed_data['processed_x']
z = processed_data['processed_z']

# %%


# Assuming these functions are already implemented:
# generate_data_mixture, regpanelmixmixturePMLE, NonParTestNoCovariates, LRTestNormal, LRTestMixture

@njit(parallel=True)  # Numba JIT compilation with parallelization
def parallel_processing(nrep, M_max, BB, alpha, tau, mu, beta, rho, sigma, gamma, mu_0, beta_0, sigma_0, gamma_0, N, T, M, K, p, q):
    # Initialize result tables
    aic_table = np.zeros((nrep, M_max))
    bic_table = np.zeros((nrep, M_max))
    lr_estim_table = np.zeros((nrep, M_max))

    rk_mean_table = np.zeros((nrep, M_max))
    rk_max_table = np.zeros((nrep, M_max))

    lr_1_table = np.zeros((nrep, M_max))
    lr_5_table = np.zeros((nrep, M_max))
    lr_10_table = np.zeros((nrep, M_max))

    lr_1_table_mixture = np.zeros((nrep, M_max))
    lr_5_table_mixture = np.zeros((nrep, M_max))
    lr_10_table_mixture = np.zeros((nrep, M_max))

    # Generate data for all repetitions
    Data = [generate_data_ar1_mixture(alpha, tau, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0,  N, T, M, K, p, q) for _ in prange(nrep)]

    # Parallel loop
    for ii in prange(nrep):  # Use prange for parallel execution
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        bootstrap_nocov = True
        bootstrap_mixture = True

        for mm in range(M_max):
            # Non-parametric test
            m = mm + 1
            n_grid = m + 1
            r_test = m
            n_bins = math.ceil((m + 1) ** (1 / (T - 1)))
            rk_stat_each = NonParTestNoCovariates(y, N, T, n_grid, n_bins, BB, r_test)

            # LR test for no-covariates model
            lr_results_nocov = LRTestNormalAR1(y, x, z, p, q, m, N, T, bootstrap=bootstrap_nocov, BB=BB)
            lr_stat_nocov, lr_90_nocov, lr_95_nocov, lr_99_nocov, aic_nocov, bic_nocov = lr_results_nocov

            # LR test for mixture model
            lr_results_mixture = LRTestAR1Mixture(y, x, z, p, q, m, 2, N, T, bootstrap=bootstrap_mixture, BB=BB)
            lr_stat_mixture, lr_90_mixture, lr_95_mixture, lr_99_mixture, aic_mixture, bic_mixture = lr_results_mixture

            # Record results
            rk_max_table[ii, mm] = rk_stat_each[0]
            rk_mean_table[ii, mm] = rk_stat_each[1]

            aic_table[ii, mm] = aic_nocov
            bic_table[ii, mm] = aic_nocov

            lr_estim_table[ii, mm] = lr_stat_nocov

            lr_1_table[ii, mm] = lr_stat_nocov > lr_99_nocov
            lr_5_table[ii, mm] = lr_stat_nocov > lr_95_nocov
            lr_10_table[ii, mm] = lr_stat_nocov > lr_90_nocov

            lr_1_table_mixture[ii, mm] = lr_stat_mixture > lr_99_mixture
            lr_5_table_mixture[ii, mm] = lr_stat_mixture > lr_95_mixture
            lr_10_table_mixture[ii, mm] = lr_stat_mixture > lr_90_mixture

            # Update bootstrap flags
            if lr_stat_nocov < lr_90_nocov:
                bootstrap_nocov = False
            if lr_stat_mixture < lr_90_mixture:
                bootstrap_mixture = False

    return (aic_table, bic_table, lr_estim_table, rk_mean_table, rk_max_table,
            lr_1_table, lr_5_table, lr_10_table,
            lr_1_table_mixture, lr_5_table_mixture, lr_10_table_mixture)

# Function to find the model where AIC stops decreasing
@njit
def find_model_stop(aic_values):
    differences = np.diff(aic_values)
    model_stop = np.where(differences > 0)[0]
    if model_stop.size == 0:
        return len(aic_values)
    return model_stop[0] + 1  # Adjust for 1-based indexing

# Function to count frequencies
@njit
def find_first_zero(sequence):
    for index, value in enumerate(sequence):
        if value == 0:
            return index + 1  # Return the first index where the value is 0
    return -1  # Return -1 if no 0 is found


def count_frequency(arr, M_max):
    """
    Counts the frequency of values in the range 1 to M_max in a NumPy array.

    Parameters:
        arr (numpy.ndarray): Input array of integers.
        M_max (int): The maximum value to consider for frequency counting.

    Returns:
        numpy.ndarray: An array of frequencies for values from 1 to M_max.
    """
    # Convert the array to integers (if not already)
    arr = arr.astype(int)
    
    # Use np.bincount with minlength to ensure range 1 to M_max is covered
    counts = np.bincount(arr, minlength=M_max + 1)[1:M_max + 1]
    return counts

# Main script
if __name__ == "__main__":
    # Parameters
    M = 3
    K = 2
    p = 0
    q = 0
    T, N = 3, 196  # Example dimensions
    M_max = 7
    nrep = 100
    BB = 199

    # Obtain DGP parameters
    out_dgp = regpanelmixAR1mixturePMLE(y, x, z, p, q, M, K)
    alpha  = out_dgp['alpha_hat'][0]
    tau  = np.ascontiguousarray(out_dgp['tau_hat'][0]).reshape(M,K)
    
    mu = np.ascontiguousarray(out_dgp['mu_hat'][0]).reshape(M,K)
    beta = out_dgp['beta_hat']
    rho = out_dgp['rho_hat'][0]
    sigma  = out_dgp['sigma_hat'][0]
    gamma  = out_dgp['gamma_hat'][0]
    
    mu_0 = np.ascontiguousarray(out_dgp['mu_0_hat']).reshape(M,K) 
    beta_0 = out_dgp['beta_0_hat']

    sigma_0  = out_dgp['sigma_0_hat'][0]
    gamma_0  = out_dgp['gamma_0_hat'][0]
    
    # Timing and execution
    start_time = time.time()
    results = parallel_processing(nrep, M_max, BB, alpha, tau, mu, beta, rho, sigma, gamma, mu_0, beta_0, sigma_0, gamma_0, N, T, M, K, p, q)
    end_time = time.time()

    print("Parallel processing completed in:", end_time - start_time, "seconds")
    
    aic_table, bic_table, lr_estim_table, rk_mean_table, rk_max_table, lr_1_table, lr_5_table, lr_10_table, lr_1_table_mixture, lr_5_table_mixture, lr_10_table_mixture = results
    result_table = np.zeros((nrep,10))
    for ii in range(nrep):
        result_table[ii,0] = find_model_stop(aic_table[ii,:])
        result_table[ii,1] = find_model_stop(bic_table[ii,:])
        result_table[ii,2] = find_first_zero(lr_1_table[ii,:])
        result_table[ii,3] = find_first_zero(lr_5_table[ii,:])
        result_table[ii,4] = find_first_zero(lr_10_table[ii,:])
        result_table[ii,5] = find_first_zero(lr_1_table_mixture[ii,:]) 
        result_table[ii,6] = find_first_zero(lr_5_table_mixture[ii,:]) 
        result_table[ii,7] = find_first_zero(lr_10_table_mixture[ii,:]) 
        result_table[ii,8] = find_first_zero(rk_mean_table[ii,:]) 
        result_table[ii,9] = find_first_zero(rk_max_table[ii,:]) 
    # Frequency table
    result_freq_table = pd.DataFrame( index =  ["aic", "bic", "lr 1%", "lr 5%", "lr 10%", "lr 1%  mixture", "lr 5%  mixture", "lr 10% mixture" , "rk mean 5%", "rk max 5%"], columns = [f"M={i}" for i in range(1, M_max + 1)])
    for kk in range(10):
        result_freq_table.iloc[kk,] = count_frequency(result_table[:,kk], M_max) / nrep
    # Set row and column names

    # Save to CSV
    result_freq_table.to_csv("empirical_test_dgp_ar1_mixture_M3.csv")

# %%
