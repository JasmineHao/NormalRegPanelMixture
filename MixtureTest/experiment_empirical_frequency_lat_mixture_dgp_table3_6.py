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
# set_num_threads(12)

# Verify the number of threads
print(f"Numba is using {get_num_threads()} threads.")


result = pyreadr.read_r('ChileanClean.rds')  # Load the RDS file
df = result[None]  # Extract the dataframe

# Call the function
processed_data = process_chilean_data(each_code=381, T=5)
# processed_data = process_chilean_data(each_code=321, T=3)

# Access the results
y = processed_data['y']
x = processed_data['x_0']
z = np.zeros((x.shape[0], 0))
z = z.astype(np.float64)  # Convert z to float64

# %%


# Assuming these functions are already implemented:
# generate_data_mixture, regpanelmixmixturePMLE, NonParTestNoCovariates, LRTestNormal, LRTestMixture

# %%
# Main script
if __name__ == "__main__":
    # Parameters
    M = 2
    K = 2
    p = 0
    q = 0
    # T, N = 5, 225 # Example dimensions
    # T, N = 3, 196  # Example dimensions
    T, N = y.shape
    M_max = 6
    nrep = 100
    BB = 199

    # Obtain DGP parameters
    out_dgp = regpanelmixmixtureAR1NoConstraintPMLE(y,x,z, p, q, m=M,k = K, alpha_bound=0.05,ninits=2)
    
    alpha  = out_dgp['alpha_hat'][0]
    tau  = np.ascontiguousarray(out_dgp['tau_hat'][0]).reshape(M, K)
    mu = np.ascontiguousarray(out_dgp['mu_hat'][0]).reshape(M,K)
    mu_0 = np.ascontiguousarray(out_dgp['mu_0_hat'][0]).reshape(M,K)
    
    mubeta = out_dgp['mubeta_hat'][0]
    mubeta = np.ascontiguousarray(mubeta)
    mubeta_mat = mubeta.reshape(((q+2),M)).T
    beta = mubeta_mat[:,1:]
    
    sigma  = out_dgp['sigma_hat'][0]
    gamma  = out_dgp['gamma_hat'][0]
    
    mubeta_0 = out_dgp['mubeta_0_hat'][0]
    mubeta_0_hat = np.ascontiguousarray(mubeta_0)
    mubeta_0_hat_mat = mubeta_0_hat.reshape((q+1, M)).T
    beta_0 = mubeta_0_hat_mat[:,1:]
    sigma_0  = out_dgp['sigma_0_hat'][0]
    gamma_0  = out_dgp['gamma_0_hat'][0]
    
    print("alpha:", alpha)
    print("tau:", tau)
    print("mu:", mu)
    print("sigma:", sigma)
    print("gamma:", gamma)
    print("mubeta:", mubeta)
    print("mu_0:", mu_0)
    print("mubeta_0:", mubeta_0)
    print("sigma_0:", sigma_0)
    print("gamma_0:", gamma_0)
    
    Data = [generate_data_ar1_mixture_noconstraint(alpha, tau, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0, N, T, M, K, p, q) for _ in prange(nrep)]
    
    # Timing and execution
    start_time = time.time()
    results = parallel_processing_empirical_test_ar1_noconstraint(nrep, M_max, BB, Data, N, T, M, K, p, q, tau_bound=0.05)
    
    end_time = time.time()

    print("Parallel processing completed in:", end_time - start_time, "seconds")
    
    aic_table, bic_table, aic_table_mixture, bic_table_mixture, lr_estim_table, rk_mean_table, rk_max_table, lr_1_table, lr_5_table, lr_10_table, lr_1_table_mixture, lr_5_table_mixture, lr_10_table_mixture = results
    result_table = np.zeros((nrep,12))
    for ii in range(nrep):
        result_table[ii,0] = find_model_stop(aic_table[ii,:])
        result_table[ii,1] = find_model_stop(bic_table[ii,:])
        result_table[ii,2] = find_model_stop(aic_table_mixture[ii,:])
        result_table[ii,3] = find_model_stop(bic_table_mixture[ii,:])
        result_table[ii,4] = find_first_zero(lr_1_table[ii,:])
        result_table[ii,5] = find_first_zero(lr_5_table[ii,:])
        result_table[ii,6] = find_first_zero(lr_10_table[ii,:])
        result_table[ii,7] = find_first_zero(lr_1_table_mixture[ii,:]) 
        result_table[ii,8] = find_first_zero(lr_5_table_mixture[ii,:]) 
        result_table[ii,9] = find_first_zero(lr_10_table_mixture[ii,:]) 
        result_table[ii,10] = find_first_zero(rk_mean_table[ii,:]) 
        result_table[ii,11] = find_first_zero(rk_max_table[ii,:]) 
    # Frequency table
    result_freq_table = pd.DataFrame( index =  ["aic", "bic", "aic mixture", "bic mixture", "lr 1%", "lr 5%", "lr 10%", "lr 1%  mixture", "lr 5%  mixture", "lr 10% mixture", "rk mean 5%", "rk max 5%"], columns = [f"M={i}" for i in range(1, M_max + 1)])
    for kk in range(12):
        result_freq_table.iloc[kk,] = count_frequency(result_table[:,kk], M_max) / nrep
    # Set row and column names
        
    # Save to CSV
    result_freq_table.to_csv("test_empirical_dgp_lat_mixture.csv")

# %%
