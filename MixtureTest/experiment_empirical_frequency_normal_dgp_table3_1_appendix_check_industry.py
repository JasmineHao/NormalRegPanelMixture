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
print(f"Numba is using {get_num_threads()} threads.")

result = pyreadr.read_r('ChileanClean.rds')  # Load the RDS file
df = result[None]  # Extract the dataframe



# Parameter tuples for one-at-a-time changes

alpha_bound, tau_bound, epsilon = 0.05, 0.05, 0.05
if __name__ == "__main__":
    # [311, 321, 381]
    for ind_code in [311, 321]:
        print(f"Processing industry code: {ind_code}")
        processed_data = process_chilean_data(each_code=ind_code, T=3)
        y = processed_data['y']
        x = processed_data['x_0']
        z = np.zeros((x.shape[0], 0))
        z = z.astype(np.float64)  # Convert z to float64
        # Parameters
        M = 3
        K = 2
        p = 0
        q = 0
        T, N = y.shape
        M_max = 6
        nrep = 100
        BB = 199

        out_dgp = regpanelmixPMLE(y, x, z, p=p, q=q, m=M)
        alpha  = out_dgp['alpha_hat'][0]
        mubeta = out_dgp['mubeta_hat'][0]
        sigma = out_dgp['sigma_hat'][0]
        gamma = out_dgp['gamma_hat'][0]
        mubeta = np.ascontiguousarray(mubeta)
        mubeta_mat = mubeta.reshape((q+1,M)).T
        beta = mubeta_mat[:,1:]
        mu =  mubeta_mat[:,0]

        print("alpha:", alpha)
        print("sigma:", sigma)
        print("gamma:", gamma)
        print("beta:", beta)
        print("mu:", mu)
        print("N", N)

        # Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]

        # print(f"Running for alpha_bound={alpha_bound}, tau_bound={tau_bound}, epsilon={epsilon}")
        # start_time = time.time()
        # results = parallel_processing_empirical_test_stationary(
        #     nrep, M_max, BB, Data, N, T, M, K, p, q,
        #     alpha_bound=alpha_bound, tau_bound=tau_bound, epsilon=epsilon
        # )
        # end_time = time.time()
        # print("Parallel processing completed in:", end_time - start_time, "seconds")

        # (aic_table, bic_table, aic_table_mixture, bic_table_mixture, lr_estim_table,
        #     rk_mean_table, rk_max_table, lr_1_table, lr_5_table, lr_10_table,
        #     lr_1_table_mixture, lr_5_table_mixture, lr_10_table_mixture) = results

        # result_table = np.zeros((nrep,12))
        # for ii in range(nrep):
        #     result_table[ii,0] = find_model_stop(aic_table[ii,:])
        #     result_table[ii,1] = find_model_stop(bic_table[ii,:])
        #     result_table[ii,2] = find_model_stop(aic_table_mixture[ii,:])
        #     result_table[ii,3] = find_model_stop(bic_table_mixture[ii,:])
        #     result_table[ii,4] = find_first_zero(lr_1_table[ii,:])
        #     result_table[ii,5] = find_first_zero(lr_5_table[ii,:])
        #     result_table[ii,6] = find_first_zero(lr_10_table[ii,:])
        #     result_table[ii,7] = find_first_zero(lr_1_table_mixture[ii,:]) 
        #     result_table[ii,8] = find_first_zero(lr_5_table_mixture[ii,:]) 
        #     result_table[ii,9] = find_first_zero(lr_10_table_mixture[ii,:]) 
        #     result_table[ii,10] = find_first_zero(rk_mean_table[ii,:]) 
        #     result_table[ii,11] = find_first_zero(rk_max_table[ii,:]) 

        # result_freq_table = pd.DataFrame(
        #     index =  ["aic", "bic", "aic mixture", "bic mixture", "lr 1%", "lr 5%", "lr 10%", "lr 1%  mixture", "lr 5%  mixture", "lr 10% mixture" , "rk mean 5%", "rk max 5%"],
        #     columns = [f"M={i}" for i in range(1, M_max + 1)]
        # )
        # for kk in range(12):
        #     result_freq_table.iloc[kk,] = count_frequency(result_table[:,kk], M_max) / nrep

        # filename = f"test_empirical_dgp_normal_M3_ind{ind_code}.csv"
        # result_freq_table.to_csv(filename)

# %%
