# %%
from MixtureTestFunctions import *
import time
import math

# Example parameter setup (adjust as needed)
N = 196
T = 3
M = 3  # Number of types
p = 0  # Number of Z
q = 0  # Number of X
M_max = 7
nrep = 100
BB = 199

# Replace with your actual parameter generation
alpha = np.array([0.1172174, 0.4086027, 0.4741799])
mu = np.array([-1.6401041, -0.3275637, 0.6876973])
sigma = np.array([1.0428541, 0.5219205, 0.6094685])
gamma = np.zeros(p)
beta = np.zeros((M,q))

Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]


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

# %%
# for m in range(M_max)
# %%
ii = 0 

start_time = time.time()
for ii in prange(nrep):
    data = Data[ii]
    y = data[0]
    x = data[1]
    z = data[2]
    bootstrap_nocov = True
    bootstrap_mixture = True
    for mm in range(M_max):
        # non-parametric test
        m = mm+1
        n_grid = m + 1
        r_test = m 
        n_bins =  math.ceil((m+1)**(1/(T - 1)))
        rk_stat_each = NonParTestNoCovariates(y, N, T, n_grid, n_bins, BB, r_test)
        # lr test
        [lr_stat_nocov, lr_90_nocov, lr_95_nocov, lr_99_nocov, aic_nocov, bic_nocov] = LRTestNormal(y, x, z, p, q, m, N, T, bootstrap = bootstrap_nocov, BB= 199)
        
        [lr_stat_mixture, lr_90_mixture, lr_95_mixture, lr_99_mixture, aic_mixture, bic_mixture] = LRTestMixture(y, x, z, p, q, m, 2, N, T, bootstrap = bootstrap_mixture, BB= 199)
        # record
        rk_max_table[ii,mm] = rk_stat_each[0]
        rk_mean_table[ii,mm] = rk_stat_each[1]
        
        aic_table[ii, mm] = aic_nocov
        bic_table[ii, mm] = aic_nocov

        lr_estim_table[ii,mm] = lr_stat_nocov

        lr_1_table[ii,mm] = lr_stat_nocov > lr_99_nocov
        lr_5_table[ii,mm] = lr_stat_nocov > lr_95_nocov
        lr_10_table[ii,mm] = lr_stat_nocov > lr_90_nocov

        lr_1_table_mixture[ii,mm] = lr_stat_mixture > lr_99_mixture
        lr_5_table_mixture[ii,mm] = lr_stat_mixture > lr_95_mixture
        lr_10_table_mixture[ii,mm] = lr_stat_mixture > lr_90_mixture

        if lr_stat_nocov < lr_90_nocov:
            bootstrap_nocov = False
        if lr_stat_mixture  < lr_90_mixture:
            bootstrap_mixture = False

end_time = time.time()
print("The time elapsed", end_time - start_time)
# %%
# Function to find the model where AIC stops decreasing
def find_model_stop(aic_values):
    differences = np.diff(aic_values)
    model_stop = np.where(differences > 0)[0]
    if model_stop.size == 0:
        return len(aic_values)
    return model_stop[0] + 1  # Adjust for 1-based indexing

# Function to count frequencies
def find_first_zero(sequence):
    for index, value in enumerate(sequence):
        if value == 0:
            return index + 1  # Return the first index where the value is 0
    return -1  # Return -1 if no 0 is found

def count_frequency(arr, M_max):
    """
    Counts the frequency of values in the range 1 to M_max in a NumPy array.

    Parameters:
        arr (numpy.ndarray): Input array of integers or floats.
        M_max (int): The maximum value to consider for frequency counting.

    Returns:
        numpy.ndarray: An array of frequencies for values from 1 to M_max.
    """
    # Convert the array to integers (floats will be truncated)
    arr = arr.astype(int)
    
    # Use np.bincount with minlength to ensure range 1 to M_max is covered
    counts = np.bincount(arr, minlength=M_max + 1)[1:M_max + 1]
    return counts

# %%
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
# result_freq_table 
    
# %%
import pandas as pd
# Frequency table
result_freq_table = pd.DataFrame( index =  ["aic", "bic", "lr 1%", "lr 5%", "lr 10%", "lr 1%  mixture", "lr 5%  mixture", "lr 10% mixture" , "rk mean 5%", "rk max 5%"], columns = [f"M={i}" for i in range(1, M_max + 1)])
for kk in range(7):
    result_freq_table.iloc[kk,] = count_frequency(result_table[:,kk], M_max) / nrep
# Set row and column names

# Save to CSV
result_freq_table.to_csv("test_empirical_M3.csv")
# %%
