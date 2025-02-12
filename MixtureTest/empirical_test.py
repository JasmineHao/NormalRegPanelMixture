# %%
from ast import match_case
from MixtureTestFunctions import *
import time
import pyreadr
from numba import jit
import sys

import warnings
warnings.filterwarnings("ignore")

from numba import set_num_threads, get_num_threads

import os

# os.environ["OPENBLAS_NUM_THREADS"] = "1"

# Your code here
# Set the number of threads you want Numba to use
# set_num_threads(1)

# Verify the number of threads
print(f"Numba is using {get_num_threads()} threads.")

# MODEL_list = ["plain", "plain_mixture", "k", "k_spline", "kl", "kl_spline", "ar1_plain", "ar1_plain_mixture", "ar1_k", "ar1_kl"]

MODEL_list = ["plain_mixture", "k_mixture", "kl_mixture", "kl_mixture_spline", "ar1_plain_mixture", "ar1_k_mixture", "ar1_kl_mixture", "ar1_kl_mixture_spline"]
MODEL_list = ['plain', 'k', 'kl', 'kl_spline']

COUNTRY = "chile"
COUNTRY_list = ['chile']

print(sys.argv)
if len(sys.argv) < 2:
    MODEL = "plain_mixture" 
else:
    MODEL = sys.argv[1]
    T = int(sys.argv[2])


BB = 199
# T = 5
p = 0
epsilon_truncation = 0.01

print("Simulation Start")
print(f"Model: {MODEL}")    
print(f"T: {T}")    

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


# %%

for COUNTRY in COUNTRY_list:
    # %%
    statistics_df = []
    result_df = []
    match COUNTRY:
        case "chile":
            chile_data = pyreadr.read_r('ChileanClean.rds')
            ind_code = [311, 381, 321]
            df = chile_data[None]
        case "japan":
            japan_data = pyreadr.read_r('JapanClean.rds')
            ind_code = [5, 13, 12]
            df = japan_data[None]

    for (each_code, ind_count) in zip(ind_code, range(len(ind_code))):
        start_time = time.time()

        # Initialize arrays for storing results
        lr_df = np.zeros(10, dtype=object)
        lr_estim = np.zeros(10, dtype=object)
        lr_crit = np.zeros(10, dtype=object)
        AIC_df = np.zeros(10)
        BIC_df = np.zeros(10)

        match COUNTRY:
            case "chile":
                ind_each = df.loc[df['ciiu_3d'] == each_code, :]
                INDNAME = ind_each['ciiu3d_descr'].iloc[0]
                ind_each['y'] = np.log(ind_each['GO'])
                ind_each['lnm'] = np.log(ind_each['WI'])
                ind_each['lnl'] = np.log(ind_each['L'])
                ind_each['lnk'] = np.log(ind_each['K'])
            case "japan":
                ind_list = ["food", "textile", "wood", "paper", "chemical",
                            "petro", "plastic", "ceramics", "steel", "othermetal",
                            "metal product", "machine", "electronics",
                            "transportation equipment", "precision instrument",
                            "other"]
                ind_each = df[df['industry_2'] == each_code]
                ind_each = ind_each[['id', 'year', 'lnmY_it', 'k_it', 'l_it']].dropna()
                ind_each['lnk'] = ind_each['k_it']
                ind_each['lnl'] = ind_each['l_it']
                ind_each['y'] = ind_each['lnmY_it']
                INDNAME = ind_list[each_code]

        year_list = sorted(ind_each['year'].unique())
        T_cap = max(year_list)

        
        t_start = T_cap - T + 1
        ind_each_t = ind_each[ind_each['year'] >= t_start].dropna()
        ind_each_y = ind_each_t.pivot(index='id', columns='year', values='y').dropna()
        id_list = ind_each_y.index
        ind_each_t = ind_each_t[ind_each_t['id'].isin(id_list)].sort_values(['id', 'year'])

        ind_each_xk = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()
        ind_each_xl = (ind_each_t['lnl'] - ind_each_t['lnl'].mean()) / ind_each_t['lnl'].std()


        y = ind_each_y.T.to_numpy()
        # y_ub = np.quantile(y, 1 - epsilon_truncation)
        # y_lb = np.quantile(y, epsilon_truncation)
        # y[y > y_ub] = y_ub
        # y[y < y_lb] = y_lb
        
        x_k = ind_each_xk.to_numpy().reshape(-1, 1)
        x_kl = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_xl.to_numpy().reshape(-1, 1)]

        N = ind_each_y.shape[0]
        x_0 = np.zeros((N * T, 0))
        z = np.zeros((N * T, p))

        bootstrap_k_cov = True
        # %%
        for m in range(1, 11):
            start_time_model_m = time.time()

            # Use precompiled LRTest function
            match MODEL:
                case "plain":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                case "k":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                case "k_spline":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= BB, spline=True)
                case "kl":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                case "kl_spline":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB, spline=True)
                    
                case "plain_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_0, z, p, 0, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "plain_mixture_3":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_0, z, p, 0, m, 3, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "k_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_k, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "k_mixture_3":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_k, z, p, 1, m, 3, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "kl_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kl, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                case "kl_mixture_spline":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kl, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB, spline=True)
                    
                
                case "ar1_plain":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                
                case "ar1_plain_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_0, z, p, 0, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                case "ar1_k_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_k, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                case "ar1_kl_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_kl, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)            
                case "ar1_kl_mixture_spline":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_kl, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB, spline=True)
                    
                case "ar1_k":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                case "ar1_kl":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                
                case _:
                    print("Model Not Recognized")
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = [0,0,0,0,0,0,0]

            
            if lr_stat_k < lr_95_k:
                bootstrap_k_cov = False

            AIC_df[m - 1] = round(aic_k, 2)
            BIC_df[m - 1] = round(bic_k, 2)
            lr_df[m - 1] = lr_stat_k > lr_95_k
            lr_estim[m - 1] = f'{lr_stat_k:.4f}'
            lr_crit[m - 1] = ", ".join(f"{x:.4f}" for x in [lr_90_k, lr_95_k, lr_99_k])

            end_time_model_m = time.time()
            print(f"Country: {COUNTRY}, Industry: {INDNAME}, Model: {MODEL}, M: {m}, Time Taken: {end_time_model_m - start_time_model_m:.2f} seconds, LR Stat: {lr_stat_k:.4f}")

        # Append results to lists
        result_each_df = [COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), find_model_stop(AIC_df), find_model_stop(BIC_df), find_first_zero(lr_df)]
        result_df.append(result_each_df)

        statistics_df.append([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'AIC'] + list(AIC_df))
        statistics_df.append([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'BIC'] + list(BIC_df))
        statistics_df.append([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'LR'] + list(lr_estim))
        statistics_df.append([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'LR Crit'] + list(lr_crit))
        end_time = time.time()
        print(COUNTRY, INDNAME, end_time - start_time)
        print(result_each_df)
        print('-' * 30)

    # Convert results to DataFrames
    result_df = pd.DataFrame(result_df, columns=['Country', 'Industry', 'Model', 'AIC', 'BIC', 'LR'])
    statistics_df = pd.DataFrame(statistics_df, columns=['Country', 'Industry', 'Model', 'Stat'] + [f'M={d}' for d in range(1, 11)])

    # Save results to CSV
    statistics_df.to_csv(f"empirical_test/statistics_{COUNTRY}_{MODEL}_{T}.csv")
    result_df.groupby(['Country', 'Industry', 'Model']).first().unstack(level="Model").T.to_csv(f"empirical_test/result_{COUNTRY}_{MODEL}_{T}.csv")

# %%
# ToDo: Create a table that record the following model index (m)
# The table should contain the following statistics for 3 industries in Chile, 
# 3 industries in Japan, different model specifications. 
# result_df_chile.to_csv("empirical_statistics_chile.csv")
# statistics_df_chile.to_csv("empirical_result_chile.csv")

