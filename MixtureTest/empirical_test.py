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


COUNTRY = "chile"
COUNTRY_list = ['chile']

# %%
print(sys.argv)
y_indicator = 'mY_share'
y_indicator = 'mL_share'
include_ciiu = True
T = 5
MODEL = "plain_mixture"
    
if len(sys.argv) < 2:
    MODEL = "plain_mixture"
    
elif len(sys.argv) <= 3:
    MODEL = sys.argv[1]
    T = int(sys.argv[2])
    
elif len(sys.argv) <= 4:
    MODEL = sys.argv[1]
    T = int(sys.argv[2])
    if sys.argv[3] == 'ciiu':
        include_ciiu = True
    else:
        include_ciiu = False
else:
    MODEL = sys.argv[1]
    T = int(sys.argv[2])
    if sys.argv[3] == 'ciiu':
        include_ciiu = True
    else:
        include_ciiu = False
    y_indicator = sys.argv[4]

# %%
BB = 199
# T = 5
p = 0
epsilon_truncation = 0.01

print("Simulation Start")
print(f"Model: {MODEL}")    
print(f"T: {T}")    
print(f"Include CIIU: {include_ciiu}")
print(f"Dependent Variable: {y_indicator}")
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
            ind_code_dict = {311.0: 'Food products',
                            322.0: 'Wearing apparel, except footwear',
                            384.0: 'Transport equipment',
                            382.0: 'Machinery, except electrical',
                            381.0: 'Fabricated metal products',
                            362.0: 'Glass and products',
                            332.0: 'Manufacture of furniture and fixtures, except primarily of metal',
                            313.0: 'Beverages',
                            371.0: 'Iron and steel',
                            342.0: 'Printing and publishing',
                            331.0: 'Wood products, except furniture',
                            372.0: 'Non-ferrous metals',
                            369.0: 'Other non-metallic mineral products',
                            383.0: 'Machinery electric',
                            390.0: 'Other manufactured products',
                            352.0: 'Other chemicals',
                            351.0: 'Industrial chemicals',
                            312.0: 'Animal feeds, etc',
                            355.0: 'Rubber products',
                            321.0: 'Textiles',
                            356.0: 'Plastic products',
                            353.0: 'Petroleum refineries',
                            354.0: 'Misc. petroleum and coal products',
                            341.0: 'Paper and products',
                            323.0: 'Leather products',
                            324.0: 'Footwear, except rubber or plastic',
                            314.0: 'Tobacco',
                            385.0: 'Professional and scientific equipment',
                            361.0: 'Manufacture of pottery, china and earthenware'}
            df = chile_data[None]
            
            # Load the export data
            df_0 = pd.read_stata("DATA_KL_EXIM5_2023.dta")
            df_0['id'] = df_0['padron']
            
            # Ensure 'id' and 'year' columns have the same data type in both DataFrames
            df_0['id'] = pd.to_numeric(df_0['padron'], errors='coerce')
            df_0['year'] = pd.to_numeric('19' + df_0['year'].astype(str), errors='coerce')
            df_0['ciiu'] = df_0['ciiu'].astype(str)
            df_0['mi_share'] = df_0['rMi'].fillna(0) / (df_0['rMi'].fillna(0) + df_0['rMd'].fillna(0) ).fillna(0)
            df_0['mY_share'] = (df_0['rMi'].fillna(0) + df_0['rMd'].fillna(0) ).fillna(0) / df_0['rY'].fillna(0)
            df_0['K'] = df_0['rkb'] + df_0['rkm'] + df_0['rkt']        
            df_0['L'] = df_0['L_wh'] * df_0['w_wh'] + df_0['L_bl'] * df_0['w_bl']    
            df_0['mL_share'] = (df_0['rMi'].fillna(0) + df_0['rMd'].fillna(0) ).fillna(0) / df_0['L'].fillna(0)
            
            df['id'] = pd.to_numeric(df['id'], errors='coerce')
            df['year'] = pd.to_numeric(df['year'], errors='coerce')
            
            # Select relevant columns from df_0
            df_0 = df_0[['id', 'year','cu' , 'ciiu', 'mi_share', 'mY_share', 'mL_share', 'rY']]

            # Perform the merge
            df = pd.merge(df[['id','ciiu_3d','year','K','L']], df_0, on=['id', 'year'])
            
            for each_code in ind_code:
                ind_each = df.loc[df['ciiu_3d'] == each_code, ['id','year','mi_share', 'mY_share', 'mL_share', 'K','L', 'rY']]
                ind_each['lnk'] = np.log(ind_each['K'])
                ind_each['lnl'] = np.log(ind_each['L'])
                ind_each['lny'] = np.log(ind_each['rY'])
                ind_each['mi_share'] = ind_each['mi_share'].fillna(0)
                year_list = sorted(ind_each['year'].unique())
                T_cap = max(year_list)
                
                t_start = T_cap - T + 1
                ind_each = ind_each[ind_each['year'] >= t_start].dropna()
                ind_each_y = ind_each.pivot(index='id', columns='year', values='mY_share').dropna()
                
                N = ind_each_y.shape[0]
                
                id_list = ind_each_y.index
                ind_each = ind_each[ind_each['id'].isin(id_list)].sort_values(['id', 'year'])
                
                # Generate a description of the data
                description = pd.DataFrame({
                    "Column Name": ind_each.columns,
                    "Mean": ind_each.mean().values,
                    "Standard Deviation": ind_each.std().values,
                    "Minimum": ind_each.min().values,
                    "Maximum": ind_each.max().values
                })
                print(ind_code_dict[each_code])
                print(f'NT={len(ind_each)},N={len(ind_each['id'].unique())}')
                description = description.set_index('Column Name')
                
                print(description.loc[['mi_share','mY_share','mL_share','lnk','lnl','lny'],'Mean'].values)
                print(description.loc[['mi_share','mY_share','mL_share','lnk','lnl','lny'],'Standard Deviation'].values)

        case "japan":
            japan_data = pyreadr.read_r('JapanClean.rds')
            ind_code = [5, 13, 12]
            df = japan_data[None]
    # %%
    for (each_code, ind_count) in zip(ind_code, range(len(ind_code))):
        # %%
        start_time = time.time()

        # Initialize arrays for storing results
        lr_df = np.zeros(10, dtype=object)
        lr_df_99 = np.zeros(10, dtype=object)
        lr_df_90 = np.zeros(10, dtype=object)
        lr_estim = np.zeros(10, dtype=object)
        lr_crit = np.zeros(10, dtype=object)
        AIC_df = np.zeros(10)
        BIC_df = np.zeros(10)

        match COUNTRY:
            case "chile":
                ind_each = df.loc[df['ciiu_3d'] == each_code, :]
                INDNAME = ind_code_dict[each_code]
                ind_each['y'] = np.log(ind_each[y_indicator])
                ind_each['lnk'] = np.log(ind_each['K'])
                ind_each['lnl'] = np.log(ind_each['L'])
                ind_each['m_share'] = ind_each['mi_share'].fillna(0)
                # if include_ciiu:
                #     ind_each = ind_each[ind_each['ciiu'] != 'nan']
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
        
        N = ind_each_y.shape[0]
        
        id_list = ind_each_y.index
        ind_each_t = ind_each_t[ind_each_t['id'].isin(id_list)].sort_values(['id', 'year'])
        
        ind_each_xk = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()
        ind_each_xl = (ind_each_t['lnl'] - ind_each_t['lnl'].mean()) / ind_each_t['lnl'].std()
        
        ind_each_m_share = ind_each_t['m_share']

        y = ind_each_y.T.to_numpy().astype(np.float64)
        # y_ub = np.quantile(y, 1 - epsilon_truncation)
        # y_lb = np.quantile(y, epsilon_truncation)
        # y[y > y_ub] = y_ub
        # y[y < y_lb] = y_lb
                
        x_k = ind_each_xk.to_numpy().reshape(-1, 1).astype(np.float64)
        x_kl = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_xl.to_numpy().reshape(-1, 1)].astype(np.float64)
        x_kmshare = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_m_share.to_numpy().reshape(-1, 1)].astype(np.float64)
        
        x_0 = np.zeros((N * T, 0))
        if include_ciiu:
            ciiu_value_count = ind_each_t.groupby('ciiu')['ciiu'].count()
            ciiu_combine = ciiu_value_count[ciiu_value_count < 0.1 * (N * T)].index
            for each_ciiu in ciiu_combine:
                ind_each_t.loc[ind_each_t['ciiu'] == each_ciiu, 'ciiu'] = 'other'
            z = 1 * pd.get_dummies(ind_each_t['ciiu'], prefix='category').values[:, 1:]
            z = z.astype(np.float64)  # Convert z to float64
            p = z.shape[1]
        else:
            p = 0
            z = np.zeros((N * T, p))
        print(N)
        bootstrap_k_cov = True
        # %%
        for m in range(1, 11):
            start_time_model_m = time.time()

            # Use precompiled LRTest function
            match MODEL:
                case "nonpar":
                    n_grid = m + 1
                    r_test = m 
                    n_bins =  math.ceil((m+1)**(1/(T - 1)))
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k]  = NonParTest(y, N, T, n_grid, n_bins, BB, r_test)
                case "plain":
                    
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                case "k":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                
                case "kmshare":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_kmshare, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                    
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
                
                case "plain_mixture_4":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_0, z, p, 0, m, 4, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "k_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_k, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "k_mixture_3":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_k, z, p, 1, m, 3, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "k_mixture_4":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_k, z, p, 1, m, 4, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                
                case "kmshare_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kmshare, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "kmshare_mixture_3":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kmshare, z, p, 2, m, 3, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                    
                case "kmshare_mixture_4":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kmshare, z, p, 2, m, 4, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                    
                case "kl_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kl, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    
                case "kl_mixture_3":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kl, z, p, 3, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                
                    
                case "kl_mixture_spline":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_kl, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB, spline=True)
                    
                
                case "ar1_plain":
                    
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                    
                case "ar1_plain_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_0, z, p, 0, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                    
                
                case "ar1_plain_mixture_3":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_0, z, p, 0, m, 3, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1MixtureNoConstraint(y, x_0, z, p, 0, m, 3, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                
                case "ar1_k":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1NormalNoConstraint(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                                        
                case "ar1_k_mixture":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_k, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1MixtureNoConstraint(y, x_k, z, p, 1, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)

                case "ar1_k_mixture_3":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_k, z, p, 1, m, 3, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1MixtureNoConstraint(y, x_k, z, p, 1, m, 3, N, T, bootstrap=bootstrap_k_cov, BB= BB)
                case "ar1_kmshare":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_kmshare, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1NormalNoConstraint(y, x_kmshare, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)                
                case "ar1_kmshare_mixture":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_kmshare, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1MixtureNoConstraint(y, x_kmshare, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                
                case "ar1_kmshare_mixture_3":
                    if y_indicator == 'mY_share':
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_kmshare, z, p, 2, m, 3, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                    else:
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1MixtureNoConstraint(y, x_kmshare, z, p, 2, m, 3, N, T, bootstrap=bootstrap_k_cov, BB=BB)
                
                    
                case "ar1_kl_mixture":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_kl, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB)            
                case "ar1_kl_mixture_spline":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_kl, z, p, 2, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= BB, spline=True)
                                        
                case "ar1_kl":
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Normal(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= BB)
                
                case _:
                    print("Model Not Recognized")
                    [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = [0,0,0,0,0,0,0]

            
            if lr_stat_k < lr_90_k:
                bootstrap_k_cov = False

            AIC_df[m - 1] = round(aic_k, 2)
            BIC_df[m - 1] = round(bic_k, 2)
            lr_df[m - 1] = lr_stat_k > lr_95_k
            lr_df_99[m - 1] = lr_stat_k > lr_95_k
            lr_df_90[m - 1] = lr_stat_k > lr_90_k
            
            lr_estim[m - 1] = f'{lr_stat_k:.4f}'
            lr_crit[m - 1] = ", ".join(f"{x:.4f}" for x in [lr_90_k, lr_95_k, lr_99_k])

            end_time_model_m = time.time()
            print(f"Country: {COUNTRY}, Industry: {INDNAME}, Model: {MODEL}, M: {m}, Time Taken: {end_time_model_m - start_time_model_m:.2f} seconds, LR Stat: {lr_stat_k:.4f}")

        # Append results to lists
        result_each_df = [COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), find_model_stop(AIC_df), find_model_stop(BIC_df), find_first_zero(lr_df), find_first_zero(lr_df_90), find_first_zero(lr_df_99)]
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
    result_df = pd.DataFrame(result_df, columns=['Country', 'Industry', 'Model', 'AIC', 'BIC', 'LR', 'LR 90', 'LR 99'])
    statistics_df = pd.DataFrame(statistics_df, columns=['Country', 'Industry', 'Model', 'Stat'] + [f'M={d}' for d in range(1, 11)])

    # Save results to CSV
    if include_ciiu:
        statistics_df.to_csv(f"empirical_test/statistics_{COUNTRY}_{''.join(y_indicator.split('_'))}_{MODEL}_ciiu_{T}.csv")
        result_df.groupby(['Country', 'Industry', 'Model']).first().unstack(level="Model").T.to_csv(f"empirical_test/result_{COUNTRY}_{''.join(y_indicator.split('_'))}_{MODEL}_ciiu_{T}.csv")
    else:
        statistics_df.to_csv(f"empirical_test/statistics_{COUNTRY}_{''.join(y_indicator.split('_'))}_{MODEL}_{T}.csv")
        result_df.groupby(['Country', 'Industry', 'Model']).first().unstack(level="Model").T.to_csv(f"empirical_test/result_{COUNTRY}_{''.join(y_indicator.split('_'))}_{MODEL}_{T}.csv")

# %%
# ToDo: Create a table that record the following model index (m)
# The table should contain the following statistics for 3 industries in Chile, 
# 3 industries in Japan, different model specifications. 
# result_df_chile.to_csv("empirical_statistics_chile.csv")
# statistics_df_chile.to_csv("empirical_result_chile.csv")

