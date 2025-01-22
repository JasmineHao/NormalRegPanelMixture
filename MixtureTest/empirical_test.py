# %%
from MixtureTestFunctions import *
import time
import pyreadr

chile_data = pyreadr.read_r('ChileanClean.rds')  # Load the RDS file
japan_data = pyreadr.read_r('JapanClean.rds')  # Load the RDS file


df = chile_data[None]  # Extract the dataframe
ind_code = [311,381,321] #chilean industry

df = japan_data[None]
ind_code = [5,13,12] #japan industry

p = 0 # No common regressor

model = "plain" 
# candidate models ["plain", "k", "k spline", "kl", "kl spline", "ar1 plain", "ar1 k", "ar1 kl"]

each_code = ind_code[0]

# %%
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


lr_df = np.zeros((3, 10), dtype=object)
AIC_df = np.zeros((3, 10))
BIC_df = np.zeros((3, 10))

# Assuming df is your DataFrame, and ind_code is a list of unique codes
for (each_code, ind_count) in zip(ind_code, range(len(ind_code))):
    
    t = time.time()  # Start timer
    # For Chilean industry
    # ---------------------
    ind_each = df.loc[df['ciiu_3d'] == each_code,:]  # Subset the DataFrame
    ind_name = ind_each['ciiu3d_descr'].iloc[0]  # Get the name
    ind_each.loc[:, 'y'] = np.log(ind_each['GO'])
    ind_each.loc[:, 'lnm'] = np.log(ind_each['WI'])
    ind_each.loc[:, 'lnl'] = np.log(ind_each['L'])
    ind_each.loc[:, 'lnk'] = np.log(ind_each['K'])

    # For Japanese industry
    # ---------------------
    ind_list = ["food","textile", "wood","paper", "chemical",
              "petro","plastic","ceramics","steel","othermetal",
              "metal product","machine","electronics",
              "transportation equipment","precision instrument",
              "other"]

    ind_each = df[df['industry_2'] == each_code]  # Subset where industry_2 matches each.code
    # Select specific columns
    ind_each = ind_each[['id', 'year', 'lnmY_it', 'k_it', 'l_it']]
    # Remove rows with missing values (complete.cases in R)
    ind_each = ind_each.dropna()
    # Standardize the 'k_it' column and store as 'lnk'
    ind_each['lnk'] = ind_each['k_it'] 
    # Standardize the 'l_it' column and store as 'lnl'
    ind_each['lnl'] = ind_each['l_it'] 
    # Standardize the 'lnmY_it' column and store as 'y'
    ind_each['y'] = ind_each['lnmY_it'] 
    # Map the industry name from ind_list using the industry code
    ind_name = ind_list[each_code]
    # ---------------------
    
    year_list = sorted(ind_each['year'].unique())
    T_cap = max(year_list)

    T = 3  
    t_start = T_cap - T + 1
    # Reshape the data
    ind_each_t = ind_each[ind_each['year'] >= t_start].dropna()
    ind_each_y = ind_each_t.pivot(index='id', columns='year', values='si')  # Reshape
    id_list = ind_each_y.dropna().index  # Get balanced panel IDs
    ind_each_t = ind_each_t[ind_each_t['id'].isin(id_list)].sort_values(['id', 'year'])

    # Reshape Y
    ind_each_y = ind_each_t.pivot(index='id', columns='year', values='y').drop(columns='id', errors='ignore')
    ind_each_y = (ind_each_y - ind_each_t['y'].mean()) / ind_each_t['y'].std()

    # Normalize X (lnk and lnl)
    ind_each_xk = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()
    ind_each_xl = (ind_each_t['lnl'] - ind_each_t['lnl'].mean()) / ind_each_t['lnl'].std()

    # Prepare data
    y = ind_each_y.T.to_numpy()  # Transpose Y
    x_k = ind_each_xk.to_numpy().reshape(-1, 1)  
    x_kl = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_xl.to_numpy().reshape(-1, 1)]
    
    N = ind_each_y.shape[0]
    x_0 = np.zeros((N*T,0)) # if only use Y 
    z = np.zeros((N*T,p))

    # %%
    bootstrap_k_cov = False    
    for m in range(1, 11):  # Loop for M = 1 to 10
        # Estimate the null model
        # M = 1        
        match model:
            case "plain":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            case "k":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            case "k spline":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= 199, spline=True)
            case "kl":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            case "kl spline":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= 199, spline=True)
            case "ar1 plain":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormalAR1(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            case "ar1 k":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormalAR1(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            case "ar1 kl":
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormalAR1(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            case _:
                print("Model Not Recognized")
                [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = [0,0,0,0,0,0,0]
        
        if lr_stat_k < lr_95_k:
                bootstrap_k_cov = False

        AIC_df[ind_count, m - 1] = round(aic_k, 2)
        BIC_df[ind_count, m - 1] = round(bic_k, 2)
        lr_df[ind_count, m -1] = lr_stat_k > lr_95_k
    

# %%
# ToDo: Create a table that record the following model index (m)
# The table should contain the following statistics for 3 industries in Chile, 
# 3 industries in Japan, different model specifications. 

find_model_stop(AIC_df[2,:])
find_model_stop(BIC_df[2,:])
find_first_zero(lr_df[2,:])
