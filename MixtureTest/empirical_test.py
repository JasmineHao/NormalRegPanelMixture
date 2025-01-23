# %%
from ast import match_case
from MixtureTestFunctions import *
import time
import pyreadr

import warnings
warnings.filterwarnings("ignore")

MODEL = "plain" 
MODEL_list = ["plain", "plain mixture", "k", "k spline", "kl", "kl spline", "ar1 plain", "ar1 plain mixture", "ar1 k", "ar1 kl"]
# candidate models ["plain", "k", "k spline", "kl", "kl spline", "ar1 plain", "ar1 k", "ar1 kl"]
COUNTRY = "chile"
COUNTRY_list = ['japan', 'chile']

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


p = 0 # No common regressor

# %%
statistics_df = pd.DataFrame()
result_df = pd.DataFrame()

for COUNTRY in COUNTRY_list:
    
    match COUNTRY:
        case "chile":
            chile_data = pyreadr.read_r('ChileanClean.rds')  # Load the RDS file
            ind_code = [311,381,321] #chilean industry
            df = chile_data[None]  # Extract the dataframe
        case "japan":
            japan_data = pyreadr.read_r('JapanClean.rds')  # Load the RDS file
            ind_code = [5,13,12] #japan industry
            df = japan_data[None]



    # Assuming df is your DataFrame, and ind_code is a list of unique codes
    for (each_code, ind_count) in zip(ind_code, range(len(ind_code))):

        start_time = time.time()
        lr_df = np.zeros(10, dtype=object)
        lr_estim = np.zeros(10, dtype=object)
        lr_crit = np.zeros(10, dtype=object)
        AIC_df = np.zeros(10)
        BIC_df = np.zeros(10)
        t = time.time()  # Start timer
        
        match COUNTRY:
            case "chile":
                # For Chilean industry
                # ---------------------
                ind_each = df.loc[df['ciiu_3d'] == each_code,:]  # Subset the DataFrame
                INDNAME = ind_each['ciiu3d_descr'].iloc[0]  # Get the name
                ind_each.loc[:, 'y'] = np.log(ind_each['GO'])
                ind_each.loc[:, 'lnm'] = np.log(ind_each['WI'])
                ind_each.loc[:, 'lnl'] = np.log(ind_each['L'])
                ind_each.loc[:, 'lnk'] = np.log(ind_each['K'])
            case "japan":
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
                # Map the industry name from ind_list using the industr code
                INDNAME = ind_list[each_code]
                # ---------------------
        
        year_list = sorted(ind_each['year'].unique())
        T_cap = max(year_list)

        T = 3  
        t_start = T_cap - T + 1
        # Reshape the data
        ind_each_t = ind_each[ind_each['year'] >= t_start].dropna()
        ind_each_y = ind_each_t.pivot(index='id', columns='year', values='y')  # Reshape
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

        for MODEL in MODEL_list:
            start_time_model = time.time()
            bootstrap_k_cov = True    
            # for m in range(1, 3):  # Loop for M = 1 to 10
            for m in range(1, 11):  # Loop for M = 1 to 10
                # Estimate the null model
                # M = 1        
                match MODEL:
                    case "plain":
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormal(y, x_0, z, p, 0, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
                    case "plain mixture":
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestMixture(y, x_0, z, p, 0, m, 2, N, T, bootstrap=bootstrap_k_cov, BB=199)
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
                    case "ar1 plain mixture":
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestAR1Mixture(y, x_0, z, p, 0, m, 2, N, T, bootstrap=bootstrap_k_cov, BB= 199)
                        
                    case "ar1 k":
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormalAR1(y, x_k, z, p, 1, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
                    case "ar1 kl":
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNormalAR1(y, x_kl, z, p, 2, m, N, T, bootstrap = bootstrap_k_cov, BB= 199)
                    
                    case _:
                        print("Model Not Recognized")
                        [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = [0,0,0,0,0,0,0]
                
                if lr_stat_k < lr_95_k:
                        bootstrap_k_cov = False

                AIC_df[m - 1] = round(aic_k, 2)
                BIC_df[m - 1] = round(bic_k, 2)
                lr_df[m -1] = lr_stat_k > lr_95_k
                lr_estim[m -1] = f'{lr_stat_k:.4f}' 
                lr_crit[m -1] = ", ".join(f"{x:.4f}" for x in [lr_90_k, lr_95_k, lr_99_k]) 
            
            result_each_df = pd.DataFrame([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), find_model_stop(AIC_df), find_model_stop(BIC_df), find_first_zero(lr_df) ]).T
            result_df = pd.concat([result_df, result_each_df],axis=0 )
            
            statistics_df = pd.concat([statistics_df,      pd.DataFrame([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'AIC' ] + list(AIC_df) ).T ])
            start_time_model = time.time()
            statistics_df = pd.concat([statistics_df,      pd.DataFrame([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'BIC' ] + list(BIC_df) ).T ])
            statistics_df = pd.concat([statistics_df,      pd.DataFrame([COUNTRY.capitalize(), INDNAME.capitalize(), MODEL.capitalize(), 'LR' ] + list(lr_estim) ).T ])
            end_time_model = time.time()
            print(COUNTRY, INDNAME, MODEL, end_time_model - start_time_model)
            
        end_time = time.time()
        print(COUNTRY, INDNAME, end_time - start_time)
        print(result_each_df)
        print('-*'*30)

result_df.columns = ['Country', 'Industry', 'Model', 'AIC', 'BIC', 'LR']
statistics_df.columns = ['Country', 'Industry', 'Model', 'Stat'] + [f'M={d}' for d in range(1,11)]

statistics_df.to_csv("empirical_statistics.csv")

print(result_df.groupby(['Country','Industry', 'Model']).first().unstack(level="Model").T)
result_df.groupby(['Country','Industry', 'Model']).first().unstack(level="Model").T.to_csv("empirical_result.csv")

# %%
# ToDo: Create a table that record the following model index (m)
# The table should contain the following statistics for 3 industries in Chile, 
# 3 industries in Japan, different model specifications. 
# result_df_chile.to_csv("empirical_statistics_chile.csv")
# statistics_df_chile.to_csv("empirical_result_chile.csv")


