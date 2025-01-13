# %%
from MixtureTestFunctions import *
import time
import pyreadr
result = pyreadr.read_r('ChileanClean.rds')  # Load the RDS file
df = result[None]  # Extract the dataframe

ind_code = [311,381,321]
p = 0 
# %%
each_code = ind_code[0]
# Assuming df is your DataFrame, and ind_code is a list of unique codes
for each_code in ind_code:
    t = time.time()  # Start timer
    ind_each = df.loc[df['ciiu_3d'] == each_code,:]  # Subset the DataFrame
    ind_name = ind_each['ciiu3d_descr'].iloc[0]  # Get the name
    ind_each.loc[:, 'y'] = np.log(ind_each['GO'])
    ind_each.loc[:, 'lnm'] = np.log(ind_each['WI'])
    ind_each.loc[:, 'lnl'] = np.log(ind_each['L'])
    ind_each.loc[:, 'lnk'] = np.log(ind_each['K'])

    ######################################################
    # Describe the data
    ######################################################
    desc_each = ind_each[ind_each['L'] != 0][['si', 'y', 'lnm', 'lnl', 'lnk']]
    year_list = sorted(ind_each['year'].unique())
    T_cap = max(year_list)

    coef_df = np.zeros((5, 10))
    lr_df = np.zeros((5, 10), dtype=object)
    AIC_df = np.zeros((5, 10))
    BIC_df = np.zeros((5, 10))

    ######################################################
    # For panel data
    ######################################################
    for T in range(3, 4):  # Loop for T = 3
        t_start = T_cap - T + 1
        # Reshape the data
        ind_each_t = ind_each[ind_each['year'] >= t_start].dropna()
        ind_each_y = ind_each_t.pivot(index='id', columns='year', values='si')  # Reshape
        id_list = ind_each_y.dropna().index  # Get balanced panel IDs
        ind_each_t = ind_each_t[ind_each_t['id'].isin(id_list)].sort_values(['id', 'year'])

        # Reshape Y
        ind_each_y = ind_each_t.pivot(index='id', columns='year', values='si').drop(columns='id', errors='ignore')
        ind_each_y = (ind_each_y - ind_each_t['si'].mean()) / ind_each_t['si'].std()

        # Normalize X
        ind_each_x = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()

        # Prepare data
        y = ind_each_y.T.to_numpy()  # Transpose Y
        x_k = ind_each_x.to_numpy().reshape(-1, 1)  # X as a matrix
        N = ind_each_y.shape[0]
        z = np.zeros((N*T,p))

        bootstrap_k_cov = True
        for M in range(1, 11):  # Loop for M = 1 to 10
            # Estimate the null model
            [lr_stat_k, lr_90_k, lr_95_k, lr_99_k, aic_k, bic_k] = LRTestNoCovariates(y, x_k, z, p, 1, M, N, T, bootstrap = bootstrap_k_cov, BB= 199)
            
            if lr_stat_k < lr_95_k:
                 bootstrap_k_cov = False

            AIC_df[T, M - 1] = round(aic_k, 2)
            BIC_df[T, M - 1] = round(bic_k, 2)
            lr_df[T, M-1] = round(lr_stat_k, 2)
