# %%

ind_code = [311, 381, 321]

each_code = 321
for (each_code, ind_count) in zip(ind_code, range(len(ind_code))):
    # %%
    ind_each = df.loc[df['ciiu_3d'] == each_code, :]
    INDNAME = ind_each['ciiu3d_descr'].iloc[0]
    ind_each['y'] = np.log(ind_each['GO'])
    ind_each['lnm'] = np.log(ind_each['WI'])
    ind_each['lnl'] = np.log(ind_each['L'])
    ind_each['lnk'] = np.log(ind_each['K'])
    
    year_list = sorted(ind_each['year'].unique())
    T_cap = max(year_list)

    T = 3
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
    print('=*'* 50)
    # %%
    for m in range(1,5):
        # out_h0 = regpanelmixAR1mixturePMLE(y, x_k, z, p, 1, m, 2)
        # out_h1 = regpanelmixAR1mixturePMLE(y, x_k, z, p, 1, m+1, 2)
        # out_h0 = regpanelmixAR1mixturePMLE(y, x_0, z, p, 0, m, 2)
        # out_h1 = regpanelmixAR1mixturePMLE(y, x_0, z, p, 0, m+1, 2)
        # out_h0 = regpanelmixAR1PMLE(y,x_k,z, p, 1, m)
        # out_h1 = regpanelmixAR1PMLE(y,x_k,z, p, 1, m+1)
        # out_h0 = regpanelmixAR1PMLE(y,x_0,z, p, 0, m)
        # out_h1 = regpanelmixAR1PMLE(y,x_0,z, p, 0, m+1)
        # out_h0 = regpanelmixmixturePMLE(y,x_k,z, p, 1, m, 2)
        # out_h1 = regpanelmixmixturePMLE(y,x_k,z, p, 1, m+1, 2)
        # out_h0 = regpanelmixmixturePMLE(y,x_0,z, p, 0, m, 2)
        # out_h1 = regpanelmixmixturePMLE(y,x_0,z, p, 0, m+1, 2)
        penloglik_h0 = out_h0['penloglik'][0,0]
        aic = out_h0['aic'][0,0]
        bic = out_h0['bic'][0,0]
        penloglik_h1 = out_h1['penloglik'][0,0]
        # print(out_h0[0])        
        lr_stat = -2 * (penloglik_h0 - penloglik_h1)
        print('--'* 30)
        print(f"m={m}")
        print(f'lr={lr_stat}, bic={bic}')
        print(f'alpha: {out_h0['alpha_hat']}')
        # print(f'tau: {out_h0['tau_hat']}')
        # print(f'beta: {out_h0['beta_hat']}')
        # print(f'mu: {out_h0['mu_hat']}')
        # print(f'mubeta: {out_h0['mubeta_hat']}')
        print(f'rho: {out_h0['rho_hat']}')
        print(f'sigma: {out_h0['sigma_hat']}')


# %%
