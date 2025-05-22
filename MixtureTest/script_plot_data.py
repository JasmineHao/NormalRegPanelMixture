# %%
from ast import match_case
from multiprocessing import process
from MixtureTestFunctions import *
import time
import pyreadr
from numba import jit
import sys
import matplotlib.pyplot as plt
import numpy as np

import warnings
import pandas as pd
warnings.filterwarnings("ignore")

# %%

T = 3
k = 2
ind_code = [311, 381, 321]
M_list = [5,4,5]
# M_list = [2,2,2]

ind_code_dict = load_ind_code_dict()
colors = load_colors()

each_code = ind_code[0]
m = M_list[0]

estimate_parameters = []

for (each_code, m) in zip(ind_code, M_list):

    processed_data = process_chilean_data(each_code=each_code, T=3)
    INDNAME = ind_code_dict[each_code]
    # Access the results
    y = processed_data['y']
    x_0 = processed_data['x_0']
    x_kmshare = processed_data['x_kmshare']
    z = processed_data['z_ciiu']
    z = z.astype(np.float64)  # Convert z to float64
    z_0 = np.zeros((z.shape[0], 0))

    # Plot the empirical distribution
    # --------------------------------- 
    fig = plt.figure(figsize=(6, 4))  # Set figure size (optional)
    num_bins = 100  # You can change this based on your data
    values = np.exp(y.T.flatten())
    plot_x_values = np.linspace(min(values) - 1, max(values) + 1, 1000)
    bin_edges = np.linspace(values.min(), values.max(), num_bins + 1)  

    # Plot the histogram   
    output_path = f"figure/empirical_distribution_{INDNAME}.png"
    plot_empirical_distribution(values, bin_edges, INDNAME, colors, output_path)
    
    # Plot the error of 3-component plain mixture 
    # ---------------------------------
    # Plain mixture model
    estimate_params = regpanelmixmixturePMLE(y, x_0, z_0, 0, 0, m, 3)
    tau_hat, mu_hat, sigma_hat = estimate_params['tau_hat'].reshape((m, 3)), estimate_params['mu_hat'].reshape((m, 3)), estimate_params['sigma_hat']
    output_path = f"figure/empirical_error_plain_{INDNAME}.png"
    plot_mixture_distribution(estimate_params, T, [y.T.flatten()] * m, np.linspace(y.min(), y.max(), 51), tau_hat, mu_hat, sigma_hat, m, 3, colors, INDNAME, np.linspace(y.min() - 1, y.max() + 1, 1000), output_path)

    # Mixture model with kmshare and ciiu
    estimate_params = regpanelmixmixturePMLE(y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], m, 3)
    values_type = [y.T.flatten() - x_kmshare @ estimate_params['beta_hat'][mm] - z @ estimate_params['gamma_hat'][0] for mm in range(m)]
    tau_hat, mu_hat, sigma_hat = estimate_params['tau_hat'].reshape((m, 3)), estimate_params['mu_hat'].reshape((m, 3)), estimate_params['sigma_hat']
    output_path = f"figure/empirical_error_kmshare_ciiu_{INDNAME}.png"
    plot_mixture_distribution(estimate_params, T, values_type, np.linspace(y.min(), y.max(), 51), tau_hat, mu_hat, sigma_hat, m, k, colors, INDNAME, np.linspace(y.min() - 1, y.max() + 1, 1000), output_path)

    
    # Plot Model Parameter Estimates
    # ---------------------------------
    
    def process_model(model_function, y, x, z, p, q, m, k, model_type, specification, output_prefix):
        if model_type in ['stationary_mixture', 'ar1_mixture']:
            model_output = model_function(y, x, z, p, q, m, k)
            if model_type == 'stationary_mixture':
                params_dict, params_array = get_params_stationary_mixture(model_output)
            elif model_type == 'ar1_mixture':
                params_dict, params_array = get_params_ar1_mixture(model_output)
            mubar = (params_dict['tau'].reshape((m, k)) * params_dict['mu'].reshape((m, k))).sum(axis=1)
        else:
            model_output = model_function(y, x, z, p, q, m)
        params_dict, params_array, standard_errors, standard_errors_dict, variable_names = compute_standard_errors(
            model_output, data=[y, x, z], model_type=model_type)
        
        # Bootstrap to compute standard error of mubar and alpha
        bootstrap_samples = 100  # Number of bootstrap samples
        mubar_bootstrap = []
        alpha_bootstrap = []
        mu_bootstrap = []
        beta_bootstrap = []
        gamma_bootstrap = []
        sigma_bootstrap = []
        rho_bootstrap = []
        tau_bootstrap = []
        T, N = y.shape
        for _ in range(bootstrap_samples):
            # Resample data with replacement
            indices = np.random.choice(N, size=N, replace=True)
            y_bootstrap = y[:, indices]
            x_bootstrap = np.zeros(x.shape)
            for qq in range(x.shape[1]):
                x_bootstrap[:, qq] = (x[:, qq].reshape((N, T)).T[:, indices]).T.flatten()
            z_bootstrap = np.zeros(z.shape)
            for pp in range(z.shape[1]):
                z_bootstrap[:, pp] = (z[:, pp].reshape((N, T)).T[:, indices]).T.flatten()

            # Refit the model on the bootstrap sample
            model_output_bootstrap = model_function(y_bootstrap, x_bootstrap, z_bootstrap, p, q, m, k)
            if model_type == 'stationary_mixture':
                params_dict_bootstrap, _ = get_params_stationary_mixture(model_output_bootstrap)
            elif model_type == 'ar1_mixture':
                params_dict_bootstrap, _ = get_params_ar1_mixture(model_output_bootstrap)
            elif model_type == 'stationary_normal':
                params_dict_bootstrap, _ = get_params_stationary_normal(model_output_bootstrap)
            elif model_type == 'ar1_normal':
                params_dict_bootstrap, _ = get_params_ar1_normal(model_output_bootstrap)
            
            # Compute mubar for the bootstrap sample
            alpha_bootstrap_sample = params_dict_bootstrap['alpha']
            alpha_bootstrap.append(alpha_bootstrap_sample)

            mu_bootstrap_sample = params_dict_bootstrap['mu']
            mu_bootstrap.append(mu_bootstrap_sample)

            beta_bootstrap_sample = params_dict_bootstrap.get('beta', None)
            if beta_bootstrap_sample is not None:
                beta_bootstrap.append(beta_bootstrap_sample)

            gamma_bootstrap_sample = params_dict_bootstrap.get('gamma', None)
            if gamma_bootstrap_sample is not None:
                gamma_bootstrap.append(gamma_bootstrap_sample)

            sigma_bootstrap_sample = params_dict_bootstrap['sigma']
            sigma_bootstrap.append(sigma_bootstrap_sample)

            rho_bootstrap_sample = params_dict_bootstrap.get('rho', None)
            if rho_bootstrap_sample is not None:
                rho_bootstrap.append(rho_bootstrap_sample)

            tau_bootstrap_sample = params_dict_bootstrap.get('tau', None)
            if tau_bootstrap_sample is not None:
                tau_bootstrap.append(tau_bootstrap_sample)
                mubar_bootstrap_sample = (params_dict_bootstrap['tau'].reshape((m, k)) * params_dict_bootstrap['mu'].reshape((m, k))).sum(axis=1)
                mubar_bootstrap.append(mubar_bootstrap_sample)
                
        # Compute standard error of mubar

        mubar_se = np.array(mubar_bootstrap).std(axis=0)
        # Compute standard error of alpha
        alpha_se = np.array(alpha_bootstrap).std(axis=0)
        if 'mubar' in locals():
            params_dict['mubar'] = mubar
        
        mu_se = np.array(mu_bootstrap).std(axis=0)
        beta_se = np.array(beta_bootstrap).std(axis=0) if beta_bootstrap else None
        gamma_se = np.array(gamma_bootstrap).std(axis=0) if gamma_bootstrap else None
        sigma_se = np.array(sigma_bootstrap).std(axis=0)
        rho_se = np.array(rho_bootstrap).std(axis=0) if rho_bootstrap else None
        tau_se = np.array(tau_bootstrap).std(axis=0) if tau_bootstrap else None
        # Store standard errors in the dictionary
        standard_errors_dict = {}
        standard_errors_dict['mu'] = mu_se
        standard_errors_dict['beta'] = beta_se
        standard_errors_dict['gamma'] = gamma_se
        standard_errors_dict['sigma'] = sigma_se
        standard_errors_dict['rho'] = rho_se
        standard_errors_dict['tau'] = tau_se
        standard_errors_dict['mubar'] = mubar_se
        standard_errors_dict['alpha'] = alpha_se

        estimate_parameters.append({
            'Industry': INDNAME,
            'Specification': specification,
            'model type': model_type,
            'params': params_dict,
            'standard errors': standard_errors_dict, 
            'variable names': variable_names
        })
        
        # PNAME = "Mean of Material Share"
        # output_path = f"figure/{output_prefix}_mu_with_error_bars_{INDNAME}.png"
        # plot_mubeta_with_error_bars(
        #     params_dict['mu'], standard_errors_dict['mu'], variable_names['mu'], INDNAME, PNAME, output_path, label=r'mu'
        # )
        
        # PNAME = "Mixing Proportion"
        # output_path = f"figure/{output_prefix}_alpha_with_error_bars_{INDNAME}.png"
        # plot_mubeta_with_error_bars(
        #     params_dict['alpha'], np.array([standard_errors_dict['alpha'][0]] * 2), ['alpha_1', 'alpha_2'], INDNAME, PNAME, output_path, ylim=(-0.1, 0.9), label=r'alpha'
        # )

        # if q > 0:
        #     PNAME = "log K coefficient"
        #     output_path = f"figure/{output_prefix}_beta_logK_with_error_bars_{INDNAME}.png"
        #     plot_mubeta_with_error_bars(
        #     params_dict['beta'][:,0], standard_errors_dict['beta'][:,0], ['beta_K_1', 'beta_K_2'], INDNAME, PNAME, output_path, label=r'beta_K'
        #     )

        #     PNAME = "import coefficient"
        #     output_path = f"figure/{output_prefix}_beta_import_with_error_bars_{INDNAME}.png"
        #     plot_mubeta_with_error_bars(
        #     params_dict['beta'][:,1], standard_errors_dict['beta'][:,1], ['beta_im_1', 'beta_im_2'], INDNAME, PNAME, output_path, label=r'beta_import'
        #     )
        # if model_type in ['ar1_mixture', 'ar1_normal']:
        #     PNAME = "AR1 coefficient"
        #     output_path = f"figure/{output_prefix}_rho_with_error_bars_{INDNAME}.png"
        #     plot_mubeta_with_error_bars(
        #         params_dict['rho'], standard_errors_dict['rho'], ['rho_1', 'rho_2'], INDNAME, PNAME, output_path, ylim=(-0.5, 1.2), label=r'rho'
        #     )

    # Stat Mixture
    # ---------------------------------
    # Process stationary mixture model
    process_model(
        regpanelmixmixturePMLE, y, x_0, z_0, 0, 0, 3, k, 'stationary_mixture', 'Stationary Mixture', 'empirical_stat_mixture'
    )
    
    # Process stationary mixture model with kmshare and ciiu
    process_model(
        regpanelmixmixturePMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 3, k, 'stationary_mixture', 'Stationary Mixture kmshare ciiu', 'empirical_stat_mixture_kmshare_ciiu'
    )

    # Stat Normal
    # ---------------------------------
    process_model(
        regpanelmixPMLE, y, x_0, z_0, 0, 0, 3, k, 'stationary_normal', 'Stationary Normal', 'empirical_stat_normal'
    )
    
    # Process stationary mixture model with kmshare and ciiu
    process_model(
        regpanelmixPMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 3, k, 'stationary_normal', 'Stationary Normal kmshare ciiu', 'empirical_stat_normal_kmshare_ciiu'
    )

    # AR1 Normal
    # ---------------------------------
    process_model(
        regpanelmixAR1PMLE, y, x_0, z_0, 0, 0, 3, k, 'ar1_normal', 'AR1 Normal', 'empirical_ar1_normal'
    )
    
    process_model(
        regpanelmixAR1PMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 3, k, 'ar1_normal', 'AR1 Normal kmshare ciiu', 'empirical_ar1_normal_kmshare_ciiu'
    )

    # AR1 Mixture
    # ---------------------------------
    process_model(
        regpanelmixAR1mixturePMLE, y, x_0, z_0, 0, 0, 3, k, 'ar1_mixture', 'AR1 Mixture', 'empirical_ar1_mixture'
    )
    
    if ind_code != 311:
        process_model(
            regpanelmixAR1mixturePMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 3, k, 'ar1_mixture', 'AR1 Mixture kmshare ciiu', 'empirical_ar1_mixture_kmshare_ciiu'
        )

# %%


def process_estimated_parameters(estimate_parameters):

    # Initialize lists to store data for the table
    industries = []
    specifications = []
    model_types = []
    alphas = []
    mus = []
    rhos = []
    betas = []
    sigmas = []
    mubars = []
    alphas_se = []
    mus_se = []
    rhos_se = []
    betas_se = []
    sigmas_se = []
    mubars_se = []

    # Iterate through the estimated parameters
    for param in estimate_parameters:
        industries.append(param['Industry'])
        specifications.append(param['Specification'])
        model_types.append(param['model type'])
        alphas.append(np.round(param['params'].get('alpha', 0.0), 4))
        mus.append(np.round(param['params'].get('mu', 0.0), 4))
        rhos.append(np.round(param['params'].get('rho', 0.0), 4))
        betas.append(np.round(param['params'].get('beta', 0.0), 4))
        sigmas.append(np.round(param['params'].get('sigma', 0.0), 4))
        mubars.append(np.round(param['params'].get('mubar', 0.0), 4))
        alphas_se.append(np.round(param['standard errors'].get('alpha', 0.0), 4))
        mus_se.append(np.round(param['standard errors'].get('mu', 0.0), 4))
        rho_se_value = param['standard errors'].get('rho', 0.0)
        rhos_se.append(np.round(rho_se_value if rho_se_value is not None else 0.0, 4))
        betas_se.append(np.round(param['standard errors'].get('beta', 0.0), 4))
        sigmas_se.append(np.round(param['standard errors'].get('sigma', 0.0), 4))
        mubars_se.append(np.round(param['standard errors'].get('mubar', 0.0), 4))

    # Create a DataFrame to tabulate the results
    df = pd.DataFrame({
        'Industry': industries,
        'Specification': specifications,
        'Model Type': model_types,
        'Alpha': alphas,
        'Alpha SE': alphas_se,
        'Mu': mus,
        'Mu SE': mus_se,
        'Rho': rhos,
        'Rho SE': rhos_se,
        'Beta': betas,
        'Beta SE': betas_se,
        'Sigma': sigmas,
        'Sigma SE': sigmas_se,
        'Mubar': mubars,
        'Mubar SE': mubars_se
    })

    # Save the table to a CSV file
    output_path = "estimated_parameters_table_bootstrap_m3.csv"
    df.to_csv(output_path, index=False)
    print(f"Estimated parameters table saved to {output_path}")

    return df

# Call the function to process and tabulate the estimated parameters
estimated_parameters_table = process_estimated_parameters(estimate_parameters)

# Filter out entries where Specification is 'AR1 Mixture kmshare ciiu' and Industry is 'Food products'
estimate_parameters = [
    entry for entry in estimate_parameters 
    if not (entry['Specification'] == 'AR1 Mixture kmshare ciiu' and entry['Industry'] == 'Food products')
]

# %%
# Restore estimated_parameters from the CSV table
def str_to_array(s):
    if pd.isnull(s):
        return None
    s = str(s).strip()
    if s.startswith('[') and s.endswith(']'):
        return np.fromstring(s[1:-1], sep=' ')
    return np.array(eval(s))

def str_to_2d_array(s):
    """
    Convert a string like '[[a b]\n [c d]]' to a numpy 2D array.
    """
    if pd.isnull(s):
        return None
    s = str(s).strip()
    if s.startswith('[[') and s.endswith(']]'):
        # Remove newlines and extra spaces
        s_clean = s.replace('\n', ' ')
        # Remove double brackets for np.fromstring
        s_clean = s_clean.replace('[', '').replace(']', '')
        arr = np.fromstring(s_clean, sep=' ')
        # Infer shape: count rows by number of '[' in original string
        n_rows = s.count('[') - 1
        n_cols = int(len(arr) / n_rows) if n_rows > 0 else len(arr)
        return arr.reshape((n_rows, n_cols))
    return np.array(eval(s))

def restore_estimated_parameters_from_table(csv_path):
    df = pd.read_csv(csv_path)
    estimate_parameters_restored = []
    for _, row in df.iterrows():
        # Convert string like '[0.143 0.857]' to numpy array
        
        entry = {
            'Industry': row['Industry'],
            'Specification': row['Specification'],
            'model type': row['Model Type'],
            'params': {
                'alpha': str_to_array(row['Alpha']),
                'mu': str_to_array(row['Mu']),
                'rho': str_to_array(row['Rho']),
                'beta': str_to_2d_array(row['Beta']),
                'sigma': str_to_array(row['Sigma']),
                'mubar': str_to_array(row['Mubar']),
            },
            'standard errors': {
                'alpha': str_to_array(row['Alpha SE']),
                'mu': str_to_array(row['Mu SE']),
                'rho': str_to_array(row['Rho SE']),
                'beta': str_to_2d_array(row['Beta SE']),
                'sigma': str_to_array(row['Sigma SE']),
                'mubar': str_to_array(row['Mubar SE']),
            },
            'variable names': {}  # Not available in CSV, can be filled if needed
        }
        estimate_parameters_restored.append(entry)
    return estimate_parameters_restored



# %%
estimate_parameters_restored = restore_estimated_parameters_from_table("estimated_parameters_table_m2.csv")

estimate_parameters_restored = [
    entry for entry in estimate_parameters_restored 
    if not (entry['Specification'] == 'AR1 Mixture kmshare ciiu' and entry['Industry'] == 'Food products')
]

# Iterate over different model types
for model_type in set([entry['model type'] for entry in estimate_parameters_restored]):
    # Get unique specifications for the current model type
    unique_specifications = list(set([entry['Specification'] for entry in estimate_parameters_restored if entry['model type'] == model_type]))
        
    for specification in unique_specifications:
        # Filter entries for the current model type and specification
        filtered_entries = [entry for entry in estimate_parameters_restored if entry['model type'] == model_type and entry['Specification'] == specification]
        # Generate a sanitized specification name
        sanitized_specification = specification.replace(" ", "_").lower()

        # Plot alpha
        output_path_alpha = f'figure/{sanitized_specification}_alpha_across_industries_m2.png'
        plot_parameter_across_industries(
            'alpha', filtered_entries, ylabel=r"$\hat{\alpha}$", 
            title=f'Alpha Across Industries ({model_type}, {specification})', 
            output_path=output_path_alpha
        )
        print(output_path_alpha)

        # Plot mu
        output_path_mu = f'figure/{sanitized_specification}_mu_across_industries_m2.png'
        plot_parameter_across_industries(
            'mu', filtered_entries, ylabel=r"$\hat{\mu}$",
            title=f'Mu Across Industries ({model_type}, {specification})', 
            output_path=output_path_mu
        )
        print(output_path_mu)

        # Plot sigma
        output_path_sigma = f'figure/{sanitized_specification}_sigma_across_industries_m2.png'
        plot_parameter_across_industries(
            'sigma', filtered_entries, ylabel=r"$\hat{\sigma}$", 
            title=f'Sigma Across Industries ({model_type}, {specification})', 
            output_path=output_path_sigma
        )
        print(output_path_sigma)

        # Plot rho if it exists
        if ('rho' in filtered_entries[0]['params']):
            if np.all(filtered_entries[0]['params']['rho'] != np.array(0.)):
                output_path_rho = f'figure/{sanitized_specification}_rho_across_industries_m2.png'
                plot_parameter_across_industries(
                'rho', filtered_entries, ylabel=r"$\hat{\rho}$", 
                title=f'Rho Across Industries ({model_type}, {specification})', 
                output_path=output_path_rho
                )
                print(output_path_rho)
            # Plot mubar if it exists
        if ('mubar' in filtered_entries[0]['params']) & (len(filtered_entries[0]['params']['mubar'].shape) > 0):
            output_path_mubar = f'figure/{sanitized_specification}_mubar_across_industries_m2.png'
            plot_parameter_across_industries(
            'mubar', filtered_entries, ylabel=r'$\hat{\mu}$', 
            title=f'Mubar Across Industries ({model_type}, {specification})', 
            output_path=output_path_mubar
            )
            print(output_path_mubar)
            
        # Plot beta
        if 'beta' in filtered_entries[0]['params'] and len(filtered_entries[0]['params']['beta']) > 0:
            # Plot beta_k
            output_path_beta_k = f'figure/{sanitized_specification}_beta_k_across_industries_m2.png'
            plot_parameter_across_industries(
            'beta_k', filtered_entries, ylabel=r"$\hat{\beta}_{\log K}$", 
            title=f'Beta K Across Industries ({model_type}, {specification})', 
            output_path=output_path_beta_k
            )
            print(output_path_beta_k)

            # Plot beta_im
            output_path_beta_im = f'figure/{sanitized_specification}_beta_im_across_industries_m2.png'
            plot_parameter_across_industries(
            'beta_im', filtered_entries, ylabel=r"$\hat{\beta}_{import}$",  
            title=f'Beta Import Across Industries ({model_type}, {specification})', 
            output_path=output_path_beta_im
            )
            print(output_path_beta_im)

# %%
len(estimate_parameters)
# %%

