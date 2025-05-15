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
warnings.filterwarnings("ignore")


# %%
T = 3
k = 2
ind_code = [311, 381, 321]
M_list = [5,4,6]
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
    estimate_params = regpanelmixmixturePMLE(y, x_0, z_0, 0, 0, m, k)
    tau_hat, mu_hat, sigma_hat = estimate_params['tau_hat'].reshape((m, k)), estimate_params['mu_hat'].reshape((m, k)), estimate_params['sigma_hat']
    output_path = f"figure/empirical_error_plain_{INDNAME}.png"
    plot_mixture_distribution(estimate_params, T, [y.T.flatten()] * m, np.linspace(y.min(), y.max(), 51), tau_hat, mu_hat, sigma_hat, m, k, colors, INDNAME, np.linspace(y.min() - 1, y.max() + 1, 1000), output_path)

    # Mixture model with kmshare and ciiu
    estimate_params = regpanelmixmixturePMLE(y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], m, k)
    values_type = [y.T.flatten() - x_kmshare @ estimate_params['beta_hat'][mm] - z @ estimate_params['gamma_hat'][0] for mm in range(m)]
    tau_hat, mu_hat, sigma_hat = estimate_params['tau_hat'].reshape((m, k)), estimate_params['mu_hat'].reshape((m, k)), estimate_params['sigma_hat']
    output_path = f"figure/empirical_error_kmshare_ciiu_{INDNAME}.png"
    plot_mixture_distribution(estimate_params, T, values_type, np.linspace(y.min(), y.max(), 51), tau_hat, mu_hat, sigma_hat, m, k, colors, INDNAME, np.linspace(y.min() - 1, y.max() + 1, 1000), output_path)


    # Plot Model Parameter Estimates
    # ---------------------------------
    
    def process_model(model_function, y, x, z, p, q, m, k, model_type, specification, output_prefix):
        if model_type in ['stationary_mixture', 'ar1_mixture']:
            model_output = model_function(y, x, z, p, q, m, k)
        else:
            model_output = model_function(y, x, z, p, q, m)
        params_dict, params_array, standard_errors, standard_errors_dict, variable_names = compute_standard_errors(
            model_output, data=[y, x, z], model_type=model_type
        )
        
        estimate_parameters.append({
            'Industry': INDNAME,
            'Specification': specification,
            'model type': model_type,
            'params': params_dict,
            'standard errors': standard_errors_dict
        })
        
        PNAME = "Mean of Material Share"
        output_path = f"figure/{output_prefix}_mu_with_error_bars_{INDNAME}.png"
        plot_mubeta_with_error_bars(
            params_dict['mu'], standard_errors_dict['mu'], variable_names['mu'], INDNAME, PNAME, output_path, label=r'mu'
        )
        
        PNAME = "Mixing Proportion"
        output_path = f"figure/{output_prefix}_alpha_with_error_bars_{INDNAME}.png"
        plot_mubeta_with_error_bars(
            params_dict['alpha'], np.array([standard_errors_dict['alpha'][0]] * 2), ['alpha_1', 'alpha_2'], INDNAME, PNAME, output_path, ylim=(-0.1, 0.9), label=r'alpha'
        )

        if q > 0:
            PNAME = "log K coefficient"
            output_path = f"figure/{output_prefix}_beta_logK_with_error_bars_{INDNAME}.png"
            plot_mubeta_with_error_bars(
            params_dict['beta'][:,0], standard_errors_dict['beta'][:,0], ['beta_K_1', 'beta_K_2'], INDNAME, PNAME, output_path, label=r'beta_K'
            )

            PNAME = "import coefficient"
            output_path = f"figure/{output_prefix}_beta_import_with_error_bars_{INDNAME}.png"
            plot_mubeta_with_error_bars(
            params_dict['beta'][:,1], standard_errors_dict['beta'][:,1], ['beta_im_1', 'beta_im_2'], INDNAME, PNAME, output_path, label=r'beta_import'
            )
        if model_type in ['ar1_mixture', 'ar1_normal']:
            PNAME = "AR1 coefficient"
            output_path = f"figure/{output_prefix}_rho_with_error_bars_{INDNAME}.png"
            plot_mubeta_with_error_bars(
                params_dict['rho'], standard_errors_dict['rho'], ['rho_1', 'rho_2'], INDNAME, PNAME, output_path, ylim=(-0.5, 1.2), label=r'rho'
            )

    # Stat Mixture
    # ---------------------------------
    # Process stationary mixture model
    process_model(
        regpanelmixmixturePMLE, y, x_0, z_0, 0, 0, 2, k, 'stationary_mixture', 'Stationary Mixture', 'empirical_stat_mixture'
    )
    
    # Process stationary mixture model with kmshare and ciiu
    process_model(
        regpanelmixmixturePMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 2, k, 'stationary_mixture', 'Stationary Mixture kmshare ciiu', 'empirical_stat_mixture_kmshare_ciiu'
    )

    # Stat Normal
    # ---------------------------------
    process_model(
        regpanelmixPMLE, y, x_0, z_0, 0, 0, 2, k, 'stationary_normal', 'Stationary Normal', 'empirical_stat_normal'
    )
    
    # Process stationary mixture model with kmshare and ciiu
    process_model(
        regpanelmixPMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 2, k, 'stationary_normal', 'Stationary Mixture kmshare ciiu', 'empirical_stat_normal_kmshare_ciiu'
    )

    # AR1 Normal
    # ---------------------------------
    process_model(
        regpanelmixAR1PMLE, y, x_0, z_0, 0, 0, 2, None, 'ar1_normal', 'AR1 Normal', 'empirical_ar1_normal'
    )
    
    process_model(
        regpanelmixAR1PMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 2, None, 'ar1_normal', 'AR1 Normal kmshare ciiu', 'empirical_ar1_normal_kmshare_ciiu'
    )

    # AR1 Mixture
    # ---------------------------------
    process_model(
        regpanelmixAR1mixturePMLE, y, x_0, z_0, 0, 0, 2, 2, 'ar1_mixture', 'AR1 Mixture', 'empirical_ar1_mixture'
    )
    
    process_model(
        regpanelmixAR1mixturePMLE, y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], 2, 2 , 'ar1_mixture', 'AR1 Mixture kmshare ciiu', 'empirical_ar1_mixture_kmshare_ciiu'
    )

# %%
