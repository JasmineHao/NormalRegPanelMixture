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
k = 3
ind_code = [311, 381, 321]
M_list = [4,5,6]

ind_code_dict = load_ind_code_dict()
colors = load_colors()

each_code = ind_code[0]
m = M_list[0]

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
    estimate_params = regpanelmixmixturePMLE(y, x_0, z_0, 0, 0, m, k)
    
    fig = plt.figure(figsize=(6, 4))  # Set figure size (optional)

    num_bins = 50  # You can change this based on your data
    
    # values = np.exp(y.T.flatten())
    values = y.T.flatten()
    plot_x_values = np.linspace(min(values) - 1, max(values) + 1, 1000)

    bin_edges = np.linspace(values.min(), values.max(), num_bins + 1)  # Define bin edges

    tau_hat = estimate_params['tau_hat'].reshape((m,k))
    mu_hat = estimate_params['mu_hat'].reshape((m,k))
    sigma_hat = estimate_params['sigma_hat']
    
    # Call the function
    values_type = [values for mm in range(m)]
    output_path = f"figure/empirical_error_plain_{INDNAME}.png"
    plot_mixture_distribution(estimate_params, T, values_type, bin_edges, tau_hat, mu_hat, sigma_hat, m, k, colors, INDNAME, plot_x_values, output_path)
    # sample of returning the estimates and standard errors

    
    # Plot the error of 3-component mixture kmshare ciiu 
    # ---------------------------------
    estimate_params = regpanelmixmixturePMLE(y, x_kmshare, z, z.shape[1], x_kmshare.shape[1], m, k)
    
    fig = plt.figure(figsize=(6, 4))  # Set figure size (optional)

    num_bins = 50  # You can change this based on your data
    
    # values = np.exp(y.T.flatten())
    values = y.T.flatten()
    plot_x_values = np.linspace(min(values) - 1, max(values) + 1, 1000)

    bin_edges = np.linspace(values.min(), values.max(), num_bins + 1)  # Define bin edges

    values_type = [values for mm in range(m)]
    for mm in range(m):
        # mm = 0
        values_type[mm] = values_type[mm] - x_kmshare @ estimate_params['beta_hat'][mm] - z @ estimate_params['gamma_hat'][0]
    tau_hat = estimate_params['tau_hat'].reshape((m,k))
    mu_hat = estimate_params['mu_hat'].reshape((m,k))
    sigma_hat = estimate_params['sigma_hat']
        
    # Call the function
    output_path = f"figure/empirical_error_kmshare_ciiu_{INDNAME}.png"
    plot_mixture_distribution(estimate_params, T, values_type, bin_edges, tau_hat, mu_hat, sigma_hat, m, k, colors, INDNAME, plot_x_values, output_path)
    # sample of returning the estimates and standard errors
    
    
    # Stat Mixture
    # ---------------------------------
    model_output = regpanelmixmixturePMLE(y, x_0, z_0, p=0, q=0, m=m, k=k)
    
    params_dict, params_array, standard_errors, standard_errors_dict, variable_names = compute_standard_errors(model_output, data=[y,x_0,z_0], model_type='stationary_mixture')
    
    
    PNAME = "Mean of Material Share"
    output_path = f"figure/empirical_mu_with_error_bars_stat_mixture_{INDNAME}.png"
    
    values=params_dict['mu']
    errors=standard_errors_dict['mu']
    types=variable_names['mu']
    plot_mubeta_with_error_bars(params_dict['mu'], standard_errors_dict['mu'], variable_names['mu'], INDNAME, PNAME, output_path)
    
    
    PARAMETER_NAME = "Mixing Proportion"
    values=params_dict['alpha']
    errors=standard_errors_dict['alpha']
    types=variable_names['alpha']
    output_path = f"figure/empirical_alpha_with_error_bars_stat_mixture_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['alpha'][:-1], standard_errors_dict['alpha'][:-1], variable_names['alpha'], INDNAME, PNAME, output_path, ylim=(-1,2))
    
    # Stat Normal
    # ---------------------------------
    
    model_output = regpanelmixPMLE(y, x_0, z_0, p=0, q=0, m=m)
    params_dict, params_array, standard_errors, standard_errors_dict, variable_names = compute_standard_errors(model_output, data=[y,x_0,z_0], model_type='stationary_normal')
    
    
    PNAME = "Mean of Material Share"
    output_path = f"figure/empirical_mu_with_error_bars_stat_normal_{INDNAME}.png"
    
    values=params_dict['mu']
    errors=standard_errors_dict['mu']
    types=variable_names['mu']
    plot_mubeta_with_error_bars(params_dict['mu'], standard_errors_dict['mu'], variable_names['mu'], INDNAME, PNAME, output_path)
    
    
    PARAMETER_NAME = "Mixing Proportion"
    values=params_dict['alpha']
    errors=standard_errors_dict['alpha']
    types=variable_names['alpha']
    output_path = f"figure/empirical_alpha_with_error_bars_stat_normal_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['alpha'][:-1], standard_errors_dict['alpha'][:-1], variable_names['alpha'], INDNAME, PNAME, output_path, ylim=(-1,2))
    
    # AR1 Normal
    # ---------------------------------
    
    model_output = regpanelmixAR1PMLE(y, x_0, z_0, p=0, q=0, m=m)
    params_dict, params_array, standard_errors, standard_errors_dict, variable_names = compute_standard_errors(model_output, data=[y,x_0,z_0], model_type='ar1_normal')   
    
    PNAME = "Mean of Material Share"
    values=params_dict['mu']
    errors=standard_errors_dict['mu']
    types=variable_names['mu']
    output_path = f"figure/empirical_mu_with_error_bars_ar1_normal_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['mu'], standard_errors_dict['mu'], variable_names['mu'], INDNAME, PNAME, output_path)

    
    PARAMETER_NAME = "Persistence of AR1 Shocks"
    values=params_dict['rho']
    errors=standard_errors_dict['rho']
    types=variable_names['rho']
    output_path = f"figure/empirical_rho_with_error_bars_ar1_normal_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['rho'], standard_errors_dict['rho'], variable_names['rho'], INDNAME, PNAME, output_path, ylim=(-1,2))

    
    PARAMETER_NAME = "Mixing Proportion"
    values=params_dict['alpha']
    errors=standard_errors_dict['alpha']
    types=variable_names['alpha']
    output_path = f"figure/empirical_alpha_with_error_bars_ar1_normal_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['alpha'][:-1], standard_errors_dict['alpha'][:-1], variable_names['alpha'], INDNAME, PNAME, output_path, ylim=(-1,2))
 
    # AR1 Mixture
    # ---------------------------------
    model_output = regpanelmixAR1mixturePMLE(y, x_0, z_0, p=0, q=0, m=2, k=k)
    
    params_dict, params_array, standard_errors, standard_errors_dict, variable_names = compute_standard_errors(model_output, data=[y,x_0,z_0], model_type='ar1_mixture')
    
    PNAME = "Mean of Material Share"
    values=params_dict['mu']
    errors=standard_errors_dict['mu']
    types=variable_names['mu']
    output_path = f"figure/empirical_mu_with_error_bars_ar1_mixture_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['mu'], standard_errors_dict['mu'], variable_names['mu'], INDNAME, PNAME, output_path)

 
    PARAMETER_NAME = "Persistence of AR1 Shocks"
    values=params_dict['rho']
    errors=standard_errors_dict['rho']
    types=variable_names['rho']
    output_path = f"figure/empirical_rho_with_error_bars_ar1_mixture_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['rho'], standard_errors_dict['rho'], variable_names['rho'], INDNAME, PNAME, output_path, ylim=(-1,2))

    
    PARAMETER_NAME = "Mixing Proportion"
    values=params_dict['alpha']
    errors=standard_errors_dict['alpha']
    types=variable_names['alpha']
    output_path = f"figure/empirical_alpha_with_error_bars_ar1_mixture_{INDNAME}.png"
    plot_mubeta_with_error_bars(params_dict['alpha'][:-1], standard_errors_dict['alpha'][:-1], variable_names['alpha'][:-1], INDNAME, PNAME, output_path, ylim=(-1,2))

    
# %%
