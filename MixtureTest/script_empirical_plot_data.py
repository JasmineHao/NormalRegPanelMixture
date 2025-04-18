# %%
from ast import match_case
from MixtureTestFunctions import *
import time
import pyreadr
from numba import jit
import sys
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import numpy as np

T = 3
def simulate_mixture(weights, means, std_devs, n_samples):
    """
    Simulate data from a mixture of normal distributions.

    Parameters:
    - weights (list or np.array): Mixture weights (must sum to 1).
    - means (list or np.array): Means of the normal components.
    - std_devs (list or np.array): Standard deviations of the normal components.
    - n_samples (int): Number of samples to generate.

    Returns:
    - data (np.array): Simulated data from the mixture distribution.
    """
    # Ensure weights sum to 1
    weights = np.array(weights)
    weights /= np.sum(weights)

    # Step 1: Sample which component each point belongs to
    components = np.random.choice(len(weights), size=n_samples, p=weights)

    # Step 2: Generate data from the selected components
    data = np.array([
        np.random.normal(loc=means[k], scale=std_devs[k]) for k in components
    ]).flatten()

    return data

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
df['id'] = pd.to_numeric(df['id'], errors='coerce')
df['year'] = pd.to_numeric(df['year'], errors='coerce')

# Select relevant columns from df_0
df_0 = df_0[['id', 'year','cu' , 'ciiu', 'mi_share', 'mY_share']]

# Perform the merge
df = pd.merge(df[['id','ciiu_3d','year','K','L']], df_0, on=['id', 'year'])
# Artistic color palette (hex codes)
colors = [
    "#FF6F61",  # Coral Red
    "#6B5B95",  # Classic Purple
    "#88B04B",  # Greenery
    "#F7CAC9",  # Rose Quartz
    "#92A8D1",  # Serenity Blue
    "#955251",  # Marsala
    "#B565A7",  # Radiant Orchid
    "#009B77",  # Teal Green
    "#DD4124",  # Fiery Red
    "#45B8AC"   # Aqua
]

# Function to darken a color by reducing its RGB values
def darken_color(hex_color, factor=0.8):
    """
    Darken a hex color by a given factor.
    Parameters:
      - hex_color (str): Hex color string (e.g., "#FF6F61")
      - factor (float): Factor to darken by (0 < factor < 1)
    Returns:
      - str: Darkened hex color
    """
    # Convert hex color to RGB
    rgb = [int(hex_color[i:i+2], 16) for i in (1, 3, 5)]
    # Apply the darkening factor
    darkened_rgb = [max(0, int(c * factor)) for c in rgb]
    # Convert back to hex
    return "#{:02x}{:02x}{:02x}".format(*darkened_rgb)

M_list = [4,5,6]
for (each_code, ind_count) in zip(ind_code, range(len(ind_code))):
    ind_each = df.loc[df['ciiu_3d'] == each_code, :]
    INDNAME = ind_code_dict[each_code]
    ind_each['y'] = np.log(ind_each['mY_share'])
    ind_each['lnk'] = np.log(ind_each['K'])
    ind_each['lnl'] = np.log(ind_each['L'])
    ind_each['m_share'] = ind_each['mi_share'].fillna(0)
    print(INDNAME)
    
    year_list = sorted(ind_each['year'].unique())
    T_cap = max(year_list)

    
    t_start = T_cap - T + 1
    ind_each_t = ind_each[ind_each['year'] >= t_start].dropna()
    ind_each_y = ind_each_t.pivot(index='id', columns='year', values='y').dropna()
    id_list = ind_each_y.index
    ind_each_t = ind_each_t[ind_each_t['id'].isin(id_list)].sort_values(['id', 'year'])
    
    ind_each_xk = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()
    ind_each_xl = (ind_each_t['lnl'] - ind_each_t['lnl'].mean()) / ind_each_t['lnl'].std()
    
    y = ind_each_y.T.to_numpy().astype(np.float64)
    
    # y_ub = np.quantile(y, 1 - epsilon_truncation)
    # y_lb = np.quantile(y, epsilon_truncation)
    # y[y > y_ub] = y_ub
    # y[y < y_lb] = y_lb
    fig = plt.figure(figsize=(6, 4))  # Set figure size (optional)

    num_bins = 100  # You can change this based on your data
    
    # values = np.exp(y.T.flatten())
    values = np.exp(y.T.flatten())
    x = np.linspace(min(values) - 1, max(values) + 1, 1000)

    bin_edges = np.linspace(values.min(), values.max(), num_bins + 1)  
    
    plt.hist(values, bins=bin_edges, density=True, alpha=0.5, label=f"Material Revenue Share", color=colors[0], edgecolor="black")

    plt.legend(loc="upper right")
    plt.title(f"Distribution of the Material Revenue Share for {INDNAME} Industry")
    plt.xlabel("Value")
    plt.ylabel("Frequency")

    # Show the plot
    plt.show()

    # Save the figure as a PNG file after showing the plot
    fig.savefig(f"figure/empirical_plot_{INDNAME}.png", dpi=300)

           
    x_k = ind_each_xk.to_numpy().reshape(-1, 1).astype(np.float64)
    x_kl = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_xl.to_numpy().reshape(-1, 1)].astype(np.float64)
    ind_each_m_share = ind_each_t['m_share']
    x_kmshare = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_m_share.to_numpy().reshape(-1, 1)].astype(np.float64)
        
    N = ind_each_y.shape[0]
    x_0 = np.zeros((N * T, 0))
    z_0 = np.zeros((N * T, 0))
    ciiu_value_count = ind_each_t.groupby('ciiu')['ciiu'].count()
    ciiu_combine = ciiu_value_count[ciiu_value_count < 0.1 * (N * T)].index
    for each_ciiu in ciiu_combine:
        ind_each_t.loc[ind_each_t['ciiu'] == each_ciiu, 'ciiu'] = 'other'
    z = 1 * pd.get_dummies(ind_each_t['ciiu'], prefix='category').values[:, 1:]
    z = z.astype(np.float64)  # Convert z to float64
    
    print(N)
    
    m = M_list[ind_count]
    k = 4
    estimate_params = regpanelmixmixturePMLE(y, x_0, z_0, 0, 0, m, k)
    
    # Create a new figure
    fig = plt.figure(figsize=(6, 4))  # Set figure size (optional)

    num_bins = 50  # You can change this based on your data
    
    # values = np.exp(y.T.flatten())
    values = y.T.flatten()
    x = np.linspace(min(values) - 1, max(values) + 1, 1000)

    bin_edges = np.linspace(values.min(), values.max(), num_bins + 1)  # Define bin edges

    tau_hat = estimate_params['tau_hat'].reshape((m,k))
    mu_hat = estimate_params['mu_hat'].reshape((m,k))
    sigma_hat = estimate_params['sigma_hat']
    for mm in range(m):
        P = np.repeat(estimate_params['post'][:,mm],T)
        P = P / np.sum(P)

        plt.hist(values, bins=bin_edges, density=True, weights= P,alpha=0.5, label=f"Component {mm+1}", color=colors[mm], edgecolor="black")

        weights = tau_hat[mm]
        means  = mu_hat[mm]
        std_devs = np.repeat(sigma_hat[:,mm],k) 
        mixture_density = np.sum([w * (1 / (np.sqrt(2 * np.pi) * s)) * np.exp(-0.5 * ((x - m) / s)**2)
                          for w, m, s in zip(weights, means, std_devs)], axis=0)
        plt.plot(x, mixture_density, color=darken_color(colors[mm], factor=0.8), label=f"Mixture Density Component {mm+1}", linewidth=2)

    # Add legend and labels
    plt.legend(loc="upper left")
    plt.title(f"Distribution of the log Material Revenue Share for {INDNAME} Industry")
    plt.xlabel("Value")
    plt.ylabel("Frequency")

    # Show the plot
    plt.show()

    # Save the figure as a PNG file after showing the plot
    fig.savefig(f"figure/distributions_plot_{INDNAME}.png", dpi=300)


# %%
# [Example] Plot weighted density 
# import numpy as np
# import matplotlib.pyplot as plt

# # Example data
# X = np.array([1, 2, 3, 4, 5])  # Variable values
# P = np.array([0.1, 0.2, 3, 0.25, 0.15])  # Weights (probabilities)

# # Normalize weights if they don't sum to 1
# P = P / np.sum(P)

# # Define the number of bins explicitly
# num_bins = 5  # You can change this based on your data
# bin_edges = np.linspace(X.min(), X.max(), num_bins + 1)  # Define bin edges

# # Plotting the histogram with explicitly defined bins
# plt.hist(X, bins=bin_edges, weights=P, edgecolor='black', alpha=0.7)

# # Add labels
# plt.title("Histogram of Probability-Weighted Distribution")
# plt.xlabel("Variable")
# plt.ylabel("Weighted Frequency")
# plt.show()
# %%
# import numpy as np
# import matplotlib.pyplot as plt

# # Generate example data for 3 distributions
# data1 = np.random.normal(0, 1, 1000)  # Distribution 1: mean=0, std=1
# data2 = np.random.normal(3, 1.5, 1000)  # Distribution 2: mean=3, std=1.5
# data3 = np.random.normal(-2, 0.5, 1000)  # Distribution 3: mean=-2, std=0.5

# # Create a new figure
# fig = plt.figure(figsize=(8, 6))  # Set figure size (optional)

# # Plot histograms for the three distributions
# bins = np.linspace(-6, 8, 40)  # Define common bin edges for all histograms
# plt.hist(data1, bins=bins, alpha=0.5, label="Distribution 1", color="blue", edgecolor="black")
# plt.hist(data2, bins=bins, alpha=0.5, label="Distribution 2", color="green", edgecolor="black")
# plt.hist(data3, bins=bins, alpha=0.5, label="Distribution 3", color="red", edgecolor="black")

# # Add legend and labels
# plt.legend(loc="upper right")
# plt.title("Multiple Distributions in One Plot")
# plt.xlabel("Value")
# plt.ylabel("Frequency")

# # Show the plot
# plt.show()

# Save the figure as a PNG file after showing the plot
# fig.savefig("distributions_plot.png", dpi=300)
# %%
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde

# # Empirical data (x) and corresponding probabilities (P)
# x = np.array([1, 1, 9, 4, 2, 3, 5])  # Observed values (empirical data)
# P = np.array([1, 0.5, 1, 0, 0, 1,0])   # Probabilities for each value

# # Normalize probabilities
# P = P / np.sum(P)

# # Repeat values in x based on probabilities to simulate empirical distribution
# scale_factor = 100  # Scale up the dataset size for better visualization
# weighted_x = np.repeat(x, (P * scale_factor).astype(int))

# # Plot the histogram of empirical data
# bins = np.arange(x.min(), x.max() + 2)  # Bin edges for integers
# plt.figure(figsize=(8, 5))
# plt.hist(weighted_x, bins=bins, density=True, color="skyblue", edgecolor="black", alpha=0.6, label="Histogram (Empirical)")

# # Fit a KDE curve to the empirical data
# kde = gaussian_kde(weighted_x, bw_method=0.3)  # Adjust bandwidth if necessary
# x_vals = np.linspace(x.min() - 1, x.max() + 1, 500)  # Smooth x values for the curve
# kde_vals = kde(x_vals)

# # Overlay the KDE curve
# plt.plot(x_vals, kde_vals, color="red", linewidth=2, label="Empirical Curve (KDE)")

# # Add titles, labels, and grid
# plt.title("Empirical Distribution with Fitted Curve", fontsize=14)
# plt.xlabel("Value", fontsize=12)
# plt.ylabel("Density", fontsize=12)
# plt.legend()
# plt.grid(alpha=0.3)

# # Save the figure (optional)
# plt.savefig("empirical_distribution.png", dpi=300)

# # Show the plot
# plt.show()
# %%

# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters of the mixture model
# weights = [0.3, 0.5, 0.2]  # Mixture weights (must sum to 1)
# means = [0, 5, 10]          # Means of the components
# std_devs = [1, 2, 1.5]      # Standard deviations of the components
# n_samples = 10000           # Number of samples to generate

# # Ensure weights sum to 1
# weights = np.array(weights)
# weights /= np.sum(weights)

# # Step 1: Sample which component each point belongs to
# components = np.random.choice(len(weights), size=n_samples, p=weights)

# # Step 2: Generate data from the selected components
# data = np.array([
#     np.random.normal(loc=means[k], scale=std_devs[k]) for k in components
# ]).flatten()

# # Plot the simulated data
# plt.figure(figsize=(10, 6))
# plt.hist(data, bins=50, density=True, alpha=0.6, color="skyblue", edgecolor="black", label="Simulated Data")

# # Overlay the true mixture density
# x = np.linspace(min(data) - 1, max(data) + 1, 1000)
# mixture_density = np.sum([w * (1 / (np.sqrt(2 * np.pi) * s)) * np.exp(-0.5 * ((x - m) / s)**2)
#                           for w, m, s in zip(weights, means, std_devs)], axis=0)
# plt.plot(x, mixture_density, color="red", label="True Mixture Density", linewidth=2)

# # Add titles and labels
# plt.title("Simulated Data from a Mixture of Normals", fontsize=14)
# plt.xlabel("Value", fontsize=12)
# plt.ylabel("Density", fontsize=12)
# plt.legend()
# plt.grid(alpha=0.3)
# plt.show()
# %%
