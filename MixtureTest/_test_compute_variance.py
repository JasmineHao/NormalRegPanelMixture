# %%
import numpy as np
from scipy.stats import norm

# import jax
# import jax.numpy as jnp
# from jax import grad, jacfwd, jacrev
from functools import partial

import numpy as np
from scipy.optimize import approx_fprime
from scipy.optimize._numdiff import approx_derivative

# def log_likelihood_array(y, mu, sigma):
#     """
#     Calculate the log-likelihood of each element in the data under a normal distribution.

#     Parameters:
#     - y: array-like, the observed data points
#     - mu: float, the mean of the normal distribution
#     - sigma: float, the standard deviation of the normal distribution

#     Returns:
#     - log_likelihoods: array, log-likelihood of each data point
#     """
#     # Precompute constants
#     constant = -0.5 * np.log(2 * np.pi)
#     variance = sigma ** 2

#     # Compute log-likelihood for each data point
#     log_likelihoods = np.array([
#         constant - np.log(sigma) - ((yi - mu) ** 2) / (2 * variance)
#         for yi in y
#     ])
    
#     return log_likelihoods

# ---'
def flatten_params_stationary_normal(params):
    """Convert parameter dictionary to flat array"""
    return np.concatenate([params['alpha'],
                          params['sigma'],
                          params['mubeta'],
                          params['gamma']])

def unflatten_params_stationary_normal(arr, original_shape):
    """Convert flat array back to parameter dictionary"""
    m = original_shape['alpha'][0]
    q_plus_1 = original_shape['mubeta'][0] // m
    p = original_shape['gamma'][0]
    
    alpha = arr[:m]
    sigma = arr[m:2*m]
    mubeta = arr[2*m:2*m+m*(q_plus_1)]
    gamma = arr[2*m+m*(q_plus_1):]
    
    return {
        'alpha': alpha,
        'sigma': sigma,
        'mubeta': mubeta,
        'gamma': gamma
    }'


# _y = np.random.normal(size=100)
# norm.logpdf(_y, 0, 1 )
# log_likelihood_array(_y, 0, 1)


def convert_out_h0_to_flatten_params_stationary_normal(out_h0):
    """
    Convert out_h0 dictionary into arrays compatible with flatten_params.
    
    Parameters:
    -----------
    out_h0 : dict
        Dictionary containing parameter arrays.
    
    Returns:
    --------
    dict
        Dictionary with flattened parameter arrays.
    """
    alpha_hat = out_h0['alpha_hat'][0]
    sigma_hat = out_h0['sigma_hat'][0]
    mubeta_hat = out_h0['mubeta_hat'][0]
    gamma_hat = out_h0['gamma_hat'][0]
    p = len(gamma_hat)
    m = len(alpha_hat)
    q = int((len(mubeta_hat) / m )- 1)
    
    mubeta_hat = np.ascontiguousarray(mubeta_hat)
    mubeta_hat_mat = mubeta_hat.reshape((q + 1, m)).T
    
    return {
        'alpha': alpha_hat,
        'sigma': sigma_hat,
        'mubeta': mubeta_hat.T.flatten(),
        'gamma': gamma_hat
    }

# %%
# Define your likelihood function




# def hessian_stationary_normal(data, params):
#     y, x, z = data
#     t, n = y.shape
#     nt = n * t
#     q = x.shape[1] if x.ndim > 1 else 1
#     p = z.shape[1] if z.ndim > 1 else 0
    
#     # Unpack parameters (alpha has m-1 free parameters)
#     alpha = params['alpha']  # First m-1 components (last is 1-sum(alpha[:-1]))
#     sigma = params['sigma']
#     mubeta = params['mubeta']
#     gamma = params['gamma']
#     m = len(alpha)  # Note: alpha includes m-1 free params + constrained last
    
    
#     # Reshape mubeta to (m, q+1) matrix
#     mubeta_mat = np.ascontiguousarray(mubeta).reshape((m, q+1))
    
#     # Compute posterior weights w using full alpha
#     _, w, _ = log_likelihood_stationary_normal(params, data)
#     w = w.T  # Shape (n, m)
    
#     # Prepare design matrix with intercept
#     x1 = np.column_stack((np.ones(nt), x)) if q > 0 else np.ones((nt, 1))
    
#     # Flatten y and compute residuals
#     y_c = y.T.flatten()
#     residuals = np.zeros((m, nt))
#     for j in range(m):
#         residuals[j] = y_c - z @ gamma - x1 @ mubeta_mat[j]
    
#     # Parameter counts (m-1 alpha, m*(q+1) mubeta, m sigma, p gamma)
#     num_alpha = m - 1
#     num_mubeta = m * (q + 1)
#     num_sigma = m
#     num_gamma = p
#     num_params = num_alpha + num_mubeta + num_sigma + num_gamma
    
#     # Initialize Hessian
#     H = np.zeros((n, num_params, num_params))
    
    
#     # ======================================================
#     # 1. Compute ∂w_j/∂θ terms (needed for Hessian), should be m x nparam
#     # ======================================================
        
#     # Compute component scores S_j = ∂logf_j/∂θ
#     S_alpha = np.zeros((m, num_alpha))  # ∂logf_j/∂α_k for k=1..m-1
#     for k in range(num_alpha):
#     S_alpha[:, k] = (np.eye(m)[:, k]/alpha[k] - 
#                     np.eye(m)[:, -1]/alpha[-1])
    
#     S_mubeta = np.zeros((n, m, num_mubeta))
#     S_sigma = np.zeros((n, m, num_sigma))
#     S_gamma = np.zeros((n, m, num_gamma))
    
#     for j in range(m):
#         # For mubeta parameters
#         for k in range(q+1):
#             idx = j*(q+1) + k
#             S_mubeta[:, j, idx] = (residuals[j] * x1[:, k] / sigma[j]**2).reshape(n, t).sum(axis=1)
        
#         # For sigma parameters
#         S_sigma[:, j, j] = (-1/sigma[j] + (residuals[j]**2 / sigma[j]**3)).reshape(n, t).sum(axis=1)
        
#         # For gamma parameters
#         for pp in range(p):
#             S_gamma[:, j, pp] = (residuals[j] * z[:, pp] / sigma[j]**2).sum()
    
#     # Combine all scores
#     S_all = np.hstack([S_alpha, S_mubeta, S_sigma, S_gamma])
    
#     # Compute ∂w_j/∂θ = w_j(S_j - ∑_k w_k S_k)
#     wS_mean = np.sum(w[:, :, None] * S_all[None, :, :], axis=1)  # (n, num_params)
#     dw_dtheta = w[:, :, None] * (S_all[None, :, :] - wS_mean[:, None, :])  # (n, m, num_params)
        
#     # ======================================================
#     # 2. Build Hessian blocks with α constraint
#     # ======================================================
    
#     # Indices for different parameter blocks
#     alpha_idx = slice(0, num_alpha)
#     mubeta_idx = slice(num_alpha, num_alpha + num_mubeta)
#     sigma_idx = slice(num_alpha + num_mubeta, num_alpha + num_mubeta + num_sigma)
#     gamma_idx = slice(num_alpha + num_mubeta + num_sigma, num_params)
    
#     # (A) α-α block: ∂²ℓ/∂α_j∂α_k 
#     for j in range(num_alpha):
#         for k in range(num_alpha):
#             term1 = -w[:, j]/alpha[j] * w[:, k]/alpha[k] - w[:, -1]/alpha[-1]**2
#             term2 = (j == k) * w[:, j]/alpha[j]**2
#             term3 = dw_dtheta[:, j, k] / alpha[j] - dw_dtheta[:, -1, k] / alpha[-1]
#             H[j, k] = np.sum(term1 + term2 + term3)
    
#     # (B) μβ-μβ block (same as before)
#     for j in range(m):
#         for k in range(m):
#             j_start = num_alpha + j*(q+1)
#             k_start = num_alpha + k*(q+1)
            
#             for jj in range(q+1):
#                 for kk in range(q+1):
#                     idx_j = j_start + jj
#                     idx_k = k_start + kk
                    
#                     if j == k:
#                         term1 = -w[:, j] * (x1[:, jj] * x1[:, kk] / sigma[j]**2).reshape(n, t).sum(axis=1)
#                     else:
#                         term1 = 0
                    
#                     term2 = dw_dtheta[:, j, num_alpha + k*(q+1) + kk] * \
#                            (residuals[j] * x1[:, jj] / sigma[j]**2).reshape(n, t).sum(axis=1)
                    
#                     H[idx_j, idx_k] = np.sum(term1 + term2)
    
#     # (C) σ-σ block (same as before)
#     for j in range(m):
#         for k in range(m):
#             if j == k:
#                 term1 = w[:, j] * (1/sigma[j]**2 - 3*(residuals[j]**2)/sigma[j]**4).reshape(n, t).sum(axis=1)
#             else:
#                 term1 = 0
                
#             term2 = dw_dtheta[:, j, sigma_idx][:, k] * \
#                    (-1/sigma[j] + (residuals[j]**2)/sigma[j]**3).reshape(n, t).sum(axis=1)
            
#             H[sigma_idx][j, k] = np.sum(term1 + term2)
#     # (D) γ-γ block (same as before)
#     for p1 in range(p):
#         for p2 in range(p):
#             term1 = 0
#             term2 = 0
#             for j in range(m):
#                 term1 += -w[:, j] * (z[:, p1] * z[:, p2] / sigma[j]**2).reshape(n, t).sum(axis=1)
#                 term2 += dw_dtheta[:, j, gamma_idx][:, p2] * \
#                         (residuals[j] * z[:, p1] / sigma[j]**2).reshape(n, t).sum(axis=1)
            
#             H[gamma_idx][p1, p2] = np.sum(term1 + term2)
    
#     # (E) Cross blocks (α-μβ, α-σ, α-γ)
#     # α-μβ block
#     for j in range(num_alpha):
#         for k in range(m):
#             for kk in range(q+1):
#                 col = num_alpha + k*(q+1) + kk
#                 term = (dw_dtheta[:, j, col] - 
#                        w[:, j] * S_all[:, k*(q+1)+kk] + 
#                        w[:, j] * wS_mean[:, col])
#                 H[j, col] = np.sum(term * (residuals[k] * x1[:, kk] / sigma[k]**2).reshape(n, t).sum(axis=1))
#                 H[col, j] = H[j, col]
    
#     # Make Hessian symmetric
#     H = (H + H.T) / 2
    
#     return H

# %%
# ----
params = convert_out_h0_to_flatten_params_stationary_normal(out_h0)
log_likelihood, w, ll = log_likelihood_stationary_normal(params, data)
score = score_stationary_normal(data, params)
