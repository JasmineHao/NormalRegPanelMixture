import numpy as np
from numba import njit, prange
from numba.typed import Dict, List
from numba.core import types
import math

# Functions for Numba
# ----------------------------------------------------------
@njit
def invert_matrix(mat, epsilon=1e-8):
    """
    Numba-compatible function to compute the inverse of a square matrix.
    If the determinant is close to zero, the matrix is regularized by adding
    epsilon to the diagonal before inversion.

    Parameters:
        mat (ndarray): Input square matrix.
        epsilon (float): Small value added to the diagonal for regularization.

    Returns:
        ndarray: Inverse of the matrix.
    """
    # Ensure the matrix is square
    if mat.shape[0] != mat.shape[1]:
        # Numba cannot raise exceptions, so we return an empty array for invalid input
        return np.zeros_like(mat)
    
    # Compute the determinant
    det_val = np.linalg.det(mat)
    
    # Regularize the matrix if the determinant is close to zero
    if abs(det_val) < epsilon:
        mat = mat + np.eye(mat.shape[0]) * epsilon
    
    # Compute and return the inverse
    return np.linalg.inv(mat)


@njit
def min_along_axis_0(r):
    # Get the shape of the array
    rows, cols = r.shape
    
    # Initialize an array to store the minimum values for each column
    min_vals = np.empty(cols)
    
    # Iterate through each column
    for j in prange(cols):
        # Initialize the minimum value for the current column
        min_val = r[0, j]
        
        # Iterate through each row in the current column
        for i in prange(1, rows):
            if r[i, j] < min_val:
                min_val = r[i, j]
        
        # Store the minimum value for the column
        min_vals[j] = min_val
    
    return min_vals


@njit
def max_along_axis_1(matrix):
    """
    Compute the maximum along axis 1 for a 2D array.
    This replaces np.max(axis=1) for Numba compatibility.
    """
    n_rows, n_cols = matrix.shape
    max_values = np.empty(n_rows)  # Array to store max values for each row
    for i in prange(n_rows):
        max_values[i] = -np.inf  # Initialize with negative infinity
        for j in prange(n_cols):
            if matrix[i, j] > max_values[i]:
                max_values[i] = matrix[i, j]
    return max_values

@njit
def mean_along_axis_1(matrix):
    """
    Compute the mean along axis 1 for a 2D array.
    This replaces np.mean(array, axis=1) for Numba compatibility.
    """
    n_rows, n_cols = matrix.shape
    mean_values = np.empty(n_rows)  # Array to store mean values for each row
    for i in prange(n_rows):
        row_sum = 0.0
        for j in prange(n_cols):
            row_sum += matrix[i, j]
        mean_values[i] = row_sum / n_cols  # Compute mean for the row
    return mean_values


@njit
def compute_quantile(data, q):
    """
    Compute the q-th quantile manually.
    This replaces np.quantile for Numba compatibility.
    """
    sorted_data = np.sort(data)  # Sort the data
    idx = int(q * (len(sorted_data) - 1))  # Find the index for the quantile
    return sorted_data[idx]



@njit
def solve_linear_system_safe(A, b):
    """
    Safely solve the linear system Ax = b.
    If A is singular or nearly singular, return a default solution (e.g., zeros).
    """
    # Check if the matrix is singular
    det = np.linalg.det(A)
    if abs(det) < 1e-12:  # Threshold for singularity
        # Handle singular matrix (e.g., return zeros or raise an error)
        return np.zeros_like(b)  # Return a vector of zeros
    else:
        # Solve the system using np.linalg.solve
        return np.linalg.solve(A, b)
    
@njit
def solve_least_squares(A, b):
    """Solve the least squares problem Ax = b using the normal equation."""
    AtA = A.T @ A  # Compute A^T * A
    Atb = A.T @ b  # Compute A^T * b
    return solve_linear_system_safe(AtA, Atb)


@njit
def generate_random_uniform(low, high, size):
    """Generate random uniform samples using Numba."""
    out = np.empty(size)
    for i in prange(size[0]):
        for j in prange(size[1]):
            out[i, j] = low + (high - low) * np.random.random()
    return out

@njit
def matrix_sqrt(A):
    """Compute the square root of a matrix using eigen-decomposition."""
    # Eigen-decomposition of the matrix
    vals, vecs = np.linalg.eigh(A)
    vals[vals < 0] = 0
    # Compute the square root of eigenvalues
    sqrt_vals = np.sqrt(vals)
    # Reconstruct the matrix square root
    sqrt_A = vecs @ np.diag(sqrt_vals) @ vecs.T
    return sqrt_A

@njit
def repeat_column(R, T, mm):
    n_rows = R.shape[0]
    R_T = np.empty((n_rows * T,), dtype=R.dtype)  # Create an empty array for the result
    for i in prange(n_rows):
        for t in prange(T):
            R_T[i * T + t] = R[i, mm]
    return R_T

@njit
def repeat_elements(arr, k):
    """
    Numba-compatible function to repeat each element in an array k times.
    
    Parameters:
        arr (array): Input 1D array.
        k (int): Number of repetitions for each element.
        
    Returns:
        array: A new array with elements repeated k times.
    """
    n = len(arr)
    result = np.empty(n * k, dtype=arr.dtype)  # Pre-allocate output array
    
    for i in prange(n):
        for j in prange(k):
            result[i * k + j] = arr[i]
    
    return result

@njit
def fill_nan(arr, fill_value):
    """
    Fill NaN values in a NumPy array with a specified value.
    
    Parameters:
        arr (np.ndarray): Input array with potential NaN values.
        fill_value (float): Value to replace NaN with.
    
    Returns:
        np.ndarray: Array with NaN values replaced.
    """
    for i in prange(arr.shape[0]):
        for j in prange(arr.shape[1]):
            # Check if the value is NaN
            if math.isnan(arr[i, j]):
                arr[i, j] = fill_value
    return arr

# %%
# Function to generate data
# ----------------------------------------------------------


@njit(parallel=True)
def generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q):
    R = np.zeros((N, M))

    # Normalize alpha if necessary
    alpha_sum = np.sum(alpha)
    if alpha_sum != 1:
        alpha = alpha / alpha_sum
    
    # Check input consistency - Numba doesn't support exceptions like Python
    if len(alpha) != M or len(mu) != M:
        raise ValueError("M must be the size of alpha and mu")
    
    # Generate prior and initialize R
    prior = np.random.random(size=N)
    alpha_cum = np.zeros(M + 1)
    for m in prange(M):
        alpha_cum[m + 1] = alpha_cum[m] + alpha[m]
    
    if M > 1:
        for m in prange(M):
            lb = alpha_cum[m]
            ub = alpha_cum[m + 1]
            for n in prange(N):
                R[n, m] = 1 if lb < prior[n] <= ub else 0
    else:
        R[:] = 1

    # Initialize output arrays
    Y = np.zeros((T, N))
    
    # Generate x and z if not provided
    if q > 0:
        x = np.empty((N * T, q))
        for i in prange(N * T):
            for j in prange(q):
                x[i, j] = np.random.normal()  # Generate one value at a time
    else:
        x = np.zeros((N * T, 1), dtype=np.float64)
        
    if p > 0:
        z = np.empty((N * T, p))
        for i in prange(N * T):
            for j in prange(p):
                z[i, j] = np.random.normal()  # Generate one value at a time
    else:
        z = np.zeros((N * T, 1), dtype=np.float64)
        
    # Precompute dot products
    mu_R = np.dot(R, mu)  # Use np.dot for matrix multiplication
    sigma_R = np.dot(R, sigma)  # Use np.dot for matrix multiplication
    beta_R = np.dot(R, beta) 
   
    # Generate u array (workaround for np.random.normal with size)
    u = np.empty((T, N))
    for t in prange(T):
        for n in prange(N):
            u[t, n] = np.random.normal()  # Generate one value at a time

    # Generate Y
    for nn in prange(N):
        y_nn = np.zeros(T)
        y_nn = mu_R[nn] + sigma_R[nn] * u[:, nn]
        
        y_nn += x[(T * nn):(T * (nn + 1)), :] @ beta_R[nn, :]
        y_nn += z[(T * nn):(T * (nn + 1)), :] @ gamma
        
        Y[:, nn] = y_nn    
    return(Y, x, z)

# %%
@njit(parallel=True)
def generate_data_ar1(alpha, rho, mu, sigma, gamma, beta, N, T, M, p, q):
    # First compute the stationary distribution 
    mu_0 = mu / (1 - rho)
    beta_0 = np.zeros(beta.shape)
    # beta 
    sigma_0_sq = sigma**2 / (1 - rho**2)
    for mm in prange(M):
        for j in prange(q):
            beta_0[mm, j] = beta[mm, j] / (1- rho[mm])
            sigma_0_sq[mm] += beta[mm, j]**2 / (1- rho[mm]**2)
    sigma_0 = np.sqrt(sigma_0_sq)
    # Generate random assignment
    R = np.zeros((N, M))

    # Normalize alpha if necessary
    alpha_sum = np.sum(alpha)
    if alpha_sum != 1:
        alpha = alpha / alpha_sum
    
    # Check input consistency - Numba doesn't support exceptions like Python
    if len(alpha) != M or len(mu) != M:
        raise ValueError("M must be the size of alpha and mu")
    
    # Generate prior and initialize R
    prior = np.random.random(size=N)
    alpha_cum = np.zeros(M + 1)
    for m in prange(M):
        alpha_cum[m + 1] = alpha_cum[m] + alpha[m]
    
    if M > 1:
        for m in prange(M):
            lb = alpha_cum[m]
            ub = alpha_cum[m + 1]
            for n in prange(N):
                R[n, m] = 1 if lb < prior[n] <= ub else 0
    else:
        R[:] = 1

    
    # Generate x and z if not provided
    if q > 0:
        x = np.empty((N * T, q))
        x_0 = np.empty((N , q))
        for i in prange(N):
            for j in prange(q):
                x_0[i, j] = np.random.normal()  # Generate one value at a time
        for i in prange(N * T):
            for j in prange(q):
                x[i, j] = np.random.normal()  # Generate one value at a time
    else:
        x = np.zeros((N * T, 1), dtype=np.float64)
        x_0 = np.zeros((N, 1), dtype=np.float64)
        
    if p > 0:
        z = np.empty((N * T, p))
        z_0 = np.empty((N, p))
        for i in prange(N):
            for j in prange(p):
                z_0[i, j] = np.random.normal()  # Generate one value at a time
        for i in prange(N * T):
            for j in prange(p):
                z[i, j] = np.random.normal()  # Generate one value at a time
    else:
        z = np.zeros((N * T, 1), dtype=np.float64)
        z_0 = np.zeros((N, 1), dtype=np.float64)
        
    # Precompute dot products
    mu_R = np.dot(R, mu)  # Use np.dot for matrix multiplication
    rho_R = np.dot(R, rho)  # Use np.dot for matrix multiplication
    mu_0_R = np.dot(R, mu_0)  # Use np.dot for matrix multiplication
    sigma_R = np.dot(R, sigma)  # Use np.dot for matrix multiplication
    sigma_0_R = np.dot(R, sigma_0)  # Use np.dot for matrix multiplication
    beta_R = np.dot(R, beta) 
    beta_0_R = np.dot(R, beta_0) 
   
    # Generate u array, including T=0 (workaround for np.random.normal with size)
    u = np.empty((T, N))
    u_0 = np.empty((N, 1))
    for n in prange(N):
        u_0[n,0] = np.random.normal()
        for t in prange(T):
            u[t, n] = np.random.normal()  # Generate one value at a time

    # Initialize output arrays
    Y = np.zeros((T+1, N))
    Y_0 = np.zeros((N, 1))

    # Generate Y0 
    for nn in prange(N):
        y_nn = mu_0_R[nn] + sigma_0_R[nn] * u_0[nn]
        y_nn += x_0[nn] @ beta_0_R[nn,:]
        Y_0[nn,:] = y_nn
    Y[0,:] = Y_0.T
    # Generate Y
    for nn in prange(N):
        for tt in prange(T):
            y_nn = 0
            y_nn = mu_R[nn] + sigma_R[nn] * u[tt, nn]
            y_nn += Y[tt, nn] * rho_R[nn]
            y_nn += x[(T * nn + tt), :] @ beta_R[nn, :]
            y_nn += z[(T * nn + tt), :] @ gamma
            Y[tt+1, nn] = y_nn    
    Y = Y[1:,:]
    return(Y, x, z, x_0, z_0)

# %%
@njit(parallel=True)
def generate_data_mixture(alpha, mu, sigma, tau, N, T, M, K, p=0, q=0):

    R = np.ones((N, M))
    R_T = np.ones((N*T, M)) 
    alpha_sum = np.sum(alpha)
    if alpha_sum != 1:
        alpha = alpha / alpha_sum
    
    prior = np.random.random(size=N)
    alpha_cum = np.zeros(M + 1)
    for mm in prange(M):
        alpha_cum[mm + 1] = alpha_cum[mm] + alpha[mm]
        lb = alpha_cum[mm]
        ub = alpha_cum[mm + 1]
        for n in prange(N):
            R[n, mm] = 1 if lb < prior[n] <= ub else 0
            R_T[:,mm] = repeat_column(R, T, mm)
                

    R_sub = np.zeros((N * T, M, K))
    mu_R_sub = np.zeros((N*T, M))
    # if q > 0:
    #     beta_R_sub = np.zeros((N*T, M, q))
    
    for mm in prange(M):
        prior_m = np.random.random(size=N*T)
        tau_cum = np.cumsum(np.concatenate((np.array([0]), tau[mm])))

        
        for kk  in prange(K):
            lb = tau_cum[kk]
            ub = tau_cum[kk + 1]
            R_sub[:,mm, kk] = ((prior_m > lb) & (prior_m <= ub)) 
        mu_R_sub[:, mm] = np.dot(R_sub[:,mm].astype(np.float64), mu[mm])
        # if q > 0:
        #     beta_R_sub[:,mm,:] = np.dot(R_sub[mm], beta[mm])[:,None]
    
    sigma_R = np.dot(R_T, sigma)    
    mu_R = (R_T * mu_R_sub).sum(axis=1)
    
    # if q > 0:
    #     beta_R = (R_T[:,:,np.newaxis] * beta_R_sub).sum(axis=1)
    u = np.empty((T, N), dtype=np.float64)
    for t in prange(T):
        for n in prange(N):
            u[t, n] = np.random.normal(0.0, 1.0)  # Mean=0, StdDev=1
    Y = np.zeros((T, N))
    
    for nn in prange(N):
        y_nn = mu_R[(T * nn):(T * (nn + 1))] + sigma_R[(T * nn):(T * (nn + 1))] * u[:, nn]
        # if q > 0:
        #     y_nn += np.dot(x[(T * nn):(T * (nn + 1)), :], beta_R[(T * nn):(T * (nn + 1)), :])
        # if p > 1:
        #     y_nn += np.dot(z[(T * nn):(T * (nn + 1)), :], gam)
        Y[:, nn] = y_nn
    x = np.zeros((N * T, 1), dtype=np.float64)
    z = np.zeros((N * T, 1), dtype=np.float64)
    return(Y, x, z)



# %%
# Nonparametric test
# ----------------------------------------------------------
@njit
def create_indicator_list(data_c, T, N, n_bins):
    """
    Create a list of indicator matrices based on quantiles for each time period.
    """
    indicator_list = List()
    
    for t in prange(T):
        # Calculate quantiles
        quantiles = np.empty(n_bins + 1)
        sorted_data = np.sort(data_c[t, :])
        
        # Define quantile boundaries
        quantiles[0] = -np.inf
        quantiles[n_bins] = np.inf
        step = N / n_bins
        for i in prange(1, n_bins):
            index = int(round(i * step)) - 1  # Use rounding for more stable indexing
            quantiles[i] = sorted_data[index]
        
        # print(quantiles)
        
        # Create indicator matrix
        indicator_matrix = np.zeros((N, n_bins))
        for n in prange(N):
            for b in prange(n_bins):
                # Assign data to bin based on quantile boundaries
                if quantiles[b] <= data_c[t, n] < quantiles[b + 1]:
                    indicator_matrix[n, b] = 1
                    break
        indicator_list.append(indicator_matrix)
    
    return indicator_list

# %%
@njit
def calculate_P_matrix(data_c, weights, n_grid=3, n_bins=2):
    """
    Calculate P matrices and Sigma matrices for triplets in a Numba-compatible way.
    """
    T = data_c.shape[0]
    N = data_c.shape[1]
    
    # Create `indicator_list_Y` with 2 bins
    indicator_list_Y = create_indicator_list(data_c, T, N, n_bins=n_bins)
    
    # Create `indicator_list_Y_ngrid` with `n_grid` bins
    indicator_list_Y_ngrid = create_indicator_list(data_c, T, N, n_bins=n_grid)
    
    # Initialize the result lists
    P_k_list = List()
    Sigma_P_k_list = List()
    
    # Iterate over the t periods
    for k in prange(T):
        # Compute the Kronecker product for each row manually
        data_partition_matrix_y = np.zeros((N, (n_bins ** (T - 1))))
        for n in prange(N):
            # Manually compute Kronecker product for the row
            kron_result = np.array([1.0])  # Start with scalar 1.0
            for t in prange(T):
                if t != k:
                    kron_result = np.kron(kron_result, indicator_list_Y[t][n, :])
            data_partition_matrix_y[n, :] = kron_result
        
        # Compute P_k
        P_k = (weights * indicator_list_Y_ngrid[k].T) @ data_partition_matrix_y
        P_k_list.append(P_k)
        
        # Compute Sigma_P_k
        P_k_vec = P_k.T.flatten()
        W_P_s = np.diag(P_k_vec) - np.outer(P_k_vec, P_k_vec)
        Sigma_P_k_list.append(W_P_s)

    
    return {
        "P_k_list": P_k_list,
        "Sigma_P_k_list": Sigma_P_k_list
    }




@njit
def compute_A_q_o(U_22, U_12):
    """Compute A_q_o."""
    sqrt_U22 = matrix_sqrt(U_22 @ U_22.T)
    inv_U22_T = invert_matrix(U_22.T)
    A_q_o = np.transpose(sqrt_U22 @ inv_U22_T @ np.hstack((U_12.T, U_22.T)))
    return A_q_o

@njit
def compute_B_q_o(V_22, V_12):
    """Compute B_q_o."""
    sqrt_V22 = matrix_sqrt(V_22 @ V_22.T)
    inv_V22_T = invert_matrix(V_22.T)
    B_q_o = sqrt_V22 @ inv_V22_T @ np.hstack((V_12.T, V_22.T))
    return B_q_o

@njit
def compute_kron_BA_o(B_q_o, A_q_o):
    """Compute the Kronecker product of B_q_o and A_q_o.T."""
    return np.kron(B_q_o, A_q_o.T)

@njit
def matrix_svd_decomposition(P, m):
    """
    Perform SVD decomposition and compute A_q_o, B_q_o, and Kronecker product.
    """
    # Perform SVD outside the Numba function
    U, S, VT = np.linalg.svd(P, full_matrices=True)
    V = VT.T
    
    # Submatrices of U and V
    U_12 = U[:m, m:]
    V_12 = V[:m, m:]
    U_22 = U[m:, m:]
    V_22 = V[m:, m:]
    
    # Compute A_q_o and B_q_o using Numba-compiled functions
    A_q_o = compute_A_q_o(U_22, U_12)
    B_q_o = compute_B_q_o(V_22, V_12)
    
    # Compute the Kronecker product
    kron_BA_o = compute_kron_BA_o(B_q_o, A_q_o)
    
    # Ensure all arrays are 2D
    S = S.reshape(-1, 1)  # Convert S to a 2D column vector
    U = np.atleast_2d(U)
    V = np.atleast_2d(V)
    U_12 = np.atleast_2d(U_12)
    V_12 = np.atleast_2d(V_12)
    U_22 = np.atleast_2d(U_22)
    V_22 = np.atleast_2d(V_22)
    A_q_o = np.atleast_2d(A_q_o)
    B_q_o = np.atleast_2d(B_q_o)
    kron_BA_o = np.atleast_2d(kron_BA_o)

    # Create a Numba-compatible dictionary
    numba_dict = Dict.empty(
        key_type=types.unicode_type,  # Keys are strings
        value_type=types.float64[:, :],  # Values are 2D arrays
    )

    # Add key-value pairs
    numba_dict["D"] = S
    numba_dict["U"] = U
    numba_dict["V"] = V
    numba_dict["U_12"] = U_12
    numba_dict["V_12"] = V_12
    numba_dict["U_22"] = U_22
    numba_dict["V_22"] = V_22
    numba_dict["A_q_o"] = A_q_o
    numba_dict["B_q_o"] = B_q_o
    numba_dict["kron_BA_o"] = kron_BA_o

    return numba_dict



    
@njit
def compute_rk_stat_given_P(P, Sigma_P, P_svd, m, N, lambda_c):
    """
    Compute statistical metrics for the "P" or "Q" transform case.
    """
    # Extract SVD components
    A_q_o = P_svd["A_q_o"]
    B_q_o = P_svd["B_q_o"]
    kron_BA_o = P_svd["kron_BA_o"]

    # Compute lambda_q
    lambda_q = A_q_o.T @ P @ B_q_o.T - lambda_c

    # Compute Omega_q
    Omega_q = kron_BA_o @ Sigma_P @ kron_BA_o.T

    # Compute rk_c
    lambda_q_flat = lambda_q.flatten()
    Omega_q_inv = invert_matrix(Omega_q)
    rk_c = N * (lambda_q_flat @ Omega_q_inv @ lambda_q_flat) 

    return lambda_q, Omega_q, rk_c

@njit
def compute_rk_information_criteria(rk_c, r, N):
    """
    Compute AIC, BIC, and HQ criteria.
    """
    AIC_c = rk_c - 2 * r
    BIC_c = rk_c - np.log(N) * r
    HQ_c = rk_c - 2 * np.log(np.log(N)) * r
    return AIC_c, BIC_c, HQ_c

@njit
def construct_stat_KP(P, Sigma_P, m, N, lambda_c=0):
    """
    Construct statistical metrics for the Kronecker Product and return 
    a Numba-compatible typed dictionary.
    """
    # Perform SVD decomposition
    P_svd = matrix_svd_decomposition(P, m)

    # Compute stats using Numba
    lambda_q, Omega_q, rk_c = compute_rk_stat_given_P(P, Sigma_P, P_svd, m, N, lambda_c)

    # Compute the rank (r)
    r = Omega_q.shape[0]

    # Compute AIC, BIC, and HQ using Numba
    AIC_c, BIC_c, HQ_c = compute_rk_information_criteria(rk_c, r, N)

    # Create a Numba-compatible dictionary to store results
    result_dict = Dict.empty(
        key_type=types.unicode_type,  # Keys are strings
        value_type=types.float64[:, :],  # Values are 2D arrays
    )

    # Add results to the dictionary
    result_dict["rk_c"] = np.array([[rk_c]])  # Scalars must be converted to 2D arrays
    result_dict["lambda_c"] = lambda_q
    result_dict["Omega_q"] = Omega_q
    result_dict["AIC_c"] = np.array([[AIC_c]])
    result_dict["BIC_c"] = np.array([[BIC_c]])
    result_dict["HQ_c"] = np.array([[HQ_c]])

    return result_dict


@njit
def NonParTestParallel(data_nopar, N, T, M, p, q, nrep, n_grid, BB, r_test):
    # Result array
    weights_equal = np.full(N, 1 / N)
    result_rk_each = np.zeros((nrep,2))
    for ii in prange(nrep):
        # Generate synthetic data (replace with actual logic)
        data_c = data_nopar[ii]  # Example: Replace this with your `generate_data` logic
        # Initialize weights
        # Compute P and Sigma matrices
        data_P_W = calculate_P_matrix(data_c, weights_equal, n_grid=n_grid, n_bins=2)
        
        # Initialize results
        rk = np.zeros(T)
        lambda_c_list = List()
        omega_c = List()
        Sigma_P_list = List()
        P_k_list = List()
        
        # Loop through T periods to compute statistics
        for k in prange(T):
            # Extract P_k and Sigma_P_k from the data_P_W object
            P_k = data_P_W["P_k_list"][k]
            Sigma_P_k = data_P_W["Sigma_P_k_list"][k]
            
            # Compute KP statistics for the k-th triplet
            stat_KP = construct_stat_KP(P_k, Sigma_P_k, r_test, N)
            
            # Store results
            rk[k] = stat_KP["rk_c"][0,0]
            lambda_c_list.append(stat_KP["lambda_c"])
            omega_c.append(stat_KP["Omega_q"])
            Sigma_P_list.append(Sigma_P_k)
            P_k_list.append(P_k)
        # Initialize result matrix
        rk_b = np.zeros((BB, T))
        
        # Smoothed Nonparametric Bootstrap
        ru = np.random.exponential(scale=1, size=(BB, N))  # Exponential random variables
        row_sums = ru.sum(axis=1).reshape(-1, 1)  # Reshape to keep dimensions
        ru /= row_sums
        
        for i in prange(BB):
            # Calculate bootstrapped P and Sigma_P matrices
            data_P_W_b = calculate_P_matrix(data_c, ru[i, :], n_grid=n_grid, n_bins=2)
            
            for k in prange(T):
                P_k = data_P_W_b['P_k_list'][k]
                Sigma_P_k = data_P_W_b['Sigma_P_k_list'][k]
                # Compute KP statistics for the k-th triplet
                rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, lambda_c_list[k])['rk_c'][0,0]
                # rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, 0)['rk_c'][0,0]
        # Compute max and mean values for rk and rk_b
        rk_b_max = max_along_axis_1(rk_b)  # Maximum of rk_b along axis 1
        rk_b_max_95 = compute_quantile(rk_b_max, 0.95)  # 95th quantile of rk_b_max

        
        # Store results
        result_rk_each[ii, 0] = 1 * (rk.max() > rk_b_max_95)
        rk_mean = np.mean(rk)  # Mean of rk (Numba supports this)
        rk_b_mean = mean_along_axis_1(rk_b)  # Mean of rk_b along axis 1
        rk_b_mean_95 = compute_quantile(rk_b_mean, 0.95)  # 95th quantile of rk_b_mean
        result_rk_each[ii, 1] = 1 * (rk_mean > rk_b_mean_95)
    return result_rk_each

@njit
def calculate_P_matrix_with_covariates(y, x, weights, n_grid=3, n_bins=2):
    
    T = y.shape[0]
    N = y.shape[1]
    indicator_list_Y = create_indicator_list(y, T, N, n_bins=n_bins)
    indicator_list_Y_ngrid = create_indicator_list(y, T, N, n_bins=n_grid)
    indicator_list_X = create_indicator_list(x.reshape((N,T)).T, T, N, n_bins=n_bins)
    
    # Initialize the result lists
    P_k_list = List()
    Sigma_P_k_list = List()
    
    n_partition = (n_bins ** (T - 1))
    # Iterate over the t periods
    for k in prange(T):
        # Compute the Kronecker product for each row manually
        data_partition_matrix_y = np.zeros((N, (n_bins ** (T - 1))))
        data_partition_matrix_x = np.zeros((N, (n_bins ** (T - 1))))    
    
        for n in prange(N):
            # Manually compute Kronecker product for the row
            kron_result_y = np.array([1.0])  # Start with scalar 1.0
            kron_result_x = np.array([1.0])  # Start with scalar 1.0
            for t in prange(T):
                if t != k:
                    kron_result_y = np.kron(kron_result_y, indicator_list_Y[t][n, :])
                    kron_result_x = np.kron(kron_result_x, indicator_list_X[t][n, :])
            data_partition_matrix_y[n, :] = kron_result_y
            data_partition_matrix_x[n, :] = kron_result_x
        
        for kk in prange(n_partition):
            data_partition_kk = data_partition_matrix_x[:,kk]
            weights_adj = data_partition_kk * weights 
            weights_adj = weights
            weights_adj_sum = max(np.sum(weights_adj), 1e-5)
            weights_adj = weights_adj / weights_adj_sum
            # Compute P_k
            P_k = (weights_adj * indicator_list_Y_ngrid[k].T) @ data_partition_matrix_y
            P_k_list.append(P_k)
            # Compute Sigma_P_k
            P_k_vec = P_k.T.flatten()
            W_P_s = np.diag(P_k_vec) - np.outer(P_k_vec, P_k_vec)
            Sigma_P_k_list.append(W_P_s)
    # print(P_k_list)
    
    return {
        "P_k_list": P_k_list,
        "Sigma_P_k_list": Sigma_P_k_list
    }
    
# %%
@njit
def NonParTestParallelCovariate(Data, N, T, M, p, q, nrep, n_grid, BB, r_test, n_bins=2):
    # Result array
    weights_equal = np.full(N, 1 / N)
    result_rk_each = np.zeros((nrep,2))
    n_partition = (n_bins ** (T - 1))
    n_stats = n_partition * T
    rk_array = np.empty((nrep,n_stats))
    for ii in prange(nrep):     
        
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        # Compute P and Sigma matrices
        data_P_W = calculate_P_matrix_with_covariates(y, x, weights_equal, n_grid=n_grid, n_bins=n_bins)
        # Initialize results
        rk = np.zeros(n_stats)
        lambda_c_list = List()
        omega_c = List()
        Sigma_P_list = List()
        P_k_list = List()
        
        for k in prange(n_stats):
            # Extract P_k and Sigma_P_k from the data_P_W object
            P_k = data_P_W["P_k_list"][k]
            Sigma_P_k = data_P_W["Sigma_P_k_list"][k]
            
            # Compute KP statistics for the k-th triplet
            stat_KP = construct_stat_KP(P_k, Sigma_P_k, r_test, N)
        
            # Store results
            rk[k] = stat_KP["rk_c"][0,0]
            lambda_c_list.append(stat_KP["lambda_c"])
            omega_c.append(stat_KP["Omega_q"])
            Sigma_P_list.append(Sigma_P_k)
            P_k_list.append(P_k)
        
        rk_array[ii,:] = rk 
        # Initialize result matrix
        rk_b = np.zeros((BB, n_stats))
        
        # Smoothed Nonparametric Bootstrap
        ru = np.random.exponential(scale=0.1, size=(BB, N))  # Exponential random variables
        row_sums = ru.sum(axis=1).reshape(-1, 1)  # Reshape to keep dimensions
        ru /= row_sums
        
        for bb in prange(BB):
            data_P_W_b = calculate_P_matrix_with_covariates(y, x, ru[bb, :], n_grid=n_grid, n_bins=n_bins)
            
            for k in prange(n_stats):
                P_k = data_P_W_b['P_k_list'][k]
                Sigma_P_k = data_P_W_b['Sigma_P_k_list'][k]
                # Compute KP statistics for the k-th triplet
                rk_b[bb, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, lambda_c_list[k])['rk_c'][0,0]
                # rk_b[bb, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, 0)['rk_c'][0,0]
        # Compute max and mean values for rk and rk_b
        rk_b_max = max_along_axis_1(rk_b)  # Maximum of rk_b along axis 1
        rk_b_max_95 = compute_quantile(rk_b_max, 0.95)  # 95th quantile of rk_b_max
        # Store results
        result_rk_each[ii, 0] = 1 * (rk.max() > rk_b_max_95)
        rk_mean = np.mean(rk)  # Mean of rk (Numba supports this)
        rk_b_mean = mean_along_axis_1(rk_b)  # Mean of rk_b along axis 1
        rk_b_mean_95 = compute_quantile(rk_b_mean, 0.95)  # 95th quantile of rk_b_mean
        result_rk_each[ii, 1] = 1 * (rk_mean > rk_b_mean_95)
        
    return result_rk_each


# %%
# LR test functions
# ----------------------------------------------------------


SINGULAR_EPS = 1e-10  # Criteria for matrix singularity
M_LN_SQRT_2PI = 0.9189385332046727  # log(sqrt(2*pi))

@njit
def process_array_or_none(arr, nt):
    if arr is None:  # Check for None
        return np.zeros((nt,0),dtype=np.float64)  # Default behavior for None
    return arr  # Process the array if it's valid


@njit
def log_likelihood_normal(y, mu, sigma):
    """
    Calculate the log-likelihood of the data under a normal distribution.

    Parameters:
    - y: array-like, the observed data points
    - mu: float, the mean of the normal distribution
    - sigma: float, the standard deviation of the normal distribution

    Returns:
    - log_likelihood: float, the log-likelihood value
    """
    n = len(y)
    term1 = -n / 2 * np.log(2 * np.pi)  # Constant term
    term2 = -n * np.log(sigma)  # Log of the standard deviation
    term3 = -1 / (2 * sigma**2) * np.sum((y - mu) ** 2)  # Data fitting term
    return term1 + term2 + term3

@njit
def log_likelihood_array(y, mu, sigma):
    """
    Calculate the log-likelihood of each element in the data under a normal distribution.

    Parameters:
    - y: array-like, the observed data points
    - mu: float, the mean of the normal distribution
    - sigma: float, the standard deviation of the normal distribution

    Returns:
    - log_likelihoods: array, log-likelihood of each data point
    """
    # Precompute constants
    constant = -0.5 * np.log(2 * np.pi)
    variance = sigma ** 2

    # Initialize result array
    log_likelihoods = np.empty(len(y))

    # Compute log-likelihood for each data point
    for i in prange(len(y)):
        log_likelihoods[i] = (
            constant 
            - np.log(sigma) 
            - ((y[i] - mu) ** 2) / (2 * variance)
        )
    
    return log_likelihoods


@njit
def compute_residual_normal_reg(m, n, t, sigma_jn, res):
    """
    Compute residuals for the EM optimization loop in a Numba-compatible way.

    Parameters:
    - m: Number of components (int)
    - n: Number of groups (int)
    - t: Number of time points per group (int)
    - sigma_jn: Array of current sigma values (1D array of floats, shape (m,))
    - res: Adjusted response variable (2D array of floats, shape (m, n * t))
    
    Returns:
    - r: Residuals array (2D array of floats, shape (m, n))
    """
    # Initialize the residuals array
    r = np.zeros((m, n), dtype=np.float64)

    # Loop over each component (m)
    for j in prange(m):
        
        # Loop over each group (n)
        for i in prange(n):
            sum_r_t = 0.0

            # Loop over each time point within the group (t)
            for k in prange(t):
                idx = i * t + k  # Compute the flattened index
                diff = res[j, idx]
                r_t = (1.0 / sigma_jn[j]) * diff
                sum_r_t += 0.5 * (r_t**2)

            # Compute residual for group i and component j
            r[j, i] = t * np.log(sigma_jn[j]) + sum_r_t
    return r



@njit
def EM_optimization(y, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=2000, tol=1e-8, tau = 0.5, epsilon=0.05):
    
    nt = len(y)
    n = nt // t

    ninits = alpha_draw.shape[1]
    # Handle x
    if q == 0:
        x1 = np.ones((nt, 1))
        q1 = 1
    else:
        x1 = np.zeros((nt, x.shape[1] + 1))
        x1[:, 0] = 1  # Add intercept
        x1[:, 1:] = x
        q1 = x1.shape[1]
    
    # Initialize variables
    lb = np.zeros(m)
    ub = np.zeros(m)
    l_j = np.zeros(m)
    w = np.zeros((m, nt))
    res = np.zeros((m, nt))
    post = np.zeros((m * n, ninits))
    notcg = np.zeros(ninits)
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    
    for jn in prange(ninits):
        alpha_jn = alpha_draw[:, jn]
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64

        mubeta_jn_mat = mubeta_jn.reshape((q+1,m)).T
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0
        
        for iter_ii in prange(maxit):
            ll = -nt * M_LN_SQRT_2PI
            
            if p > 0:
                ytilde = y - np.dot(z, gamma_jn)
            else:
                ytilde = y
            
            for j in prange(m):
                res[j] = ytilde - x1 @ mubeta_jn_mat[j]
            r = compute_residual_normal_reg(m, n, t, sigma_jn, res)
            
            minr = min_along_axis_0(r)
            
            # Initialize arrays
            l_j = np.zeros((m,n))  # Same shape as `r`
            sum_l_j = np.zeros(n)   # Sum along axis 0
            w = np.zeros((m,n))    # Weights
            ll = 0.0                # Log-likelihood accumulator

            # Compute l_j = alpha_jn[:, None] * exp(minr - r)
            for i in prange(n):
                for j in prange(m):
                    l_j[j, i] = alpha_jn[j] * np.exp(minr[i] - r[j,i])
            
            # Compute sum_l_j = np.sum(l_j, axis=0)
            for j in prange(m):
                for i in prange(n):
                    sum_l_j[i] += l_j[j, i]
            
            # Compute w = l_j / sum_l_j
            for i in prange(n):
                for j in prange(m):
                    w[j, i] = l_j[j, i] / sum_l_j[i]
            
            # Compute ll += np.sum(np.log(sum_l_j) - minr)
            for i in prange(n):
                ll += np.log(sum_l_j[i]) - minr[i]
            
            penloglik = ll + np.log(2.0) + min(np.log(tau), np.log(1 - tau))
            
            for j in prange(m):
                s0j = sigma_0[j] / sigma_jn[j]
                penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
                # penloglik += min(np.log(alpha_jn[j]), np.log(1 - alpha_jn[j]))
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            emit += 1
            
            if np.max(np.abs(diff)) < tol:
                break
            
            # Update parameters
            mubeta_jn_mat = np.zeros((m,q1),dtype=np.float64)
            wtilde = np.zeros(nt)
            for j in prange(m):
                alpha_jn[j] = np.mean(w[j, :])
                wtilde = w[j, :].T
                w_j = np.zeros(nt)
                for i in prange(n):
                    w_j[i * t : (i + 1) * t] = wtilde[i]
                xtilde = np.zeros((nt, q1))
                for ii in prange(q1):
                    xtilde[:, ii] = w_j * x1[:, ii]
                # design_matrix = xtilde.T @ x1
                # solve_linear_system_safe(xtilde.T @ x1, xtilde.T @ ytilde)
                # xtilde.T @ ytilde
                mubeta_jn_mat[j,:] = solve_linear_system_safe(xtilde.T @ x1, xtilde.T @ ytilde)
                ssr_j = np.sum(w_j * (ytilde - x1 @ mubeta_jn_mat[j,:])**2)
                sigma_jn[j] = np.sqrt((ssr_j + 2.0 * an * sigma_0[j]**2) / (np.sum(w_j) + 2.0 * an))
                sigma_jn[j] = max(sigma_jn[j], epsilon * sigma_0[j])
            
            # update alpha
            total_alpha = np.sum(alpha_jn)
            for j in prange(m):
                alpha_jn[j] = max(0.01, alpha_jn[j] / total_alpha)
            
            # update gamma
            if p > 0:
                ztilde = np.zeros((nt, p), dtype=np.float64) 
                zz = np.zeros((p, p), dtype=np.float64) 
                ze = np.zeros((p, 1), dtype=np.float64) 
                for j in prange(m):
                    wtilde = w[j, :]
                    w_j = np.zeros(nt)
                    for i in prange(n):
                        w_j[i * t : (i + 1) * t] = wtilde[i]
                    for ii in prange(p):
                        ztilde[:, ii] = w_j * z[:, ii]
                    zz += ztilde.T @ z / (sigma_jn[j]**2)
                    ze += ztilde.T @( y - x1 @ mubeta_jn_mat[j,:]) / (sigma_jn[j]**2)
                gamma_jn = solve_linear_system_safe(zz,ze).flatten()
            
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mubeta_draw[:, jn] = mubeta_jn_mat.T.flatten()
        sigma_draw[:, jn] = sigma_jn
        if p > 0:
            gamma_draw[:, jn] = gamma_jn
    return(alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset ,post)
                
# %%   
@njit
def regpanelmixPMLE(y,x,z, p, q, m, ninits=10, epsilon=1e-6, maxit=2000, epsilon_short=1e-2, maxit_short=200)  : 
    
    t,n = y.shape
    nt = n * t
    y = y.T.flatten()
    
    # y.reshape((n,t)).T - data_lr[0][0] # check equivalence
    # Handle x
    
    if q > 0:
        x1 = np.hstack((np.ones((nt, 1)), x))
    else:
        x1 = np.ones((nt, 1))
    q1 = q + 1
    
    if p > 0:
        xz = np.hstack((x1, z))
    else:
        xz = x1 
        
    
    out_coef = solve_least_squares(xz, y)  # Replace np.linalg.lstsq
    residuals = y - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q1 + 1) * m + p
    # ninits_short = ninits * 10 * (q1 + p) * m
    ninits_short = ninits * 10 
    
    if (m == 1) :
        mubeta = out_coef[:q1]
        if p > 0:
            gamma = out_coef[q1:(q1 + p)]
        else:
            gamma = np.array([0.0])
        res = y - xz @ out_coef
        sigma = np.sqrt(np.mean(res**2))
        loglik = log_likelihood_normal(res,0,sigma)

        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        penloglik = loglik
        alpha = np.array([1])
        postprobs = np.ones(n)
    else: 
        # First draw random start point
        if p > 0:
            gamma = out_coef[q1:(q1 + p)]
            # Perform least squares regression with both x and z
            gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) * gamma
            mubeta_hat = out_coef[:q1]
            y = y - z @ gamma
        else:
            # Perform least squares regression with x only
            
            gamma = np.array([0.0])
            mubeta_hat = out_coef[:-1]
            gamma_draw = np.zeros((1,ninits_short), dtype=np.float64)

        # Initialize alpha
        alpha_draw = generate_random_uniform(0, 1, (m, ninits_short))
        alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))

        # Initialize mubeta
        if q > 0:
            minMU = np.min(y - x @ mubeta_hat[1:])
            maxMU = np.max(y - x @ mubeta_hat[1:])
            mubeta_draw = np.zeros((q1 * m, ninits_short))
            for j in prange(m):
                mubeta_draw[q1 * j, :] = np.random.uniform(minMU, maxMU, size=ninits_short)
                for i in prange(1, q1):
                    mubeta_draw[q1 * j + i, :] = mubeta_hat[i] * np.random.uniform(-2, 2, size=ninits_short)
        else:
            minMU = np.min(y)
            maxMU = np.max(y)
            mubeta_draw = np.zeros((q1 * m, ninits_short))
            for mm in prange(m):
                mubeta_draw[q1 * mm, :] = np.random.uniform(minMU, maxMU, size=ninits_short)
        
        an = 1 / n    
        sigma_0 = np.full(m, stdR)
    
        # Initialize sigma
        sigma_draw = generate_random_uniform(0.01, 1, (m, ninits_short)) * stdR
        
        alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=maxit_short, tol=epsilon_short)

        components = np.argsort(penloglikset)[::-1][:ninits]
        alpha_draw = alpha_draw[:,components]
        mubeta_draw = mubeta_draw[:,components]
        sigma_draw = sigma_draw[:,components]
        gamma_draw = gamma_draw[:,components]
        
        
        alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=maxit, tol=epsilon)
        
        index = np.argmax(penloglikset)
        alpha_hat = alpha_draw[:,index]
        mubeta_hat = mubeta_draw[:,index]
        sigma_hat = sigma_draw[:,index]
        gamma_hat = gamma_draw[:,index]
        post = post[:, index]
        penloglik = penloglikset[index]
        loglik = loglikset[index]
        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        
        # Create a Numba-compatible dictionary to store results
        result_dict = Dict.empty(
            key_type=types.unicode_type,  # Keys are strings
            value_type=types.float64[:, :],  # Values are 2D arrays
        )
        
        result_dict['penloglik'] = np.array([[penloglik]])
        result_dict['loglik'] = np.array([[loglik]])
        result_dict['aic'] = np.array([[aic]])
        result_dict['bic'] = np.array([[bic]])
        
        result_dict['alpha_hat']  = alpha_hat[np.newaxis,:]
        result_dict['sigma_hat']  = sigma_hat[np.newaxis,:]
        result_dict['mubeta_hat'] = mubeta_hat[np.newaxis,:]
        result_dict['gamma_hat'] = gamma_hat[np.newaxis,:]
        
        return result_dict

# %%@njit
def EM_optimization_AR1(y_c, xz, y_0, xz_0, p, q, sigma_0, alpha_draw, rho_draw, mubeta_draw, sigma_draw, gamma_draw, m, n, t, an, maxit=2000, tol=1e-8, tau = 0.5, epsilon=0.05):
    nt = n * t
    ninits = alpha_draw.shape[1]
    # Handle x
    q1 = 2 * q + 1
    # Initialize variables    
    l_j = np.zeros(m)
    w = np.zeros((m, n * (t-1)))
    res = np.zeros((m, n * (t-1)))
    res_0 = np.zeros((m, n ))
    post = np.zeros((m * n, ninits))
    notcg = np.zeros(ninits)
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    for jn in prange(ninits):
        alpha_jn = alpha_draw[:, jn]
        rho_jn = rho_draw[:, jn]
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64

        mubeta_jn_mat = mubeta_jn.reshape((q+1,m)).T
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0
        
        for iter_ii in prange(maxit):
            ll = -nt * M_LN_SQRT_2PI
            
            mubeta_0_jn_mat = np.zeros((m,q+1),dtype=np.float64)
            mubeta_jn_mat = np.zeros((m,q1+1),dtype=np.float64)
            for mm in prange(m):
                mubeta_jn_mat[mm,-1] = rho_jn[mm]
                for j in prange(q+1):
                    mubeta_0_jn_mat[mm, j] = mubeta_jn[mm*q1+j,]/ (1 - rho_jn[mm])
                    if j == 0:
                        mubeta_jn_mat[mm,j] = mubeta_jn[mm*q1+j,]
                    else:
                        mubeta_jn_mat[mm,j] = mubeta_jn[mm*q1+j,] # coefficient on x
                        mubeta_jn_mat[mm,q1+j] = - mubeta_jn[mm*q1+j,] * rho_jn[mm] # coefficient on x lag
                        
            # beta 
            sigma_0_sq = sigma_jn**2 / (1 - rho_jn**2)
            for mm in prange(m):
                for j in prange(q):
                    sigma_0_sq[mm] += mubeta_jn_mat[mm, j+1]**2 / (1- rho_jn[mm]**2)
            sigma_0_jn = np.sqrt(sigma_0_sq)
            
            r = np.zeros((m, n*(t-1)))
            r_0 = np.zeros((m, n))
            for mm in prange(m):
                res[mm] = y_c - xz @ mubeta_jn_mat[mm]
                res_0[mm] = y_0 - xz_0 @ mubeta_0_jn_mat[mm]
                r[mm] = log_likelihood_array(res[mm], 0.0, sigma_jn[mm])
                r_0[mm] = log_likelihood_array(res_0[mm], 0.0, sigma_0_jn[mm])
            
            # --------------
            r_sum = np.zeros((m,n))
            for mm in prange(m):
                for nn in prange(n):
                    r_sum[mm,nn] += res_0[mm, nn]
                    for tt in prange(t-1):
                        r_sum[mm,nn] += res[mm, nn * (t-1) + tt ]
            minr = min_along_axis_0(r_sum)
            
            # Initialize arrays
            l_j = np.zeros((m,n))  # Same shape as `r`
            sum_l_j = np.zeros(n)   # Sum along axis 0
            w = np.zeros((m,n))    # Weights
            ll = 0.0                # Log-likelihood accumulator

            # Compute l_j = alpha_jn[:, None] * exp(minr - r)
            for nn in prange(n):
                for mm in prange(m):
                    l_j[mm, nn] = alpha_jn[mm] * np.exp(r_sum[mm,nn] - minr[nn])
                    sum_l_j[nn] += l_j[mm, nn]
                    
            # Compute w = l_j / sum_l_j
            for nn in prange(n):
                for mm in prange(m):
                    w[mm, nn] = l_j[mm, nn] / sum_l_j[nn]
            
            # Compute ll += np.sum(np.log(sum_l_j) - minr)
            for nn in prange(n):
                ll += np.log(sum_l_j[nn]) + minr[nn]
            
            penloglik = ll + np.log(2.0)
            
            for j in prange(m):
                s0j = sigma_0[j] / sigma_jn[j]
                penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
                
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            emit += 1
            
            if np.max(np.abs(diff)) < tol:
                break
            
            # Update parameters
            wtilde = np.zeros(nt)
            for j in prange(m):
                alpha_jn[j] = np.mean(w[j, :])
                wtilde = w[j, :].T
                w_j = np.zeros(n*(t-1))
                for i in prange(n):
                    w_j[i * (t-1) : (i + 1) * (t-1)] = wtilde[i]
                xtilde = np.zeros(xz.shape)
                for ii in prange(xz.shape[1]):
                    xtilde[:, ii] = w_j * xz[:, ii]
                
                coef_jn = solve_linear_system_safe(xtilde.T @ xz, xtilde.T @ y_c)
                mubeta_jn_mat[j,:] = coef_jn[:q1]
                rho_jn[j] = coef_jn[-1]
                if rho_jn[j] > 0.99:
                    rho_jn[j] = 0.99
                ssr_j = np.sum(w_j * (y_c - xz @ coef_jn)**2)
                sigma_jn[j] = np.sqrt((ssr_j + 2.0 * an * sigma_0[j]**2) / (np.sum(w_j) + 2.0 * an))
                sigma_jn[j] = max(sigma_jn[j], epsilon * sigma_0[j])
            
            # update alpha
            for j in prange(m):
                alpha_jn[j] = max(0.01, alpha_jn[j] )
            
            total_alpha = np.sum(alpha_jn)
            for j in prange(m):
                alpha_jn[j] = alpha_jn[j] / total_alpha
            
            # update gamma
            # if p > 0:
            #     ztilde = np.zeros((nt, p), dtype=np.float64) 
            #     zz = np.zeros((p, p), dtype=np.float64) 
            #     ze = np.zeros((p, 1), dtype=np.float64) 
            #     for j in prange(m):
            #         wtilde = w[j, :]
            #         w_j = np.zeros(nt)
            #         for i in prange(n):
            #             w_j[i * t : (i + 1) * t] = wtilde[i]
            #         for ii in prange(p):
            #             ztilde[:, ii] = w_j * z[:, ii]
            #         zz += ztilde.T @ z / (sigma_jn[j]**2)
            #         ze += ztilde.T @( y - x1 @ mubeta_jn_mat[j,:]) / (sigma_jn[j]**2)
            #     gamma_jn = solve_linear_system_safe(zz,ze).flatten()
            
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mubeta_draw[:, jn] = mubeta_jn_mat.T.flatten()
        sigma_draw[:, jn] = sigma_jn
        if p > 0:
            gamma_draw[:, jn] = gamma_jn
    return(alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset ,post)

# %%

@njit
def regpanelmixPMLEAR1(y, x, z, p, q, m, ninits=10, epsilon=1e-6, maxit=2000, epsilon_short=1e-2, maxit_short=200)  : 
    t,n = y.shape
    nt = n * t
    
    y_l = y[:-1,:]
    y_0 = y[0,:]
    y_c = y[1:,:]
    
    y_c = y_c.T.flatten()
    y_l = y_l.T.flatten()
    
    # y.reshape((n,t)).T - data_lr[0][0] # check equivalence
    # Handle x
    x_0 = np.zeros((n,q))
    x_l = np.zeros((n*(t-1),q))
    x_c = np.zeros((n*(t-1),q))
    for j in prange(q):
        x_j_mat = x[:,j].reshape((n,t)).T
        x_0[:,j] = x_j_mat[0,:]
        x_l[:,j] = x_j_mat[:-1,:].T.flatten()
        x_c[:,j] = x_j_mat[1:,:].T.flatten()
        
    # Handle x
    z_0 = np.zeros((n,p))
    z_l = np.zeros((n*(t-1),p))
    z_c = np.zeros((n*(t-1),p))
    for j in prange(p):
        z_j_mat = z[:,j].reshape((n,t)).T
        z_0[:,j] = z_j_mat[0,:]
        z_l[:,j] = z_j_mat[:-1,:].T.flatten()
        z_c[:,j] = z_j_mat[1:,:].T.flatten()
        
    x1 = np.hstack((np.ones((n*(t-1), 1)), x_c, x_l, y_l[:,np.newaxis]))
    xz = np.hstack((x1, z_l, z_c))
    q1 = 2*q + 1
    
    xz_0 = np.hstack((np.ones((n, 1)), x_0, z_0) )
        
    out_coef = solve_least_squares(xz, y_c)  # Replace np.linalg.lstsq
    residuals = y_c - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q1 + 1) * m + p + 1
    # ninits_short = ninits * 10 * (q1 + p) * m
    ninits_short = ninits * 10 
    
    if (m == 1) :
        mubeta = out_coef[:q1]
        if p > 0:
            gamma = out_coef[q1:(q1 + p)]
        else:
            gamma = np.array([0.0])
        res = y - xz @ out_coef
        sigma = np.sqrt(np.mean(res**2))
        loglik = log_likelihood_normal(res,0,sigma)

        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        penloglik = loglik
        alpha = np.array([1])
        postprobs = np.ones(n)
    else: 
        # First draw random start point
        if p > 0:
            gamma = out_coef[(q1+1):(q1 + p + 1)]
            # Perform least squares regression with both x and z
            gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) * gamma
            mubeta_hat = out_coef[:q1]
            y = y - z @ gamma
        else:
            # Perform least squares regression with x only
            gamma = np.array([0.0])
            gamma_draw = np.zeros((1,ninits_short), dtype=np.float64)
            
        # Initialize alpha
        alpha_draw = generate_random_uniform(0, 1, (m, ninits_short))
        alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))
        
        # Initialize rho
        rho_draw = generate_random_uniform(0, 1, (m, ninits_short))
        
        # Initialize mubeta
        minMU = np.min(y_c - xz[:,1:] @ out_coef[1:])
        maxMU = np.max(y_c - xz[:,1:] @ out_coef[1:])
        mubeta_draw = np.zeros((q1 * m, ninits_short))
        
        for mm in prange(m):
            mubeta_draw[mm, :] = np.random.uniform(minMU, maxMU, size=ninits_short)
            if q > 0:
                for j in prange(q):
                    mubeta_draw[(q1+1) * mm + j + 1, :] = out_coef[1+j] * np.random.uniform(-2, 2, size=ninits_short)
                    mubeta_draw[(q1+1) * mm + j + q + 1, :] = mubeta_draw[(q1+1) * mm + j + 1, :] * rho_draw[mm,:]    
        an = 1 / n    
        sigma_0 = np.full(m, stdR)
    
        # Initialize sigma
        sigma_draw = generate_random_uniform(0.01, 1, (m, ninits_short)) * stdR
        
        alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization_AR1(y_c, xz, y_0, xz_0, p, q, sigma_0, alpha_draw, rho_draw, mubeta_draw, sigma_draw, gamma_draw, m, n, t, an, maxit=maxit_short, tol=epsilon_short)

        components = np.argsort(penloglikset)[::-1][:ninits]
        alpha_draw = alpha_draw[:,components]
        mubeta_draw = mubeta_draw[:,components]
        sigma_draw = sigma_draw[:,components]
        gamma_draw = gamma_draw[:,components]
        
        
        alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=maxit, tol=epsilon)
        
        index = np.argmax(penloglikset)
        alpha_hat = alpha_draw[:,index]
        mubeta_hat = mubeta_draw[:,index]
        sigma_hat = sigma_draw[:,index]
        gamma_hat = gamma_draw[:,index]
        post = post[:, index]
        penloglik = penloglikset[index]
        loglik = loglikset[index]
        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        
        # Create a Numba-compatible dictionary to store results
        result_dict = Dict.empty(
            key_type=types.unicode_type,  # Keys are strings
            value_type=types.float64[:, :],  # Values are 2D arrays
        )
        
        result_dict['penloglik'] = np.array([[penloglik]])
        result_dict['loglik'] = np.array([[loglik]])
        result_dict['aic'] = np.array([[aic]])
        result_dict['bic'] = np.array([[bic]])
        
        result_dict['alpha_hat']  = alpha_hat[np.newaxis,:]
        result_dict['sigma_hat']  = sigma_hat[np.newaxis,:]
        result_dict['mubeta_hat'] = mubeta_hat[np.newaxis,:]
        result_dict['gamma_hat'] = gamma_hat[np.newaxis,:]
        
        return result_dict
    
# %%
def EM_optimization_mixture(y, x, z, p, q, sigma_0, alpha_draw, tau_draw, mubeta_draw, sigma_draw, gamma_draw, m, k, t, an, maxit=1000, tol= 1e-8, epsilon =0.05):
    
    nt = len(y)
    n = nt // t
    mk = int(m*k)
    ninits = alpha_draw.shape[1]
    if q == 0:
        x1 = np.ones((nt, 1))
        q1 = 1
    else:
        x1 = np.zeros((nt, x.shape[1] + 1))
        x1[:, 0] = 1  # Add intercept
        x1[:, 1:] = x
        q1 = x1.shape[1]

    post = np.zeros((m * n, ninits))
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    for jn in prange(ninits):        
        
        
        alpha_jn = alpha_draw[:, jn]
        tau_jn   = tau_draw[:, jn]
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64

        mubeta_jn_mat = mubeta_jn.reshape((q+1,m*k)).T
        
        oldpenloglik = -np.inf
        diff = 1.0
        
        for iter_ii in prange(maxit):
            
            ll = - nt * M_LN_SQRT_2PI
            
            if p > 0:
                ytilde = y - np.dot(z, gamma_jn)
            else:
                ytilde = y
            
            r = np.zeros((mk,nt)) 
            sigma_jn_rep = repeat_elements(sigma_jn, k)
            res = np.zeros((mk, nt))
            
            for j in prange(mk):
                res[j] = ytilde - x1 @ mubeta_jn_mat[j]
                r[j,:] = log_likelihood_array(res[j], 0.0, sigma_jn_rep[j])
                # r[j,:] = log_likelihood_array(res[j], 0.0, sigma_jn_rep[0])
                
            
            l_m = np.zeros((m,n))  # Same shape as `r`
            w_mk = np.zeros((mk,nt))    # Weights
                
            for mm in prange(m):
                r_m = r[ mm*k: (mm+1)*k, :]
                tau_m = tau_jn[ mm*k: (mm+1)*k] 
                minr = min_along_axis_0(r_m) 
                w_m = np.zeros((k,nt))    # Weights
                
                l_m_k = np.zeros((k,nt))  # Same shape as `r`                            
                for i in prange(nt):
                    for kk in prange(k):
                        l_m_k[kk,i] = max(min(tau_m[kk] * np.exp( r_m[kk,i] - minr[i] ), 1e10),1e-10)

                sum_l_m_k = np.zeros(nt)   # Sum along axis 0
                for i in prange(nt):
                    for kk in prange(k):
                        sum_l_m_k[i] += l_m_k[kk, i]
                
                for i in prange(nt):
                    for kk in prange(k):
                        w_m[kk, i] = l_m_k[kk, i] / sum_l_m_k[i]
                
                w_mk[mm*k: (mm+1)*k, :] = w_m
                
                # compute l_m 
                for nn in prange(n):
                    sum_l_m_k_nn = 0.0
                    for tt in prange(t):
                        idx = nn * t + tt  # Compute the flattened index
                        sum_l_m_k_nn += np.log(sum_l_m_k[idx]) + minr[idx]
                    l_m[mm,nn] = sum_l_m_k_nn
            
            sum_l_m = np.zeros(n)   # Sum along axis 0
            l_m_weighted = np.zeros((m,n))
            w = np.zeros((m, n))
            min_l_m = min_along_axis_0(l_m)
            for i in prange(n):
                for mm in prange(m):
                    l_m_weighted[mm, i] = max(min(alpha_jn[mm] * np.exp(l_m[mm,i] - min_l_m[i]), 1e10),1e-10)
                    sum_l_m[i] += l_m_weighted[mm, i]
                        
            # update weight
            # Compute w = l_j / sum_l_j
            for i in prange(n):
                for mm in prange(m):
                    w[mm, i] = l_m_weighted[mm, i] / sum_l_m[i]
                    w_mk[mm*k: (mm+1)*k, i * t: (i+1)*t] = w_mk[mm*k: (mm+1)*k, i * t: (i+1)*t] * w[mm, i]
            
            for i in prange(n):
                ll += np.log(sum_l_m[i]) + min_l_m[i]
            
            penloglik = ll + np.log(2.0) 
            
            # for mm in prange(m):
            #     s0j = sigma_0[mm] / sigma_jn[mm]
            #     penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
            
            # for j in prange(mk):
            #     penloglik += min(np.log(tau_jn[j]), np.log(1 - tau_jn[j]))
            
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            
            if np.max(np.abs(diff)) < tol:
                break
            
            # Fill w_mk nan values, for the program to run
            w_mk = fill_nan(w_mk, 1)
            # Make sure w_mk each row add up to 1
            for i in prange(nt):
                w_mk[:,i] = w_mk[:,i] / np.sum(w_mk[:,i])
            
            # Update parameters
            for mm in prange(m):
                
                alpha_jn[mm] = np.mean(w[mm])
                res_mm_sq = np.zeros(nt)
                w_mk_sum = np.zeros(nt)
                
                sum_tau_jn = 0 
                for kk in prange(k):
                    idx_type = mm * k + kk
                    xtilde = np.zeros((nt, q1))
                    for ii in prange(q1):
                        xtilde[:, ii] = w_mk[idx_type] * x1[:, ii]    
                    mubeta_jn_mat[idx_type,:] = solve_linear_system_safe(xtilde.T @ x1, xtilde.T @ ytilde)
                    res_mm_sq +=  w_mk[idx_type] * (ytilde - x1 @ mubeta_jn_mat[idx_type,:])**2
                    w_mk_sum += w_mk[idx_type]
                    tau_jn[idx_type] = np.mean(w_mk[idx_type] )
                    sum_tau_jn += tau_jn[idx_type]
                
                for kk in prange(k):
                    idx_type = mm * k + kk
                    tau_jn[idx_type] = tau_jn[idx_type] / sum_tau_jn
                
                sigma_jn[mm] = np.sqrt(( np.sum(res_mm_sq) ) / (np.sum(w_mk_sum) )  )
                sigma_jn[mm] = max(sigma_jn[mm], epsilon * sigma_0[mm])
            
            # update alpha
            for mm in prange(m):
                alpha_jn[mm] = min(max(0.01, alpha_jn[mm]), 0.99)
            
            total_alpha = np.sum(alpha_jn)
            for mm in prange(m):
                alpha_jn[mm] = alpha_jn[mm] / total_alpha
            # print(diff)
            # print(tau_jn)
            # print(mubeta_jn)
            # print(sigma_jn)
            
        # print(iter_ii)
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mubeta_draw[:, jn] = mubeta_jn_mat.T.flatten()
        sigma_draw[:, jn] = sigma_jn
        tau_draw[:, jn] = tau_jn
        
    return(alpha_draw,tau_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset ,post)

# m = M
# k = K
def regpanelmixmixturePMLE(y,x,z, p, q, m, k, ninits=2, epsilon=1e-6, maxit=2000, epsilon_short=1e-2, maxit_short=200):    # Extract the generated data
    
    t,n = y.shape
    nt = n * t
    y = y.T.flatten()
    
    if q > 0:
        x1 = np.hstack((np.ones((nt, 1)), x))
    else:
        x1 = np.ones((nt, 1))
    q1 = q + 1
    
    if p > 0:
        xz = np.hstack((x1, z))
    else:
        xz = x1 
        
    out_coef = solve_least_squares(xz, y)  # Replace np.linalg.lstsq
    residuals = y - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q1 + 1) * m + p
    # ninits_short = ninits * 10 * (q1 + p) * m
    ninits_short = ninits * 10 
    
    if (m == 1) & (k == 1):
        mubeta = out_coef[:q1]
        if p > 0:
            gamma = out_coef[q1:(q1 + p)]
        else:
            gamma = np.array([0.0])
        res = y - xz @ out_coef
        sigma = np.sqrt(np.mean(res**2))
        loglik = log_likelihood_normal(res,0,sigma)

        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        penloglik = loglik
        alpha = np.array([1])
        postprobs = np.ones(n)
        
    else: #either m >=2 or k >= 2
        
        # disregard p > 0  case
        # if p > 0:
        #     gamma = out_coef[q1:(q1 + p)]
        #     # Perform least squares regression with both x and z
        #     gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) * gamma
        #     mubeta_hat = out_coef[:q1]
        #     y = y - z @ gamma
        # else:
        gamma = np.array([0.0])
        mubeta_hat = out_coef[:-1]

        alpha_draw = generate_random_uniform(0, 1, (m, ninits_short))
        alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))

        # Draw mubeta and tau
        # minMU = np.min(y)
        # maxMU = np.max(y)
        mk = m * k 
        mubeta_draw = np.zeros((q1 * m * k, ninits_short))
        tau_draw = np.zeros((m * k, ninits_short))
        for mm in prange(m):
            tau_draw_m = generate_random_uniform(0, 1, (k, ninits_short))
            tau_draw_m = (tau_draw_m / np.sum(tau_draw_m, axis=0))
            tau_draw[(mm * k):((mm +1) * k),:] = tau_draw_m
            
            for kk in prange(k):
                idx_type = mm * k + kk 
                lb = compute_quantile(y, idx_type/(mk+1))
                ub = compute_quantile(y, (idx_type+1)/(mk+1))
                
                mubeta_draw[q1 * ( mm * k + kk ), :] = np.random.uniform(lb, ub, size=ninits_short)
        
        sigma_draw = generate_random_uniform(0.5, 2, (m, ninits_short)) * stdR
        gamma_draw = np.zeros((1,ninits_short), dtype=np.float64)
        an = 1 / n    
        sigma_0 = np.full(m, stdR)
        alpha_draw,tau_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization_mixture(y, x, z, p, q, sigma_0, alpha_draw, tau_draw, mubeta_draw, sigma_draw, gamma_draw, m, k, t, an, maxit=maxit_short, tol=epsilon_short)

        components = np.argsort(penloglikset)[::-1][:ninits]
        alpha_draw = alpha_draw[:,components]
        mubeta_draw = mubeta_draw[:,components]
        sigma_draw = sigma_draw[:,components]
        gamma_draw = gamma_draw[:,components]
        tau_draw = tau_draw[:,components]
        
        alpha_draw,tau_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization_mixture(y, x, z, p, q, sigma_0, alpha_draw, tau_draw, mubeta_draw, sigma_draw, gamma_draw, m, k, t, an,  maxit=maxit, tol=epsilon)
        
        index = np.argmax(penloglikset)
        alpha_hat = alpha_draw[:,index]
        tau_hat = tau_draw[:,index]
        mubeta_hat = mubeta_draw[:,index]
        sigma_hat = sigma_draw[:,index]
        gamma_hat = gamma_draw[:,index]
        post = post[:, index]
        penloglik = penloglikset[index]
        loglik = loglikset[index]
        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        
    # ----
    # Create a Numba-compatible dictionary to store results
    result_dict = Dict.empty(
        key_type=types.unicode_type,  # Keys are strings
        value_type=types.float64[:, :],  # Values are 2D arrays
    )
    result_dict['penloglik'] = np.array([[penloglik]])
    result_dict['loglik'] = np.array([[loglik]])
    result_dict['aic'] = np.array([[aic]])
    result_dict['bic'] = np.array([[bic]])
    
    result_dict['alpha_hat']  = alpha_hat[np.newaxis,:]
    result_dict['tau_hat']  = tau_hat[np.newaxis,:]
    result_dict['sigma_hat']  = sigma_hat[np.newaxis,:]
    result_dict['mubeta_hat'] = mubeta_hat[np.newaxis,:]
    result_dict['gamma_hat'] = gamma_hat[np.newaxis,:]
    return result_dict

# %%

@njit(parallel=True)
def compute_lr_BB(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T, m, p, q, BB):
    # Preallocate a matrix and initialize it with values
    mubeta_hat_mat = np.empty((m, q+1)) 
    mubeta_hat_mat[:, 0] = mu_hat
    mubeta_hat_mat[:, 1:] = beta_hat
    mubeta_hat = mubeta_hat_mat.T.flatten()

    # Generate data for bootstrap
    Data = [generate_data(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T, m, p, q) for _ in prange(BB)]
    
    # Preallocate lr_stat as a 1D array (Numba-compatible)
    lr_stat_bb = np.zeros(BB, dtype=np.float64)
    
    for bb in prange(BB):
        data = Data[bb]
        y_bb = data[0]
        x_bb = data[1]
        z_bb = data[2]
        
        # Call regpanelmixPMLE for m components
        out_h0 = regpanelmixPMLE(y_bb, x_bb, z_bb, p, q, m)
        out_h1 = regpanelmixPMLE(y_bb, x_bb, z_bb, p, q, m + 1)
        penloglik_h0 = out_h0['penloglik'][0, 0]
        penloglik_h1 = out_h1['penloglik'][0, 0]
        
        # Compute likelihood ratio statistic
        lr_stat_bb[bb] = -2 * (penloglik_h0 - penloglik_h1)
    
    return lr_stat_bb

@njit(parallel=True)
def LRTestParallel(Data, N, T, M, p, q, nrep, BB = 199):
    result_lr_each = np.zeros((nrep,1))
    for ii in prange(nrep): 
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        out_h0 = regpanelmixPMLE(y,x,z, p, q, M)
        out_h1 = regpanelmixPMLE(y,x,z, p, q, M+1)
        alpha_hat  = out_h0['alpha_hat'][0]
        mubeta_hat = out_h0['mubeta_hat'][0]
        sigma_hat  = out_h0['sigma_hat'][0]
        gamma_hat  = out_h0['gamma_hat'][0]
        penloglik_h0 = out_h0['penloglik'][0,0]
        penloglik_h1 = out_h1['penloglik'][0,0]
        # print(out_h0[0])        
        lr_stat = -2 * (penloglik_h0 - penloglik_h1)
        
        mubeta_hat = np.ascontiguousarray(mubeta_hat)
        mubeta_hat_mat = mubeta_hat.reshape((q+1,M)).T
        beta_hat = mubeta_hat_mat[:,1:]
        mu_hat =  mubeta_hat_mat[:,0]

        lr_stat_bb =  compute_lr_BB(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T, M, p, q, BB)
        lr_95 = compute_quantile(lr_stat_bb, 0.95)
        result_lr_each[ii,0] = 1 * ( lr_stat >  lr_95)
    return(result_lr_each)

# %%

# @njit
# def reshape_example(mubeta_hat, q, M):
#     # Ensure the array is contiguous
#     mubeta_hat = np.ascontiguousarray(mubeta_hat)
#     # Perform the reshape
#     mubeta_hat_mat = mubeta_hat.reshape((q + 1, M)).T
#     return mubeta_hat_mat