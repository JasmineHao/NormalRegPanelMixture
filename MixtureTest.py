# %%
import numpy as np
from numba import njit, prange
from numba.typed import Dict, List
from numba.core import types

# %%
# Functions for Numba
# ----------------------------------------------------------
@njit
def invert_matrix(mat, epsilon=1e-6):
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
    for j in range(cols):
        # Initialize the minimum value for the current column
        min_val = r[0, j]
        
        # Iterate through each row in the current column
        for i in range(1, rows):
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
    for i in range(n_rows):
        max_values[i] = -np.inf  # Initialize with negative infinity
        for j in range(n_cols):
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
    for i in range(n_rows):
        row_sum = 0.0
        for j in range(n_cols):
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
    for i in range(size[0]):
        for j in range(size[1]):
            out[i, j] = low + (high - low) * np.random.random()
    return out


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
    for m in range(M):
        alpha_cum[m + 1] = alpha_cum[m] + alpha[m]
    
    if M > 1:
        for m in range(M):
            lb = alpha_cum[m]
            ub = alpha_cum[m + 1]
            for n in range(N):
                R[n, m] = 1 if lb < prior[n] <= ub else 0
    else:
        R[:] = 1

    # Initialize output arrays
    Y = np.zeros((T, N))
    
    # Generate x and z if not provided
    if q > 0:
        x = np.empty((N * T, q))
        for i in range(N * T):
            for j in range(q):
                x[i, j] = np.random.normal()  # Generate one value at a time
    else:
        x = np.zeros((N * T, 1), dtype=np.float64)
        
    if p > 0:
        z = np.empty((N * T, p))
        for i in range(N * T):
            for j in range(p):
                z[i, j] = np.random.normal()  # Generate one value at a time
    else:
        z = np.zeros((N * T, 1), dtype=np.float64)
        
    # Precompute dot products
    mu_R = np.dot(R, mu)  # Use np.dot for matrix multiplication
    sigma_R = np.dot(R, sigma)  # Use np.dot for matrix multiplication
    beta_R = np.dot(R, beta) 
   
    # Generate u array (workaround for np.random.normal with size)
    u = np.empty((T, N))
    for t in range(T):
        for n in range(N):
            u[t, n] = np.random.normal()  # Generate one value at a time

    # Generate Y
    for nn in range(N):
        y_nn = np.zeros(T)
        y_nn = mu_R[nn] + sigma_R[nn] * u[:, nn]
        
        y_nn += x[(T * nn):(T * (nn + 1)), :] @ beta_R[nn, :]
        y_nn += z[(T * nn):(T * (nn + 1)), :] @ gamma
        
        Y[:, nn] = y_nn    
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
    for t in range(T):
        # Calculate quantiles manually
        quantiles = np.empty(n_bins + 1)
        sorted_data = np.sort(data_c[t, :])
        for i in range(n_bins + 1):
            if i == 0:
                quantiles[i] = -np.inf
            elif i == n_bins:
                quantiles[i] = np.inf
            else:
                quantiles[i] = sorted_data[int(i * N / n_bins)]

        # Create indicator matrix
        indicator_matrix = np.zeros((N, n_bins))
        for n in range(N):
            for b in range(n_bins):
                if quantiles[b] <= data_c[t, n] < quantiles[b + 1]:
                    indicator_matrix[n, b] = 1
                    break
        indicator_list.append(indicator_matrix)
    return indicator_list


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
    for k in range(T):
        # Compute the Kronecker product for each row manually
        data_partition_matrix_y = np.zeros((N, (n_bins ** (T - 1))))
        for n in range(N):
            # Manually compute Kronecker product for the row
            kron_result = np.array([1.0])  # Start with scalar 1.0
            for t in range(T):
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
    for ii in range(nrep):
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
        for k in range(T):
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
        
        for i in range(BB):
            # Calculate bootstrapped P and Sigma_P matrices
            data_P_W_b = calculate_P_matrix(data_c, ru[i, :], n_grid=n_grid, n_bins=2)
            
            for k in range(T):
                P_k = data_P_W_b['P_k_list'][k]
                Sigma_P_k = data_P_W_b['Sigma_P_k_list'][k]
                # Compute KP statistics for the k-th triplet
                rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, lambda_c_list[k])['rk_c'][0,0]
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

@njit
def calculate_P_matrix_with_covariates(y, x, weights, ngrid=3, n_bins=2):
    T = y.shape[0]
    N = y.shape[1]
    indicator_list_Y = create_indicator_list(y, T, N, n_bins=n_bins)
    indicator_list_Y_ngrid = create_indicator_list(y, T, N, n_bins=n_grid)
    
    indicator_list_X = create_indicator_list(x, T, N, n_bins=n_bins)
    
    # Initialize the result lists
    P_k_list = List()
    Sigma_P_k_list = List()
    
    n_partition = (n_bins ** (T - 1))
    # Iterate over the t periods
    for k in range(T):
        # Compute the Kronecker product for each row manually
        data_partition_matrix_y = np.zeros((N, (n_bins ** (T - 1))))
        data_partition_matrix_x = np.zeros((N, (n_bins ** (T - 1))))    
    
        for n in range(N):
            # Manually compute Kronecker product for the row
            kron_result_y = np.array([1.0])  # Start with scalar 1.0
            kron_result_x = np.array([1.0])  # Start with scalar 1.0
            for t in range(T):
                if t != k:
                    kron_result_y = np.kron(kron_result_y, indicator_list_Y[t][n, :])
                    kron_result_x = np.kron(kron_result_x, indicator_list_X[t][n, :])
            data_partition_matrix_y[n, :] = kron_result_y
            data_partition_matrix_x[n, :] = kron_result_x
        
        for kk in range(n_partition):
            data_partition_kk = data_partition_matrix_x[:,kk]
            weights_adj = data_partition_kk * weights 
            weights_adj_sum = max(np.sum(weights_adj), 1e-6)
            weights_adj = weights_adj / np.sum(weights_adj)
            # Compute P_k
            P_k = (weights_adj * indicator_list_Y_ngrid[k].T) @ data_partition_matrix_y
            P_k_list.append(P_k)
            # Compute Sigma_P_k
            P_k_vec = P_k.T.flatten()
            W_P_s = np.diag(P_k_vec) - np.outer(P_k_vec, P_k_vec)
            Sigma_P_k_list.append(W_P_s)
    
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
    for ii in range(nrep):        
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        # Compute P and Sigma matrices
        data_P_W = calculate_P_matrix_with_covariates(y, x, weights_equal, ngrid=3, n_bins=2)
        
        # Initialize results
        rk = np.zeros(n_stats)
        lambda_c_list = List()
        omega_c = List()
        Sigma_P_list = List()
        P_k_list = List()
        
        for k in range(n_stats):
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
        
        for bb in range(BB):
            data_P_W_b = calculate_P_matrix_with_covariates(y, x, ru[bb, :], ngrid=3, n_bins=2)
            
            for k in range(n_stats):
                P_k = data_P_W_b['P_k_list'][k]
                Sigma_P_k = data_P_W_b['Sigma_P_k_list'][k]
                # Compute KP statistics for the k-th triplet
                # rk_b[bb, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, lambda_c_list[k])['rk_c'][0,0]
                rk_b[bb, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, 0)['rk_c'][0,0]
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
    for i in range(len(y)):
        log_likelihoods[i] = (
            constant 
            - np.log(sigma) 
            - ((y[i] - mu) ** 2) / (2 * variance)
        )
    
    return log_likelihoods


@njit
def compute_residual_normal_reg(m, n, t, sigma_jn, ytilde, mubeta_jn):
    """
    Compute residuals for the EM optimization loop in a Numba-compatible way.

    Parameters:
    - m: Number of components (int)
    - n: Number of groups (int)
    - t: Number of time points per group (int)
    - sigma_jn: Array of current sigma values (1D array of floats, shape (m,))
    - ytilde: Adjusted response variable (1D array of floats, shape (n * t,))
    - mubeta_jn: Array of current beta means (1D array of floats, shape (m,))

    Returns:
    - r: Residuals array (2D array of floats, shape (m, n))
    """
    # Initialize the residuals array
    r = np.zeros((m, n), dtype=np.float64)

    # Loop over each component (m)
    for j in range(m):
        
        # Loop over each group (n)
        for i in range(n):
            sum_r_t = 0.0

            # Loop over each time point within the group (t)
            for k in range(t):
                idx = i * t + k  # Compute the flattened index
                diff = ytilde[idx] - mubeta_jn[j]
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
    post = np.zeros((m * n, ninits))
    notcg = np.zeros(ninits)
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    
    for jn in range(ninits):
        alpha_jn = alpha_draw[:, jn]
        mubeta_jn = mubeta_draw[:, jn]
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64
    
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0
        
        for iter_ii in range(maxit):
            ll = -nt * M_LN_SQRT_2PI
            
            if p > 0:
                ytilde = y - np.dot(z, gamma_jn)
            else:
                ytilde = y
            
            r = compute_residual_normal_reg(m, n, t, sigma_jn, ytilde, mubeta_jn)
            
            minr = min_along_axis_0(r)
            
            # Initialize arrays
            l_j = np.zeros((m,n))  # Same shape as `r`
            sum_l_j = np.zeros(n)   # Sum along axis 0
            w = np.zeros((m,n))    # Weights
            ll = 0.0                # Log-likelihood accumulator

            # Compute l_j = alpha_jn[:, None] * exp(minr - r)
            for i in range(n):
                for j in range(m):
                    l_j[j, i] = alpha_jn[j] * np.exp(minr[i] - r[j,i])
            
            # Compute sum_l_j = np.sum(l_j, axis=0)
            for j in range(m):
                for i in range(n):
                    sum_l_j[i] += l_j[j, i]
            
            # Compute w = l_j / sum_l_j
            for i in range(n):
                for j in range(m):
                    w[j, i] = l_j[j, i] / sum_l_j[i]
            
            # Compute ll += np.sum(np.log(sum_l_j) - minr)
            for i in range(n):
                ll += np.log(sum_l_j[i]) - minr[i]
            
            penloglik = ll + np.log(2.0) + min(np.log(tau), np.log(1 - tau))
            
            for j in range(m):
                s0j = sigma_0[j] / sigma_jn[j]
                penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
                penloglik += min(np.log(alpha_jn[j]), np.log(1 - alpha_jn[j]))
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            emit += 1
            
            # Update parameters
            mubeta_jn_mat = np.zeros((m,q1),dtype=np.float64)
            wtilde = np.zeros(nt)
            for j in range(m):
                alpha_jn[j] = np.mean(w[j, :])
                wtilde = w[j, :].T
                w_j = np.zeros(nt)
                for i in range(n):
                    w_j[i * t : (i + 1) * t] = wtilde[i]
                xtilde = np.zeros((nt, q1))
                for ii in range(q1):
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
            for j in range(m):
                alpha_jn[j] = max(0.01, alpha_jn[j] / total_alpha)
            
            # update gamma
            if p > 0:
                ztilde = np.zeros((nt, p), dtype=np.float64) 
                zz = np.zeros((p, p), dtype=np.float64) 
                ze = np.zeros((p, 1), dtype=np.float64) 
                for j in range(m):
                    wtilde = w[j, :]
                    w_j = np.zeros(nt)
                    for i in range(n):
                        w_j[i * t : (i + 1) * t] = wtilde[i]
                    for ii in range(p):
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
def regpanelmixPMLE(y,x,z, p, q, m, ninits=10, epsilon=1e-8, maxit=2000, epsilon_short=1e-2, maxit_short=500)  : 
    
    t,n = y.shape
    nt = n * t
    y = y.T.flatten()
    
    # y.reshape((n,t)).T - data_lr[0][0] # check equivalence
    # Handle x
    
    x1 = np.hstack((np.ones((nt, 1)), x))
    q1 = q + 1
    

    xz = np.hstack((x1, z))
    
    out_coef = solve_least_squares(xz, y)  # Replace np.linalg.lstsq
    residuals = y - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q1 + 1) * m + p
    ninits_short = ninits * 10 * (q1 + p) * m
    
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
            for j in range(m):
                mubeta_draw[q1 * j, :] = np.random.uniform(minMU, maxMU, size=ninits_short)
                for i in range(1, q1):
                    mubeta_draw[q1 * j + i, :] = mubeta_hat[i] * np.random.uniform(-2, 2, size=ninits_short)
        else:
            minMU = np.min(y)
            maxMU = np.max(y)
            mubeta_draw = np.zeros((q1 * m, ninits_short))
            for j in range(m):
                mubeta_draw[q1 * j, :] = np.random.uniform(minMU, maxMU, size=ninits_short)
        
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

# %%
@njit
def regpanelmixPMLE_Bootstrap(y,x,z, p, q, m, alpha_hat, mubeta_hat, sigma_hat, gamma_hat, test_h0 = 1, ninits=10, epsilon=1e-8, maxit=2000)  : 
    # The bootstrap computation should be faster
    t,n = y.shape
    nt = n * t
    y = y.T.flatten()
    
    # Handle x
    
    x1 = np.hstack((np.ones((nt, 1)), x))
    q1 = q + 1
    
    xz = np.hstack((x1, z))
    out_coef = solve_least_squares(xz, y)  # Replace np.linalg.lstsq
    residuals = y - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q1 + 1) * m + p
    
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
            # Perform least squares regression with both x and z
            gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits)) * gamma_hat
            y = y - z @ gamma
        else:
            # Perform least squares regression with x only
            gamma_draw = np.zeros((1,ninits), dtype=np.float64)

        if test_h0:
            # Initialize alpha
            alpha_draw = generate_random_uniform(0, 1, (m, ninits))
            alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))

            # Initialize mubeta
            
            mubeta_draw = np.empty((q1 * m, ninits))
            for jj in range((q1 * m)):
                mubeta_draw[jj, :] = np.random.uniform(0.5, 2, size=ninits) * mubeta_hat[jj]
            
            an = 1 / n    
            sigma_0 = np.full(m, stdR)
        
            # Initialize sigma
            sigma_draw = generate_random_uniform(0.01, 1, (m, ninits)) * stdR
            
            alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=maxit, tol=epsilon)
        
            
        else:
            # Initialize alpha
            alpha_draw = generate_random_uniform(0, 1, (m+1, ninits))
            alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))
            # Initialize mubeta
            
            mubeta_draw = np.empty((q1 * (m+1), ninits))
            for jj in range((q1 * (m+1))):
                mubeta_draw[jj, :] = np.random.uniform(0.5, 2, size=ninits) * mubeta_hat[jj]
            
            an = 1 / n    
            sigma_0 = np.full((m+1), stdR)
        
            # Initialize sigma
            sigma_draw = generate_random_uniform(0.01, 1, (m+1, ninits)) * stdR
            
            alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m+1, t, an, maxit=maxit, tol=epsilon)
        penloglik = np.max(penloglikset)
    return(penloglik)

# %%

@njit(parallel=True)
def compute_lr_BB(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T,  m, p, q, BB):
    
    mubeta_hat_mat = np.empty((m,q+1)) 
    mubeta_hat_mat[:,0] = mu_hat
    mubeta_hat_mat[:,1:] = beta_hat
    mubeta_hat = mubeta_hat_mat.T.flatten()
    Data = [generate_data(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T, m, p, q) for _ in range(BB)]
    
    # Preallocate lr_stat as a 1D array (Numba-compatible)
    lr_stat_bb = np.zeros(BB, dtype=np.float64)
    
    for bb in prange(BB):
        data = Data[bb]
        y_bb = data[0]
        x_bb = data[1]
        z_bb = data[2]
        
        # Call regpanelmixPMLE for m components
        penloglik = regpanelmixPMLE_Bootstrap(y_bb, x_bb, z_bb, p, q, m, alpha_hat, mubeta_hat, sigma_hat, gamma_hat)
        # Call regpanelmixPMLE for m+1 components
        penloglik_m1 = regpanelmixPMLE_Bootstrap(y_bb, x_bb, z_bb, p, q, m, alpha_hat, mubeta_hat, sigma_hat, gamma_hat,test_h0=0)
        # Compute likelihood ratio statistic
        lr_stat_bb[bb] = -2 * (penloglik_m1 - penloglik)
    return lr_stat_bb


@njit(parallel=True)
def LRTestParallel(Data, N, T, M, p, q, nrep, BB = 199):
    result_lr_each = np.zeros((nrep,1))
    for ii in range(nrep): 
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
        lr_stat_bb =  compute_lr_BB(alpha, mu, sigma, gamma, beta, N, T, M, p, q, BB)
        lr_95 = compute_quantile(lr_stat_bb, 0.95)
        result_lr_each[ii,0] = 1 * ( lr_stat >  lr_95)
    return(result_lr_each)
# q1 = q + 1yu709tr
# mubeta_hat_mat = mubeta_hat.reshape((q1,m)).T
# beta_hat = mubeta_hat_mat[:,1:]
# mu_hat =  mubeta_hat_mat[:,0]

# %%
# Simulation
# ------------------------------------------
   
import time
# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
alphaset = [np.array([0.5, 0.5]), np.array([0.2, 0.8])]
muset = [np.array([-1.0, 1.0]), np.array([-0.5, 0.5])]
sigmaset = [np.array([0.8, 1.2])]
beta = np.array([[1.0],[1.0]])
gamma = np.array([0.0])

# Test panel mixture
N = 200
T = 3
M = 2
p = 0
q = 0
nrep = 100
BB = 199

alpha = alphaset[0]
mu = muset[0]
sigma = sigmaset[0]
n_grid=3
r_test=2
# Generate data
weights_equal = np.full(N, 1 / N)

# Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]
# %%
 # Determine the total number of parameter combinations
total_combinations = len(alphaset) * len(muset)

# Initialize the result matrix to store mean results for each parameter combination
simulation_result_matrix = np.zeros((total_combinations, 2))

# Optional: Track parameter combinations (for debugging or analysis)
parameter_combinations = []
count = 0
# Loop over parameters
for alpha in alphaset:
    for mu in muset:       
        start_time = time.time()
        result_rk_each = np.zeros((nrep,2))
        Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]
        data_nopar = [data[0] for data in Data]
        # Nonparametric test
        result_rk_each = NonParTestParallel(data_nopar, N, T, M, p, q, nrep, n_grid, BB, r_test)
        
        # Compute the mean results across replications and store them
        simulation_result_matrix[count, :] = result_rk_each.mean(axis=0)
        
         # Print execution time for this parameter combination
        print(f"Execution time for alpha={alpha}, mu={mu}: {time.time() - start_time:.2f} seconds")

        # Increment the counter
        count += 1

# Print the final result matrix
print("Simulation Result Matrix:")
print(simulation_result_matrix)

np.savetxt("NonparTestNoCovariates.txt", simulation_result_matrix, delimiter=',', fmt='%d')


# %%

# Test panel mixture
N = 800
T = 3
M = 2
p = 0
q = 1
nrep = 100
BB = 199

alpha = alphaset[0]
mu = np.array([-5.0, 5.0])
sigma = np.array([0.5, 0.5])
beta = np.array([[0.0],[0.0]])
n_grid=3
r_test=2
# Generate data
weights_equal = np.full(N, 1 / N)

Data = [generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q) for _ in range(nrep)]

result_rk_each = NonParTestParallelCovariate(Data, N, T, M, p, q, nrep, n_grid=3, BB=BB, r_test=2, n_bins=2)

result_rk_each.mean(axis=0)

# %%
file = open('LRTestOutput.txt', 'w')


file.close() 
