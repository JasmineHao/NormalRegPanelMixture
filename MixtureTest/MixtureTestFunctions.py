import numpy as np
from numba import njit, prange
from numba.typed import Dict, List
from numba.core import types
import math
import pandas as pd
import time
import pyreadr

from numpy._typing._array_like import NDArray

np.set_printoptions(
    precision=3,
    threshold=None,
    edgeitems=None,
    linewidth=100,
    suppress=True,
    nanstr=None,
    infstr=None,
    formatter=None,
    sign=None,
    floatmode=None,
    legacy=None
)


# Functions for Numba
# ----------------------------------------------------------

@njit
def compute_quantiles(data, probs):
    """
    Compute quantiles manually in a Numba-compatible way.

    Parameters:
        data: np.ndarray
            Input data to compute quantiles.
        probs: list or np.ndarray
            List of probabilities for which quantiles are computed.

    Returns:
        np.ndarray:
            Quantiles corresponding to the given probabilities.
    """
    sorted_data = np.sort(data)
    n = len(sorted_data)
    quantiles = np.zeros(len(probs))
    
    for i, p in enumerate(probs):
        pos = p * (n - 1)  # Position in sorted data
        lower = int(np.floor(pos))
        upper = int(np.ceil(pos))
        
        if lower == upper:
            quantiles[i] = sorted_data[lower]
        else:
            # Linear interpolation
            quantiles[i] = sorted_data[lower] + (pos - lower) * (sorted_data[upper] - sorted_data[lower])
    
    return quantiles

@njit
def b_spline_basis(x, knots, degree, i):
    """
    Compute the B-spline basis function recursively in a Numba-compatible way.

    Parameters:
        x: np.ndarray
            Input values where the basis function is evaluated.
        knots: np.ndarray
            Array of knot positions.
        degree: int
            Degree of the B-spline.
        i: int
            Index of the basis function.

    Returns:
        np.ndarray:
            The B-spline basis function evaluated at x.
    """
    n = len(x)
    basis = np.zeros(n)

    if degree == 0:  # Base case: degree 0 basis functions
        for j in range(n):
            if knots[i] <= x[j] < knots[i + 1]:
                basis[j] = 1.0
        if i == len(knots) - 2:  # Include right boundary for the last knot interval
            basis[x == knots[-1]] = 1.0
    else:  # Recursive case
        for j in range(n):
            left = 0.0
            if knots[i + degree] != knots[i]:  # Avoid division by zero
                left = ((x[j] - knots[i]) / (knots[i + degree] - knots[i])) * b_spline_basis(np.array([x[j]]), knots, degree - 1, i)[0]

            right = 0.0
            if knots[i + degree + 1] != knots[i + 1]:  # Avoid division by zero
                right = ((knots[i + degree + 1] - x[j]) / (knots[i + degree + 1] - knots[i + 1])) * b_spline_basis(np.array([x[j]]), knots, degree - 1, i + 1)[0]

            basis[j] = left + right

    return basis

@njit
def generate_b_spline_basis(x, knots, degree):
    """
    Generate all B-spline basis functions for a given degree and knot sequence.

    Parameters:
        x: np.ndarray
            The input values where the basis functions are evaluated.
        knots: np.ndarray
            Array of knot positions.
        degree: int
            Degree of the B-spline.

    Returns:
        np.ndarray:
            A 2D array where each column is a B-spline basis function evaluated at x.
    """
    n_basis = len(knots) - degree - 1
    n = len(x)
    basis_matrix = np.zeros((n, n_basis))

    for i in range(n_basis):
        basis_matrix[:, i] = b_spline_basis(x, knots, degree, i)

    return basis_matrix


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
def max_along_axis_0(r):
    # Get the shape of the array
    rows, cols = r.shape
    
    # Initialize an array to store the minimum values for each column
    max_vals = np.empty(cols)
    
    # Iterate through each column
    for j in prange(cols):
        # Initialize the minimum value for the current column
        max_val = r[0, j]
        
        # Iterate through each row in the current column
        for i in prange(1, rows):
            if r[i, j] > max_val:
                max_val = r[i, j]
        
        # Store the minimum value for the column
        max_vals[j] = max_val
    
    return max_vals

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
        for t in range(T):
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
    
    for i in range(N):
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


@njit(parallel=False)
def generate_data(alpha, mu, sigma, gamma, beta, N, T, M, p, q, spline=False, z_input = np.zeros((0,0)), x_input = np.zeros((0,0))):
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
    
    # Generate x and z if not provided
    if q > 0:
        if (len(x_input) == 0):
            x = np.empty((N * T, q))
            for i in prange(N * T):
                for j in range(q):
                    x[i, j] = np.random.normal()  # Generate one value at a time
        else:
            x = x_input
    else:
        x = np.zeros((N * T, 0), dtype=np.float64)
    
    if spline:
        bs_degree = 2    
        x_spline = np.zeros((N*T, 0))
        for qq in range(q):
            quantiles = np.quantile(x[:,qq], [0.33, 0.67])
            knots = np.concatenate((
                np.array([quantiles[0]] * (bs_degree)),  # This is fine
                quantiles,
                np.array([quantiles[1]] * (bs_degree))  # This is fine
            ))
            # Generate B-spline basis columns
            bs_columns = generate_b_spline_basis(x[:,qq], knots, bs_degree)
            x_spline = np.concatenate((x_spline, bs_columns),axis=1)
    
    if p > 0:
        if (len(z_input) == 0):
            z = np.empty((N * T, p))
            for i in prange(N * T):
                for j in prange(p):
                    z[i, j] = np.random.normal()  # Generate one value at a time
        else:
            z = z_input
    else:
        z = np.zeros((N * T, 0), dtype=np.float64)
    
    
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
    Y = np.zeros((T, N))
    
    for nn in range(N):
        Y[:, nn] += mu_R[nn] + sigma_R[nn] * u[:, nn]
        if spline:
            Y[:, nn] += x_spline[(T * nn):(T * (nn + 1)), :] @ beta_R[nn, :]
        else:
            Y[:, nn] += x[(T * nn):(T * (nn + 1)), :] @ beta_R[nn, :]
        Y[:, nn] += z[(T * nn):(T * (nn + 1)), :] @ gamma
        
    return(Y, x, z)


# %%
@njit(parallel=False)
def generate_data_ar1(alpha, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0,  N, T, M, p, q, z_input = np.zeros((0,0)), x_input = np.zeros((0,0))):
    # 
    # First compute the stationary distribution 
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

    
    # Generate x and z if not provided
    if q > 0:
        if len(x_input)==0:
            x = np.empty((N * T, q))
            for i in range(N * T):
                for j in range(q):
                    x[i, j] = np.random.normal()  # Generate one value at a time
        else:
            x = x_input
    else:
        x = np.zeros((N * T, 0), dtype=np.float64)
        
    if p > 0:
        if len(z_input) == 0:
            z = np.empty((N * T, p))
            for i in range(N * T):
                for j in range(p):
                    z[i, j] = np.random.normal()  # Generate one value at a time
        else:
            z = z_input
    else:
        z = np.zeros((N * T, 0), dtype=np.float64)
        
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
    for tt in range(T):
        for nn in range(N):
            u[tt, nn] = np.random.normal()  # Generate one value at a time
    
    # Initialize output arrays
    Y = np.zeros((T, N))
    # Generate Y
    for nn in range(N):
        Y[0,nn] += mu_0_R[nn] + sigma_0_R[nn] * u[0,nn] 
        Y[0,nn] += x[(T * nn)] @ beta_0_R[nn,:]
        Y[0,nn] += z[(T * nn)] @ gamma
        for tt in range(T-1):
            Y[tt+1, nn] += mu_R[nn] * (1 - rho_R[nn]) + sigma_R[nn] * u[tt+1, nn]
            Y[tt+1, nn] += Y[tt, nn] * rho_R[nn]
            Y[tt+1, nn] += x[(T * nn + tt+1), :] @ beta_R[nn, :]
            Y[tt+1, nn] -= x[(T * nn + tt), :] @ (beta_R[nn, :] * rho_R[nn])            
            Y[tt+1, nn] += z[(T * nn + tt+1), :] @ gamma
            Y[tt+1, nn] -= z[(T * nn + tt), :] @ (gamma * rho_R[nn]) 
    return(Y, x, z)


# %%
@njit(parallel=False)
def generate_data_ar1_mixture(alpha, tau, rho, mu, sigma, beta, gamma, mu_0, sigma_0, beta_0, gamma_0,  N, T, M, K, p, q, spline=False, z_input = np.zeros((0,0)), x_input = np.zeros((0,0))):
    R = np.ones((N, M))
    R_T = np.ones((N*(T-1), M)) 
    
    alpha_sum = np.sum(alpha)
    if alpha_sum != 1:
        alpha = alpha / alpha_sum
    
    prior = np.random.random(size=N)
    alpha_cum = np.zeros(M + 1)
    for mm in range(M):
        alpha_cum[mm + 1] = alpha_cum[mm] + alpha[mm]
        lb = alpha_cum[mm]
        ub = alpha_cum[mm + 1]
        for n in range(N):
            R[n, mm] = 1 if lb < prior[n] <= ub else 0
            R_T[:,mm] = repeat_column(R, T-1, mm)
                

    R_sub = np.zeros((N * (T-1), M, K))
    R_sub_0 = np.ones((N, M, K))
    mu_R_sub = np.zeros((N * (T-1), M))
    mu_R_sub_0 = np.zeros((N, M))
    
    for mm in range(M):
        prior_m = np.random.random(size=N*(T-1))
        prior_m_0 = np.random.random(size=N)
        
        tau_cum = np.cumsum(np.concatenate((np.array([0]), tau[mm])))
        for kk  in range(K):
            lb = tau_cum[kk]
            ub = tau_cum[kk + 1]
            R_sub[:,mm, kk] = ((prior_m > lb) & (prior_m <= ub)) 
            R_sub_0[:, mm, kk] = ((prior_m_0 > lb) & (prior_m_0 <= ub)) 
        mu_R_sub[:, mm] = np.dot(R_sub[:,mm].astype(np.float64), mu[mm])
        mu_R_sub_0[:, mm] = np.dot(R_sub_0[:,mm].astype(np.float64), mu_0[mm])
    
    mu_R = (R_T * mu_R_sub).sum(axis=1)
    mu_0_R = (R * mu_R_sub_0).sum(axis=1)  # Use np.dot for matrix multiplication
    rho_R = np.dot(R, rho)  # Use np.dot for matrix multiplication
    sigma_R = np.dot(R, sigma)  # Use np.dot for matrix multiplication
    sigma_0_R = np.dot(R, sigma_0)  # Use np.dot for matrix multiplication
    beta_R = np.dot(R, beta) 
    beta_0_R = np.dot(R, beta_0) 
    
    # Generate x and z if not provided
    if q > 0:
        if len(x_input) == 0:
            x = np.empty((N * T, q))
            for i in range(N * T):
                for j in range(q):
                    x[i, j] = np.random.normal()  # Generate one value at a time
        else:
            x = x_input
    else:
        x = np.zeros((N * T, 0), dtype=np.float64)
        
    if p > 0:
        if len(z_input) == 0:
            z = np.empty((N * T, p))
            for i in range(N * T):
                for j in range(p):
                    z[i, j] = np.random.normal()  # Generate one value at a time
        else:
            z = z_input 
    else:
        z = np.zeros((N * T, 0), dtype=np.float64)
    
    if spline:
        bs_degree = 2    
        x_spline = np.zeros((N*T, 0))
        for qq in range(q):
            quantiles = np.quantile(x[:,qq], [0.33, 0.67])
            knots = np.concatenate((
                np.array([quantiles[0]] * (bs_degree)),  # This is fine
                quantiles,
                np.array([quantiles[1]] * (bs_degree))  # This is fine
            ))
            # Generate B-spline basis columns
            bs_columns = generate_b_spline_basis(x[:,qq], knots, bs_degree)
            x_spline = np.concatenate((x_spline, bs_columns),axis=1)
        
    # Generate u array, including T=0 (workaround for np.random.normal with size)
    u = np.empty((T, N))
    for tt in range(T):
        for nn in range(N):
            u[tt, nn] = np.random.normal()  # Generate one value at a time
    
    # Initialize output arrays
    Y = np.zeros((T, N))
    # Generate Y
    for nn in range(N):
        Y[0,nn] += mu_0_R[nn] + sigma_0_R[nn] * u[0,nn] 
        Y[0,nn] += x[(T * nn)] @ beta_0_R[nn,:]
        Y[0,nn] += z[(T * nn)] @ gamma
        for tt in range(T-1):
            Y[tt+1, nn] += mu_R[ ( (T -1) * nn + tt)] * (1 - rho_R[nn]) + sigma_R[nn] * u[tt+1, nn]
            Y[tt+1, nn] += Y[tt, nn] * rho_R[nn] + x[(T * nn + tt+1), :] @ beta_R[nn, :] - x[(T * nn + tt), :] @ (beta_R[nn, :] * rho_R[nn])      
            Y[tt+1, nn] += z[(T * nn + tt+1), :] @ gamma - z[(T * nn + tt), :] @ (gamma * rho_R[nn]) 
    return(Y, x, z)

# %%
@njit(parallel=False)
def generate_data_mixture(alpha, mu, beta, sigma, tau, gamma, N, T, M, K, p=0, q=0, spline=False, z_input = np.zeros((0,0)), x_input = np.zeros((0,0))):

    R = np.ones((N, M))
    R_T = np.ones((N*T, M)) 
    alpha_sum = np.sum(alpha)
    if alpha_sum != 1:
        alpha = alpha / alpha_sum
    
    prior = np.random.random(size=N)
    alpha_cum = np.zeros(M + 1)
    for mm in range(M):
        alpha_cum[mm + 1] = alpha_cum[mm] + alpha[mm]
        lb = alpha_cum[mm]
        ub = alpha_cum[mm + 1]
        for n in range(N):
            R[n, mm] = 1 if lb < prior[n] <= ub else 0
            R_T[:,mm] = repeat_column(R, T, mm)
    
    R_sub = np.zeros((N * T, M, K))
    mu_R_sub = np.zeros((N*T, M))
    
    for mm in range(M):
        prior_m = np.random.random(size=N*T)
        
        tau_cum = np.cumsum(np.concatenate((np.array([0]), tau[mm])))

        for kk  in range(K):
            lb = tau_cum[kk]
            ub = tau_cum[kk + 1]
            R_sub[:,mm, kk] = ((prior_m > lb) & (prior_m <= ub)) 
        mu_R_sub[:, mm] = np.dot(R_sub[:,mm].astype(np.float64), mu[mm])
    
    sigma_R = np.dot(R, sigma)  # Use np.dot for matrix multiplication
    mu_R = (R_T * mu_R_sub).sum(axis=1)
    beta_R = np.dot(R, beta) 
    
    if q > 0:
        if len(x_input)==0:
            x = np.empty((N * T, q))
            for i in range(N * T):
                for j in range(q):
                    x[i, j] = np.random.normal()  # Generate one value at a time
        else:
            x = x_input
    else:
        x = np.zeros((N * T, 0), dtype=np.float64)
    
    if spline:
        bs_degree = 2    
        x_spline = np.zeros((N*T, 0))
        for qq in range(q):
            quantiles = np.quantile(x[:,qq], [0.33, 0.67])
            knots = np.concatenate((
                np.array([quantiles[0]] * (bs_degree)),  # This is fine
                quantiles,
                np.array([quantiles[1]] * (bs_degree))  # This is fine
            ))
            # Generate B-spline basis columns
            bs_columns = generate_b_spline_basis(x[:,qq], knots, bs_degree)
            x_spline = np.concatenate((x_spline, bs_columns),axis=1)
    
    if p > 0:
        if len(z_input) == 0:
            z = np.empty((N * T, p))
            for i in range(N * T):
                for j in range(p):
                    z[i, j] = np.random.normal()  # Generate one value at a time
        else:
            z = z_input
    else:
        z = np.zeros((N * T, 0), dtype=np.float64)
    
    u = np.empty((T, N))
    for nn in range(N):
        for tt in range(T):
            u[tt, nn] = np.random.normal()  # Generate one value at a time
    
    Y = np.zeros((T, N))
    
    for nn in range(N):
        for tt in range(T):
            Y[tt, nn] += mu_R[(T * nn + tt)] + sigma_R[nn] * u[tt, nn]
            if spline:
                Y[tt, nn] += x_spline[(T * nn + tt), :] @ beta_R[nn, :]
            else:
                Y[tt, nn] += x[(T * nn + tt), :] @ beta_R[nn, :]
            Y[tt, nn] += z[(T * nn + tt), :] @ gamma
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
        for n in range(N):
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
        
        for i in prange(BB):
            # Calculate bootstrapped P and Sigma_P matrices
            data_P_W_b = calculate_P_matrix(data_c, ru[i, :], n_grid=n_grid, n_bins=2)
            
            for k in range(T):
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




# %%
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
        for i in range(N):
            sum_r_t = 0.0

            # Loop over each time point within the group (t)
            for k in range(T):
                idx = i * t + k  # Compute the flattened index
                diff = res[j, idx]
                r_t = (1.0 / sigma_jn[j]) * diff
                sum_r_t += 0.5 * (r_t**2)

            # Compute residual for group i and component j
            r[j, i] = t * np.log(sigma_jn[j]) + sum_r_t
    return r


# %%
@njit
def EM_optimization(y_c, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=2000, tol=1e-8, tau = 0.5, epsilon=0.1, alpha_bound = 0.05):
    
    nt = len(y_c)
    n = nt // t

    ninits = alpha_draw.shape[1]
    # Handle x
    x1 = np.hstack((np.ones((nt, 1)), x))
    xz = np.hstack((x1, z))
    q1 = q + 1
        
    # Initialize variables
    
    l_j = np.zeros(m)
    
    
    post = np.zeros((m * n, ninits))
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    for jn in range(ninits):
        alpha_jn = alpha_draw[:, jn]
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64

        mubeta_jn_mat = mubeta_jn.reshape((q+1,m)).T
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        
        for iter_ii in range(maxit):
            
            if p > 0:
                ytilde = y_c - np.dot(z, gamma_jn)
            else:
                ytilde = y_c
            r = np.zeros((m,nt))
            for mm in range(m):
                res = ytilde - x1 @ mubeta_jn_mat[mm]
                r[mm] = log_likelihood_array(res, 0.0, sigma_jn[mm])
            
            
            r_sum = np.zeros((m,n))
            for mm in range(m):
                for nn in range(n):
                    for tt in range(t):
                        r_sum[mm,nn] += r[mm, nn * t + tt ]
            
            minr = max_along_axis_0(r_sum)
            
            # Initialize arrays
            l_j = np.zeros((m,n))  # Same shape as `r`
            sum_l_j = np.zeros(n)   # Sum along axis 0
            w = np.zeros((m,n))    # Weights
            ll = 0.0                # Log-likelihood accumulator

            # Compute l_j = alpha_jn[:, None] * exp(minr - r)
            for nn in range(n):
                for mm in range(m):
                    l_j[mm, nn] = alpha_jn[mm] * np.exp(r_sum[mm,nn] - minr[nn])
                    sum_l_j[nn] += l_j[mm, nn]
                    
            
            # Compute w = l_j / sum_l_j
            for nn in range(n):
                for mm in range(m):
                    w[mm, nn] = l_j[mm, nn]  / sum_l_j[nn]
            
            # Compute ll += np.sum(np.log(sum_l_j) - minr)
            for nn in range(n):
                ll += np.log(sum_l_j[nn]) + minr[nn]
            
            penloglik = ll + np.log(2.0) 
            
            for mm in range(m):
                s0j = sigma_0[mm] / sigma_jn[mm]
                penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
                # penloglik += min(np.log(alpha_jn[j]), np.log(1 - alpha_jn[j]))
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            emit += 1
            
            if abs(diff) < tol:
                break
            
            # Update parameters
            # mubeta_jn_mat = np.zeros((m,q1),dtype=np.float64)
            
            
            for mm in range(m):
                alpha_jn[mm] = np.mean(w[mm, :])
                wtilde = w[mm, :].T
                w_j = np.zeros(nt)
                for nn in range(n):
                    w_j[nn * t : (nn + 1) * t] = wtilde[nn]
                xtilde = np.zeros((nt, q1))
                for ii in range(q1):
                    xtilde[:, ii] = w_j * x1[:, ii]
                
                mubeta_jn_mat[mm,:] = solve_linear_system_safe(xtilde.T @ x1, xtilde.T @ ytilde)
                ssr_j = np.sum(w_j * (ytilde - x1 @ mubeta_jn_mat[mm,:])**2)
                sigma_jn[mm] = np.sqrt((ssr_j + 2.0 * an * sigma_0[mm]**2) / (np.sum(w_j) + 2.0 * an))
                sigma_jn[mm] = max(sigma_jn[mm], epsilon * sigma_0[mm])
            
            # update alpha
            for j in range(m):
                alpha_jn[j] = max(alpha_bound, alpha_jn[j] )
            
            total_alpha = np.sum(alpha_jn)
            for j in range(m):
                alpha_jn[j] = alpha_jn[j] / total_alpha
            
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
                    ze += ztilde.T @( y_c - x1 @ mubeta_jn_mat[j,:]) / max(sigma_jn[j]**2,0.01)
                gamma_jn = solve_linear_system_safe(zz,ze).flatten()
        
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mubeta_draw[:, jn] = mubeta_jn_mat.T.flatten()
        sigma_draw[:, jn] = sigma_jn
        
        if p > 0:
            gamma_draw[:, jn] = gamma_jn
    return(alpha_draw,mubeta_draw,sigma_draw,gamma_draw, penloglikset, loglikset ,post)


@njit
def regpanelmixPMLE(y,x,z, p, q, m, ninits=10, tol_long=1e-6, maxit=2000, tol_short=1e-2, maxit_short=200, alpha_bound=0.05): 
    
    t,n = y.shape
    nt = n * t
    y_c = y.T.flatten()
    
    # y.reshape((n,t)).T - data_lr[0][0] # check equivalence
    # Handle x
    
    x1 = np.hstack((np.ones((nt, 1)), x))
    xz = np.hstack((x1, z))
    q1 = q + 1
    
    out_coef = solve_least_squares(xz, y_c)  # Replace np.linalg.lstsq
    residuals = y_c - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q1 + 1) * m + p
    # ninits_short = ninits * 10 * (q1 + p) * m
    ninits_short = ninits * 10 
        
    # First draw random start point
    gamma_hat = out_coef[q1:(q1 + p)]
    mubeta_hat = out_coef[:q1]
    if p > 0:
        # Perform least squares regression with both x and z
        gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) 
        for pp in range(p):
            gamma_draw[:,pp] = gamma_draw[:,pp] * gamma_hat[pp]
    else:
        # Perform least squares regression with x only
        gamma_draw = np.zeros((0,ninits_short), dtype=np.float64)

    # Initialize alpha
    alpha_draw = generate_random_uniform(0, 1, (m, ninits_short))
    alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))

    mubeta_draw = np.zeros((q1 * m, ninits_short))
    y_mean = y_c - x @ mubeta_hat[1:]
    for mm in range(m):
        lb = compute_quantile(y_mean, mm/(m))
        ub = compute_quantile(y_mean, (mm+1)/(m))
        
        mubeta_draw[mm, :] = np.random.uniform(lb, ub, size=ninits_short)
        for qq in range(1, q1):
            mubeta_draw[qq * m + mm, :] = mubeta_hat[qq] * np.random.uniform(-2, 2, size=ninits_short)
    
    an = 1 / n    
    sigma_0 = np.full(m, stdR)

    # Initialize sigma
    sigma_draw = generate_random_uniform(0.01, 1, (m, ninits_short)) * stdR
    
    alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y_c, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=maxit_short, tol=tol_short,alpha_bound=alpha_bound)

    components = np.argsort(penloglikset)[::-1][:ninits]
    alpha_draw = alpha_draw[:,components]
    mubeta_draw = mubeta_draw[:,components]
    sigma_draw = sigma_draw[:,components]
    gamma_draw = gamma_draw[:,components]
    
    
    alpha_draw,mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization(y_c, x, z, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, m, t, an, maxit=maxit, tol=tol_long, alpha_bound=alpha_bound)
    
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
    result_dict['post'] = np.ascontiguousarray(post).reshape((n,m))
    
    return result_dict

# %%

@njit
def EM_optimization_AR1(y_c, xz, y_l, x_c, x_l, z_c, z_l, y_0, xz_0, p, q,  sigma_0, alpha_draw,  mubeta_draw, sigma_draw, gamma_draw, mubeta_0_draw, sigma_0_draw, gamma_0_draw, m, n, t, an, maxit=2000, tol=1e-8, epsilon=0.05):
    nt = n * (t-1)
    ninits = alpha_draw.shape[1]
    # Handle x
    q1 = 2 * q + 2
    # Initialize variables    
    l_j = np.zeros(m)
    w = np.zeros((m, n * (t-1)))
    post = np.zeros((m * n, ninits))
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)

    x1 = np.hstack((np.ones((nt, 1)), x_c, x_l, y_l[:,np.newaxis]))
    z1 = np.hstack((z_c, z_l))
    
    x1_0 = xz_0[:, :(q+1)]
    z_0 = xz_0[:, (q+1):]
    
    for jn in prange(ninits):    
        alpha_jn = alpha_draw[:, jn]    
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64
        
        mubeta_0_jn = np.ascontiguousarray(mubeta_0_draw[:, jn])
        sigma_0_jn = sigma_0_draw[:, jn]
        gamma_0_jn = gamma_0_draw[:, jn]  # Likely float64
        
        mubeta_0_jn_mat = mubeta_0_jn.reshape((q+1,m)).T             
        mubeta_jn_mat = mubeta_jn.reshape((q1,m)).T
        
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0
        
        for iter_ii in range(maxit):           
            r = np.zeros((m, nt))
            r_0 = np.zeros((m, n))
            for mm in range(m):
                res = y_c - x1 @ mubeta_jn_mat[mm] - z_c @ gamma_jn + z_l @ gamma_jn * mubeta_jn_mat[mm,-1]
                res_0 = y_0 - x1_0 @ mubeta_0_jn_mat[mm] - z_0 @ gamma_0_jn
                r[mm] = log_likelihood_array(res, 0.0, sigma_jn[mm])
                r_0[mm] = log_likelihood_array(res_0, 0.0, sigma_0_jn[mm])
            
            r_sum = np.zeros((m,n))
            for mm in range(m):
                for nn in range(n):
                    r_sum[mm,nn] += r_0[mm, nn]
                    for tt in range(t-1):
                        r_sum[mm,nn] += r[mm, nn * (t-1) + tt ]
            maxr = max_along_axis_0(r_sum)
            
            # Initialize arrays
            l_j = np.zeros((m,n))  # Same shape as `r`
            sum_l_j = np.zeros(n)   # Sum along axis 0
            w = np.zeros((m,n))    # Weights
            ll = 0.0                # Log-likelihood accumulator
            
            # Compute l_j = alpha_jn[:, None] * exp(minr - r)
            for nn in range(n):
                for mm in range(m):
                    l_j[mm, nn] = alpha_jn[mm] * np.exp( r_sum[mm,nn] - maxr[nn]) # add a cap to avoid overflow
                    sum_l_j[nn] += l_j[mm, nn]
                    
            # Compute w = l_j / sum_l_j
            for nn in range(n):
                for mm in range(m):
                    w[mm, nn] = l_j[mm, nn] / sum_l_j[nn]
            
            # Compute ll += np.sum(np.log(sum_l_j) - minr)
            for nn in range(n):
                ll += np.log(sum_l_j[nn]) + maxr[nn]
            
            penloglik = ll + np.log(2.0)
            
            # for j in range(m):
            #     s0j = sigma_0[j] / sigma_jn[j]
            #     penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
                
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            emit += 1
            
            if abs(diff) < tol:
                break
            # Update parameters
            
            for j in range(m):
                alpha_jn[j] = np.mean(w[j, :])
                wtilde = w[j, :].T
                
                # estimate param
                w_j = np.zeros(nt)
                for i in range(n):
                    w_j[i * (t-1) : (i + 1) * (t-1)] = wtilde[i]
                
                # xtilde = np.zeros(xz.shape)
                # for ii in range(xtilde.shape[1]):
                #     xtilde[:, ii] = wtilde * xz[:, ii]
            
                rho_estim = mubeta_jn_mat[j, -1]
                    
                # # Step 1: estimate beta given rho 
                y_tmp = y_c - rho_estim * y_l - z_c @ gamma_jn + z_l @ gamma_jn * rho_estim
                x_tmp = np.hstack(( (1 - rho_estim) * np.ones((nt, 1)), x_c - rho_estim * x_l)) 

                x_tmp_w = np.zeros(x_tmp.shape)
                for ii in range(x_tmp.shape[1]):
                    x_tmp_w[:, ii] = w_j * x_tmp[:, ii]
                
                mubeta_estim = solve_linear_system_safe(x_tmp_w.T @ x_tmp, x_tmp_w.T @ y_tmp)
                
                # # Step 2: estimate rho given beta
                y_tmp = y_c - np.hstack((np.ones((nt, 1)), x_c)) @ mubeta_estim - z_c @ gamma_jn
                x_tmp = y_l - np.hstack((np.ones((nt, 1)), x_l)) @ mubeta_estim - z_l @ gamma_jn
                
                if np.abs(np.sum(w_j)) < 1e-3:
                    w_j[w_j < 1e-5] = 1e-5
                
                rho_estim = np.sum(y_tmp * x_tmp * w_j) / np.sum(x_tmp * x_tmp * w_j)
                rho_estim = min(0.99, max(0.01, rho_estim))
                # end mubeta and rho update
                
                mubeta_jn_mat[j,0] = mubeta_estim[0] * (1 - rho_estim)
                mubeta_jn_mat[j,1:(q+1)] = mubeta_estim[1:(q+1)]
                mubeta_jn_mat[j,(q+1):(2*q+1)] = - mubeta_estim[1:(q+1)] * rho_estim
                mubeta_jn_mat[j,-1] = rho_estim
                
                                
                # estimate mubeta_0 for initial period
                xtilde_0 = np.zeros(x1_0.shape)
                for ii in range(x1_0.shape[1]):
                    xtilde_0[:, ii] = wtilde * x1_0[:, ii]
                
                coef_0_jn = solve_linear_system_safe(xtilde_0.T @ x1_0, xtilde_0.T @ (y_0 - z_0 @ gamma_0_jn))
                
                mubeta_0_jn_mat[j,:(q+1)] = coef_0_jn[:(q+1)]

                
                # print(coef_jn)
                ssr_j = np.sum(w_j * (y_c - x1 @ mubeta_jn_mat[j] - z_c @ gamma_jn + z_l @ gamma_jn * mubeta_jn_mat[j,-1])**2)
                ssr_j_0 = np.sum(wtilde * (y_0 - x1_0 @ mubeta_0_jn_mat[j] - z_0 @ gamma_0_jn)**2)
                sigma_jn[j] = np.sqrt((ssr_j + 2.0 * an * sigma_0[j]**2) / (np.sum(w_j) + 2.0 * an))
                sigma_jn[j] = max(sigma_jn[j], epsilon * sigma_0[j])
                
                sigma_0_jn[j] = np.sqrt((ssr_j_0 + 2.0 * an * sigma_0[j]**2) / (np.sum(wtilde) + 2.0 * an))
                sigma_0_jn[j] = max(sigma_0_jn[j], epsilon * sigma_0[j])
                
                
            mubeta_jn = mubeta_jn_mat.T.flatten() 
            mubeta_0_jn = mubeta_0_jn_mat[:,:q+1].T.flatten()
            
            for j in range(m):
                alpha_jn[j] = max(0.05, alpha_jn[j] )
            
            total_alpha = np.sum(alpha_jn)
            for j in range(m):
                alpha_jn[j] = alpha_jn[j] / total_alpha

            # update gamma
            if p > 0:
                ztilde = np.zeros((nt, p), dtype=np.float64) 
                zz = np.zeros((p, p), dtype=np.float64) 
                ze = np.zeros((p, 1), dtype=np.float64) 
                for j in range(m):
                    wtilde = w[j, :]
                    w_j = np.zeros(nt)
                    for i in range(n):
                        w_j[i * (t - 1) : (i + 1) * (t - 1)] = wtilde[i]
                    for ii in range(p):
                        ztilde[:, ii] = w_j * (z_c[:, ii] - z_l[:, ii] * mubeta_jn_mat[j,-1])
                    zz += ztilde.T @ (z_c[:, ii] - z_l[:, ii] * mubeta_jn_mat[j,-1]) / (sigma_jn[j]**2)
                    ze += ztilde.T @( y_c - x1 @ mubeta_jn_mat[j,:]) / max(sigma_jn[j]**2,0.01)
                gamma_jn = solve_linear_system_safe(zz,ze).flatten()
                   
        
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mubeta_draw[:, jn] = mubeta_jn_mat[:,:q1].T.flatten()
        sigma_draw[:, jn] = sigma_jn
        
        mubeta_0_draw[:, jn] = mubeta_0_jn_mat[:,:(q+1)].T.flatten()
        sigma_0_draw[:, jn] = sigma_0_jn
        if p > 0:
            gamma_draw[:, jn] = gamma_jn
            gamma_0_draw[:, jn] = gamma_0_jn
    return(alpha_draw,  mubeta_draw,sigma_draw,gamma_draw,mubeta_0_draw,sigma_0_draw,gamma_0_draw, penloglikset, loglikset ,post)

# %%

@njit
def regpanelmixAR1PMLE(y, x, z, p, q, m, ninits=10, tol_long=1e-6, maxit=2000, tol_short=1e-2, maxit_short=200)  : 
    t,n = y.shape
    nt = n * (t-1)
    
    y_l = y[:-1,:]
    y_0 = y[0,:]
    y_c = y[1:,:]
    
    y_c = y_c.T.flatten()
    y_l = y_l.T.flatten()
    
    # y.reshape((n,t)).T - data_lr[0][0] # check equivalence
    # Handle x
    x_0 = np.zeros((n,q))
    x_l = np.zeros((nt,q))
    x_c = np.zeros((nt,q))
    for j in range(q):
        x_j_mat = np.ascontiguousarray(x[:,j]).reshape((n,t)).T
        x_0[:,j] = x_j_mat[0,:]
        x_l[:,j] = x_j_mat[:-1,:].T.flatten()
        x_c[:,j] = x_j_mat[1:,:].T.flatten()
        
    # Handle x
    z_0 = np.zeros((n,p))
    z_l = np.zeros((nt,p))
    z_c = np.zeros((nt,p))
    for j in range(p):
        z_j_mat = np.ascontiguousarray(z[:,j]).reshape((n,t)).T
        z_0[:,j] = z_j_mat[0,:]
        z_l[:,j] = z_j_mat[:-1,:].T.flatten()
        z_c[:,j] = z_j_mat[1:,:].T.flatten()
        
    x1 = np.hstack((np.ones((nt, 1)), x_c, x_l, y_l[:,np.newaxis]))
    xz = np.hstack((x1, z_c, z_l))
    q1 = 2*q + 2
    
    xz_0 = np.hstack((np.ones((n, 1)), x_0, z_0) )
        
    out_coef = solve_least_squares(xz, y_c)  # Replace np.linalg.lstsq
    out_coef_0 = solve_least_squares(xz_0, y_0)  # Replace np.linalg.lstsq
    residuals = y_c - xz @ out_coef
    residuals_0 = y_0 - xz_0 @ out_coef_0
    stdR = np.std(residuals)
    stdR_0 = np.std(residuals_0)
    npar = m - 1 + (q1) * m + p 
    npar_0 =  (q + 1) * m + p 
    
    ninits_short = ninits * 5 * (q1 + p) * m
    
    # First draw random start point
    if p > 0:
        gamma_hat = out_coef[q1:(q1 + p)]
        gamma_0_hat = out_coef_0[(q+1):]
        # Perform least squares regression with both x and z
        gamma_0_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) 
        gamma_draw = generate_random_uniform(-1, 1, (p, ninits_short))
        
        for pp in range(p):
            gamma_draw[:,pp] = gamma_draw[:,pp] * gamma_hat[pp]
            gamma_0_draw[:,pp] = gamma_0_draw[:,pp] * gamma_0_hat[pp]
       
    else:
        # Perform least squares regression with x only
        gamma_draw = np.zeros((0,ninits_short), dtype=np.float64)
        gamma_0_draw = np.zeros((0,ninits_short), dtype=np.float64)
        
    # Initialize alpha
    alpha_draw = generate_random_uniform(0, 1, (m, ninits_short))
    alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))
    
    # Initialize rho and mu
    rho_draw = generate_random_uniform(0, 1, (m, ninits_short))
    mu_draw  = np.zeros((m, ninits_short))
    
    # Initialize mubeta
    minMU = np.min(y_c - xz[:,1:] @ out_coef[1:])
    maxMU = np.max(y_c - xz[:,1:] @ out_coef[1:])
    
    minMU_0 = np.min(y_0 - xz_0[:,1:] @ out_coef_0[1:])
    maxMU_0 = np.max(y_0 - xz_0[:,1:] @ out_coef_0[1:])
    
    mubeta_draw = np.zeros(( (q1) * m, ninits_short))
    mubeta_0_draw = np.zeros(( (q + 1) * m, ninits_short))
    
    for mm in range(m):
        mu_draw[mm, :] = np.random.uniform(minMU, maxMU, size=ninits_short) #assign mu
        mubeta_draw[mm, :] = mu_draw[mm,:] * (1 - rho_draw[mm,:])
        #assign rho
        mubeta_draw[m*(q1 - 1)+mm, :] = rho_draw[mm,:]
        
        mubeta_0_draw[mm, :] = np.random.uniform(minMU_0, maxMU_0, size=ninits_short)
        if q > 0:
            for j in range(q):
                mubeta_draw[ m*(j+1)+mm , :] = out_coef[1+j] * np.random.uniform(-0.5, 0.5, size=ninits_short)
                mubeta_draw[ m*(j+q+1)+mm , :] = mubeta_draw[ m*(j+1)+mm , :] * mubeta_draw[m*(q1 - 1)+mm, :]  
                mubeta_0_draw[ m*(j+1)+mm , :] = out_coef_0[1+j] * np.random.uniform(-0.5, 0.5, size=ninits_short)
                                
    an = 1 / n    
    sigma_0 = np.full(m, stdR)

    # Initialize sigma
    sigma_draw = generate_random_uniform(0.5, 1.5, (m, ninits_short)) * stdR
    sigma_0_draw = generate_random_uniform(0.5, 1.5, (m, ninits_short)) * stdR_0
    
    alpha_draw,  mubeta_draw,sigma_draw,gamma_draw,mubeta_0_draw,sigma_0_draw,gamma_0_draw ,penloglikset, loglikset, post = EM_optimization_AR1(y_c, xz, y_l, x_c, x_l, z_c, z_l, y_0, xz_0, p, q, sigma_0, alpha_draw, mubeta_draw, sigma_draw, gamma_draw, mubeta_0_draw, sigma_0_draw, gamma_0_draw, m, n, t, an, maxit=maxit_short, tol=tol_short)

    components = np.argsort(penloglikset)[::-1][:ninits]
    alpha_draw = alpha_draw[:,components]
    mubeta_draw = mubeta_draw[:,components]
    sigma_draw = sigma_draw[:,components]
    gamma_draw = gamma_draw[:,components]
    mubeta_0_draw = mubeta_0_draw[:,components]
    sigma_0_draw = sigma_0_draw[:,components]
    gamma_0_draw = gamma_0_draw[:,components]
    
    alpha_draw,  mubeta_draw,sigma_draw,gamma_draw,mubeta_0_draw,sigma_0_draw,gamma_0_draw ,penloglikset, loglikset, post = EM_optimization_AR1(y_c, xz, y_l, x_c, x_l, z_c, z_l, y_0, xz_0, p, q,  sigma_0, alpha_draw,  mubeta_draw, sigma_draw, gamma_draw, mubeta_0_draw, sigma_0_draw, gamma_0_draw, m, n, t, an, maxit=maxit, tol=tol_long)

    index = np.argmax(penloglikset)
    alpha_hat = alpha_draw[:,index]
    mubeta_mat =  np.ascontiguousarray(mubeta_draw[:,index]).reshape((q1,m)).T
    rho_hat = mubeta_mat[:,-1]
    mu_hat = mubeta_mat[:,0] / (1 - rho_hat)
    
    sigma_hat = sigma_draw[:,index]
    gamma_hat = gamma_draw[:,index]        
    mubeta_0_hat = mubeta_0_draw[:,index]
    
    sigma_0_hat = sigma_0_draw[:,index]
    gamma_0_hat = gamma_0_draw[:,index]
    post = post[:, index]
    penloglik = penloglikset[index]
    loglik = loglikset[index]
    aic = -2 * loglik + 2 * (npar + npar_0)
    bic = -2 * loglik + np.log(n) * (npar + npar_0)
    
    # Create a Numba-compatible dictionary to store results
    result_dict = Dict.empty(
        key_type=types.unicode_type,  # Keys are strings
        value_type=types.float64[:, :],  # Values are 2D arrays
    )
    
    result_dict['penloglik'] = np.array([[penloglik]])
    result_dict['loglik'] = np.array([[loglik]])
    result_dict['aic'] = np.array([[aic]])
    result_dict['bic'] = np.array([[bic]])
    result_dict['post'] = np.ascontiguousarray(post).reshape((n,m))
    
    sorted_indices = np.argsort(mu_hat)

    mubeta_hat = mubeta_mat[sorted_indices,:(q+1)].T.flatten() 
    mubeta_0_mat =  np.ascontiguousarray(mubeta_0_hat).reshape((q+1,m)).T
    
    mubeta_0_hat = mubeta_0_mat[sorted_indices,:(q+1)].T.flatten() 
    
    result_dict['alpha_hat']  = alpha_hat[np.newaxis,sorted_indices]
    result_dict['rho_hat']  = rho_hat[np.newaxis,sorted_indices]
    result_dict['mu_hat']  = mu_hat[np.newaxis,sorted_indices]
    result_dict['sigma_hat']  = sigma_hat[np.newaxis,sorted_indices]
    
    result_dict['mubeta_hat'] = mubeta_hat[np.newaxis,:]
    result_dict['gamma_hat'] = gamma_hat[np.newaxis,:]
    
    result_dict['sigma_0_hat']  = sigma_0_hat[np.newaxis,sorted_indices]
    result_dict['mubeta_0_hat'] = mubeta_0_hat[np.newaxis,:]
    result_dict['gamma_0_hat'] = gamma_0_hat[np.newaxis,:]
    return result_dict


# %%

@njit
def EM_optimization_AR1_mixture(y_c,  y_l, x_c, x_l, z_c, z_l,  xz, y_0, xz_0, p, q, sigma_0, alpha_draw, tau_draw, mu_draw, mubeta_draw,sigma_draw,gamma_draw,mu_0_draw, mubeta_0_draw,sigma_0_draw,gamma_0_draw, m, k, n, t, an, maxit=2000, tol=1e-8, epsilon=0.1):
    
    nt = n * (t-1)
    mk = int(m*k)
    ninits = alpha_draw.shape[1]
    # Handle x
    q1 = 2 * q + 2
    # Initialize variables    
    w = np.zeros((m, n * (t-1)))
    post = np.zeros((m * n, ninits))
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    x1 = np.hstack((np.ones((nt, 1)), x_c, x_l, y_l[:,np.newaxis]))
    z1 = np.hstack((z_c, z_l))
    
    x1_0 = xz_0[:, :(q+1)]
    z_0 = xz_0[:, (q+1):]
    for jn in prange(ninits):
        
        alpha_jn = alpha_draw[:, jn]
        tau_jn = tau_draw[:, jn]
        
        mu_jn = mu_draw[:, jn]
        mu_0_jn = mu_0_draw[:, jn]
        
        
        sigma_jn = sigma_draw[:, jn]
        sigma_0_jn = sigma_0_draw[:, jn]
        
        gamma_jn = gamma_draw[:, jn]  # Likely float64
        gamma_0_jn = gamma_0_draw[:, jn]  # Likely float64
        
        
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        mubeta_0_jn = np.ascontiguousarray(mubeta_0_draw[:, jn])
        
        mubeta_0_jn_mat = mubeta_0_jn.reshape((q+1,m)).T             
        mubeta_jn_mat = mubeta_jn.reshape((q1,m)).T
        
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0
        
        for iter_ii in range(maxit):
            
            r = np.zeros((mk, nt))
            r_0 = np.zeros((mk, n))
            for mm in range(m):
                for kk in range(k):
                    
                    res =  y_c - x1 @ mubeta_jn_mat[mm] - z_c @ gamma_jn + z_l @ gamma_jn * mubeta_jn_mat[mm,-1] - mu_jn[mm*k + kk] * ( 1 - mubeta_jn_mat[mm,-1])
                    res_0 = y_0 - x1_0 @ mubeta_0_jn_mat[mm] - z_0 @ gamma_0_jn - mu_0_jn[mm*k + kk]
                    r[mm*k + kk] = log_likelihood_array(res, 0.0, sigma_jn[mm])
                    r_0[mm*k + kk] = log_likelihood_array(res_0, 0.0, sigma_0_jn[mm])
            
            l_m = np.zeros((m,n))  # Same shape as `r`
            w_mk = np.zeros((mk,nt))    # Weights
            w_mk_0 = np.zeros((mk,n))    # Weights
            ll = 0.0
            
            
            for mm in range(m):
                r_m = r[ mm*k: (mm+1)*k, :]
                r_m_0 = r_0[ mm*k: (mm+1)*k, :]
                tau_m = tau_jn[ mm*k: (mm+1)*k] 
                minr = max_along_axis_0(r_m)
                minr_0 = max_along_axis_0(r_m_0)
                 
                w_m = np.zeros((k,nt))    # Weights
                l_m_k = np.zeros((k,nt))  # Same shape as `r`
                
                w_m_0 = np.zeros((k,n))    # Weights
                l_m_k_0 = np.zeros((k,n))  # Same shape as `r`
                                            
                for i in range(nt):
                    for kk in range(k):
                        l_m_k[kk,i] = tau_m[kk] * np.exp( r_m[kk,i] - minr[i] )
                
                for i in range(n):
                    for kk in range(k):
                        l_m_k_0[kk,i] = tau_m[kk] * np.exp( r_m_0[kk,i] - minr_0[i] )
                        
                sum_l_m_k = np.zeros(nt)   # Sum along axis 0
                for i in range(nt):
                    for kk in range(k):
                        sum_l_m_k[i] += l_m_k[kk, i]
                
                sum_l_m_k_0 = np.zeros(n)   # Sum along axis 0
                for i in range(n):
                    for kk in range(k):
                        sum_l_m_k_0[i] += l_m_k_0[kk, i]
                
                for i in range(nt):
                    for kk in range(k):
                        w_m[kk, i] = l_m_k[kk, i] / max(sum_l_m_k[i], 1e-6)
                
                for i in range(n):
                    for kk in range(k):
                        w_m_0[kk, i] = l_m_k_0[kk, i] / max(sum_l_m_k_0[i], 1e-6)
                
                w_mk[mm*k: (mm+1)*k, :] = w_m
                w_mk_0[mm*k: (mm+1)*k, :] = w_m_0
                
                # compute l_m 
                for nn in range(n):
                    sum_l_m_k_nn = 0.0
                    sum_l_m_k_nn += np.log(sum_l_m_k_0[nn]) + minr_0[nn]
                    for tt in range(t-1):
                        idx = nn * (t-1) + tt  # Compute the flattened index
                        sum_l_m_k_nn += np.log(sum_l_m_k[idx]) + minr[idx]
                    l_m[mm,nn] = sum_l_m_k_nn
            
            # construct 
            sum_l_m = np.zeros(n)   # Sum along axis 0
            l_m_weighted = np.zeros((m,n))
            w = np.zeros((m, n))
            min_l_m = max_along_axis_0(l_m)
            for i in range(n):
                for mm in range(m):
                    l_m_weighted[mm, i] = alpha_jn[mm] * np.exp(l_m[mm,i] - min_l_m[i])
                    sum_l_m[i] += l_m_weighted[mm, i]
                        
            # update weight
            # Compute w = l_j / sum_l_j
            for i in range(n):
                for mm in range(m):
                    w[mm, i] = l_m_weighted[mm, i] /  max(sum_l_m[i], 0.01)
                    w_mk[mm*k: (mm+1)*k, i*(t-1): (i+1)*(t-1)] = w_mk[mm*k: (mm+1)*k, i*(t-1):(i+1)*(t-1)] * w[mm, i]
                    w_mk_0[mm*k: (mm+1)*k, i] = w_mk_0[mm*k: (mm+1)*k, i] * w[mm, i]
                    
            for i in range(n):
                ll += np.log(sum_l_m[i]) + min_l_m[i]
            
            penloglik = ll + np.log(2.0)
            
            # for j in range(m):
            #     s0j = sigma_0[j] / sigma_jn[j]
            #     penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
                
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            emit += 1
            
            if abs(diff) < tol:
                break
            
            mu_mk_weighted = np.zeros((m, nt))
            mu_mk_0_weighted = np.zeros((m, n))
                
            # 
            # Update parameters
            for mm in range(m):
                alpha_jn[mm] = np.mean(w[mm, :])
                # First estimate beta, rho (mm-specific)
                wtilde = w[mm, :].T
                # estimate param
                
                w_mk_sum = np.zeros(nt)
                w_mk_0_sum = np.zeros(n)
                
                # First compute mu_mk_weighted
                
                for kk in range(k):
                    idx_type = mm * k + kk
                    mu_mk_weighted[mm] += w_mk[idx_type] * mu_jn[idx_type]
                    mu_mk_0_weighted[mm] += w_mk_0[idx_type] * mu_0_jn[idx_type]
                    w_mk_sum += w_mk[idx_type]
                    w_mk_0_sum += w_mk_0[idx_type]
                    
                w_mk_sum[w_mk_sum < 1e-12] = 1e-12
                w_mk_0_sum[w_mk_0_sum < 1e-12] = 1e-12
                mu_mk_weighted[mm] = mu_mk_weighted[mm] / w_mk_sum
                mu_mk_0_weighted[mm] = mu_mk_0_weighted[mm] / w_mk_0_sum
                
                # Step 1: estimate beta given rho 
                
                rho_estim = mubeta_jn_mat[mm, -1]
                y_tmp = y_c - rho_estim * y_l - mu_mk_weighted[mm] * (1 - rho_estim) - z_c @ gamma_jn + z_l @ gamma_jn * rho_estim
                x_tmp = ( x_c - rho_estim * x_l) 
                
                w_m = np.zeros(nt)
                for i in range(n):
                    w_m[i * (t-1) : (i + 1) * (t-1)] = wtilde[i]
                
                x_tmp_w = np.zeros(x_tmp.shape)
                for ii in range(x_tmp.shape[1]):
                    x_tmp_w[:, ii] = w_m * x_tmp[:, ii]
                
                mubeta_estim = solve_linear_system_safe(x_tmp_w.T @ x_tmp, x_tmp_w.T @ y_tmp)
                
                # # Step 2: estimate rho given beta
                y_tmp = y_c - mu_mk_weighted[mm] - x_c @ mubeta_estim - z_c @ gamma_jn
                x_tmp = y_l - mu_mk_weighted[mm] - x_l @ mubeta_estim - z_l @ gamma_jn
                if np.abs(np.sum(w_m)) < 1e-6:
                    w_m[w_m < 1e-12] = 1e-12
                # np.sum(y_tmp * w_m) / np.sum(x_tmp * w_m)                
                rho_estim = np.sum(y_tmp * x_tmp * w_m) / np.sum(x_tmp * x_tmp * w_m)
                rho_estim = min(0.99, max(0.01, rho_estim))
                
                # end mubeta and rho update
                mubeta_jn_mat[mm,1:(q+1)] = mubeta_estim[:(q)]
                mubeta_jn_mat[mm,(q+1):(2*q+1)] = - mubeta_estim[:(q)] * rho_estim
                mubeta_jn_mat[mm,-1] = rho_estim
                
                # estimate param_0 for initial period
                
                # Estimate mu: mk-specific
                ytilde = y_c - x1[:,1:] @ mubeta_jn_mat[mm,1:] - z_c @ gamma_jn + z_l @ gamma_jn * mubeta_jn_mat[mm,-1] 
                if q > 0:
                    xtilde_0 = np.zeros((n,q))
                    for ii in range(q):
                        xtilde_0[:, ii] = wtilde * x1_0[:, (ii+1)]
                    
                    coef_0_jn = solve_linear_system_safe(xtilde_0.T @ x1_0[:, 1:], xtilde_0.T @ (y_0 - z_0 @ gamma_0_jn - mu_mk_0_weighted[mm]) )
                    
                    mubeta_0_jn_mat[mm,1:(q+1)] = coef_0_jn[:(q)]
                    
                    ytilde_0 = y_0 - z_0 @ gamma_0_jn  - x1_0[:, 1:] @ coef_0_jn
                    
                else:
                    ytilde_0 = y_0 - z_0 @ gamma_0_jn
                    
                w_mk_sum = np.zeros(nt)
                w_mk_0_sum = np.zeros(n)                
                res_mm_sq = np.zeros(nt)
                res_mm_0_sq = np.zeros(n)                
                sum_tau_jn = 0                 
                
                
                for kk in range(k):
                    idx_type = mm * k + kk
                    
                    mu_jn[idx_type] = (np.mean(ytilde * w_mk[idx_type,:]) / max(np.mean(w_mk[idx_type,:]), 1e-6) ) / (1 - mubeta_jn_mat[mm,-1])
                    
                    mu_0_jn[idx_type] = np.mean(ytilde_0 * w_mk_0[idx_type,:]) / max(np.mean(w_mk_0[idx_type,:]), 1e-6)
                
                    res_mm_sq += w_mk[idx_type,:] * (ytilde - mu_jn[idx_type] *  (1 - mubeta_jn_mat[mm,-1]) )**2
                    res_mm_0_sq += w_mk_0[idx_type,:] * (ytilde_0 - mu_0_jn[idx_type])**2
                    w_mk_sum += w_mk[idx_type,:]
                    w_mk_0_sum += w_mk_0[idx_type,:]
                    
                    tau_jn[idx_type] = min(max(np.mean(np.concatenate((w_mk[idx_type,:], w_mk_0[idx_type,:]))),  0.05),0.95)
                    
                    sum_tau_jn+=tau_jn[idx_type]
                
                for kk in range(k):
                    idx_type = mm * k + kk
                    tau_jn[idx_type] = tau_jn[idx_type] / sum_tau_jn

                sigma_jn[mm] = np.sqrt(( np.sum(res_mm_sq) ) / max(np.sum(w_mk_sum), 1e-6)  )
                sigma_jn[mm] = max(sigma_jn[mm], epsilon * sigma_0[mm])

                sigma_0_jn[mm] = np.sqrt(( np.sum(res_mm_0_sq) ) / max(np.sum(w_mk_0_sum), 1e-6)  )
                sigma_0_jn[mm] = max(sigma_0_jn[mm], epsilon * sigma_0[mm])
           
                
            mubeta_jn = mubeta_jn_mat.T.flatten()
            mubeta_0_jn = mubeta_0_jn_mat.T.flatten()
            # update alpha
            for j in range(m):
                alpha_jn[j] = max(0.05, alpha_jn[j] )
            
            total_alpha = np.sum(alpha_jn)
            for j in range(m):
                alpha_jn[j] = alpha_jn[j] / total_alpha
             
            # update gamma
            if p > 0:
                ztilde = np.zeros((nt, p), dtype=np.float64) 
                ztilde_0 = np.zeros((n, p), dtype=np.float64) 
                zz = np.zeros((p, p), dtype=np.float64) 
                ze = np.zeros((p, 1), dtype=np.float64) 
                
                zz_0 = np.zeros((p, p), dtype=np.float64) 
                ze_0 = np.zeros((p, 1), dtype=np.float64) 
                for j in range(m):
                    wtilde = w[j, :]
                    w_j = np.zeros(nt)
                    for i in range(n):
                        w_j[i * (t - 1) : (i + 1) * (t - 1)] = wtilde[i]
                    
                    for ii in range(p):
                        ztilde[:, ii] = w_j * (z_c[:, ii] - z_l[:, ii] * mubeta_jn_mat[j,-1])
                        ztilde_0[:, ii] = wtilde * z_0[:, ii]
                        
                    zz += ztilde.T @ (z_c[:, ii] - z_l[:, ii] * mubeta_jn_mat[j,-1]) / (sigma_jn[j]**2)
                    ze += ztilde.T @( y_c - x1 @ mubeta_jn_mat[j,:] - mu_mk_weighted[j]) / max(sigma_jn[j]**2, 1e-6)
                    
                    zz_0 += ztilde_0.T @ z_0 / (sigma_0_jn[j]**2)
                    ze_0 += ztilde_0.T @( y_0 - x1_0 @ mubeta_0_jn_mat[j,:] - mu_mk_0_weighted[j]) / max(sigma_0_jn[j]**2, 1e-6)
                gamma_jn = solve_linear_system_safe(zz,ze).flatten()
                gamma_0_jn = solve_linear_system_safe(zz_0,ze_0).flatten()
                
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mubeta_draw[:, jn] = mubeta_jn_mat[:,:q1].T.flatten()
        sigma_draw[:, jn] = sigma_jn
        mu_draw[:,jn] = mu_jn
        mu_0_draw[:,jn] = mu_0_jn
        tau_draw[:,jn] = tau_jn
        mubeta_0_draw[:, jn] = mubeta_0_jn_mat[:,:(q+1)].T.flatten()
        sigma_0_draw[:, jn] = sigma_0_jn
        
        if p > 0:
            gamma_draw[:, jn] = gamma_jn
            gamma_0_draw[:, jn] = gamma_0_jn
    
    return(alpha_draw,tau_draw, mu_draw, mubeta_draw,sigma_draw,gamma_draw,mu_0_draw, mubeta_0_draw,sigma_0_draw,gamma_0_draw, penloglikset, loglikset ,post)

# %%  
@njit
def regpanelmixAR1mixturePMLE(y, x, z, p, q, m, k, ninits=10, tol_long=1e-6, maxit=2000, tol_short=1e-2, maxit_short=200)  : 
    t,n = y.shape
    nt = n * (t-1)
    y_l = y[:-1,:]
    y_0 = y[0,:]
    y_c = y[1:,:]
    
    y_c = y_c.T.flatten()
    y_l = y_l.T.flatten()
    
    # y.reshape((n,t)).T - data_lr[0][0] # check equivalence
    # Handle x
    x_0 = np.zeros((n,q))
    x_l = np.zeros((nt,q))
    x_c = np.zeros((nt,q))
    for j in range(q):
        x_j_mat = np.ascontiguousarray(x[:,j]).reshape((n,t)).T
        x_0[:,j] = x_j_mat[0,:]
        x_l[:,j] = x_j_mat[:-1,:].T.flatten()
        x_c[:,j] = x_j_mat[1:,:].T.flatten()
        
    # Handle x
    z_0 = np.zeros((n,p))
    z_l = np.zeros((nt,p))
    z_c = np.zeros((nt,p))
    for j in range(p):
        z_j_mat = np.ascontiguousarray(z[:,j]).reshape((n,t)).T
        z_0[:,j] = z_j_mat[0,:]
        z_l[:,j] = z_j_mat[:-1,:].T.flatten()
        z_c[:,j] = z_j_mat[1:,:].T.flatten()
        
    x1 = np.hstack((np.ones((nt, 1)), x_c, x_l, y_l[:,np.newaxis]))
    xz = np.hstack((x1, z_c, z_l))
    q1 = 2*q + 2
    
    xz_0 = np.hstack((np.ones((n, 1)), x_0, z_0) )
    
    out_coef = solve_least_squares(xz, y_c)  # Replace np.linalg.lstsq
    out_coef_0 = solve_least_squares(xz_0, y_0)  # Replace np.linalg.lstsq
    residuals = y_c - xz @ out_coef
    residuals_0 = y_0 - xz_0 @ out_coef_0
    stdR = np.std(residuals)
    stdR_0 = np.std(residuals_0)
    npar = m - 1 + (q + 2 * k - 1) * m + p 
    npar_0 =  (q + k) * m + p 
    # ninits_short = ninits * 10 * (q1 + p) * m
    ninits_short = ninits * 10 
    
    # -----
    estim_init = regpanelmixAR1PMLE(y, x, z, p, q, m)
    alpha_hat  = estim_init['alpha_hat'][0]
    mubeta_hat = estim_init['mubeta_hat'][0]
    mubeta_hat = np.ascontiguousarray(mubeta_hat)
    mubeta_hat_mat = mubeta_hat.reshape((q+1,m)).T
    mu_hat = estim_init['mu_hat'][0]
    rho_hat = estim_init['rho_hat'][0]
    sigma_hat  = estim_init['sigma_hat'][0]
    gamma_hat  = estim_init['gamma_hat'][0]
    mubeta_0_hat= estim_init['mubeta_0_hat'][0]
    mubeta_0_hat = np.ascontiguousarray(mubeta_0_hat)
    mubeta_0_hat_mat = mubeta_0_hat.reshape((q+1,m)).T
    mu_0_hat = mubeta_0_hat_mat[:,0]
    sigma_0_hat = estim_init['sigma_0_hat'][0]
    gamma_0_hat = estim_init['gamma_0_hat'][0]
    # -----
    # First draw random start point
    if p > 0:
        # Perform least squares regression with both x and z
        gamma_0_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) 
        gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short)) 
        for pp in range(p):
            gamma_draw[:,pp] = gamma_draw[:,pp] * gamma_hat[pp]
            gamma_0_draw[:,pp] = gamma_0_draw[:,pp] * gamma_0_hat[pp]
       
        
    else:
        # Perform least squares regression with x only
        gamma_draw = np.zeros((0,ninits_short), dtype=np.float64)
        gamma_0_draw = np.zeros((0,ninits_short), dtype=np.float64)
    
    # Initialize alpha
    alpha_draw = (alpha_hat +  generate_random_uniform(0, 0.1, (m, ninits_short)).T).T
    alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))

    # Initialize rho and mu
    rho_draw = generate_random_uniform(0.1, 0.9, (m, ninits_short))
    
    # Draw mu and tau
    mk = m * k 
    mu_0_draw = np.zeros((m * k, ninits_short))
    mu_draw = np.zeros((m * k, ninits_short))
    tau_draw = np.zeros((m * k, ninits_short))
    
    for mm in range(m):
        tau_draw_m = generate_random_uniform(0, 1, (k, ninits_short))
        tau_draw_m = (tau_draw_m / np.sum(tau_draw_m, axis=0))
        tau_draw[(mm * k):((mm +1) * k),:] = tau_draw_m
        
        quantile_cutoff = np.linspace(0, 1, k+1)
        for kk in range(k):
            idx_type = mm * k + kk 
            mu_draw[( idx_type), :] = mu_hat[mm] + np.random.uniform(quantile_cutoff[kk], quantile_cutoff[kk+1], size=ninits_short) * np.abs(mu_hat[mm])
            mu_0_draw[( idx_type), :] = mu_hat[mm] + np.random.uniform(quantile_cutoff[kk], quantile_cutoff[kk+1], size=ninits_short) * np.abs(mu_hat[mm])
    
    # Initialize mubeta
    mubeta_draw = np.zeros((q1 * m, ninits_short))
    mubeta_0_draw = np.zeros(( (q + 1) * m, ninits_short))
    
    for mm in range(m):
        # Leave mu blank
        mubeta_draw[mm, :] = 0 #leave mu blank for mubeta draw
        mubeta_0_draw[mm, :] = 0
                
        #assign rho
        mubeta_draw[m*(q1 - 1)+mm, :] = rho_draw[mm,:]
        
        if q > 0:
            for j in range(q):
                mubeta_draw[ m*(j+1)+mm , :] = mubeta_hat_mat[mm,j+1] * np.random.uniform(0.9, 1.1, size=ninits_short)
                mubeta_draw[ m*(j+q+1)+mm , :] = mubeta_draw[ m*(j+1)+mm , :] * mubeta_draw[m*(q1 - 1)+mm, :]  
                mubeta_0_draw[ m*(j+1)+mm , :] = mubeta_0_hat_mat[mm,j+1] * np.random.uniform(0.9, 1.1, size=ninits_short)
                
    an = 1 / n    
    sigma_0 = np.full(m, stdR)

    # Initialize sigma
    sigma_draw = (generate_random_uniform(0.1,0.2, (m, ninits_short)).T * sigma_hat).T
    sigma_0_draw = (generate_random_uniform(0.1, 0.2, (m, ninits_short)).T * sigma_0_hat).T        
    
    alpha_draw,tau_draw, mu_draw, mubeta_draw,sigma_draw,gamma_draw,mu_0_draw, mubeta_0_draw,sigma_0_draw,gamma_0_draw, penloglikset, loglikset, post = EM_optimization_AR1_mixture(y_c,  y_l,  x_c, x_l, z_c, z_l, xz, y_0, xz_0, p, q, sigma_0, alpha_draw, tau_draw, mu_draw, mubeta_draw,sigma_draw,gamma_draw,mu_0_draw, mubeta_0_draw,sigma_0_draw,gamma_0_draw, m, k, n, t, an, maxit=maxit_short, tol=tol_short)

    components = np.argsort(penloglikset)[::-1][:ninits]
    alpha_draw = alpha_draw[:,components]
    tau_draw = tau_draw[:,components]
    mu_draw = mu_draw[:,components]
    mubeta_draw = mubeta_draw[:,components]
    sigma_draw = sigma_draw[:,components]
    gamma_draw = gamma_draw[:,components]
    
    mu_0_draw = mu_0_draw[:,components]
    mubeta_0_draw = mubeta_0_draw[:,components]
    sigma_0_draw = sigma_0_draw[:,components]
    gamma_0_draw = gamma_0_draw[:,components]
    
    alpha_draw,tau_draw, mu_draw, mubeta_draw, sigma_draw,gamma_draw,mu_0_draw, mubeta_0_draw,sigma_0_draw,gamma_0_draw, penloglikset, loglikset, post = EM_optimization_AR1_mixture(y_c,  y_l,  x_c, x_l, z_c, z_l, xz, y_0, xz_0, p, q, sigma_0, alpha_draw, tau_draw, mu_draw, mubeta_draw,sigma_draw,gamma_draw,mu_0_draw, mubeta_0_draw,sigma_0_draw,gamma_0_draw, m, k, n, t, an, maxit=maxit, tol=tol_long)

    
    index = np.argmax(penloglikset)
    alpha_hat = alpha_draw[:,index]
    tau_hat = tau_draw[:,index]
    mu_hat = mu_draw[:,index]
    mubeta_mat =  np.ascontiguousarray(mubeta_draw[:,index]).reshape((q1,m)).T
    rho_hat = mubeta_mat[:,-1]
    
    beta_hat = mubeta_mat[:,1:(q+1)]
    sigma_hat = sigma_draw[:,index]
    gamma_hat = gamma_draw[:,index]        
    
    mu_0_hat = mu_0_draw[:,index]
    mubeta_0_hat = mubeta_0_draw[:,index]
    mu_0_hat = mu_0_draw[:,index]
    mubeta_0_mat =  np.ascontiguousarray(mubeta_0_draw[:,index]).reshape((q+1,m)).T
    # beta_0_hat = mubeta_0_mat[:, 1:(q+1)].T.flatten()
    beta_0_hat = mubeta_0_mat[:, 1:(q+1)]
    
    sigma_0_hat = sigma_0_draw[:,index]
    gamma_0_hat = gamma_0_draw[:,index]
    
    post = post[:, index]
    penloglik = penloglikset[index]
    loglik = loglikset[index]
    aic = -2 * loglik + 2 * (npar + npar_0)
    bic = -2 * loglik + np.log(n) * (npar + npar_0)
    
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
    result_dict['rho_hat']  = rho_hat[np.newaxis,:]
    result_dict['sigma_hat']  = sigma_hat[np.newaxis,:]
    result_dict['mu_hat'] = mu_hat[np.newaxis,:]
    result_dict['beta_hat'] = beta_hat
    result_dict['gamma_hat'] = gamma_hat[np.newaxis,:]
    result_dict['post'] = np.ascontiguousarray(post).reshape((n,m))
    
    result_dict['mu_0_hat'] = mu_0_hat[np.newaxis,:]
    result_dict['beta_0_hat'] = beta_0_hat
    result_dict['sigma_0_hat']  = sigma_0_hat[np.newaxis,:]
    result_dict['gamma_0_hat'] = gamma_0_hat[np.newaxis,:]
        
    return result_dict 
# %%
@njit
def EM_optimization_mixture(y_c, x, z, p, q, sigma_0, alpha_draw, tau_draw, mu_draw, mubeta_draw, sigma_draw, gamma_draw, m, k, t, an, maxit=1000, tol= 1e-8, epsilon = 1e-6):
    
    nt = len(y_c)
    n = nt // t
    mk = int(m*k)
    ninits = alpha_draw.shape[1]
    x1 = np.hstack((np.ones((nt, 1)), x))
    q1 = q + 1
    
    post = np.zeros((m * n, ninits))
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
              
    for jn in prange(ninits):            
        alpha_jn = alpha_draw[:, jn]
        tau_jn   = tau_draw[:, jn]
        mu_jn = mu_draw[:, jn]
        mubeta_jn = np.ascontiguousarray(mubeta_draw[:, jn])
        sigma_jn = sigma_draw[:, jn]
        gamma_jn = gamma_draw[:, jn]  # Likely float64
        mubeta_jn_mat = mubeta_jn.reshape((q+1,m)).T
        
        oldpenloglik = -np.inf
        diff = 1.0
        
        for iter_ii in range(maxit):
            
            ll = 0
            
            if p > 0:
                ytilde = y_c - np.dot(z, gamma_jn)
            else:
                ytilde = y_c
            
            r = np.zeros((mk, nt))
            for mm in range(m):
                for kk in range(k):     
                    res = ytilde - x1 @ mubeta_jn_mat[mm] - mu_jn[mm*k + kk]
                    r[mm*k + kk] = log_likelihood_array(res, 0.0, sigma_jn[mm])
            
            l_m = np.zeros((m,n))  # Same shape as `r`
            w_mk = np.zeros((mk,nt))    # Weights
                
            for mm in range(m):
                r_m = r[ mm*k: (mm+1)*k, :]
                tau_m = tau_jn[ mm*k: (mm+1)*k] 
                minr = max_along_axis_0(r_m)
                 
                w_m = np.zeros((k,nt))    # Weights
                l_m_k = np.zeros((k,nt))  # Same shape as `r`
                
                                            
                for i in range(nt):
                    for kk in range(k):
                        l_m_k[kk,i] = tau_m[kk] * np.exp( r_m[kk,i] - minr[i] )
                        
                sum_l_m_k = np.zeros(nt)   # Sum along axis 0
                for i in range(nt):
                    for kk in range(k):
                        sum_l_m_k[i] += l_m_k[kk, i]
                
                
                for i in range(nt):
                    for kk in range(k):
                        w_m[kk, i] = l_m_k[kk, i] / max(sum_l_m_k[i], 1e-6)
                                
                w_mk[mm*k: (mm+1)*k, :] = w_m
                
                # compute l_m 
                for nn in range(n):
                    sum_l_m_k_nn = 0.0
                    for tt in range(t):
                        idx = nn * (t) + tt  # Compute the flattened index
                        sum_l_m_k_nn += np.log(sum_l_m_k[idx]) + minr[idx]
                    l_m[mm,nn] = sum_l_m_k_nn
            
            # construct 
            sum_l_m = np.zeros(n)   # Sum along axis 0
            l_m_weighted = np.zeros((m,n))
            w = np.zeros((m, n))
            min_l_m = max_along_axis_0(l_m)
            for i in range(n):
                for mm in range(m):
                    l_m_weighted[mm, i] = alpha_jn[mm] * np.exp(l_m[mm,i] - min_l_m[i])
                    sum_l_m[i] += l_m_weighted[mm, i]
                        
            # update weight
            # Compute w = l_j / sum_l_j
            for i in range(n):
                for mm in range(m):
                    w[mm, i] = l_m_weighted[mm, i] /  max(sum_l_m[i], 1e-6)
                    w_mk[mm*k: (mm+1)*k, i*(t): (i+1)*(t)] = w_mk[mm*k: (mm+1)*k, i*(t):(i+1)*(t)] * w[mm, i]
                    
            for i in range(n):
                ll += np.log(sum_l_m[i]) + min_l_m[i]
            
            penloglik = ll + np.log(2.0) 
            
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik
            
            if abs(diff) < tol:
                break
            
            # Update parameters
            mu_mk_weighted = np.zeros((m, nt))
            for mm in range(m):
                alpha_jn[mm] = np.mean(w[mm, :])
                
                # First estimate beta, rho (mm-specific)
                wtilde = w[mm, :].T
                # estimate param
                w_m = np.zeros(nt)
                for i in range(n):
                    w_m[i * (t) : (i + 1) * (t)] = wtilde[i]
                xtilde = np.zeros((nt,q1-1)) # remove the constant term
                for ii in range(1, x1.shape[1]):
                    xtilde[:, (ii-1)] = w_m * x1[:, ii]
                
                w_mk_sum = np.zeros(nt)
                               
                mu_mk_weighted[mm] = np.zeros(nt)
                for kk in range(k):
                    idx_type = mm * k + kk
                    mu_mk_weighted[mm] += w_mk[idx_type] * mu_jn[idx_type]
                    w_mk_sum += w_mk[idx_type]
                w_mk_sum[w_mk_sum < 1e-6] = 1e-6
                mu_mk_weighted[mm] = mu_mk_weighted[mm] / w_mk_sum
                
                coef_jn = solve_linear_system_safe(xtilde.T @ x1[:,1:], xtilde.T @ (y_c - mu_mk_weighted[mm] - z @ gamma_jn ) )
                
                
                # Estimate mu: mk-specific
                ytilde = y_c - x1[:,1:] @ coef_jn - z @ gamma_jn
                mubeta_jn_mat[mm,1:q1] = coef_jn[:q1-1]
                    
                res_mm_sq = np.zeros(nt)
                
                sum_tau_jn = 0 
                
                for kk in range(k):
                    idx_type = mm * k + kk
                    mu_jn[idx_type] = np.mean(ytilde * w_mk[idx_type,:]) / max(np.mean(w_mk[idx_type,:]), 1e-6)
                
                    res_mm_sq += w_mk[idx_type,:] * (ytilde - mu_jn[idx_type])**2
                    
                    tau_jn[idx_type] = min(max(np.mean(w_mk[idx_type,:]),  0.05),0.95)
                    sum_tau_jn+=tau_jn[idx_type]
                
                for kk in range(k):
                    idx_type = mm * k + kk
                    tau_jn[idx_type] = tau_jn[idx_type] / max(sum_tau_jn, 1e-6)

                sigma_jn[mm] = np.sqrt( np.sum(res_mm_sq) / np.sum(w_mk_sum)  )
                sigma_jn[mm] = max(sigma_jn[mm], epsilon * sigma_0[mm])

            mubeta_jn = mubeta_jn_mat.T.flatten()
            # update alpha
            for j in range(m):
                alpha_jn[j] = max(0.05, alpha_jn[j] )
            
            total_alpha = np.sum(alpha_jn)
            for j in range(m):
                alpha_jn[j] = alpha_jn[j] / total_alpha
            
            # update gamma
            if p > 0:
                ztilde = np.zeros((nt, p), dtype=np.float64) 
                zz = np.zeros((p, p), dtype=np.float64) 
                ze = np.zeros((p, 1), dtype=np.float64) 
                for j in range(m):
                    wtilde = w[j, :]
                    w_j = np.zeros(nt)
                    for i in range(n):
                        w_j[i * (t) : (i + 1) * (t )] = wtilde[i]
                    for ii in range(p):
                        ztilde[:, ii] = w_j * z[:, ii]
                         
                    zz += ztilde.T @ (z) / (sigma_jn[j]**2)
                    ze += ztilde.T @( y_c - x1 @ mubeta_jn_mat[j,:] - mu_mk_weighted[j]) / max(sigma_jn[j]**2,1e-6)
                gamma_jn = solve_linear_system_safe(zz,ze).flatten()
                
        # print(iter_ii)
        penloglikset[jn] = penloglik
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        alpha_draw[:, jn] = alpha_jn
        mu_draw[:, jn] = mu_jn
        mubeta_draw[:, jn] = mubeta_jn_mat.T.flatten()
        sigma_draw[:, jn] = sigma_jn
        tau_draw[:, jn] = tau_jn
        gamma_draw[:, jn] = gamma_jn
    
    return(alpha_draw,tau_draw,mu_draw, mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset ,post)



@njit
def regpanelmixmixturePMLE(y, x, z, p, q, m, k, ninits=2, tol=1e-6, maxit=2000, tol_short=1e-2, maxit_short=50,alpha_bound=0.05):    # Extract the generated data
    t,n = y.shape
    nt = n * t
    y_c = y.T.flatten()
    
    x1 = np.hstack((np.ones((nt, 1)), x))
    xz = np.hstack((x1, z))
    q1 = q + 1
    
    if p > 0:
        xz = np.hstack((x1, z))
    else:
        xz = x1 
        
    out_coef = solve_least_squares(xz, y_c)  # Replace np.linalg.lstsq
    residuals = y_c - xz @ out_coef
    stdR = np.std(residuals)
    npar = m - 1 + (q + 2 * k - 1) * m + p 
    # ninits_short = ninits * 10 * (q1 + p) * m
    ninits_short = ninits * 10 
    # -----
    estim_init = regpanelmixPMLE(y, x, z, p, q, m, alpha_bound=alpha_bound)
    alpha_hat  = estim_init['alpha_hat'][0]
    mubeta_hat = estim_init['mubeta_hat'][0]
    mubeta_hat = np.ascontiguousarray(mubeta_hat)
    mubeta_hat_mat = mubeta_hat.reshape((q+1,m)).T
    mu_hat = mubeta_hat_mat[:,0]
    sigma_hat  = estim_init['sigma_hat'][0]
    gamma_hat  = estim_init['gamma_hat'][0]
    # -----

    if p > 0:
        gamma_draw = generate_random_uniform(0.5, 1.5, (p, ninits_short))
        for pp in range(p):
            gamma_draw[:,pp] = gamma_draw[:,pp] * gamma_hat[pp]
       
    else:
        # Perform least squares regression with x only
        gamma_draw = np.zeros((0,ninits_short), dtype=np.float64)
        
    # Initialize alpha
    alpha_draw = (alpha_hat +  generate_random_uniform(0, 0.1, (m, ninits_short)).T).T
    alpha_draw = (alpha_draw / np.sum(alpha_draw, axis=0))

    mk = m * k 
    mu_draw = np.zeros((m * k, ninits_short))
    tau_draw = np.zeros((m * k, ninits_short))
    for mm in range(m):
        tau_draw_m = generate_random_uniform(0, 1, (k, ninits_short))
        tau_draw_m = (tau_draw_m / np.sum(tau_draw_m, axis=0))
        tau_draw[(mm * k):((mm +1) * k),:] = tau_draw_m
        
        quantile_cutoff = np.linspace(0, 1, k+1)
        for kk in range(k):
            idx_type = mm * k + kk 
            mu_draw[( idx_type), :] = mu_hat[mm] + np.random.uniform(quantile_cutoff[kk], quantile_cutoff[kk+1], size=ninits_short) * np.abs(mu_hat[mm])
    
    mubeta_draw = np.zeros((q1 * m, ninits_short))
    
    for mm in range(m):
        mubeta_draw[mm, :] = 0 #leave mu blank for mubeta draw    
        if q > 0:
            for j in range(q):
                mubeta_draw[ m*(j+1)+mm , :] = mubeta_hat_mat[mm,j+1] * np.random.uniform(0.9, 1.1, size=ninits_short)
                
    sigma_draw = (generate_random_uniform(0.9, 1.1, (m, ninits_short)).T * sigma_hat).T
    
    an = 1 / n    
    sigma_0 = np.full(m, stdR)
    
    alpha_draw,tau_draw,mu_draw, mubeta_draw,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization_mixture(y_c, x, z, p, q, sigma_0, alpha_draw, tau_draw, mu_draw, mubeta_draw, sigma_draw, gamma_draw, m, k, t, an, maxit=maxit_short, tol=tol_short)

    components = np.argsort(penloglikset)[::-1][:ninits]
    alpha_draw = alpha_draw[:,components]
    mu_draw = mu_draw[:,components]
    mubeta_draw = mubeta_draw[:,components]
    sigma_draw = sigma_draw[:,components]
    gamma_draw = gamma_draw[:,components]
    tau_draw = tau_draw[:,components]
    
    alpha_draw,tau_draw, mu_draw, mubeta_draw ,sigma_draw,gamma_draw,penloglikset, loglikset, post = EM_optimization_mixture(y_c, x, z, p, q, sigma_0, alpha_draw, tau_draw, mu_draw, mubeta_draw, sigma_draw, gamma_draw, m, k, t, an,  maxit=maxit, tol=tol)
    
    index = np.argmax(penloglikset)
    alpha_hat = alpha_draw[:,index]
    tau_hat = tau_draw[:,index]
    mu_hat = mu_draw[:,index]
    
    mubeta_mat =  np.ascontiguousarray(mubeta_draw[:,index]).reshape((q1,m)).T
    beta_hat = mubeta_mat[:,1:(q+1)]
    
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
    result_dict['post'] = np.ascontiguousarray(post).reshape((n,m))
    result_dict['alpha_hat']  = alpha_hat[np.newaxis,:]
    result_dict['tau_hat']  = tau_hat[np.newaxis,:]
    result_dict['sigma_hat']  = sigma_hat[np.newaxis,:]
    result_dict['mu_hat'] = mu_hat[np.newaxis,:]
    result_dict['beta_hat'] = beta_hat    
    result_dict['gamma_hat'] = gamma_hat[np.newaxis,:]
    return result_dict

# %%

@njit(parallel=True)
def LRTestParallel(Data, N, T, M, p, q, nrep, BB = 199, alpha_bound=0.05):
    result_lr_each = np.zeros((nrep,6))
    for ii in prange(nrep): 
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        out_h0 = regpanelmixPMLE(y,x,z, p, q, M, alpha_bound=alpha_bound)
        out_h1 = regpanelmixPMLE(y,x,z, p, q, M+1, alpha_bound=alpha_bound)
        alpha_hat  = out_h0['alpha_hat'][0]
        mubeta_hat = out_h0['mubeta_hat'][0]
        sigma_hat  = out_h0['sigma_hat'][0]
        gamma_hat  = out_h0['gamma_hat'][0]
        penloglik_h0 = out_h0['penloglik'][0,0]
        penloglik_h1 = out_h1['penloglik'][0,0]
        # print(out_h0[0])        
        lr_stat = -2 * (penloglik_h0 - penloglik_h1)
        aic_crit = out_h1['aic'][0,0] < out_h0['aic'][0,0]
        bic_crit = out_h1['bic'][0,0] < out_h0['bic'][0,0]
        mubeta_hat = np.ascontiguousarray(mubeta_hat)
        mubeta_hat_mat = mubeta_hat.reshape((q+1,M)).T
        beta_hat = mubeta_hat_mat[:,1:]
        mu_hat =  mubeta_hat_mat[:,0]

        # Generate data for bootstrap
        Data_bb = [generate_data(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T, M, p, q) for _ in prange(BB)]
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data = Data_bb[bb]
            y_bb = data[0]
            x_bb = data[1]
            z_bb = data[2]
            
            # Call regpanelmixPMLE for m components
            out_h0 = regpanelmixPMLE(y_bb, x_bb, z_bb, p, q, M, alpha_bound=alpha_bound)
            out_h1 = regpanelmixPMLE(y_bb, x_bb, z_bb, p, q, M + 1, alpha_bound=alpha_bound)
            penloglik_h0 = out_h0['penloglik'][0, 0]
            penloglik_h1 = out_h1['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0 - penloglik_h1)
        
        result_lr_each[ii,0] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.90))
        result_lr_each[ii,1] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.95))
        result_lr_each[ii,2] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.99))
        result_lr_each[ii,3] = aic_crit
        result_lr_each[ii,4] = bic_crit
        result_lr_each[ii,5] = lr_stat
        
    return(result_lr_each)
# %%

@njit(parallel=True)
def LRTestAR1Parallel(Data, N, T, M, p, q, nrep, BB = 199):
    result_lr_each = np.zeros((nrep,6))
    for ii in prange(nrep): 
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        out_h0 = regpanelmixAR1PMLE(y, x, z, p, q, M)
        out_h1 = regpanelmixAR1PMLE(y, x, z, p, q, M+1)
        alpha_hat  = out_h0['alpha_hat'][0]
        mubeta_hat = out_h0['mubeta_hat'][0]
        rho_hat = out_h0['rho_hat'][0]
        mu_hat =  out_h0['mu_hat'][0]

        sigma_hat  = out_h0['sigma_hat'][0]
        gamma_hat  = out_h0['gamma_hat'][0]
        
        mubeta_0_hat = out_h0['mubeta_0_hat'][0]
        
        sigma_0_hat  = out_h0['sigma_0_hat'][0]
        gamma_0_hat  = out_h0['gamma_0_hat'][0]
        
        penloglik_h0 = out_h0['penloglik'][0,0]
        penloglik_h1 = out_h1['penloglik'][0,0]
        # print(out_h0[0])        
        lr_stat = -2 * (penloglik_h0 - penloglik_h1)
        aic_crit = out_h1['aic'][0,0] < out_h0['aic'][0,0]
        bic_crit = out_h1['bic'][0,0] < out_h0['bic'][0,0]
        
        mubeta_hat = np.ascontiguousarray(mubeta_hat)
        mubeta_hat_mat = mubeta_hat.reshape((q+1,M)).T
        beta_hat = mubeta_hat_mat[:,1:]
        

        mubeta_0_hat = np.ascontiguousarray(mubeta_0_hat)
        mubeta_0_hat_mat = mubeta_0_hat.reshape((q+1,M)).T
        beta_0_hat = mubeta_0_hat_mat[:,1:]
        mu_0_hat =  mubeta_0_hat_mat[:,0]
        
        # mubeta_hat_mat = np.empty((M, q+1)) 
        # mubeta_hat_mat[:, 0] = mu_hat
        # mubeta_hat_mat[:, 1:] = beta_hat
        # mubeta_hat_mat.T.flatten()

        # Generate data for bootstrap        
        Data_bb = [generate_data_ar1(alpha_hat, rho_hat, mu_hat, sigma_hat, beta_hat, gamma_hat, mu_0_hat, sigma_0_hat, beta_0_hat, gamma_0_hat, N, T, M, p, q) for _ in prange(BB)]
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data_bb = Data_bb[bb]
            y_bb = data_bb[0]
            x_bb = data_bb[1]
            z_bb = data_bb[2]
            
            # Call regpanelmixPMLE for m components
            out_h0_bb = regpanelmixAR1PMLE(y_bb, x_bb, z_bb, p, q, M)
            out_h1_bb = regpanelmixAR1PMLE(y_bb, x_bb, z_bb, p, q, M + 1)
            penloglik_h0_bb = out_h0_bb['penloglik'][0, 0]
            penloglik_h1_bb = out_h1_bb['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0_bb - penloglik_h1_bb)
        
        result_lr_each[ii,0] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.90))
        result_lr_each[ii,1] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.95))
        result_lr_each[ii,2] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.99))
        result_lr_each[ii,3] = aic_crit
        result_lr_each[ii,4] = bic_crit
        result_lr_each[ii,5] = lr_stat
    return(result_lr_each)


# %%
@njit(parallel=True)
def LRTestAR1MixtureParallel(Data, N, T, M, K, p, q, nrep, BB = 199):
    result_lr_each = np.zeros((nrep,6))
    for ii in prange(nrep): 
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        out_h0 = regpanelmixAR1mixturePMLE(y, x, z, p, q, M, K)
        out_h1 = regpanelmixAR1mixturePMLE(y, x, z, p, q, M+1, K)
        alpha_hat  = out_h0['alpha_hat'][0]
        tau_hat  = np.ascontiguousarray(out_h0['tau_hat'][0]).reshape(M,K)
        
        mu_hat = np.ascontiguousarray(out_h0['mu_hat'][0]).reshape(M,K)
        beta_hat = out_h0['beta_hat']
        rho_hat = out_h0['rho_hat'][0]
        sigma_hat  = out_h0['sigma_hat'][0]
        gamma_hat  = out_h0['gamma_hat'][0]
        
        mu_0_hat = np.ascontiguousarray(out_h0['mu_0_hat']).reshape(M,K) 
        beta_0_hat = out_h0['beta_0_hat']

        sigma_0_hat  = out_h0['sigma_0_hat'][0]
        gamma_0_hat  = out_h0['gamma_0_hat'][0]
        
        
        penloglik_h0 = out_h0['penloglik'][0,0]
        penloglik_h1 = out_h1['penloglik'][0,0]
        # print(out_h0[0])        
        lr_stat = -2 * (penloglik_h0 - penloglik_h1)
        aic_crit = out_h1['aic'][0,0] < out_h0['aic'][0,0]
        bic_crit = out_h1['bic'][0,0] < out_h0['bic'][0,0]

        # Generate data for bootstrap        
        Data_bb = [generate_data_ar1_mixture(alpha_hat, tau_hat, rho_hat, mu_hat, sigma_hat, beta_hat, gamma_hat, mu_0_hat, sigma_0_hat, beta_0_hat, gamma_0_hat, N, T, M, K, p, q) for _ in prange(BB)]
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data_bb = Data_bb[bb]
            y_bb = data_bb[0]
            x_bb = data_bb[1]
            z_bb = data_bb[2]
            
            # Call regpanelmixPMLE for m components
            out_h0_bb = regpanelmixAR1mixturePMLE(y_bb, x_bb, z_bb, p, q, M, K) 
            out_h1_bb = regpanelmixAR1mixturePMLE(y_bb, x_bb, z_bb, p, q, M+1, K) 
            penloglik_h0_bb = out_h0_bb['penloglik'][0, 0]
            penloglik_h1_bb = out_h1_bb['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0_bb - penloglik_h1_bb)
        result_lr_each[ii,0] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.90))
        result_lr_each[ii,1] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.95))
        result_lr_each[ii,2] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.99))
        result_lr_each[ii,3] = aic_crit
        result_lr_each[ii,4] = bic_crit
        result_lr_each[ii,5] = lr_stat
    return(result_lr_each)

# %%
@njit(parallel=True)
def LRTestMixtureParallel(Data, N, T, M, K, p, q, nrep, BB = 199):
    result_lr_each = np.zeros((nrep,6))
    for ii in prange(nrep): 
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        out_h0 = regpanelmixmixturePMLE(y,x,z, p, q, M, K)
        out_h1 = regpanelmixmixturePMLE(y,x,z, p, q, M+1, K)
        alpha_hat  = out_h0['alpha_hat'][0]
        tau_hat  = np.ascontiguousarray(out_h0['tau_hat'][0]).reshape(M,K)
        mubeta_hat = out_h0['mubeta_hat'][0]
        sigma_hat  = out_h0['sigma_hat'][0]
        gamma_hat  = out_h0['gamma_hat'][0]
        
        penloglik_h0 = out_h0['penloglik'][0,0]
        penloglik_h1 = out_h1['penloglik'][0,0]
        # print(out_h0[0])        
        lr_stat = -2 * (penloglik_h0 - penloglik_h1)
        aic_crit = out_h1['aic'][0,0] < out_h0['aic'][0,0]
        bic_crit = out_h1['bic'][0,0] < out_h0['bic'][0,0]
        mubeta_hat = np.ascontiguousarray(mubeta_hat)
        mubeta_hat_mat = mubeta_hat.reshape((q+1,M*K)).T
        beta_hat = mubeta_hat_mat[:,1:]
        mu_hat =  np.ascontiguousarray(mubeta_hat_mat[:,0]).reshape(M,K)



        # Generate data for bootstrap        
        Data_bb = [generate_data_mixture(alpha_hat, mu_hat, beta_hat, sigma_hat, tau_hat, gamma_hat, N, T, M, K, p, q) for _ in prange(BB)]
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data_bb = Data_bb[bb]
            y_bb = data_bb[0]
            x_bb = data_bb[1]
            z_bb = data_bb[2]
            
            # Call regpanelmixPMLE for m components
            out_h0_bb = regpanelmixmixturePMLE(y_bb,x_bb,z_bb, p, q, M, K)
            out_h1_bb = regpanelmixmixturePMLE(y_bb,x_bb,z_bb, p, q, M+1, K)
            penloglik_h0_bb = out_h0_bb['penloglik'][0, 0]
            penloglik_h1_bb = out_h1_bb['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0_bb - penloglik_h1_bb)
        result_lr_each[ii,0] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.90))
        result_lr_each[ii,1] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.95))
        result_lr_each[ii,2] = 1 * ( lr_stat >  compute_quantile(lr_stat_bb, 0.99))
        result_lr_each[ii,3] = aic_crit
        result_lr_each[ii,4] = bic_crit
        result_lr_each[ii,5] = lr_stat
        
    return(result_lr_each)

# %%
# For sequential test
# %%

@njit
def NonParTestNoCovariates(y, N, T, n_grid, n_bins, BB, r_test):
    weights_equal = np.full(N, 1 / N)
    
    data_P_W = calculate_P_matrix(y, weights_equal, n_grid=n_grid, n_bins=n_bins)
        
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
    
    for i in prange(BB):
        # Calculate bootstrapped P and Sigma_P matrices
        data_P_W_b = calculate_P_matrix(y, ru[i, :], n_grid=n_grid, n_bins=n_bins)
        
        for k in range(T):
            P_k = data_P_W_b['P_k_list'][k]
            Sigma_P_k = data_P_W_b['Sigma_P_k_list'][k]
            # Compute KP statistics for the k-th triplet
            rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, lambda_c_list[k])['rk_c'][0,0]
            # rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, 0)['rk_c'][0,0]
    # Compute max and mean values for rk and rk_b
    rk_b_max = max_along_axis_1(rk_b)  # Maximum of rk_b along axis 1
    rk_b_max_95 = compute_quantile(rk_b_max, 0.95)  # 95th quantile of rk_b_max

    
    # Store results

    rk_mean = np.mean(rk)  # Mean of rk (Numba supports this)
    rk_b_mean = mean_along_axis_1(rk_b)  # Mean of rk_b along axis 1
    rk_b_mean_95 = compute_quantile(rk_b_mean, 0.95)  # 95th quantile of 
    return np.array([(rk.max() > rk_b_max_95), 1 * (rk_mean > rk_b_mean_95)])


# %%
@njit(parallel=False) 
def NonParTest(y, N, T, n_grid, n_bins, BB, r_test):
    weights_equal = np.full(N, 1 / N)
    
    data_P_W = calculate_P_matrix(y, weights_equal, n_grid=n_grid, n_bins=n_bins)
        
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
    
    for i in prange(BB):
        # Calculate bootstrapped P and Sigma_P matrices
        data_P_W_b = calculate_P_matrix(y, ru[i, :], n_grid=n_grid, n_bins=n_bins)
        
        for k in range(T):
            P_k = data_P_W_b['P_k_list'][k]
            Sigma_P_k = data_P_W_b['Sigma_P_k_list'][k]
            # Compute KP statistics for the k-th triplet
            rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, lambda_c_list[k])['rk_c'][0,0]
            # rk_b[i, k] = construct_stat_KP(P_k, Sigma_P_k, r_test, N, 0)['rk_c'][0,0]
    # Compute max and mean values for rk and rk_b
    rk_b_max = max_along_axis_1(rk_b)  # Maximum of rk_b along axis 1
    lr_95 = compute_quantile(rk_b_max, 0.95)  # 95th quantile of rk_b_max
    lr_90 = compute_quantile(rk_b_max, 0.90)  # 95th quantile of rk_b_max
    lr_99 = compute_quantile(rk_b_max, 0.99)  # 95th quantile of rk_b_max
    lr_stat = rk.max()      
    aic = 0
    bic = 0                   
    
    # Store results

    rk_mean = np.mean(rk)  # Mean of rk (Numba supports this)
    rk_b_mean = mean_along_axis_1(rk_b)  # Mean of rk_b along axis 1
    rk_b_mean_95 = compute_quantile(rk_b_mean, 0.95)  # 95th quantile of 
    # return np.array([(rk.max() > rk_b_max_95), 1 * (rk_mean > rk_b_mean_95)])
    return np.array([lr_stat, lr_90, lr_95, lr_99, aic, bic])
   
# %%

@njit(parallel=False) 
def LRTestNormal(y, x, z, p, q, m, N, T, bootstrap = True, BB= 199, spline=False):
    if (spline) & (q > 0):
        bs_degree = 2                # Replace with the degree of the B-spline
        # Compute quantiles (equivalent to probs = 0.33 and probs = 0.67)
        x_spline = np.zeros((N*T, 0))
        for qq in range(q):
            quantiles = np.quantile(x[:,qq], [0.33, 0.67])
            knots = np.concatenate((
                np.array([quantiles[0]] * (bs_degree)),  # This is fine
                quantiles,
                np.array([quantiles[1]] * (bs_degree))  # This is fine
            ))
            # Generate B-spline basis columns
            bs_columns = generate_b_spline_basis(x[:,qq], knots, bs_degree)
            x_spline = np.concatenate((x_spline, bs_columns),axis=1)
        q_spline = x_spline.shape[1]
        out_h0 = regpanelmixPMLE(y,x_spline,z, p, q_spline, m)
        out_h1 = regpanelmixPMLE(y,x_spline,z, p, q_spline, m+1)
    else:
        out_h0 = regpanelmixPMLE(y,x,z, p, q, m)
        out_h1 = regpanelmixPMLE(y,x,z, p, q, m+1)
        # print(out_h0)
        # print(out_h1)
        # -2 * (out_h0['penloglik'] - out_h1['penloglik'])
    alpha_hat  = out_h0['alpha_hat'][0]
    mubeta_hat = out_h0['mubeta_hat'][0]
    sigma_hat  = out_h0['sigma_hat'][0]
    gamma_hat  = out_h0['gamma_hat'][0]
    penloglik_h0 = out_h0['penloglik'][0,0]
    aic = out_h0['aic'][0,0]
    bic = out_h0['bic'][0,0]
    penloglik_h1 = out_h1['penloglik'][0,0]
    # print(out_h0[0])        
    lr_stat = -2 * (penloglik_h0 - penloglik_h1)
    
    if (spline) & (q > 0):
        mubeta_hat = np.ascontiguousarray(mubeta_hat)
        mubeta_hat_mat = mubeta_hat.reshape((q_spline+1,m)).T
        beta_hat = mubeta_hat_mat[:,1:]
        mu_hat =  mubeta_hat_mat[:,0]
    else:
        mubeta_hat = np.ascontiguousarray(mubeta_hat)
        mubeta_hat_mat = mubeta_hat.reshape((q+1,m)).T
        beta_hat = mubeta_hat_mat[:,1:]
        mu_hat =  mubeta_hat_mat[:,0]
    
    if bootstrap:
        # Generate data for bootstrap
        Data = [generate_data(alpha_hat, mu_hat, sigma_hat, gamma_hat, beta_hat, N, T, m, p, q, spline=spline, z_input=z, x_input=x) for _ in prange(BB)]
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data = Data[bb]
            y_bb = data[0]
            x_bb = data[1]
            z_bb = data[2]
            
            if (spline) & (q > 0):
                # Compute quantiles (equivalent to probs = 0.33 and probs = 0.67)
                x_spline_bb = np.zeros((N*T, 0))
                for qq in range(q):
                    quantiles = np.quantile(x_bb[:,qq], [0.33, 0.67])
                    knots = np.concatenate((
                        np.array([quantiles[0]] * (bs_degree)),  # This is fine
                        quantiles,
                        np.array([quantiles[1]] * (bs_degree))  # This is fine
                    ))
                    # Generate B-spline basis columns
                    bs_columns = generate_b_spline_basis(x_bb[:,qq], knots, bs_degree)
                    x_spline_bb = np.concatenate((x_spline_bb, bs_columns),axis=1)
                
                out_h0 = regpanelmixPMLE(y_bb,x_spline_bb,z_bb, p, q_spline, m)
                out_h1 = regpanelmixPMLE(y_bb,x_spline_bb,z_bb, p, q_spline, m+1)
            else:
                # Call regpanelmixPMLE for m components
                out_h0 = regpanelmixPMLE(y_bb, x_bb, z_bb, p, q, m)
                out_h1 = regpanelmixPMLE(y_bb, x_bb, z_bb, p, q, m + 1)
            penloglik_h0 = out_h0['penloglik'][0, 0]
            penloglik_h1 = out_h1['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0 - penloglik_h1)
        
        lr_99 = compute_quantile(lr_stat_bb, 0.99)
        lr_95 = compute_quantile(lr_stat_bb, 0.95)
        lr_90 = compute_quantile(lr_stat_bb, 0.90)
    else:
        lr_99 = np.inf
        lr_95 = np.inf
        lr_90 = np.inf
    return np.array([lr_stat, lr_90, lr_95, lr_99, aic, bic])
   
# %%
@njit(parallel=False)
def LRTestMixture(y, x, z, p, q, m, k, N, T, bootstrap = True, BB= 199, spline=False):
    if (spline) & (q > 0):    
        bs_degree = 2                # Replace with the degree of the B-spline
        # Compute quantiles (equivalent to probs = 0.33 and probs = 0.67)
        x_spline = np.zeros((N*T, 0))
        for qq in range(q):
            quantiles = np.quantile(x[:,qq], [0.33, 0.67])
            knots = np.concatenate((
                np.array([quantiles[0]] * (bs_degree)),  # This is fine
                quantiles,
                np.array([quantiles[1]] * (bs_degree))  # This is fine
            ))
            # Generate B-spline basis columns
            bs_columns = generate_b_spline_basis(x[:,qq], knots, bs_degree)
            x_spline = np.concatenate((x_spline, bs_columns),axis=1)
        q_spline = x_spline.shape[1]
        out_h0 = regpanelmixmixturePMLE(y,x_spline,z, p, q_spline, m,k)
        out_h1 = regpanelmixmixturePMLE(y,x_spline,z, p, q_spline, m+1,k)
    else:
        out_h0 = regpanelmixmixturePMLE(y,x,z, p, q, m, k)
        out_h1 = regpanelmixmixturePMLE(y,x,z, p, q, m+1, k)
    alpha_hat  = out_h0['alpha_hat'][0]
    tau_hat  = np.ascontiguousarray(out_h0['tau_hat'][0]).reshape(m,k)
    mu_hat = np.ascontiguousarray(out_h0['mu_hat'][0]).reshape(m, k)
    beta_hat = out_h0['beta_hat']
    sigma_hat  = out_h0['sigma_hat'][0]
    gamma_hat  = out_h0['gamma_hat'][0]
    penloglik_h0 = out_h0['penloglik'][0,0]
    aic = out_h0['aic'][0,0]
    bic = out_h0['bic'][0,0]
    penloglik_h1 = out_h1['penloglik'][0,0]
    # print(out_h0[0])        
    lr_stat = -2 * (penloglik_h0 - penloglik_h1)
    
    if bootstrap:
        # Generate data for bootstrap
        Data = [generate_data_mixture(alpha_hat, mu_hat, beta_hat, sigma_hat, tau_hat, gamma_hat, N, T, m, k, p, q, spline=spline, z_input=z, x_input=x) for _ in prange(BB)]
        
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data = Data[bb]
            y_bb = data[0]
            x_bb = data[1]
            z_bb = data[2]
            
            if (spline) & (q > 0):
                # Compute quantiles (equivalent to probs = 0.33 and probs = 0.67)
                x_spline_bb = np.zeros((N*T, 0))
                for qq in range(q):
                    quantiles = np.quantile(x_bb[:,qq], [0.33, 0.67])
                    knots = np.concatenate((
                        np.array([quantiles[0]] * (bs_degree)),  # This is fine
                        quantiles,
                        np.array([quantiles[1]] * (bs_degree))  # This is fine
                    ))
                    # Generate B-spline basis columns
                    bs_columns = generate_b_spline_basis(x_bb[:,qq], knots, bs_degree)
                    x_spline_bb = np.concatenate((x_spline_bb, bs_columns),axis=1)
                out_h0 = regpanelmixmixturePMLE(y_bb,x_spline_bb,z_bb, p, q_spline, m, k)
                out_h1 = regpanelmixmixturePMLE(y_bb,x_spline_bb,z_bb, p, q_spline, m+1, k)
            else:               # Call regpanelmixPMLE for m components
                out_h0 = regpanelmixmixturePMLE(y_bb, x_bb, z_bb, p, q, m, k)
                out_h1 = regpanelmixmixturePMLE(y_bb, x_bb, z_bb, p, q, m+1, k)
            penloglik_h0 = out_h0['penloglik'][0, 0]
            penloglik_h1 = out_h1['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0 - penloglik_h1)
        
        lr_99 = compute_quantile(lr_stat_bb, 0.99)
        lr_95 = compute_quantile(lr_stat_bb, 0.95)
        lr_90 = compute_quantile(lr_stat_bb, 0.90)
    else:
        lr_99 = np.inf
        lr_95 = np.inf
        lr_90 = np.inf
    return np.array([lr_stat, lr_90, lr_95, lr_99, aic, bic])

# %%
@njit(parallel=False) 
def LRTestAR1Mixture(y, x, z, p, q, m, k, N, T, bootstrap = True, BB= 199, spline=False):
    if (spline) & (q > 0):    
        bs_degree = 2                # Replace with the degree of the B-spline
        # Compute quantiles (equivalent to probs = 0.33 and probs = 0.67)
        x_spline = np.zeros((N*T, 0))
        for qq in range(q):
            quantiles = np.quantile(x[:,qq], [0.33, 0.67])
            knots = np.concatenate((
                np.array([quantiles[0]] * (bs_degree)),  # This is fine
                quantiles,
                np.array([quantiles[1]] * (bs_degree))  # This is fine
            ))
            # Generate B-spline basis columns
            bs_columns = generate_b_spline_basis(x[:,qq], knots, bs_degree)
            x_spline = np.concatenate((x_spline, bs_columns),axis=1)
        q_spline = x_spline.shape[1]
        out_h0 = regpanelmixAR1mixturePMLE(y,x_spline,z, p, q_spline, m,k)
        out_h1 = regpanelmixAR1mixturePMLE(y,x_spline,z, p, q_spline, m+1,k)
    else:
        out_h0 = regpanelmixAR1mixturePMLE(y, x, z, p, q, m, k)
        out_h1 = regpanelmixAR1mixturePMLE(y, x, z, p, q, m+1, k)
    
    alpha_hat  = out_h0['alpha_hat'][0]
    tau_hat  = np.ascontiguousarray(out_h0['tau_hat'][0]).reshape(m, k)
    
    mu_hat = np.ascontiguousarray(out_h0['mu_hat'][0]).reshape(m, k)
    beta_hat = out_h0['beta_hat']
    rho_hat = out_h0['rho_hat'][0]
    sigma_hat  = out_h0['sigma_hat'][0]
    gamma_hat  = out_h0['gamma_hat'][0]
    
    mu_0_hat = np.ascontiguousarray(out_h0['mu_0_hat']).reshape(m, k) 
    beta_0_hat = out_h0['beta_0_hat']

    sigma_0_hat  = out_h0['sigma_0_hat'][0]
    gamma_0_hat  = out_h0['gamma_0_hat'][0]
    
    penloglik_h0 = out_h0['penloglik'][0,0]
    aic = out_h0['aic'][0,0]
    bic = out_h0['bic'][0,0]
    penloglik_h1 = out_h1['penloglik'][0,0]
    # print(out_h0[0])        
    lr_stat = -2 * (penloglik_h0 - penloglik_h1)
    
    
    if bootstrap:
        # Generate data for bootstrap
        Data = [generate_data_ar1_mixture(alpha_hat, tau_hat, rho_hat, mu_hat, sigma_hat, beta_hat, gamma_hat, mu_0_hat, sigma_0_hat, beta_0_hat, gamma_0_hat, N, T, m, k, p, q, spline=spline, z_input=z, x_input=x) for _ in prange(BB)]
        
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data = Data[bb]
            y_bb = data[0]
            x_bb = data[1]
            z_bb = data[2]
            
            if (spline) & (q > 0):
                # Compute quantiles (equivalent to probs = 0.33 and probs = 0.67)
                x_spline_bb = np.zeros((N*T, 0))
                for qq in range(q):
                    quantiles = np.quantile(x_bb[:,qq], [0.33, 0.67])
                    knots = np.concatenate((
                        np.array([quantiles[0]] * (bs_degree)),  # This is fine
                        quantiles,
                        np.array([quantiles[1]] * (bs_degree))  # This is fine
                    ))
                    # Generate B-spline basis columns
                    bs_columns = generate_b_spline_basis(x_bb[:,qq], knots, bs_degree)
                    x_spline_bb = np.concatenate((x_spline_bb, bs_columns),axis=1)
                out_h0 = regpanelmixAR1mixturePMLE(y_bb,x_spline_bb,z_bb, p, q_spline, m, k)
                out_h1 = regpanelmixAR1mixturePMLE(y_bb,x_spline_bb,z_bb, p, q_spline, m+1, k)
            else:
                # Call regpanelmixPMLE for m components
                out_h0 =regpanelmixAR1mixturePMLE(y_bb, x_bb, z_bb, p, q, m, k) 
                out_h1 =regpanelmixAR1mixturePMLE(y_bb, x_bb, z_bb, p, q, m+1, k) 
            penloglik_h0 = out_h0['penloglik'][0, 0]
            penloglik_h1 = out_h1['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0 - penloglik_h1)
        
        lr_99 = compute_quantile(lr_stat_bb, 0.99)
        lr_95 = compute_quantile(lr_stat_bb, 0.95)
        lr_90 = compute_quantile(lr_stat_bb, 0.90)
    else:
        lr_99 = np.inf
        lr_95 = np.inf
        lr_90 = np.inf
    return np.array([lr_stat, lr_90, lr_95, lr_99, aic, bic])

# %%
@njit(parallel=False)
def LRTestAR1Normal(y, x, z, p, q, m, N, T, bootstrap = True, BB= 199, spline=False):
    out_h0 = regpanelmixAR1PMLE(y,x,z, p, q, m)
    out_h1 = regpanelmixAR1PMLE(y,x,z, p, q, m+1)
    alpha_hat  = out_h0['alpha_hat'][0]
    mubeta_hat = out_h0['mubeta_hat'][0]
    mu_hat = out_h0['mu_hat'][0]
    rho_hat = out_h0['rho_hat'][0]
    sigma_hat  = out_h0['sigma_hat'][0]
    gamma_hat  = out_h0['gamma_hat'][0]
    
    mubeta_0_hat = out_h0['mubeta_0_hat'][0]
    
    sigma_0_hat  = out_h0['sigma_0_hat'][0]
    gamma_0_hat  = out_h0['gamma_0_hat'][0]
    
    
    penloglik_h0 = out_h0['penloglik'][0,0]
    penloglik_h1 = out_h1['penloglik'][0,0]
    
    penloglik_h0 = out_h0['penloglik'][0,0]
    aic = out_h0['aic'][0,0]
    bic = out_h0['bic'][0,0]
    penloglik_h1 = out_h1['penloglik'][0,0]
    # print(out_h0[0])        
    lr_stat = -2 * (penloglik_h0 - penloglik_h1)
    
    mubeta_hat = np.ascontiguousarray(mubeta_hat)
    mubeta_hat_mat = mubeta_hat.reshape((q+1,m)).T
    beta_hat = mubeta_hat_mat[:,1:]

    mubeta_0_hat = np.ascontiguousarray(mubeta_0_hat)
    mubeta_0_hat_mat = mubeta_0_hat.reshape((q+1,m)).T
    beta_0_hat = mubeta_0_hat_mat[:,1:]
    mu_0_hat =  mubeta_0_hat_mat[:,0]
    
    if bootstrap:
        # Generate data for bootstrap
        Data = [generate_data_ar1(alpha_hat, rho_hat, mu_hat, sigma_hat, beta_hat, gamma_hat, mu_0_hat, sigma_0_hat, beta_0_hat, gamma_0_hat, N, T, m, p, q, z_input=z,x_input=x) for _ in prange(BB)]
        
        
        # Preallocate lr_stat as a 1D array (Numba-compatible)
        lr_stat_bb = np.zeros(BB, dtype=np.float64)
        
        for bb in prange(BB):
            data = Data[bb]
            y_bb = data[0]
            x_bb = data[1]
            z_bb = data[2]
            
            # Call regpanelmixPMLE for m components
            
            out_h0 = regpanelmixAR1PMLE(y_bb, x_bb, z_bb, p, q, m)
            out_h1 = regpanelmixAR1PMLE(y_bb, x_bb, z_bb, p, q, m+1)
            penloglik_h0 = out_h0['penloglik'][0, 0]
            penloglik_h1 = out_h1['penloglik'][0, 0]
            
            # Compute likelihood ratio statistic
            lr_stat_bb[bb] = -2 * (penloglik_h0 - penloglik_h1)
        
        lr_99 = compute_quantile(lr_stat_bb, 0.99)
        lr_95 = compute_quantile(lr_stat_bb, 0.95)
        lr_90 = compute_quantile(lr_stat_bb, 0.90)
    else:
        lr_99 = np.inf
        lr_95 = np.inf
        lr_90 = np.inf
    return np.array([lr_stat, lr_90, lr_95, lr_99, aic, bic])
          
# %%

# %%

@njit(parallel=True)  # Numba JIT compilation with parallelization
def parallel_processing_empirical_test_stationary(nrep, M_max, BB, Data, N, T, M, K, p, q):
    # Initialize result tables
    aic_table = np.zeros((nrep, M_max))
    bic_table = np.zeros((nrep, M_max))
    aic_table_mixture = np.zeros((nrep, M_max))
    bic_table_mixture = np.zeros((nrep, M_max))

    lr_estim_table = np.zeros((nrep, M_max))

    rk_mean_table = np.zeros((nrep, M_max))
    rk_max_table = np.zeros((nrep, M_max))

    lr_1_table = np.zeros((nrep, M_max))
    lr_5_table = np.zeros((nrep, M_max))
    lr_10_table = np.zeros((nrep, M_max))

    lr_1_table_mixture = np.zeros((nrep, M_max))
    lr_5_table_mixture = np.zeros((nrep, M_max))
    lr_10_table_mixture = np.zeros((nrep, M_max))

    # Generate data for all repetitions
    

    # Parallel loop
    for ii in prange(nrep):  # Use prange for parallel execution
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        bootstrap_nocov = True
        bootstrap_mixture = True

        for mm in range(M_max):
            # Non-parametric test
            m = mm + 1
            n_grid = m + 1
            r_test = m
            n_bins = math.ceil((m + 1) ** (1 / (T - 1)))
            rk_stat_each = NonParTestNoCovariates(y, N, T, n_grid, n_bins, BB, r_test)

            # LR test for no-covariates model
            lr_results_nocov = LRTestNormal(y, x, z, p, q, m, N, T, bootstrap=bootstrap_nocov, BB=BB)
            lr_stat_nocov, lr_90_nocov, lr_95_nocov, lr_99_nocov, aic_nocov, bic_nocov = lr_results_nocov

            # LR test for mixture model
            lr_results_mixture = LRTestMixture(y, x, z, p, q, m, 2, N, T, bootstrap=bootstrap_mixture, BB=BB)
            lr_stat_mixture, lr_90_mixture, lr_95_mixture, lr_99_mixture, aic_mixture, bic_mixture = lr_results_mixture

            # Record results
            rk_max_table[ii, mm] = rk_stat_each[0]
            rk_mean_table[ii, mm] = rk_stat_each[1]

            aic_table[ii, mm] = aic_nocov
            bic_table[ii, mm] = bic_nocov
            aic_table_mixture[ii, mm] = aic_mixture
            bic_table_mixture[ii, mm] = bic_mixture

            lr_estim_table[ii, mm] = lr_stat_nocov

            lr_1_table[ii, mm] = lr_stat_nocov > lr_99_nocov
            lr_5_table[ii, mm] = lr_stat_nocov > lr_95_nocov
            lr_10_table[ii, mm] = lr_stat_nocov > lr_90_nocov

            lr_1_table_mixture[ii, mm] = lr_stat_mixture > lr_99_mixture
            lr_5_table_mixture[ii, mm] = lr_stat_mixture > lr_95_mixture
            lr_10_table_mixture[ii, mm] = lr_stat_mixture > lr_90_mixture

            # Update bootstrap flags
            if lr_stat_nocov < lr_90_nocov:
                bootstrap_nocov = False
            if lr_stat_mixture < lr_90_mixture:
                bootstrap_mixture = False

    return (aic_table, bic_table, aic_table_mixture, bic_table_mixture, lr_estim_table, rk_mean_table, rk_max_table,
            lr_1_table, lr_5_table, lr_10_table,
            lr_1_table_mixture, lr_5_table_mixture, lr_10_table_mixture)



@njit(parallel=True)  # Numba JIT compilation with parallelization
def parallel_processing_empirical_test_ar1(nrep, M_max, BB, Data, N, T, M, K, p, q):
    # Initialize result tables
    aic_table = np.zeros((nrep, M_max))
    bic_table = np.zeros((nrep, M_max))
    aic_table_mixture = np.zeros((nrep, M_max))
    bic_table_mixture = np.zeros((nrep, M_max))

    lr_estim_table = np.zeros((nrep, M_max))

    rk_mean_table = np.zeros((nrep, M_max))
    rk_max_table = np.zeros((nrep, M_max))

    lr_1_table = np.zeros((nrep, M_max))
    lr_5_table = np.zeros((nrep, M_max))
    lr_10_table = np.zeros((nrep, M_max))

    lr_1_table_mixture = np.zeros((nrep, M_max))
    lr_5_table_mixture = np.zeros((nrep, M_max))
    lr_10_table_mixture = np.zeros((nrep, M_max))

    # Generate data for all repetitions
    

    # Parallel loop
    for ii in prange(nrep):  # Use prange for parallel execution
        data = Data[ii]
        y = data[0]
        x = data[1]
        z = data[2]
        bootstrap_nocov = True
        bootstrap_mixture = True

        for mm in range(M_max):
            # Non-parametric test
            m = mm + 1
            n_grid = m + 1
            r_test = m
            n_bins = math.ceil((m + 1) ** (1 / (T - 1)))
            rk_stat_each = NonParTestNoCovariates(y, N, T, n_grid, n_bins, BB, r_test)

            # LR test for no-covariates model
            lr_results_nocov = LRTestAR1Normal(y, x, z, p, q, m, N, T, bootstrap=bootstrap_nocov, BB=BB)
            lr_stat_nocov, lr_90_nocov, lr_95_nocov, lr_99_nocov, aic_nocov, bic_nocov = lr_results_nocov

            # LR test for mixture model
            lr_results_mixture = LRTestAR1Mixture(y, x, z, p, q, m, 2, N, T, bootstrap=bootstrap_mixture, BB=BB)
            lr_stat_mixture, lr_90_mixture, lr_95_mixture, lr_99_mixture, aic_mixture, bic_mixture = lr_results_mixture

            # Record results
            rk_max_table[ii, mm] = rk_stat_each[0]
            rk_mean_table[ii, mm] = rk_stat_each[1]

            aic_table[ii, mm] = aic_nocov
            bic_table[ii, mm] = bic_nocov
            aic_table_mixture[ii, mm] = aic_mixture
            bic_table_mixture[ii, mm] = bic_mixture

            lr_estim_table[ii, mm] = lr_stat_nocov

            lr_1_table[ii, mm] = lr_stat_nocov > lr_99_nocov
            lr_5_table[ii, mm] = lr_stat_nocov > lr_95_nocov
            lr_10_table[ii, mm] = lr_stat_nocov > lr_90_nocov

            lr_1_table_mixture[ii, mm] = lr_stat_mixture > lr_99_mixture
            lr_5_table_mixture[ii, mm] = lr_stat_mixture > lr_95_mixture
            lr_10_table_mixture[ii, mm] = lr_stat_mixture > lr_90_mixture

            # Update bootstrap flags
            if lr_stat_nocov < lr_90_nocov:
                bootstrap_nocov = False
            if lr_stat_mixture < lr_90_mixture:
                bootstrap_mixture = False

    return (aic_table, bic_table, aic_table_mixture, bic_table_mixture, lr_estim_table, rk_mean_table, rk_max_table,
            lr_1_table, lr_5_table, lr_10_table,
            lr_1_table_mixture, lr_5_table_mixture, lr_10_table_mixture)

# Function to find the model where AIC stops decreasing
@njit
def find_model_stop(aic_values):
    differences = np.diff(aic_values)
    model_stop = np.where(differences > 0)[0]
    if model_stop.size == 0:
        return len(aic_values)
    return model_stop[0] + 1  # Adjust for 1-based indexing

# Function to count frequencies
@njit
def find_first_zero(sequence):
    for index, value in enumerate(sequence):
        if value == 0:
            return index + 1  # Return the first index where the value is 0
    return -1  # Return -1 if no 0 is found


def count_frequency(arr, M_max):
    """
    Counts the frequency of values in the range 1 to M_max in a NumPy array.

    Parameters:
        arr (numpy.ndarray): Input array of integers.
        M_max (int): The maximum value to consider for frequency counting.

    Returns:
        numpy.ndarray: An array of frequencies for values from 1 to M_max.
    """
    # Convert the array to integers (if not already)
    arr = arr.astype(int)
    
    # Use np.bincount with minlength to ensure range 1 to M_max is covered
    counts = np.bincount(arr, minlength=M_max + 1)[1:M_max + 1]
    return counts
# %%
# for computing variance
# ----------------------------------------------------------

# %%
def numerical_derivative_score_i(y_it, x_it, z_it, params_array, m, k, q, p, score_function, epsilon=1e-6):
    """
    Compute the numerical derivative of the score function with respect to parameters.

    Parameters:
    -----------
    y_it : ndarray
        Array of dependent variable values for all t for a single observation i.
    x_it : ndarray
        Array of covariates in x for all t for a single observation i.
    z_it : ndarray
        Array of covariates in z for all t for a single observation i.
    params_array : ndarray
        Flattened array of parameters.
    m : int
        Number of components.
    k : int
        Number of sub-components per component.
    q : int
        Number of covariates in x.
    p : int
        Number of covariates in z.
    score_function : callable
        Function to compute the score for a single observation.
    epsilon : float, optional
        Small perturbation for numerical differentiation (default is 1e-6).

    Returns:
    --------
    ndarray
        Jacobian matrix of the score function with respect to parameters.
    """
    n_params = len(params_array)
    score_base = score_function(y_it, x_it, z_it, params_array, m, k, q, p)
    jacobian = np.zeros((len(score_base), n_params))

    for i in range(n_params):
        params_perturbed = params_array.copy()
        params_perturbed[i] += epsilon
        score_perturbed = score_function(y_it, x_it, z_it, params_perturbed, m, k, q, p)
        jacobian[:, i] = (score_perturbed - score_base) / epsilon

    return jacobian


def get_params_stationary_normal(out_h0):
    """
    Convert out_h0 dictionary into arrays compatible with flatten_params.

    Parameters:
    -----------
    out_h0 : dict
        Dictionary containing parameter arrays.

    Returns:
    --------
    tuple
        - params_dict: Dictionary with flattened parameter arrays.
        - params_array: Flattened array of parameters.
    """
    alpha_hat = out_h0['alpha_hat'][0]
    sigma_hat = out_h0['sigma_hat'][0]
    mubeta_hat = out_h0['mubeta_hat'][0]
    gamma_hat = out_h0['gamma_hat'][0]
    m = len(alpha_hat)
    q = int((len(mubeta_hat) / m) - 1)

    mubeta_hat = np.ascontiguousarray(mubeta_hat)
    mubeta_hat_mat = mubeta_hat.reshape((q + 1, m)).T

    params_dict = {
        'alpha': alpha_hat,
        'sigma': sigma_hat,
        'mubeta': mubeta_hat.T.flatten(),
        'gamma': gamma_hat
    }
    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['mubeta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma']
    ])
    return params_dict, params_array

def get_params_dict_from_array_stationary_normal(p_array, m, k, q, p):
    """
    Convert a flattened parameter array back into a params_dict for the stationary normal model.

    Parameters:
    -----------
    p_array : ndarray
        Flattened array of parameters.
    m : int
        Number of components.
    q : int
        Number of covariates in x.
    p : int
        Number of covariates in z.

    Returns:
    --------
    dict
        Dictionary containing the parameters.
    """
    alpha_f = p_array[:m - 1]
    mubeta = p_array[m - 1:m - 1 + m * (q + 1)]
    sigma = p_array[m - 1 + m * (q + 1):m - 1 + m * (q + 1) + m]
    gamma = p_array[m - 1 + m * (q + 1) + m:]

    # Reconstruct alpha with the last component as 1 - sum of others
    alpha = np.append(alpha_f, 1 - np.sum(alpha_f))

    return {
        'alpha': alpha,
        'mubeta': mubeta,
        'sigma': sigma,
        'gamma': gamma
}


def log_likelihood_stationary_normal(params, data):
    """
    Compute the log-likelihood for the stationary normal model.

    Parameters:
    -----------
    params : dict
        Dictionary containing model parameters.
    data : tuple
        Tuple containing (y, x, z) data arrays.

    Returns:
    --------
    tuple
        - log_likelihood: Log-likelihood matrix for each component and unit.
        - w: Weights for each component and unit.
        - ll: Total log-likelihood value.
    """
    y, x, z = data
    q = x.shape[1]
    p = z.shape[1]

    # Unpack parameters
    alpha = params['alpha']
    sigma = params['sigma']
    mubeta = params['mubeta']
    gamma = params['gamma']
    m = len(alpha)

    mubeta = np.ascontiguousarray(mubeta)
    mubeta_mat = mubeta.reshape((q + 1, m)).T

    # Get dimensions
    t, n = y.shape
    nt = n * t

    # Flatten y and adjust for z if needed
    y_c = y.T.flatten()  # (n*t,)

    # Prepare x1 with intercept column
    x1 = np.column_stack((np.ones(nt), x)) if q > 0 else np.ones((nt, 1))

    # Compute residuals and log-likelihoods
    r = np.zeros((m, nt))
    for mm in range(m):
        res = y_c - x1 @ mubeta_mat[mm] - z @ gamma
        r[mm] = log_likelihood_array(res, 0.0, sigma[mm])

    # Sum over time periods for each unit
    r_sum = np.zeros((m, n))
    for mm in range(m):
        for nn in range(n):
            r_sum[mm, nn] = np.sum(r[mm, nn * t: (nn + 1) * t])

    # Normalization for numerical stability
    minr = np.max(r_sum, axis=0)  # (n,) max per unit

    # Compute exponentiated terms
    log_likelihood = np.zeros((m, n))
    sum_log_likelihood = np.zeros(n)
    for nn in range(n):
        for mm in range(m):
            log_likelihood[mm, nn] = alpha[mm] * np.exp(r_sum[mm, nn] - minr[nn])
        sum_log_likelihood[nn] = np.sum(log_likelihood[:, nn])

    # Compute weights and final log-likelihood
    w = np.zeros((m, n))
    ll = 0.0
    for nn in range(n):
        w[:, nn] = log_likelihood[:, nn] / sum_log_likelihood[nn]
        ll += np.log(sum_log_likelihood[nn]) + minr[nn]

    return log_likelihood, w, ll


def score_i_stationary_normal(y_it, x_it, z_it, params_array, m, k, q, p):
    """
    Compute the score for a single observation (i) in the stationary normal model.

    Parameters:
    -----------
    y_it : ndarray
        Array of dependent variable values for all t for a single observation i.
    x_it : ndarray
        Array of covariates in x for all t for a single observation i.
    z_it : ndarray
        Array of covariates in z for all t for a single observation i.
    params_array : ndarray
        Flattened array of parameters.

    Returns:
    --------
    ndarray
        Score vector for the single observation.
    """
    t = len(y_it)
    params_dict = get_params_dict_from_array_stationary_normal(params_array, m, k, q, p)
    x1_it = np.column_stack((np.ones(t), x_it)) if q > 0 else np.ones((t, 1))
    alpha = params_dict['alpha']
    sigma = params_dict['sigma']
    mubeta = params_dict['mubeta']
    gamma = params_dict['gamma']

    mubeta_mat = mubeta.reshape((q + 1, m)).T

    residuals_matrix = np.zeros((m, t))
    w = np.zeros(m)

    for j in range(m):
        residuals = y_it - z_it @ gamma - np.dot(x1_it, mubeta_mat[j]) 
        residuals_matrix[j] = residuals

        log_likelihood = (
            -0.5 * t * np.log(2 * np.pi * sigma[j] ** 2)
            - 0.5 * np.sum(residuals ** 2) / sigma[j] ** 2
        )
        w[j] = alpha[j] * np.exp(log_likelihood)

    w /= np.sum(w)

    score_alpha = np.zeros(m - 1)
    score_mubeta = np.zeros(m * (q + 1))
    score_sigma = np.zeros(m)
    score_gamma = np.zeros(p)

    for j in range(m - 1):
        score_alpha[j] = w[j] / alpha[j] - w[-1] / alpha[-1]

    for j in range(m):
        residuals = residuals_matrix[j]

        for qq in range(q + 1):
            idx = j + qq * m
            score_mubeta[idx] += w[j] * np.sum(residuals * x1_it[:, qq]) / sigma[j] ** 2

        term1 = -t / sigma[j]
        term2 = np.sum(residuals ** 2) / sigma[j] ** 3
        score_sigma[j] += w[j] * (term1 + term2)

        for pp in range(p):
            score_gamma[pp] += w[j] * np.sum(residuals * z_it[:, pp]) / sigma[j] ** 2

    score_i = np.concatenate([score_alpha, score_mubeta, score_sigma, score_gamma])
    return score_i

def score_stationary_normal(data, params_dict):
    y, x, z = data
    t, n = y.shape
    q = x.shape[1]
    p = z.shape[1]

    scores = []
    hessians = []
    alpha = params_dict['alpha']
    sigma = params_dict['sigma']
    mubeta = params_dict['mubeta']
    gamma = params_dict['gamma']
    m = len(alpha)

    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['mubeta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma']
    ])

    for i in range(n):
        y_it = y[:, i]
        z_it = z[(i * t):((i + 1) * t), :]
        x_it = x[(i * t):((i + 1) * t), :]
        score_i = score_i_stationary_normal(y_it, x_it, z_it, params_array, m, 1, q, p)
        hessian_i = numerical_derivative_score_i(y_it, x_it, z_it, params_array, m, 0, q, p, score_i_stationary_normal)
        
        scores.append(score_i)
        hessians.append(hessian_i)

    score = np.vstack(scores)
    hessian = np.array(hessians)
    return score, hessian


def loglik_obj_stationary_normal(p_array, data, M, p, q):
    """
    Objective function for log-likelihood minimization.

    Parameters:
    -----------
    p_array : ndarray
        Flattened array of parameters.
    data : tuple
        Tuple containing (y, x, z) data arrays.
    M : int
        Number of components.
    p : int
        Number of covariates in z.
    q : int
        Number of covariates in x.

    Returns:
    --------
    float
        Negative log-likelihood value for minimization.
    """
    alpha_f = p_array[:M - 1]
    mubeta = p_array[M - 1:M - 1 + M * (q + 1)]
    sigma = p_array[M - 1 + M * (q + 1):M - 1 + M * (q + 1) + M]
    gamma = p_array[M - 1 + M * (q + 1) + M:]

    # Construct parameter dictionary
    params_dict = {
        "alpha": np.append(alpha_f, np.clip(1 - np.sum(alpha_f), 0, 1)),
        "mubeta": mubeta,
        "sigma": sigma,
        "gamma": gamma
    }

    # Compute log-likelihood
    _, _, ll = log_likelihood_stationary_normal(params_dict, data)
    return -ll  # Negative for minimization

# %%
def get_params_stationary_mixture(out_h0):
    """
    Convert out_h0 dictionary into arrays compatible with flatten_params for the stationary mixture model.

    Parameters:
    -----------
    out_h0 : dict
        Dictionary containing parameter arrays.

    Returns:
    --------
    tuple
        - params_dict: Dictionary with flattened parameter arrays.
        - params_array: Flattened array of parameters.
    """
    alpha_hat = out_h0['alpha_hat'][0]
    tau_hat = out_h0['tau_hat'][0]
    sigma_hat = out_h0['sigma_hat'][0]
    mu_hat = out_h0['mu_hat'][0]
    beta_hat = out_h0['beta_hat'].T.flatten()
    gamma_hat = out_h0['gamma_hat'][0]
    m = len(alpha_hat)
    k = len(tau_hat) // m
    q = int((len(beta_hat) / m))

    # mubeta_hat = np.ascontiguousarray(mubeta_hat)
    # beta_hat_mat = mubeta_hat.reshape((q , m)).T
    
    params_dict = {
        'alpha': alpha_hat,
        'tau': tau_hat,
        'mu': mu_hat, 
        'beta': beta_hat, 
        'sigma': sigma_hat,
        'gamma': gamma_hat
    }
    
    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['tau'].reshape((m, k))[:,:-1].flatten(),
        params_dict['mu'].flatten(),
        params_dict['beta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma']
    ])
    return params_dict, params_array

def get_params_dict_from_array_stationary_mixture(params_array, m, k, q, p):
    """
    Convert a flattened parameter array back into a params_dict.

    Parameters:
    -----------
    params_array : ndarray
        Flattened array of parameters.
    m : int
        Number of components.
    k : int
        Number of sub-components per component.
    q : int
        Number of covariates in x.
    p : int
        Number of covariates in z.

    Returns:
    --------
    dict
        Dictionary containing the parameters.
    """
    alpha_f = params_array[:m - 1]
    tau_f = params_array[m - 1:m - 1 + m * (k - 1)]
    mu = params_array[m - 1 + m * (k - 1):m - 1 + m * (k - 1) + m * k]
    beta = params_array[m - 1 + m * (k - 1) + m * k:m - 1 + m * (k - 1) + m * k + m * q]
    sigma = params_array[m - 1 + m * (k - 1) + m * k + m * q:m - 1 + m * (k - 1) + m * k + m * q + m]
    gamma = params_array[m - 1 + m * (k - 1) + m * k + m * q + m:]

    # Reconstruct tau with the last component as 1 - sum of others
    tau = np.zeros(m * k)
    for mm in range(m):
        tau[mm * k:(mm + 1) * k - 1] = tau_f[mm * (k - 1):(mm + 1) * (k - 1)]
        tau[mm * k + k - 1] = 1 - np.sum(tau[mm * k:(mm + 1) * k - 1])

    # Reconstruct alpha with the last component as 1 - sum of others
    alpha = np.append(alpha_f, 1 - np.sum(alpha_f))

    return {
        'alpha': alpha,
        'tau': tau,
        'mu': mu,
        'beta': beta,
        'sigma': sigma,
        'gamma': gamma
    }
    return params_dict

def log_likelihood_stationary_mixture(params_dict, data):
    """
    Compute the log-likelihood for the stationary mixture model.

    Parameters:
    -----------
    params_dict : dict
        Dictionary containing model parameters.
    data : tuple
        Tuple containing (y, x, z) data arrays.

    Returns:
    --------
    tuple
        - log_likelihood: Log-likelihood matrix for each component and unit.
        - w: Weights for each component and unit.
        - ll: Total log-likelihood value.
    """
    y, x, z = data
    q = x.shape[1]
    p = z.shape[1]

    # Unpack parameters
    alpha = params_dict['alpha']
    tau = params_dict['tau']
    sigma = params_dict['sigma']
    mu = params_dict['mu']
    beta = params_dict['beta']
    gamma = params_dict['gamma']
    m = len(alpha)
    k = len(tau) // m
    
    beta_mat = beta.reshape((q, m)).T
    mu_mat = mu.reshape((m, k))
    
    # Get dimensions
    t, n = y.shape
    nt = n * t

    # Flatten y and adjust for z if needed
    y_c = y.T.flatten()  # (n*t,)

    
    # Compute residuals and log-likelihoods
    r = np.zeros((m * k, nt))
    for mm in range(m):
        for kk in range(k):
            res = y_c - x @ beta_mat[mm] - z @ gamma - mu_mat[mm, kk]
            r[mm * k + kk] = log_likelihood_array(res, 0.0, sigma[mm])

    # Sum over time periods for each unit
    r_sum = np.zeros((m * k, n))
    for mm in range(m * k):
        for nn in range(n):
            r_sum[mm, nn] = np.sum(r[mm, nn * t: (nn + 1) * t])

    # Normalization for numerical stability
    minr = np.max(r_sum, axis=0)  # (n,) max per unit

    # Compute exponentiated terms
    log_likelihood = np.zeros((m * k, n))
    sum_log_likelihood = np.zeros(n)
    for nn in range(n):
        for mm in range(m):
            for kk in range(k):
                idx = mm * k + kk
                log_likelihood[idx, nn] = alpha[mm] * tau[idx] * np.exp(r_sum[idx, nn] - minr[nn])
        sum_log_likelihood[nn] = np.sum(log_likelihood[:, nn])

    # Compute weights and final log-likelihood
    w = np.zeros((m * k, n))
    ll = 0.0
    for nn in range(n):
        w[:, nn] = log_likelihood[:, nn] / sum_log_likelihood[nn]
        ll += np.log(sum_log_likelihood[nn]) + minr[nn]

    return log_likelihood, w, ll


# %%
def score_i_stationary_mixture(y_it, x_it, z_it, p_array, m, k, q, p):
    """
    Compute the score for a single observation (i) in the stationary mixture model.

    Parameters:
    -----------
    y_it : ndarray
        Array of dependent variable values for all t for a single observation i.
    x_it : ndarray
        Array of covariates in x for all t for a single observation i.
    z_it : ndarray
        Array of covariates in z for all t for a single observation i.
    params_dict : dict
        Dictionary containing model parameters.

    Returns:
    --------
    ndarray
        Score vector for the single observation.
    """
    t = len(y_it)
    q = x_it.shape[1]
    p = z_it.shape[1]

    params_dict = get_params_dict_from_array_stationary_mixture(p_array, m, k, q, p)
    
    alpha = params_dict['alpha']  # Last component is 1-sum(alpha[:-1])
    tau = params_dict['tau']
    sigma = params_dict['sigma']
    mu = params_dict['mu']
    beta = params_dict['beta']
    gamma = params_dict['gamma']
    m = len(alpha)
    k = len(tau) // m

    # Compute weights using full alpha (including constrained last component)
    residuals_matrix = np.zeros((m * k, t))
    w = np.zeros(m * k)

    beta_mat = beta.reshape((q, m)).T
    mu_mat = mu.reshape((m, k))

    # Compute residuals and weights
    for j in range(m):
        for kk in range(k):
            idx = j * k + kk
            residuals = y_it - z_it @ gamma - x_it @ beta_mat[j] - mu_mat[j, kk]
            residuals_matrix[idx] = residuals
            log_likelihood = -0.5 * t * np.log(2 * np.pi * sigma[j] ** 2) - 0.5 * np.sum(residuals ** 2) / sigma[j] ** 2
            w[idx] = alpha[j] * tau[idx] * np.exp(log_likelihood)

    # Normalize weights
    w /= np.sum(w)

    # Derivatives - now with m-1 alpha components
    score_alpha = np.zeros(m - 1)
    score_tau = np.zeros(m * (k - 1))
    score_beta = np.zeros(m * q)
    score_sigma = np.zeros(m)
    score_gamma = np.zeros(p)
    score_mu = np.zeros(m * k)

    # For alpha (first m-1 components)
    for j in range(m - 1):
        score_alpha[j] = np.sum(w[j * k:(j + 1) * k]) / alpha[j] - np.sum(w[(m - 1) * k:]) / alpha[-1]

    # For tau (leave-one-out for each component in m)
    for j in range(m):
        for kk in range(k - 1):  # Exclude the last tau in each component
            idx_tau = j * (k - 1) + kk
            idx_full = j * k + kk
            score_tau[idx_tau] = w[idx_full] / tau[idx_full] - w[j * k + (k - 1)] / tau[j * k + (k - 1)]

    # For other parameters
    for j in range(m):
        for kk in range(k):
            idx = j * k + kk
            residuals = residuals_matrix[idx]

            # beta derivatives
            for qq in range(q):
                beta_id = j + m * qq
                score_beta[beta_id] += w[idx] * np.sum(residuals * x_it[:, qq]) / sigma[j] ** 2

            # sigma derivatives
            term1 = -t / sigma[j]
            term2 = np.sum(residuals ** 2) / sigma[j] ** 3
            score_sigma[j] += w[idx] * (term1 + term2)

            # gamma derivatives
            for pp in range(p):
                score_gamma[pp] += w[idx] * np.sum(residuals * z_it[:, pp]) / sigma[j] ** 2

            # mu derivatives
            score_mu[idx] = w[idx] * np.sum(residuals) / sigma[j] ** 2

    # Combine all scores
    score_i = np.concatenate([score_alpha, score_tau, score_mu, score_beta, score_sigma, score_gamma])
    return score_i


def score_stationary_mixture(data, params_dict):
    y, x, z = data
    t, n = y.shape
    q = x.shape[1]
    p = z.shape[1]

    scores = []  # Initialize a list to store scores for each observation
    hessians = []
    alpha = params_dict['alpha']  # Last component is 1-sum(alpha[:-1])
    tau = params_dict['tau']
    sigma = params_dict['sigma']
    mu = params_dict['mu']
    beta = params_dict['beta']
    gamma = params_dict['gamma']
    m = len(alpha)
    k = len(tau) // m

    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['tau'].reshape((m, k))[:,:-1].flatten(),
        params_dict['mu'].flatten(),
        params_dict['beta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma']
    ])
    
    for i in range(n):
        y_it = y[:, i]
        z_it = z[(i * t):((i + 1) * t), :]
        x_it = x[(i * t):((i + 1) * t), :]
        score_i = score_i_stationary_mixture(y_it, x_it, z_it, params_array,  m, k, q, p)
        hessian_i = numerical_derivative_score_i(y_it, x_it, z_it, params_array, m, k, q, p, score_i_stationary_mixture)
        
        
        scores.append(score_i)  # Append the score for the current observation
        hessians.append(hessian_i)
    # Ensure consistency with score_stationary_mixture by summing scores across observations
    score = np.vstack(scores)  # Concatenate all scores into a single array
    hessian = np.array(hessians)
    return score, hessian

# %%



# %%
# Variance AR1 Mixture
# ----------------------------------------------------------


def get_params_ar1_mixture(out_h0):
    """
    Convert out_h0 dictionary into arrays compatible with flatten_params for the AR(1) mixture model.

    Parameters:
    -----------
    out_h0 : dict
        Dictionary containing parameter arrays.

    Returns:
    --------
    tuple
        - params_dict: Dictionary with flattened parameter arrays.
        - params_array: Flattened array of parameters.
    """
    alpha_hat = out_h0['alpha_hat'][0]
    tau_hat = out_h0['tau_hat'][0]
    rho_hat = out_h0['rho_hat'][0]
    sigma_hat = out_h0['sigma_hat'][0]
    mu_hat = out_h0['mu_hat'][0]
    beta_hat = out_h0['beta_hat'].T.flatten()
    gamma_hat = out_h0['gamma_hat'][0]
    mu_0_hat = out_h0['mu_0_hat'][0]
    beta_0_hat = out_h0['beta_0_hat'].T.flatten()
    sigma_0_hat = out_h0['sigma_0_hat'][0]
    gamma_0_hat = out_h0['gamma_0_hat'][0]

    m = len(alpha_hat)
    k = len(tau_hat) // m
    q = int(len(beta_hat) / m)
    p = len(gamma_hat)

    params_dict = {
        'alpha': alpha_hat,
        'tau': tau_hat,
        'rho': rho_hat,
        'sigma': sigma_hat,
        'mu': mu_hat,
        'beta': beta_hat,
        'gamma': gamma_hat,
        'mu_0': mu_0_hat,
        'beta_0': beta_0_hat,
        'sigma_0': sigma_0_hat,
        'gamma_0': gamma_0_hat
    }

    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['tau'].reshape((m, k))[:, :-1].flatten(),
        params_dict['mu'].flatten(),
        params_dict['beta'].T.flatten(),
        params_dict['sigma'],
        params_dict['gamma'],
        params_dict['rho'],
        params_dict['mu_0'].flatten(),
        params_dict['beta_0'].flatten(),
        params_dict['sigma_0'],
        params_dict['gamma_0']
    ])
    return params_dict, params_array


def get_params_dict_from_array_ar1_mixture(p_array, m, k, q, p):
    """
    Convert a flattened parameter array back into a params_dict for the AR(1) mixture model.

    Parameters:
    -----------
    p_array : ndarray
        Flattened array of parameters.
    m : int
        Number of components.
    k : int
        Number of sub-components per component.
    q : int
        Number of covariates in x.
    p : int
        Number of covariates in z.

    Returns:
    --------
    dict
        Dictionary containing the parameters.
    """
    alpha_f = p_array[:m - 1]
    tau_f = p_array[m - 1:m - 1 + m * (k - 1)]
    mu = p_array[m - 1 + m * (k - 1):m - 1 + m * (k - 1) + m * k]
    beta = p_array[m - 1 + m * (k - 1) + m * k:m - 1 + m * (k - 1) + m * k + m * q]
    sigma = p_array[m - 1 + m * (k - 1) + m * k + m * q:m - 1 + m * (k - 1) + m * k + m * q + m]
    gamma = p_array[m - 1 + m * (k - 1) + m * k + m * q + m:m - 1 + m * (k - 1) + m * k + m * q + m + p]
    rho = p_array[m - 1 + m * (k - 1) + m * k + m * q + m + p:m - 1 + m * (k - 1) + m * k + m * q + m + p + m]
    mu_0 = p_array[m - 1 + m * (k - 1) + m * k + m * q + m + p + m:m - 1 + m * (k - 1) + m * k + m * q + m + p + m + m * k]
    beta_0 = p_array[m - 1 + m * (k - 1) + m * k + m * q + m + p + m + m * k:m - 1 + m * (k - 1) + m * k + m * q + m + p + m + m * k + m * q]
    sigma_0 = p_array[m - 1 + m * (k - 1) + m * k + m * q + m + p + m + m * k + m * q:m - 1 + m * (k - 1) + m * k + m * q + m + p + m + m * k + m * q + m]
    gamma_0 = p_array[m - 1 + m * (k - 1) + m * k + m * q + m + p + m + m * k + m * q + m:]

    # Reconstruct tau with the last component as 1 - sum of others
    tau = np.zeros(m * k)
    for mm in range(m):
        tau[mm * k:(mm + 1) * k - 1] = tau_f[mm * (k - 1):(mm + 1) * (k - 1)]
        tau[mm * k + k - 1] = 1 - np.sum(tau[mm * k:(mm + 1) * k - 1])

    # Reconstruct alpha with the last component as 1 - sum of others
    alpha = np.append(alpha_f, 1 - np.sum(alpha_f))

    return {
        'alpha': alpha,
        'tau': tau,
        'mu': mu,
        'beta': beta,
        'sigma': sigma,
        'gamma': gamma,
        'rho': rho,
        'mu_0': mu_0,
        'beta_0': beta_0,
        'sigma_0': sigma_0,
        'gamma_0': gamma_0
    }

def score_i_ar1_mixture(y_it, x_it, z_it, params_array, m, k, q, p):
    """
    Compute the score for a single observation (i) in the AR(1) mixture model.

    Parameters:
    -----------
    y_it : ndarray
        Array of dependent variable values for all t for a single observation i.
    x_it : ndarray
        Array of covariates in x for all t for a single observation i.
    z_it : ndarray
        Array of covariates in z for all t for a single observation i.
    params_array : ndarray
        Flattened array of parameters.

    Returns:
    --------
    ndarray
        Score vector for the single observation.
    """
    t = len(y_it)
    params_dict = get_params_dict_from_array_ar1_mixture(params_array, m, k, q, p)
    
    alpha = params_dict['alpha']
    tau = params_dict['tau']
    sigma = params_dict['sigma']
    mu = params_dict['mu']
    beta = params_dict['beta']
    gamma = params_dict['gamma']
    rho = params_dict['rho']
    mu_0 = params_dict['mu_0']
    beta_0 = params_dict['beta_0']
    sigma_0 = params_dict['sigma_0']
    gamma_0 = params_dict['gamma_0']

    beta_mat = beta.reshape((q, m)).T
    mu_mat = mu.reshape((m, k))
    beta_0_mat = beta_0.reshape((q, m)).T
    mu_0_mat = mu_0.reshape((m, k))

    residuals_matrix = np.zeros((m * k, t - 1))
    residuals_0_matrix = np.zeros((m * k, 1))
    w = np.zeros(m * k)
    for j in range(m):
        for kk in range(k):
            idx = j * k + kk
            residuals = (y_it[1:] - x_it[1:] @ beta_mat[j] - z_it[1:] @ gamma - mu_mat[j, kk]) - rho[j] * (y_it[:-1] - x_it[:-1] @ beta_mat[j] - z_it[:-1] @ gamma - mu_mat[j, kk])
            residuals_0 = y_it[0] - x_it[0] @ beta_0_mat[j] - z_it[0] @ gamma_0 - mu_0_mat[j, kk]
            residuals_matrix[idx] = residuals
            residuals_0_matrix[idx] = residuals_0

            log_likelihood = (
                -0.5 * (t - 1) * np.log(2 * np.pi * sigma[j] ** 2)
                - 0.5 * np.sum(residuals ** 2) / sigma[j] ** 2
                - 0.5 * np.log(2 * np.pi * sigma_0[j] ** 2)
                - 0.5 * residuals_0 ** 2 / sigma_0[j] ** 2
            )
            
            w[idx] = alpha[j] * tau[idx] * np.exp(log_likelihood)

    w /= np.sum(w)

    score_alpha = np.zeros(m - 1)
    score_tau = np.zeros(m * (k - 1))
    score_beta = np.zeros(m * q)
    score_sigma = np.zeros(m)
    score_gamma = np.zeros(p)
    score_mu = np.zeros(m * k)
    score_rho = np.zeros(m)
    score_beta_0 = np.zeros(m * q)
    score_sigma_0 = np.zeros(m)
    score_gamma_0 = np.zeros(p)
    score_mu_0 = np.zeros(m * k)

    for j in range(m - 1):
        score_alpha[j] = np.sum(w[j * k:(j + 1) * k]) / alpha[j] - np.sum(w[(m - 1) * k:]) / alpha[-1]

    for j in range(m):
        for kk in range(k - 1):
            idx_tau = j * (k - 1) + kk
            idx_full = j * k + kk
            score_tau[idx_tau] = w[idx_full] / tau[idx_full] - w[j * k + (k - 1)] / tau[j * k + (k - 1)]

    for j in range(m):
        for kk in range(k):
            idx = j * k + kk
            residuals = residuals_matrix[idx]
            residuals_0 = residuals_0_matrix[idx]

            for qq in range(q):
                beta_id = j + m * qq
                score_beta[beta_id] += w[idx] * np.sum(residuals * (x_it[1:, qq] - rho[j] * x_it[:-1, qq])) / sigma[j] ** 2
                score_beta_0[beta_id] += w[idx] * residuals_0 * x_it[0, qq] / sigma_0[j] ** 2

            term1 = -(t - 1) / sigma[j]
            term2 = np.sum(residuals ** 2) / sigma[j] ** 3
            score_sigma[j] += w[idx] * (term1 + term2)

            term1_0 = -1 / sigma_0[j]
            term2_0 = residuals_0 ** 2 / sigma_0[j] ** 3
            score_sigma_0[j] += w[idx] * (term1_0 + term2_0)

            for pp in range(p):
                score_gamma[pp] += w[idx] * np.sum(residuals * (z_it[1:, pp] - rho[j] * z_it[:-1, pp])) / sigma[j] ** 2
                score_gamma_0[pp] += w[idx] * residuals_0 * z_it[0, pp] / sigma_0[j] ** 2

            score_mu[idx] = w[idx] * np.sum(residuals) / sigma[j] ** 2
            score_mu_0[idx] = w[idx] * residuals_0 / sigma_0[j] ** 2

            score_rho[j] += w[idx] * np.sum(residuals * (y_it[:-1] - x_it[:-1] @ beta_mat[j] - z_it[:-1] @ gamma)) / sigma[j] ** 2
            
    score_i = np.concatenate([
        score_alpha, score_tau, score_mu, score_beta, score_sigma, score_gamma, score_rho,
        score_mu_0, score_beta_0, score_sigma_0, score_gamma_0
    ])
    return score_i
    

# %%

def score_ar1_mixture(data, params_dict):
    y, x, z = data
    t, n = y.shape
    q = x.shape[1]
    p = z.shape[1]

    scores = []
    hessians = []
    alpha = params_dict['alpha']
    tau = params_dict['tau']
    sigma = params_dict['sigma']
    mu = params_dict['mu']
    beta = params_dict['beta']
    gamma = params_dict['gamma']
    rho = params_dict['rho']
    mu_0 = params_dict['mu_0']
    beta_0 = params_dict['beta_0']
    sigma_0 = params_dict['sigma_0']
    gamma_0 = params_dict['gamma_0']
    m = len(alpha)
    k = len(tau) // m

    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['tau'].reshape((m, k))[:, :-1].flatten(),
        params_dict['mu'].flatten(),
        params_dict['beta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma'],
        params_dict['rho'],
        params_dict['mu_0'].flatten(),
        params_dict['beta_0'].flatten(),
        params_dict['sigma_0'],
        params_dict['gamma_0']
    ])

    for i in range(n):
        y_it = y[:, i]
        z_it = z[(i * t):((i + 1) * t), :]
        x_it = x[(i * t):((i + 1) * t), :]
        score_i = score_i_ar1_mixture(y_it, x_it, z_it, params_array, m, k, q, p)
        hessian_i = numerical_derivative_score_i(y_it, x_it, z_it, params_array, m, k, q, p, score_i_ar1_mixture)
        
        scores.append(score_i)
        hessians.append(hessian_i)

    score = np.vstack(scores)
    hessian = np.array(hessians)
    return score, hessian

# %%
# Variance AR1 Normal
# ----------------------------------------------------------

def get_params_ar1_normal(out_h0):
    """
    Convert out_h0 dictionary into arrays compatible with flatten_params for the AR(1) normal model.

    Parameters:
    -----------
    out_h0 : dict
        Dictionary containing parameter arrays.

    Returns:
    --------
    tuple
        - params_dict: Dictionary with flattened parameter arrays.
        - params_array: Flattened array of parameters.
    """
    alpha_hat = out_h0['alpha_hat'][0]
    mu_hat = out_h0['mu_hat'][0]
    mubeta_hat = out_h0['mubeta_hat'][0]
    
    sigma_hat = out_h0['sigma_hat'][0]
    gamma_hat = out_h0['gamma_hat'][0]
    rho_hat = out_h0['rho_hat'][0]
    mu_0_hat = out_h0['mubeta_0_hat'][0]
    mubeta_0_hat = out_h0['mubeta_0_hat'][0]
    
    sigma_0_hat = out_h0['sigma_0_hat'][0]
    gamma_0_hat = out_h0['gamma_0_hat'][0]

    m = len(alpha_hat)
    q = (len(mubeta_hat) // m) - 1
    p = len(gamma_hat)

    mubeta_hat = np.ascontiguousarray(mubeta_hat)
    mubeta_hat_mat = mubeta_hat.reshape((q+1,m)).T
    beta_hat = mubeta_hat_mat[:,1:]
    mubeta_0_hat = np.ascontiguousarray(mubeta_0_hat)
    mubeta_0_hat_mat = mubeta_0_hat.reshape((q+1,m)).T
    beta_0_hat = mubeta_0_hat_mat[:,1:]
    mu_0_hat =  mubeta_0_hat_mat[:,0]

    params_dict = {
        'alpha': alpha_hat,
        'mu': mu_hat,
        'beta': beta_hat,
        'sigma': sigma_hat,
        'gamma': gamma_hat,
        'rho': rho_hat,
        'mu_0': mu_0_hat,
        'beta_0': beta_0_hat,
        'sigma_0': sigma_0_hat,
        'gamma_0': gamma_0_hat
    }

    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['mu'],
        params_dict['beta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma'],
        params_dict['rho'],
        params_dict['mu_0'],
        params_dict['beta_0'].flatten(),
        params_dict['sigma_0'],
        params_dict['gamma_0']
    ])
    return params_dict, params_array


def get_params_dict_from_array_ar1_normal(p_array, m, q, p):
    """
    Convert a flattened parameter array back into a params_dict for the AR(1) normal model.

    Parameters:
    -----------
    p_array : ndarray
        Flattened array of parameters.
    m : int
        Number of components.
    q : int
        Number of covariates in x.
    p : int
        Number of covariates in z.

    Returns:
    --------
    dict
        Dictionary containing the parameters.
    """
    alpha_f = p_array[:m - 1]
    mu = p_array[m - 1:m - 1 + m]
    beta = p_array[m - 1 + m:m - 1 + m + m * q]
    sigma = p_array[m - 1 + m + m * q:m - 1 + m + m * q + m]
    gamma = p_array[m - 1 + m + m * q + m:m - 1 + m + m * q + m + p]
    rho = p_array[m - 1 + m + m * q + m + p:m - 1 + m + m * q + m + p + m]
    mu_0 = p_array[m - 1 + m + m * q + m + p + m:m - 1 + m + m * q + m + p + m + m]
    beta_0 = p_array[m - 1 + m + m * q + m + p + m + m:m - 1 + m + m * q + m + p + m + m + m * q]
    sigma_0 = p_array[m - 1 + m + m * q + m + p + m + m + m * q:m - 1 + m + m * q + m + p + m + m + m * q + m]
    gamma_0 = p_array[m - 1 + m + m * q + m + p + m + m + m * q + m:]

    # Reconstruct alpha with the last component as 1 - sum of others
    alpha = np.append(alpha_f, 1 - np.sum(alpha_f))

    return {
        'alpha': alpha,
        'mu': mu,
        'beta': beta.reshape((m, q)),
        'sigma': sigma,
        'gamma': gamma,
        'rho': rho,
        'mu_0': mu_0,
        'beta_0': beta_0.reshape((m, q)),
        'sigma_0': sigma_0,
        'gamma_0': gamma_0
    }
    
def score_i_ar1_normal(y_it, x_it, z_it, params_array, m, k, q, p):
    """
    Compute the score for a single observation (i) in the AR(1) normal model.

    Parameters:
    -----------
    y_it : ndarray
        Array of dependent variable values for all t for a single observation i.
    x_it : ndarray
        Array of covariates in x for all t for a single observation i.
    z_it : ndarray
        Array of covariates in z for all t for a single observation i.
    params_array : ndarray
        Flattened array of parameters.

    Returns:
    --------
    ndarray
        Score vector for the single observation.
    """
    t = len(y_it)
    params_dict = get_params_dict_from_array_ar1_normal(params_array, m, q, p)

    alpha = params_dict['alpha']
    mu = params_dict['mu']
    beta = params_dict['beta']
    sigma = params_dict['sigma']
    gamma = params_dict['gamma']
    rho = params_dict['rho']
    mu_0 = params_dict['mu_0']
    beta_0 = params_dict['beta_0']
    sigma_0 = params_dict['sigma_0']
    gamma_0 = params_dict['gamma_0']

    residuals_matrix = np.zeros((m, t - 1))
    residuals_0_matrix = np.zeros((m, 1))
    w = np.zeros(m)

    for j in range(m):
        residuals = (y_it[1:] - x_it[1:] @ beta[j] - z_it[1:] @ gamma - mu[j]) - rho[j] * (y_it[:-1] - x_it[:-1] @ beta[j] - z_it[:-1] @ gamma - mu[j])
        residuals_0 = y_it[0] - x_it[0] @ beta_0[j] - z_it[0] @ gamma_0 - mu_0[j]
        residuals_matrix[j] = residuals
        residuals_0_matrix[j] = residuals_0

        log_likelihood = (
            -0.5 * (t - 1) * np.log(2 * np.pi * sigma[j] ** 2)
            - 0.5 * np.sum(residuals ** 2) / sigma[j] ** 2
            - 0.5 * np.log(2 * np.pi * sigma_0[j] ** 2)
            - 0.5 * residuals_0 ** 2 / sigma_0[j] ** 2
        )

        w[j] = alpha[j] * np.exp(log_likelihood)

    w /= np.sum(w)

    score_alpha = np.zeros(m - 1)
    score_mu = np.zeros(m)
    score_beta = np.zeros(m * q)
    score_sigma = np.zeros(m)
    score_gamma = np.zeros(p)
    score_rho = np.zeros(m)
    score_mu_0 = np.zeros(m)
    score_beta_0 = np.zeros(m * q)
    score_sigma_0 = np.zeros(m)
    score_gamma_0 = np.zeros(p)

    for j in range(m - 1):
        score_alpha[j] = w[j] / alpha[j] - w[-1] / alpha[-1]

    for j in range(m):
        residuals = residuals_matrix[j]
        residuals_0 = residuals_0_matrix[j]

        score_mu[j] = w[j] * np.sum(residuals * (1 - rho[j])) / sigma[j] ** 2
        score_mu_0[j] = w[j] * residuals_0 / sigma_0[j] ** 2

        for qq in range(q):
            score_beta[j * q + qq] = w[j] * np.sum(residuals * (x_it[1:, qq] - rho[j] * x_it[:-1, qq])) / sigma[j] ** 2
            score_beta_0[j * q + qq] = w[j] * residuals_0 * x_it[0, qq] / sigma_0[j] ** 2

        term1 = -(t - 1) / sigma[j]
        term2 = np.sum(residuals ** 2) / sigma[j] ** 3
        score_sigma[j] = w[j] * (term1 + term2)

        term1_0 = -1 / sigma_0[j]
        term2_0 = residuals_0 ** 2 / sigma_0[j] ** 3
        score_sigma_0[j] = w[j] * (term1_0 + term2_0)

        for pp in range(p):
            score_gamma[pp] += w[j] * np.sum(residuals * (z_it[1:, pp] - rho[j] * z_it[:-1, pp])) / sigma[j] ** 2
            score_gamma_0[pp] += w[j] * residuals_0 * z_it[0, pp] / sigma_0[j] ** 2

        score_rho[j] = w[j] * np.sum(residuals * (y_it[:-1] - x_it[:-1] @ beta[j] - z_it[:-1] @ gamma)) / sigma[j] ** 2

    score_i = np.concatenate([
        score_alpha, score_mu, score_beta, score_sigma, score_gamma, score_rho,
        score_mu_0, score_beta_0, score_sigma_0, score_gamma_0
    ])
    return score_i

def score_ar1_normal(data, params_dict):
    y, x, z = data
    t, n = y.shape
    q = x.shape[1]
    p = z.shape[1]

    scores = []
    hessians = []
    alpha = params_dict['alpha']
    mu = params_dict['mu']
    beta = params_dict['beta']
    sigma = params_dict['sigma']
    gamma = params_dict['gamma']
    rho = params_dict['rho']
    mu_0 = params_dict['mu_0']
    beta_0 = params_dict['beta_0']
    sigma_0 = params_dict['sigma_0']
    gamma_0 = params_dict['gamma_0']
    m = len(alpha)

    params_array = np.concatenate([
        params_dict['alpha'][:-1],
        params_dict['mu'],
        params_dict['beta'].flatten(),
        params_dict['sigma'],
        params_dict['gamma'],
        params_dict['rho'],
        params_dict['mu_0'],
        params_dict['beta_0'].flatten(),
        params_dict['sigma_0'],
        params_dict['gamma_0']
    ])

    for i in range(n):
        y_it = y[:, i]
        z_it = z[(i * t):((i + 1) * t), :]
        x_it = x[(i * t):((i + 1) * t), :]
        score_i = score_i_ar1_normal(y_it, x_it, z_it, params_array, m, 1, q, p)
        hessian_i = numerical_derivative_score_i(y_it, x_it, z_it, params_array, m, 1, q, p, score_i_ar1_normal)
        
        scores.append(score_i)
        hessians.append(hessian_i)

    score = np.vstack(scores)
    hessian = np.array(hessians)
    return score, hessian
# %%
# ----------------------------------------------------------
# Plotting Data Functions
# ----------------------------------------------------------

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

# Store the dictionary
ind_code_dict = {
    311.0: 'Food products',
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
    361.0: 'Manufacture of pottery, china and earthenware'
}

# Ensure the dictionary is loaded when the function is loaded
def load_ind_code_dict():
    return ind_code_dict


def process_chilean_data(each_code, T=3, p=0, y_indicator='mY_share'):
    """
    Process Chilean dataset given a specific `each_code` and `T`.
    
    Parameters:
    - df (pd.DataFrame): The input DataFrame containing the Chilean dataset.
    - each_code (str): The specific 'ciiu_3d' code to filter the dataset.
    - T (int): The time window to reshape the panel data (default is 3).
    - p (int): Dimensionality of the Z matrix (default is 1).
    
    Returns:
    - dict: A dictionary containing processed data and results.
    """
    t = time.time()  # Start timer
    chile_data = pyreadr.read_r('ChileanClean.rds')
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
    df_0['mL_share'] = (df_0['rMi'].fillna(0) + df_0['rMd'].fillna(0) ).fillna(0) / df_0['L'].fillna(0)
    df['id'] = pd.to_numeric(df['id'], errors='coerce')
    df['year'] = pd.to_numeric(df['year'], errors='coerce')

    # Select relevant columns from df_0
    df_0 = df_0[['id', 'year','cu' , 'ciiu', 'mi_share', 'mY_share', 'mL_share', 'rY']]

    # Perform the merge
    df = pd.merge(df[['id','ciiu_3d','year','K','L']], df_0, on=['id', 'year'])
    
    # Filter dataset for the code
    ind_each = df.loc[df['ciiu_3d'] == each_code, ['id','year','mi_share', 'mY_share', 'mL_share', 'K','L', 'rY', 'ciiu']]
    ind_each['lnk'] = np.log(ind_each['K'])
    ind_each['lnl'] = np.log(ind_each['L'])
    ind_each['lny'] = np.log(ind_each['rY'])
    ind_each['y'] = np.log(ind_each[y_indicator])
    ind_each['m_share'] = ind_each['mi_share'].fillna(0)
    year_list = sorted(ind_each['year'].unique())
    T_cap = max(year_list)
    
    t_start = T_cap - T + 1
    ind_each_t = ind_each[ind_each['year'] >= t_start].dropna()
    
    
    ind_each_y = ind_each_t.pivot(index='id', columns='year', values='y').dropna()
    
    N = ind_each_y.shape[0]
        
    id_list = ind_each_y.index
    ind_each_t = ind_each_t[ind_each_t['id'].isin(id_list)].sort_values(['id', 'year'])
    
    ind_each_xk = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()
    ind_each_xl = (ind_each_t['lnl'] - ind_each_t['lnl'].mean()) / ind_each_t['lnl'].std()
    
    ind_each_m_share = ind_each_t['m_share']

    y = ind_each_y.T.to_numpy().astype(np.float64)
    # y_ub = np.quantile(y, 1 - epsilon_truncation)
    # y_lb = np.quantile(y, epsilon_truncation)
    # y[y > y_ub] = y_ub
    # y[y < y_lb] = y_lb
            
    x_k = ind_each_xk.to_numpy().reshape(-1, 1).astype(np.float64)
    x_kl = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_xl.to_numpy().reshape(-1, 1)].astype(np.float64)
    x_kmshare = np.c_[ind_each_xk.to_numpy().reshape(-1, 1), ind_each_m_share.to_numpy().reshape(-1, 1)].astype(np.float64)
    
    x_0 = np.zeros((N * T, 0))
    x_0 = x_0.astype(np.float64)  # Convert x_0 to float64
    # Normalize X
    ind_each_x = (ind_each_t['lnk'] - ind_each_t['lnk'].mean()) / ind_each_t['lnk'].std()

    # Prepare data
    y = ind_each_y.T.to_numpy().astype(np.float64)
    x_k = ind_each_x.to_numpy().reshape(-1, 1)  # X as a matrix
    N = ind_each_y.shape[0]
    
    
    ciiu_value_count = ind_each_t.groupby('ciiu')['ciiu'].count()
    ciiu_combine = ciiu_value_count[ciiu_value_count < 0.1 * (N * T)].index
    for each_ciiu in ciiu_combine:
        ind_each_t.loc[ind_each_t['ciiu'] == each_ciiu, 'ciiu'] = 'other'
    z_ciiu = 1 * pd.get_dummies(ind_each_t['ciiu'], prefix='category').values[:, 1:]
    z_ciiu = z_ciiu.astype(np.float64)  # Convert z to float64
    
    # Return processed data and results
    results = {
        'year_list': year_list,
        'y': y,
        'x_0': x_0,
        'x_k': x_k,
        'x_kmshare': x_kmshare,
        'z_ciiu': z_ciiu,
        'execution_time': time.time() - t
    }

    return results
# %%
