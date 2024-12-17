
# %%
import numpy as np


# Function to generate data
def generate_data(alpha, mu, sigma, gam, beta, N, T, M, p, q, x=None, z=None):
    # print(f"N = {N}")
    # print(f"T = {T}")

    R = np.zeros((N, M))
    if sum(alpha) != 1:
        alpha = np.array(alpha) / sum(alpha)

    if len(alpha) != M or len(mu) != M:
        raise ValueError("M must be the size of alpha and mu")
    
    prior = np.random.uniform(size=N)
    alpha_cum = np.cumsum([0] + list(alpha))
    
    if M > 1:
        for m in range(M):
            lb = alpha_cum[m]
            ub = alpha_cum[m + 1]
            R[:, m] = ((prior > lb) & (prior <= ub)).astype(int)
    else:
        R = np.ones((N, M))

    Y = np.zeros((T, N))
    
    if q != 0 and x is None:
        x = np.random.normal(size=(N * T, q))
    if p != 0 and z is None:
        z = np.random.normal(size=(N * T, p))
    
    mu_R = np.dot(R, mu)
    sigma_R = np.dot(R, sigma)
    u = np.random.normal(size=(T, N))
    
    for nn in range(N):
        y_nn = np.zeros(T)
        y_nn = mu_R[nn] + sigma_R[nn] * u[:, nn]
        
        if q > 1:
            beta_R = np.dot(R, beta)
            y_nn += np.dot(x[(T * nn):(T * (nn + 1)), :], beta_R[nn, :])
        elif q == 1:
            beta_R = np.dot(R, np.ravel(beta))
            y_nn += x[(T * nn):(T * (nn + 1)), 0] * beta_R[nn]
        
        if p > 1:
            y_nn += np.dot(z[(T * nn):(T * (nn + 1)), :], gam)
        elif p == 1:
            y_nn += z[(T * nn):(T * (nn + 1)), 0] * gam
        
        Y[:, nn] = y_nn
    
    if p == 0:
        z = None
    if q == 0:
        x = None
    
    return {"Y": Y, "Z": z, "X": x}

# %%
import numpy as np
from scipy.linalg import svd, sqrtm, inv, det, block_diag
from scipy.stats import rankdata
from functools import reduce
from numpy.random import default_rng
from scipy.stats import expon

# Function to compute the square root of a matrix using SVD
def matrix_sqrt_svd(mat):
    if not isinstance(mat, np.ndarray):
        raise ValueError("Input must be a matrix (NumPy array).")
    
    # Check if the matrix is square
    if mat.shape[0] != mat.shape[1]:
        raise ValueError("Input must be a square matrix.")
    
    # Perform SVD decomposition
    U, S, VT = svd(mat)
    
    # Compute the square root of singular values
    sqrt_singular_values = np.sqrt(S)
    
    # Reconstruct the square root matrix
    sqrt_mat = U @ np.diag(sqrt_singular_values) @ VT
    return sqrt_mat


# Function to invert a matrix, with regularization for singular matrices
def invert_matrix(mat, epsilon=1e-8):
    if not isinstance(mat, np.ndarray) or mat.shape[0] != mat.shape[1]:
        raise ValueError("The input must be a square matrix.")
    
    # Check if the determinant is close to zero
    det_val = np.linalg.det(mat)
    if abs(det_val) < epsilon:
        # Regularize the matrix by adding epsilon to the diagonal
        mat = mat + np.eye(mat.shape[0]) * epsilon
    
    # Return the inverse
    return np.linalg.inv(mat)

# Function to perform SVD decomposition and compute various matrices
def matrix_svd_decomposition(P, m):
    # Perform SVD
    U, S, VT = svd(P, full_matrices=True)
    
    # Submatrices of U and VT
    U_12 = U[:m, m:]
    V_12 = VT[:m, m:]
    
    U_22 = U[m:, m:]
    V_22 = VT[m:, m:]
    
    # Construct A_q_o and B_q_o matrices
    A_q_o = np.transpose(
        sqrtm(U_22 @ U_22.T) @ invert_matrix(U_22.T) @ np.hstack([U_12.T, U_22.T])
    )
    B_q_o = sqrtm(V_22 @ V_22.T) @ invert_matrix(V_22.T) @ np.hstack([V_12.T, V_22.T])
    
    # Compute the Kronecker product of B_q_o and A_q_o.T
    kron_BA_o = np.kron(B_q_o, A_q_o.T)
    
    return {
        "D": S,
        "U": U,
        "V": VT.T,
        "U_12": U_12,
        "V_12": V_12,
        "U_22": U_22,
        "V_22": V_22,
        "A_q_o": A_q_o,
        "B_q_o": B_q_o,
        "kron_BA_o": kron_BA_o,
    }

# Function to construct the statistic KP
def construct_stat_KP(P, Sigma_P, m, n_size, lambda_c=0, transform="P"):
    allowed_strings = ["P", "Q"]
    
    if transform not in allowed_strings:
        raise ValueError(f"Invalid input! The transform must be one of: {', '.join(allowed_strings)}")
    
    if transform == "P":
        # Perform SVD decomposition
        P_svd = matrix_svd_decomposition(P, m)
        lambda_q = (
            P_svd["A_q_o"].T @ P @ P_svd["B_q_o"].T - lambda_c
        )
        Omega_q = P_svd["kron_BA_o"] @ Sigma_P @ P_svd["kron_BA_o"].T
        
        r = Omega_q.shape[0]
        rk_c = n_size * np.sum(
            np.array(lambda_q).flatten() @ invert_matrix(Omega_q) @ np.array(lambda_q).flatten()
        )
    else:
        Q = P @ P.T
        J_P = np.kron(P, np.eye(n_size))
        Sigma_Q = J_P @ Sigma_P @ J_P.T
        
        Q_svd = matrix_svd_decomposition(Q, m)
        lambda_q = (
            Q_svd["A_q_o"].T @ Q @ Q_svd["B_q_o"].T - lambda_c
        )
        Omega_q = Q_svd["kron_BA_o"] @ Sigma_Q @ Q_svd["kron_BA_o"].T
        
        r = Omega_q.shape[0]
        rk_c = n_size * np.sum(
            np.array(lambda_q).flatten() @ invert_matrix(Omega_q) @ np.array(lambda_q).flatten()
        )
    
    AIC_c = rk_c - 2 * r
    BIC_c = rk_c - np.log(n_size) * r
    HQ_c = rk_c - 2 * np.log(np.log(n_size)) * r
    
    return {
        "rk_c": rk_c,
        "lambda_c": lambda_q,
        "Omega_q": Omega_q,
        "AIC_c": AIC_c,
        "BIC_c": BIC_c,
        "HQ_c": HQ_c,
    }
    

def create_indicator_list(data_c, T, N, n_bins):
    indicator_list = []
    for t in range(T):
        # Calculate quantiles for the given number of bins
        quantiles = np.quantile(data_c[t, :], q=np.linspace(0, 1, n_bins + 1))
        quantiles[0] = -np.inf
        quantiles[-1] = np.inf
        
        # Create indicator matrix using quantile bins
        cut_indices = np.digitize(data_c[t, :], bins=quantiles, right=False) - 1
        indicator_matrix = np.zeros((N, n_bins))
        indicator_matrix[np.arange(N), cut_indices] = 1
        indicator_list.append(indicator_matrix)
    return indicator_list

# Function to calculate P matrices and Sigma matrices for triplets
def calculate_P_matrix(data_c, weights, n_grid=3):
    T = data_c.shape[0]
    N = data_c.shape[1]
    
    # Create `indicator_list_Y` with 2 bins
    indicator_list_Y = create_indicator_list(data_c, T, N, n_bins=2)

    # Create `indicator_list_Y_ngrid` with `n_grid` bins
    indicator_list_Y_ngrid = create_indicator_list(data_c, T, N, n_bins=n_grid)
    
    # Initialize the result lists
    P_k_list = []
    Sigma_P_k_list = []
    
    # Iterate over the t periods
    for k in range(T):
        # Compute the Kronecker product for each row
        result_matrix = np.array([
            reduce(
                np.kron,
                [indicator_list_Y[t][n, :] for t in [i for i in range(T) if i != k]]
            )
            for n in range(N)
        ])
        
        # Compute P_k
        P_k = (weights * indicator_list_Y_ngrid[k].T) @ result_matrix
        P_k_list.append(P_k)
        
        # Compute Sigma_P_k
        P_k_vec = P_k.flatten()
        W_P_s = np.diag(P_k_vec) - np.outer(P_k_vec, P_k_vec)
        Sigma_P_k_list.append(W_P_s)
    
    return {
        "P_k_list": P_k_list,
        "Sigma_P_k_list": Sigma_P_k_list
    }

# Function to compute rk statistics for triplet T
def compute_rk_statistics(data, N, m, n_grid=3):
    # Initialize weights
    weights = np.full(N, 1 / N)
    
    # Compute P and Sigma matrices based on T_triplet_list
    data_P_W = calculate_P_matrix(data["Y"], weights, n_grid=n_grid)
    
    # Initialize results
    rk = np.zeros(T)
    lambda_c = []
    omega_c = []
    Sigma_P_list = []
    P_k_list = []
    
    # Loop through T_triplet_list to compute statistics
    for k in range(T):
        # Extract P_k and Sigma_P_k from the data_P_W object
        P_k = data_P_W["P_k_list"][k]
        Sigma_P_k = data_P_W["Sigma_P_k_list"][k]
        
        # Compute KP statistics for the k-th triplet
        stat_KP = construct_stat_KP(P_k, Sigma_P_k, m, N, transform="P")
        
        # Store results
        rk[k] = stat_KP["rk_c"]
        lambda_c.append(stat_KP["lambda_c"])
        omega_c.append(stat_KP["Omega_q"])
        Sigma_P_list.append(Sigma_P_k)
        P_k_list.append(P_k)
    
    return {
        "rk": rk,
        "lambda_c": lambda_c,
        "omega_c": omega_c,
        "Sigma_P_list": Sigma_P_list,
        "P_k_list": P_k_list
    }


def construct_stat_KP_P_triplet_bootstrap_combined(
    data=None,
    P_k_list=None,
    Sigma_P_list=None,
    N=None,
    T=None,
    BB=199,
    r_test=2,
    lambda_c=0,
    n_grid=3,
    transform="P",
    method="parametric"
):
    """
    Combined function for both parametric and smoothed nonparametric bootstrap.

    Parameters:
        data (dict): Dataset containing "Y" for nonparametric bootstrap.
        P_k_list (list): List of P matrices for parametric bootstrap.
        Sigma_P_list (list): List of Sigma matrices for parametric bootstrap.
        T_triplet_list (list): List of time triplets.
        N (int): Number of observations.
        BB (int): Number of bootstrap repetitions.
        r_test (int): Test parameter for `construct_stat_KP`.
        lambda_c (list): List of lambda_c values for each triplet.
        n_grid (int): Number of bins for quantiles in grid.
        transform (str): Transformation type ("P" or "Q").
        method (str): Bootstrap method ("parametric" or "nonparametric").

    Returns:
        dict: Dictionary containing bootstrap results for rk_b.
    """
    # Initialize result matrices
    rk_b = np.zeros((BB, T))
    
    if method == "parametric":
        # Parametric Bootstrap
        rng = default_rng()  # Random number generator
        for k in range(T):
            Sigma_k_0 = Sigma_P_list[k]
            vec_k_0 = P_k_list[k]
            
            # Generate bootstrap samples from multivariate normal
            vec_BB = rng.multivariate_normal(
                mean=vec_k_0.flatten(),
                cov=Sigma_k_0 / N,
                size=BB
            )
            
            for i in range(BB):
                # Reshape bootstrap sample to matrix form
                vec_b = vec_BB[i].reshape(vec_k_0.shape)
                
                # Compute KP statistics for the k-th triplet
                stat_KP = construct_stat_KP(
                    vec_b,
                    Sigma_k_0,
                    r_test,
                    N,
                    lambda_c=lambda_c[k],
                    transform=transform
                )
                
                # Update result matrix
                rk_b[i, k] = stat_KP["rk_c"]
    
    elif method == "nonparametric":
        # Smoothed Nonparametric Bootstrap
        # Generate ru matrix (random weights normalized by row sum)
        ru = expon.rvs(scale=1, size=(BB, N))  # Exponential random variables
        ru /= ru.sum(axis=1, keepdims=True)   # Normalize by row sums
        
        for i in range(BB):
            # Calculate bootstrapped P and Q matrices
            data_P_W_b = calculate_P_matrix(
                data["Y"], ru[i, :], n_grid=n_grid
            )
            
            for k in range(T):
                P_k = data_P_W_b["P_k_list"][k]
                Sigma_P_k = data_P_W_b["Sigma_P_k_list"][k]
                
                # Compute KP statistics for the k-th triplet
                stat_KP = construct_stat_KP(
                    P_k,
                    Sigma_P_k,
                    r_test,
                    N,
                    lambda_c=lambda_c[k],
                    transform=transform
                )
                
                # Update result matrix
                rk_b[i, k] = stat_KP["rk_c"]
    
    else:
        raise ValueError("Invalid method. Choose 'parametric' or 'nonparametric'.")
    
    # Return results as a dictionary
    return {
        "rk_b": rk_b
    }


def compute_statistics_for_rep(ii, Data, N, M, BB, r_test=2, n_grid=3):
    """
    Function to compute statistics for a single replication (similar to one iteration of R's foreach loop).
    """
    # Extract the data for the current replication
    data = Data[ii]
    
    # Compute pairwise statistics
    result_rk = compute_rk_statistics(data, N, m=2, n_grid=3)

    stats_KP_boot = construct_stat_KP_P_triplet_bootstrap_combined(
    data=data,
    P_k_list=result_rk['P_k_list'],
    Sigma_P_list=result_rk['Sigma_P_list'],
    N=N,
    T=T,
    BB=199,
    r_test=2,
    lambda_c=result_rk['lambda_c'],
    n_grid=3,
    transform="P",
    method="nonparametric") 
    
    # Return results as a dictionary
    return {
        "rk": result_rk["rk"],
        "rk_b": stats_KP_boot["rk_b"],
        "omega_c": result_rk["omega_c"],
        
    }
# Define a regular function to replace the lambda
def execute_compute_statistics(ii, Data, N, M, BB, r_test=2, n_grid=3):
    return compute_statistics_for_rep(ii, Data, N, M, BB, r_test, n_grid)

# %%
import concurrent.futures
from functools import partial

# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
alphaset = [[0.5, 0.5], [0.2, 0.8]]
muset = [[-1, 1], [-0.5, 0.5]]
sigmaset = [[0.8, 1.2]]

# Test panel mixture
alpha = [0.5, 0.5]
mu = [-1, 1]
sigma = [0.8, 1.2]
gam = None
beta = None
N = 200
T = 3
M = 2
p = 0
q = 0
nrep = 100
BB = 199

import numpy as np
import concurrent.futures
import time

# Initialize result matrix
result_matrix = np.zeros((len(alphaset) * len(muset), 2))  # Adjust dimensions based on `alphaset` and `muset`
count = 0

# Loop over parameters
for alpha in alphaset:
    for mu in muset:
        start_time = time.time()  # Start timer for this iteration
        
        # Generate data
        Data = [generate_data(alpha, mu, sigma, gam, beta, N, T, M, p, q) for _ in range(nrep)]
        
        # Use ProcessPoolExecutor for parallel execution
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(
                execute_compute_statistics,
                range(nrep),
                [Data] * nrep,
                [N] * nrep,
                [M] * nrep,
                [BB] * nrep,
                [2] * nrep,  # r_test
                [3] * nrep   # n_grid
            ))

        # Extract results and compute metrics
        rk_max = np.array([np.max(res['rk']) for res in results])
        rk_mean = np.array([np.mean(res['rk']) for res in results])
        rk_max_crit = np.array([
            np.quantile(np.apply_along_axis(np.max, 1, res['rk_b']), 0.95) for res in results
        ])
        rk_mean_crit = np.array([
            np.quantile(np.mean(res['rk_b'], axis=1), 0.95) for res in results
        ])

        # Update result matrix
        result_matrix[count, :] = [
            np.mean(rk_mean > rk_mean_crit),
            np.mean(rk_max > rk_max_crit),
        ]

        count += 1  # Increment count
        
        # Print time taken for this iteration
        end_time = time.time()
        print(f"Completed count {count} in {end_time - start_time:.2f} seconds.")

print(result_matrix)

# %%