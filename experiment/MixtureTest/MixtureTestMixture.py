# %%
import numpy as np


def generate_data(alpha, mu, sigma, tau, N, T, M, K, p, q, x=None, z=None):

    prior = np.random.uniform(size=N)
    alpha_cum = np.cumsum([0] + list(alpha))
    if M > 1:
        for m in range(M):
            lb = alpha_cum[m]
            ub = alpha_cum[m + 1]
            R[:, m] = ((prior > lb) & (prior <= ub)).astype(int)
    else:
        R = np.ones((N, M))
    
    R_sub = {}
    m = 1
    for m in range(M):
        prior_m = np.random.uniform(size=N)
        tau_cum = np.cumsum([0] + list(tau[m]))
        R_sub[m] = np.zeros((N,M))
        for k in range(K):
            lb = tau_cum[k]
            ub = tau_cum[k + 1]
            R_sub[m][:, k] = ((prior_m > lb) & (prior_m <= ub)).astype(int)
        
    return {"Y": Y, "Z": z, "X": x}

# %%
# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
# alphaset = [[0.5, 0.5], [0.2, 0.8]]

# Test panel mixture
alpha = [0.5, 0.5]
mu = [[-1, 1],[-0.5, 0.5]]
tau = [[0.5, 0.5],[0.5, 0.5]]
sigma = [0.8, 1.2]
gam = None
beta = None

N = 200
T = 3
M = 2
K = 2 # conditional probability is k-components model
p = 0
q = 0


data = generate_data(alpha, mu, sigma, gam, beta, N, T, M, p, q)
y = data["Y"]
z = data["Z"]
x = data["X"]

regpanelmixPMLE(y,x, m = 2)