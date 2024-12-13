# %%
import numpy as np

# Input Parameters
Nset = [200, 400]
Tset = [3, 5, 8]
alphaset = [[0.5, 0.5], [0.2, 0.8]]
muset = [[-1, 1], [-0.5, 0.5]]
sigmaset = [[0.8, 1.2]]

# Test panel mixture
alpha = [0.5, 0.5]
mu = [-1, 1]
mu = [-0.5, 0.5]
sigma = [0.8, 1.2]
gam = None
beta = None
N = 200
T = 3
M = 2
p = 0
q = 0

data = generate_data(alpha, mu, sigma, gam, beta, N, T, M, p, q)
y = data["Y"]
z = data["Z"]
x = data["X"]

regpanelmixPMLE(y,x, m = 2)