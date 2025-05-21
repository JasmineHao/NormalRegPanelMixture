import numpy as np

# Example parameter values

alpha_hat = np.array([0.5, 0.5])  # mixing proportions
tau_hat = np.array([[0.5, 0.5], [0.5, 0.5]])  # mixing proportions
mu_hat = np.array([[-2, -1.5], [-1, -0.5]])
sigma_hat = np.array([0.173, 0.403])
beta_hat = np.array([
    [0.052,  0.293, -0.141, -0.022,  0.047],
    [0.104, -0.466, -0.056,  0.111,  0.016]
])
gamma_hat = np.array([0.,0.])
mu_0_hat = np.array([[-2, -1.5], [-1, -0.5]])
sigma_0_hat = np.array([0.161, 0.351])
beta_0_hat = np.array([
    [-0.04,  0.027],
    [-0.007, 0.099]
])
gamma_0_hat = np.array([0.])


data = generate_data_ar1_mixture_noconstraint(alpha_hat, tau_hat, mu_hat, sigma_hat, beta_hat, gamma_hat, mu_0_hat, sigma_0_hat, beta_0_hat, gamma_0_hat, N, T, m, k, p, q)

regpanelmixmixtureAR1NoConstraintPMLE(data[0], x, z, p, q, m, k, alpha_bound=alpha_bound,ninits=2)

