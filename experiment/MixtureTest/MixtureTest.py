# %%
import numpy as np


# Function to generate data
def generate_data(alpha, mu, sigma, gam, beta, N, T, M, p, q, x=None, z=None):
    print(f"N = {N}")
    print(f"T = {T}")

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
from scipy.linalg import solve, LinAlgError
from scipy.optimize import minimize
from numpy.linalg import inv, cond

SINGULAR_EPS = 1e-10  # Criteria for matrix singularity
M_LN_SQRT_2PI = 0.9189385332046727  # log(sqrt(2*pi))

def EM_optimization(y, x, z, mu_0, sigma_0, initial_pars, m, t, an, maxit=2000, ninits=10, tol=1e-8, tau=0.5, h=0, k=0, epsilon=0.05):
    nt = len(y)
    n = nt // t
    
    if x is not None:
        if x.shape[0] != nt:
            raise ValueError("y and x must have the same number of rows.")
        x1 = np.hstack((np.ones((x.shape[0], 1)), x))
        q1 = x1.shape[1]
        q = q1 - 1
    else:
        x1 = np.ones((nt,1))
        q = 0
        q1 = 1
    
    # Handle z
    if z is not None:
        z = np.array(z)
        p = z.shape[1]
        if z.shape[0] != nt:
            raise ValueError("y and z must have the same number of rows.")
    else:
        p = 0
        gam = None
    
    lb = np.zeros(m)
    ub = np.zeros(m)
    
    alp_sig = np.zeros(m)
    r = np.zeros((m,n))

    l_j = np.zeros(m)
    w = np.zeros((m, nt))
    post = np.zeros((m * n, ninits))
    notcg = np.zeros(ninits)
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    gamma = np.zeros(p)
    ytilde = np.zeros(nt)
    xtilde = np.zeros((nt, q1))
    wtilde = np.zeros(nt)
    ztilde = np.zeros((nt, p))
    zz = np.zeros((p, p))
    ze = np.zeros((p, 1))
    
    # if k == 1:
    #     mu0[0] = -np.inf
    #     mu0[m] = np.inf
    #     for j in range(h):
    #         lb[j] = (mu0[j] + mu0[j + 1]) / 2.0
    #         ub[j] = (mu0[j + 1] + mu0[j + 2]) / 2.0
    #     for j in range(h, m):
    #         lb[j] = (mu0[j - 1] + mu0[j]) / 2.0
    #         ub[j] = (mu0[j] + mu0[j + 1]) / 2.0
    
    for jn in range(ninits):
        
        alpha = initial_pars['alpha'][:,jn]
        mubeta = initial_pars['mubeta'][:,jn] 
        sigma = initial_pars['sigma'][:,jn]
        
        if p > 0:
            gamma = initial_pars['gam'][:,jn]
        
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0

        for iter_ii in range(maxit):
            
            ll = -nt * M_LN_SQRT_2PI
            if p == 0:
                ytilde = y
            else:
                ytilde = y - z @ gamma
            # compute residual
            for j in range(m):
                if q == 0:
                    r_t = (1.0 / sigma[j]) * (ytilde - mubeta[j])
                    r_t = 0.5 * (r_t**2)
                    r[j] = t * np.log(sigma[j]) + r_t.reshape((n,t)).T.sum(axis=0)
            minr = np.min(r,axis=0)
            l_j = alpha[:, None] * np.exp(minr - r)
            sum_l_j = np.sum(l_j,axis=0)
               
            w = l_j / sum_l_j            
            ll += np.sum(np.log(sum_l_j) - minr)
            penloglik = ll + np.log(2.0) + min(np.log(tau), np.log(1 - tau))
            
            for j in range(m):
                s0j = sigma_0[j] / sigma[j]
                penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik

            # if np.max(np.abs(diff)) < tol:
            #     break
            emit += 1

            for j in range(m):
                alpha[j] = np.mean(w[j, :]) 
                wtilde = w[j, :].T
                w_j = np.repeat(w[j, :], t).flatten()
                if q == 0:
                    mubeta[j] = np.sum(w_j * ytilde) / np.sum(w_j)                    
                    ssr_j = np.sum( w_j * ((ytilde - mubeta[j])**2 ) )
                else:
                    for ii in range(q1):
                        xtilde[:, ii] = w_j * x1[:, ii]
                    design_matrix = xtilde.T @ x1
                    if np.linalg.cond(design_matrix) < SINGULAR_EPS:
                        sing = 1
                        break
                    mubeta[:, j] = np.linalg.solve(design_matrix, xtilde.T @ ytilde)
                    ssr_j = np.sum( np.repeat(w[j, :], t)  * (ytilde - x1 @ mubeta[:, j])**2)
                sigma[j] = np.sqrt((ssr_j + 2.0 * an * sigma_0[j]**2) / ( np.sum(w_j) + 2.0 * an) )
                sigma[j] = max(sigma[j], epsilon * sigma_0[j])
                # if k == 1:
                #     mubeta[0, j] = min(max(mubeta[0, j], lb[j]), ub[j])

            alpha = np.clip(alpha, 0.01, 0.99)
            alpha /= np.sum(alpha)

            if k == 1:
                alphah = alpha[h - 1] + alpha[h]
                alpha[h - 1] = alphah * tau
                alpha[h] = alphah * (1 - tau)
            elif k > 1:
                alphah = alpha[h - 1] + alpha[h]
                tauhat = alpha[h - 1] / (alpha[h - 1] + alpha[h])
                if tauhat <= 0.5:
                    tau = min((alpha[h - 1] * n + 1.0) / (alpha[h - 1] * n + alpha[h] * n + 1.0), 0.5)
                else:
                    tau = max(alpha[h - 1] * n / (alpha[h - 1] * n + alpha[h] * n + 1.0), 0.5)
                alpha[h - 1] = alphah * tau
                alpha[h] = alphah * (1 - tau)

            if p > 0:
                zz.fill(0)
                ze.fill(0)
                for j in range(m):
                    wtilde = w[j, :]
                    for ii in range(p):
                        ztilde[:, ii] = wtilde * z[:, ii]
                    zz += (ztilde.T @ z) / (sigma[j]**2)
                    ze += (ztilde.T @ (y - x1 @ mubeta[:, j])) / (sigma[j]**2)
                # if np.linalg.cond(zz) < SINGULAR_EPS or np.isnan(ze).any():
                #     sing = 1
                #     break
                gamma = solve(zz, ze)

            for j in range(m):
                if alpha[j] < 1e-8 or np.isnan(alpha[j]) or sigma[j] < 1e-8:
                    sing = 1

            # if sing:
            #     notcg[jn] = 1
            #     break
        
        penloglikset[jn] = penloglik.sum()
        loglikset[jn] = ll
        post[:, jn] = w.T.flatten()
        
        initial_pars['alpha'][:,jn] = alpha
        initial_pars['mubeta'][:,jn] = mubeta
        initial_pars['sigma'][:,jn] = sigma
        
        if p > 0:
            initial_pars['gam'][:,jn] = gamma
        
    
    return {
        "parameters" : initial_pars,
        "penloglikset": penloglikset,
        "loglikset": loglikset,
        "notcg": notcg,
        "post": post
    }
    
# %%
import numpy as np

from itertools import product

def regpanelmixPMLEinit(y, x, z=None, ninits=1, m=2):
    """
    Initialize parameters for regpanelmixPMLE.

    Parameters:
    - y: Response vector (1D array-like).
    - x: Predictor matrix (2D array-like).
    - z: Optional covariate matrix (2D array-like). Default is None.
    - ninits: Number of initializations. Default is 1.
    - m: Number of components in the mixture. Default is 2.
    - model_ar1: Whether the model is AR1. Default is False.
    - z_init: Initialization for AR1 model. Default is None.

    Returns:
    - A dictionary containing initialized `alpha`, `mubeta`, `sigma`, `gam`, and `gamma0`.
    """
    y = y.flatten()
    if x is not None:
        x1 = np.hstack((np.ones((x.shape[0], 1)), x))
        q1 = x1.shape[1]
    else:
        q = 0
        q1 = 1
        x1 = np.ones((len(y),1))
    # Handle z
    if z is not None:
        z = np.array(z)
        p = z.shape[1]
    else:
        p = 0
        gam = None

    
    if x is not None and z is not None:
        # Case 1: Both x and z are not None
        xz = np.hstack((x1, z))
        
    elif x is None and z is not None:
        # Case 2: x is None, z is not None
        xz = z
    elif x is not None and z is None:
        # Case 3: x is not None, z is None
        xz = x1
    else:
        # Case 4: Both x and z are None
        xz = np.ones((len(y),1))
    
    out_coef = np.linalg.lstsq(xz, y, rcond=None)[0]  # Linear regression
    residuals = y - xz @ out_coef
    stdR = np.std(residuals)
    
    if z is not None:
        # Perform least squares regression with both x and z
        gam0 = out_coef[q1:(q1 + p)]
        gam = np.random.uniform(0.5, 1.5, size=(p, ninits)) * gam0[:, None]
        mubeta_hat = out_coef[:q1]
        y = y - z @ gam0
        
    else:
        # Perform least squares regression with x only
        mubeta_hat = out_coef
        gam = None

    # Initialize alpha
    alpha = np.random.uniform(size=(m, ninits))
    alpha = (alpha / np.sum(alpha, axis=0))

    # Initialize mubeta
    if x is not None:
        minMU = np.min(y - x @ mubeta_hat[1:])
        maxMU = np.max(y - x @ mubeta_hat[1:])
        mubeta = np.zeros((q1 * m, ninits))
        for j in range(m):
            mubeta[q1 * j, :] = np.random.uniform(minMU, maxMU, size=ninits)
            for i in range(1, q1):
                mubeta[q1 * j + i, :] = mubeta_hat[i] * np.random.uniform(-2, 2, size=ninits)
    else:
        minMU = np.min(y)
        maxMU = np.max(y)
        mubeta = np.zeros((q1 * m, ninits))
        for j in range(m):
            mubeta[q1 * j, :] = np.random.uniform(minMU, maxMU, size=ninits)
        
    # Initialize sigma
    sigma = np.random.uniform(0.01, 1, size=(m, ninits)) * stdR

    # Return the initialized values
    return {
        "alpha": alpha,
        "mubeta": mubeta,
        "sigma": sigma,
        "gam": gam,
    }
    
# %%
def regpanelmixPMLE(y, x, m=2, z=None, ninits=10, epsilon=1e-8, maxit=2000, epsilon_short=1e-2, maxit_short=500, binit=None):
    
    
    t = y.shape[0]
    n = y.shape[1]
    nt = n * t
    y_mat = y.copy()
    y = y.T.flatten()
    
    
    if x is not None:
        if x.shape[0] != nt:
            raise ValueError("y and x must have the same number of rows.")
        x1 = np.hstack((np.ones((x.shape[0], 1)), x))
        q1 = x1.shape[1]
    else:
        q = 0
        q1 = 1
    
    # Handle z
    if z is not None:
        z = np.array(z)
        p = z.shape[1]
        if z.shape[0] != nt:
            raise ValueError("y and z must have the same number of rows.")
    else:
        p = 0
        gam = None
    
    ninits_short = ninits * 10 * (q1 + p) * m
    
    
    npar = m - 1 + (q1 + 1) * m + p
    
    if x is not None and z is not None:
        # Case 1: Both x and z are not None
        xz = np.hstack((x, z))
    elif x is None and z is not None:
        # Case 2: x is None, z is not None
        xz = z
    elif x is not None and z is None:
        # Case 3: x is not None, z is None
        xz = x
    else:
        # Case 4: Both x and z are None
        xz = np.ones((nt,1))
    
    ls_out = np.linalg.lstsq(xz, y, rcond=None)
    sd0 = np.sqrt(ls_out[1]/nt)
    
    if m == 1:
        mubeta = ls_out[0][:q1]
        if z is not None:
            gam = ls_out[0][q1:(q1 + p)]
        res = y - xz @ ls_out[0]
        sigma = np.sqrt(np.mean(res**2))
        
        loglik = - (nt / 2) * (1 + np.log(2 * np.pi) + 2 * np.log(sigma))
        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        penloglik = loglik
        
        parlist = {"alpha": 1, "mubeta": mubeta, "sigma": sigma, "gam": gam}
        coefficients = {"alpha": 1, "mubeta": mubeta, "sigma": sigma, "gam": gam}
        postprobs = np.ones(n)
    else:  # m >= 2
        # Generate initial values.
        initial_pars = regpanelmixPMLEinit(y=y, x=x, z=z, ninits=ninits_short, m=m)
        
        h = 0  # Setting h=0 gives PMLE
        tau = 0.5  # Setting tau=0.5 gives PMLE
        k = 0  # Setting k=0 gives PMLE
        an = 1 / n  # Penalty term for variance
        an_0 = 0.3  # Default value from KS 15
        
        sigma_0 = np.full(m, sd0)
        mu_0 = np.zeros(m + 1)  # Dummy
    
        ztilde = np.zeros_like(z) if z is not None else None
        
        out_short = EM_optimization(y,x,z, mu_0, sigma_0, initial_pars, m, t, an )
        
        
        components = np.argsort(out_short["penloglikset"])[::-1][:ninits]
        if z is not None:
            long_pars = {'alpha': initial_pars['alpha'][:,components], 'mubeta': initial_pars['mubeta'][:,components], 'sigma': initial_pars['sigma'][:,components], 'gam': initial_pars['gam'][:,components]}
        else:
            long_pars = {'alpha': initial_pars['alpha'][:,components], 'mubeta': initial_pars['mubeta'][:,components], 'sigma': initial_pars['sigma'][:,components], 'gam': None}
            
        out = EM_optimization(y,x,z, mu_0, sigma_0, long_pars, m, t, an)
    
        index = np.argmax(out["penloglikset"])
        alpha = long_pars['alpha'][:,index]
        mubeta = long_pars['mubeta'][:,index]
        sigma = long_pars['sigma'][:,index]
        
        
        penloglik = out["penloglikset"][index]
        loglik = out["loglikset"][index]
        postprobs = out["post"][:, index].reshape(n, -1)
        
        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        
        mu_order = np.argsort(mubeta)
        alpha = alpha[mu_order]
        mubeta = mubeta[mu_order]
        sigma = sigma[mu_order]
        
        
        postprobs = postprobs[:, mu_order]
        
        parlist = {
            "alpha": alpha,
            "mubeta": mubeta,
            "sigma": sigma,
            "gam": gam
        }
        
        
    result = {
        "parlist": parlist,
        "loglik": loglik,
        "penloglik": penloglik,
        "aic": aic,
        "bic": bic,
        "postprobs": postprobs,
        "m": m,
    }
    
    return result
# %%
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
# %%
