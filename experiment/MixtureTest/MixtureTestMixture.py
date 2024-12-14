# %%
from types import NoneType
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import copy

def generate_data(alpha, mu, sigma, tau, N, T, M, K, p, q, x=None, z=None):
    if q != 0 and x is None:
        x = np.random.normal(size=(N * T, q))
    if p != 0 and z is None:
        z = np.random.normal(size=(N * T, p))
    
    prior = np.random.uniform(size=N)
    alpha_cum = np.cumsum([0] + list(alpha))
    R = np.ones((N, M))
    if M > 1:
        for mm in range(M):
            lb = alpha_cum[mm]
            ub = alpha_cum[mm + 1]
            R[:, mm] = ((prior > lb) & (prior <= ub)).astype(int)
    
    R_sub = {}
    mu_R_sub = np.zeros((N, M))
    if q > 0:
        beta_R_sub = np.zeros((N, M, q))
    
    for mm in range(M):
        prior_m = np.random.uniform(size=N)
        tau_cum = np.cumsum([0] + list(tau[mm]))
        R_sub[mm] = np.zeros((N,M))
        for kk in range(K):
            lb = tau_cum[kk]
            ub = tau_cum[kk + 1]
            R_sub[mm][:, kk] = ((prior_m > lb) & (prior_m <= ub)).astype(int)
        mu_R_sub[:, mm] = np.dot(R_sub[mm], mu[mm])
        if q > 0:
            beta_R_sub[:,mm,:] = np.dot(R_sub[mm], beta[mm])[:,None]
    
    sigma_R = np.dot(R, sigma)    
    mu_R = (R * mu_R_sub).sum(axis=1)
    
    if q > 0:
        beta_R = (R[:,:,np.newaxis] * beta_R_sub).sum(axis=1)
    u = np.random.normal(size=(T, N))
    Y = np.zeros((T, N))
    
    for nn in range(N):
        y_nn = np.zeros(T)
        y_nn = mu_R[nn] + sigma_R[nn] * u[:, nn]
        if q > 0:
            y_nn += np.dot(x[(T * nn):(T * (nn + 1)), :], beta_R[nn, :])
        if p > 1:
            y_nn += np.dot(z[(T * nn):(T * (nn + 1)), :], gam)
        Y[:, nn] = y_nn
    return {"Y": Y, "Z": z, "X": x}

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
    term1 = -n / 2 * np.log(2 * np.pi)
    term2 = -n * np.log(sigma)
    term3 = -1 / (2 * sigma**2) * np.sum((y - mu)**2)
    return term1 + term2 + term3

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
    # Constant term: -0.5 * log(2 * pi)
    constant = -0.5 * np.log(2 * np.pi)
    
    # Log-likelihood for each element
    log_likelihoods = (
        constant 
        - np.log(sigma) 
        - ((y - mu) ** 2) / (2 * sigma ** 2)
    )
    return log_likelihoods

# %%
def reg_panelmix_pmle(data,initial_pars,m,k, sigma_0, maxit=1000, tol = 1e-8, epsilon=0.05):

    tau_penalty = 0.5
    
    tol = 1e-8
    my_inf = 1e30
    t,n = data["Y"].shape
    nt = n * t
    an = 1/n
    y = data["Y"].T.flatten()
    x = data["X"]
    z = data["Z"]

    if x is not None:
        x1 = np.hstack((np.ones((x.shape[0], 1)), x))
        q1 = x1.shape[1]
    else:
        q = 0
        q1 = 1
    ninits = initial_pars['alpha'].shape[1]
    lb = np.zeros(m)
    ub = np.zeros(m)
    
    l_m = np.zeros((m,n))
    w   = np.zeros((m,n))

    l_m_k = [np.zeros((k,nt)) for mm in range(m)]
    w_m = [np.zeros((k, nt)) for mm in range(m)]
    res = [np.zeros((m, nt)) for mm in range(m)]
    
    zz = np.zeros((p, p))
    ze = np.zeros(p)
    
    post = np.zeros((m * n, ninits))
    notcg = np.zeros(ninits)
    penloglikset = np.zeros(ninits)
    loglikset = np.zeros(ninits)
    
    for jn  in range(ninits):        
        
        alpha = initial_pars['alpha'][:,jn]
        sigma = initial_pars['sigma'][:,jn]
        
        if p > 0:
            gamma = initial_pars['gam'][:,jn]
        
        mubeta = [initial_pars['mubeta'][mm][:,jn] for mm in range(m)] 
        tau = [initial_pars['tau'][mm][:,jn] for mm in range(m)] 
        
        oldpenloglik = -np.inf
        emit = 0
        diff = 1.0
        sing = 0
        
        
        for iter_ii in range(maxit):
            
            # ll = -nt * M_LN_SQRT_2PI
            if p == 0:
                ytilde = y
            else:
                ytilde = y - z @ gamma
            # compute residual
            for mm in range(m):
                mubeta_m = mubeta[mm].reshape((q1,k))
                for kk in range(k):
                    res[mm][kk] = (ytilde - x1 @ mubeta_m[:,kk])
                    # r_t.reshape((n,t)).T.sum(axis=0)
                l_m_k[mm] = log_likelihood_array(res[mm], 0, sigma[mm])
                min_l_m_k = l_m_k[mm].min(axis=0)
                # update weight conditional on type mm 
                l_m_k_weighted = tau[mm][:, None] * np.exp(l_m_k[mm] - min_l_m_k)
                # l_m_k_weighted_2 = tau[mm][:, None] * np.exp(l_m_k[mm])
                w_m[mm] = l_m_k_weighted / np.sum(l_m_k_weighted,axis=0)    
                w_m[mm] = np.nan_to_num(w_m[mm], nan=0.5)
            
                # compute f_i conditional on alpha
                l_m[mm] = (np.log(np.sum(l_m_k_weighted,axis=0)) + min_l_m_k ).reshape((n,t)).T.sum(axis=0)
                # (np.log(np.sum(l_m_k_weighted_2,axis=0)) ).reshape((n,t)).T.sum(axis=0)
            
            min_l_m = l_m.min(axis=0)
            l_m_weighted = alpha[:,None] * np.exp(l_m - min_l_m)
            
            # update weight
            w = l_m_weighted / np.sum(l_m_weighted,axis=0)
            w = np.nan_to_num(w, nan=0.5)
            penloglik = (np.log(np.sum(l_m_weighted,axis=0)) + min_l_m).sum()
            
            for mm in range(m):
                s0j = sigma_0[mm] / sigma[mm]
                penloglik += -an * (s0j**2 - 2.0 * np.log(s0j) - 1.0)
            diff = penloglik - oldpenloglik
            oldpenloglik = penloglik

            if np.max(np.abs(diff)) < tol:
                break

            # update parameters
            for mm in range(m):
                
                mubeta_m = np.zeros((q1,k))
                alpha[mm] = np.mean(w[mm])
                tau[mm]   = np.mean(w_m[mm],axis=1)
                w_j = np.repeat(w[mm, :], t).flatten()
                tau[mm] = np.clip(tau[mm], 0.01, 0.99)
                tau[mm] /= np.sum(tau[mm])
                
                
                for kk in range(k):
                    w_jk = w_j * w_m[mm][kk]
                    xtilde = w_jk[:,np.newaxis] * x1 
                    design_matrix = xtilde.T @ x1
                    det_A = np.linalg.det(design_matrix)
                    if np.abs(det_A) < 1e-5:
                         sing = 1
                         break
                    mubeta_m[:,kk] = np.linalg.solve(design_matrix, xtilde.T @ ytilde)
                    res[mm][kk] = ytilde - x1 @ mubeta_m[:,kk]
                    
                sigma[mm] = np.sqrt( ((((res[mm] **2 ) * w_m[mm]) * w_j[np.newaxis,:]).sum() + 2.0 * an * sigma_0[mm]**2)  / (w_j[np.newaxis,:].sum() + 2.0 * an))
                sigma[mm] = max(sigma[mm], epsilon * sigma_0[mm])
                
                mubeta[mm] = mubeta_m.flatten()
            alpha = np.clip(alpha, 0.01, 0.99)
            alpha /= np.sum(alpha)
            
            if p > 0:
                zz.fill(0)
                ze.fill(0)
                for mm in range(m):
                    mubeta_m = mubeta[mm].reshape((q1,k))
                    for kk in range(k):
                        w_jk = w_j * w_m[mm][kk]
                        ztilde = w_jk[:,np.newaxis] * z
                        zz += (ztilde.T @ z) / (sigma[mm]**2)
                        ze += (ztilde.T @ (y - x1 @ mubeta_m[:,kk])) / (sigma[mm]**2)
                det_A = np.linalg.det(zz)
                if np.abs(det_A) < 1e-5:
                    sing = 1
                    break
                # if np.linalg.cond(zz) < SINGULAR_EPS or np.isnan(ze).any():
                #     sing = 1
                #     break
                gamma = np.linalg.solve(zz, ze)
                
        # if np.any(np.isnan(alpha)):
        #     print(jn)
            
            
        penloglikset[jn] = penloglik.sum()
        loglikset[jn] = penloglik.sum()
        post[:, jn] = w.T.flatten()
        
        initial_pars['alpha'][:,jn] = alpha
        initial_pars['sigma'][:,jn] = sigma
        
        for mm in range(m):
            initial_pars['mubeta'][mm][:,jn] = mubeta[mm]
            initial_pars['tau'][mm][:,jn] = tau[mm]
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
def regpanelmixPMLE(data, m, k,  ninits=10, epsilon=1e-8, maxit=2000, epsilon_short=1e-2, maxit_short=500):
    # Extract the generated data
    t,n = data["Y"].shape
    nt = n * t
    
    y = data["Y"].T.flatten()
    x = data["X"]
    z = data["Z"]
    
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
        
    npar = m - 1 + (k * q1 + 1) * m + p
    ninits_short = ninits * 10 * npar
    
    if x is not None and z is not None:
        # Case 1: Both x and z are not None
        xz = np.hstack((np.ones((nt,1)), x, z))
        
    elif x is None and z is not None:
        # Case 2: x is None, z is not None
        xz = np.hstack((np.ones((nt,1)),z))
    elif x is not None and z is None:
        # Case 3: x is not None, z is None
        xz = np.hstack((np.ones((nt,1)),x))
    else:
        # Case 4: Both x and z are None
        xz = np.ones((nt,1))

    ls_out = np.linalg.lstsq(xz, y, rcond=None)
    sd0 = np.sqrt(ls_out[1]/nt)
    
    if (m == 1) & (k == 1):
        mubeta = ls_out[0][:q1]
        if z is not None:
            gam = ls_out[0][q1:(q1 + p)]
        else:
            gam = None
        res = y - xz @ ls_out[0]
        sigma = np.sqrt(np.mean(res**2))
        loglik = log_likelihood_normal(res,0,sigma)
        # loglik = - (nt / 2) * (1 + np.log(2 * np.pi) + 2 * np.log(sigma))
        aic = -2 * loglik + 2 * npar
        bic = -2 * loglik + np.log(n) * npar
        penloglik = loglik
        coefficients = {"alpha": 1, "mubeta": mubeta, "sigma": sigma, "gam": gam}
    else: #either m >=2 or k >= 2
        
        initial_pars = {}
        # Generate initial guesses
        mubeta0 = ls_out[0][:q1]
        if p > 0:
            gam0 = ls_out[0][q1:(q1 + p)]
        
        # Generate initial parameters
        initial_pars['alpha'] = np.random.uniform(size=(m, ninits_short))
        initial_pars['alpha'] = initial_pars['alpha'] / initial_pars['alpha'].sum(axis=0)
        
        initial_pars['sigma'] = np.random.uniform(0.01, 1, size=(m, ninits_short)) * sd0
        
        initial_pars['tau'] = [[] for mm in range(m)]
        for mm in range(m):
             initial_pars['tau'][mm] = np.random.uniform(size=(k, ninits_short))
             initial_pars['tau'][mm] = initial_pars['tau'][mm] / initial_pars['tau'][mm].sum(axis=0)
        if p > 0:
            initial_pars['gam'] = np.random.uniform(0.5, 1.5, size=(p, ninits_short)) * gam0[:, None]
        minMU = np.min(y - xz @ ls_out[0]) + mubeta0[0]
        maxMU = np.max(y - xz @ ls_out[0]) + mubeta0[0]
        initial_pars['mubeta'] = [np.zeros((q1 * k, ninits_short)) for mm in range(m)]
        for mm in range(m):
            mubeta_m = np.zeros((q1 * k, ninits_short)) 
            # mubeta_m[:,0].reshape((q1,k)) This will get you the mubeta for mk-th component
            mubeta_m[:k,] = np.random.uniform(minMU, maxMU, size=(k, ninits_short))
            for kk in range(k):
                for i in range(1, q1):
                    mubeta_m[i * k + kk , :] = mubeta0[i] * np.random.uniform(-2, 2, size=ninits_short)
            initial_pars['mubeta'][mm] = mubeta_m        
        sigma_0 = np.full(m, sd0)
        
        out_short = reg_panelmix_pmle(data,initial_pars,m,k, sigma_0, maxit=1000, tol = 1e-8, epsilon=0.05)

        components = np.argsort(out_short["penloglikset"])[::-1][:ninits]
        if z is not None:
            long_pars = {'alpha': initial_pars['alpha'][:,components], 'mubeta': initial_pars['mubeta'][:,components], 'sigma': initial_pars['sigma'][:,components], 'gam': initial_pars['gam'][:,components]}
        else:
            long_pars = {'alpha': initial_pars['alpha'][:,components], 'mubeta': initial_pars['mubeta'][:,components], 'sigma': initial_pars['sigma'][:,components], 'gam': None}
        
# %%
# Define parameters
alpha = [0.5, 0.5]  # Category probabilities (for M categories)
mu = [[2, 4], [1, 3]]  # Mean values for each subcategory (M x K)
sigma = [1, 2]  # Standard deviation for each category (length M)
tau = [[0.6, 0.4], [0.7, 0.3]]  # Subcategory probabilities for each M (M x K)
beta = [[1, 0.5], [0.2, 0.3]]  # Coefficients for q covariates (M x q)
gam = None  # Coefficients for p covariates (length p)

N = 100  # Number of individuals
T = 5  # Number of time periods
M = 2  # Number of categories
K = 2  # Number of subcategories per category
p = 0  # Number of covariates in z
q = 2  # Number of covariates in x

# Call the function
data = generate_data(alpha, mu, sigma, tau, N, T, M, K, p, q)

# %%
# Input Parameters


regpanelmixPMLE(data, M, K)