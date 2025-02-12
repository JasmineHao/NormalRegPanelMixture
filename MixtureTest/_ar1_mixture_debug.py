

from math import gamma
from readline import read_history_file

# Define the true weights, need to compare whether the computed weight using the true parameters corresponds to the true weights
weight_true = np.array([ each.flatten() for each in R_sub]) * np.repeat(np.repeat(R, m, axis=1), t-1, axis=0)

weight_true_0 = np.array([ each.flatten() for each in R_sub_0]) * np.repeat(R, m, axis=1)


# %%
m = M
k = 2

y = Y
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
xz = np.hstack((x1, z_l, z_c))
q1 = 2*q + 2

xz_0 = np.hstack((np.ones((n, 1)), x_0, z_0) )
    

# %%
# deal with parameters
alpha_jn = alpha
mubeta_jn_mat = np.zeros((2,4))
mubeta_jn_mat[:,-1] = rho 
mubeta_jn_mat[:,1] = beta[:,0]
mubeta_jn_mat[:,2] = - beta[:,0] * rho
 
mubeta_0_jn_mat = np.zeros((2,2))
mubeta_0_jn_mat[:,1] = beta_0[:,0]

mu_jn = mu.flatten()
mu_0_jn = mu_0.flatten()
tau_jn = tau.flatten()
sigma_jn = sigma 
sigma_0_jn = sigma_0
gamma_jn = gamma
gamma_0_jn = gamma_0

# %%
# Given true parameters, get the estimates of weights, check weight magnitude

np.abs(w_mk_0.T - weight_true_0).max()

np.abs(w_mk.T - weight_true).max()

tmp = np.abs(w_mk.T - weight_true)

indices = np.where(tmp > 0.1)

# (array([  38,   38,  537,  537,  700,  700, 1003, 1003, 1444, 1444]),
#  array([2, 3, 2, 3, 2, 3, 2, 3, 2, 3]))