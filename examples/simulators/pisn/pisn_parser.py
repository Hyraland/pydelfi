import numpy as np
from astropy.cosmology import Planck15, z_at_value
import h5py
import astropy.units as u

def b14_parse(z_min=None, z_max=None, qual_cut=False, \
        pisn_path='/Users/hyraland/Projects/GW/Cosmology/PISNLineCosmography/'):

    filename = 'm90.txt'
    datar = np.loadtxt(pisn_path + filename)

    # keys: etaobs,m1s,m2s,mcobs,posterior,rhoobs,sigma_eta,sigma_mc,sigma_rho,sigma_t,thetaobs,thetas
    dl = datar[0]
    m90 = datar[1]
    m90std = datar[2]
    dlstd = datar[3]

    # take m1 as data for now
    data = {'dl':dl, 'm90':m90,'m90std':m90std,'dlstd':dlstd} 
    n_sn = len(data['dl'])

    cmats = {'dl': np.zeros((n_sn, n_sn))}
    for cmat in cmats:
        d = data['dlstd']
        for i in range(n_sn):
            cmats[cmat][i, i] = d[i]
    return data, cmats

# def b14_covariance(data, cmats, alpha, beta):
#     n_sn = len(data)
#     c_mat = cmats['v0'] + alpha ** 2 * cmats['va'] + \
#             beta ** 2 * cmats['vb'] + 2 * alpha * cmats['v0a'] - \
#             2 * beta * cmats['v0b'] - \
#             2 * alpha * beta * cmats['vab']
#     d_mat = data['dmb'] ** 2 + (alpha * data['dx1']) ** 2 + \
#             (beta * data['dcolor']) ** 2 + \
#             2 * alpha * data['cov_m_s'] - \
#             2 * beta * data['cov_m_c'] - \
#             2 * alpha * beta * data['cov_s_c']
#     return c_mat + np.diag(d_mat)

# def b14_covariance_derivative(data, cmats, alpha, beta):

#     n_sn = len(data)

#     dC = np.zeros((2, n_sn, n_sn))

#     # Derivatives w.r.t. alpha
#     c_mat = 2*alpha*cmats['va'] + 2*cmats['v0a'] - 2*beta*cmats['vab']
#     d_mat = 2*alpha*data['dx1']**2 + 2*data['cov_m_s'] - 2*beta*data['cov_s_c']
#     dC[0,:,:] = c_mat + np.diag(d_mat)

#     # Derivatives w.r.t. beta
#     c_mat = 2*beta*cmats['vb'] - 2*cmats['v0b'] - 2*alpha*cmats['vab']
#     d_mat = 2*beta*data['dcolor']**2 - 2*data['cov_m_c'] - 2*alpha*data['cov_s_c']
#     dC[1,:,:] = c_mat + np.diag(d_mat)

#     return dC

# def b14_chi_sq(pars, data, cmats, h_0 = 70.0, delta_m_cut = 10.0):
#     alpha, beta, abs_mag, delta_m, om_m = pars
#     mu = z2mu(data['zcmb'], om_m, h_0)
#     res = data['mb'] - (abs_mag - alpha * data['x1'] + \
#                         beta * data['color'] + \
#                         delta_m * (data['3rdvar'] > delta_m_cut)) - mu
#     cov_mat = b14_cov_mat(data, cmats, alpha, beta)
#     cov_mat_chol = np.linalg.cholesky(cov_mat)
#     white_res = np.linalg.solve(cov_mat_chol, res)
#     chi_sq = np.dot(white_res, white_res)
#     return chi_sq

# def emcee_b14_ln_p(pars, data, cmats, h_0 = 70.0, delta_m_cut = 10.0):
#     if 0.0 < pars[0] < 0.5 and 0.0 < pars[1] < 6.0 and -25.0 < pars[2] < -15.0 and -0.5 < pars[3] < 0.5 and 0.2 < pars[4] < 0.4:
#         return -b14_chi_sq(pars, data, cmats, h_0, delta_m_cut) / 2.0
#     else:
#         return -np.inf


