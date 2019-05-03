import numpy as np
import scipy.integrate as integrate
from .moped import *
import h5py
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u

# luminosity distance
def dlofz(theta, auxiliary_data):

    # Cosmological parameters
    h0 = theta[0]
    omegam = theta[1]
    om = theta[2]

    # Pull out the relevant things from the data
    z = auxiliary_data[:,0]

    c =  3*10**5
    omegal = 1- omegam
    tau = lambda x, omegalf, omegamf, omf : 1/(np.sqrt(omegalf * (1+x)**(3*(1+omf)) + omegamf * (1+x)**3))
    dl = []
    for i in z:
        inte, err = integrate.quad(tau, 0, i, args=(omegal, omegam, om, ))
        dl.append((1+i) * c/h0/1000 * inte)
        #dl.append((1+i) * c * inte)
    return np.array(dl)

# Generate realisation of \mu
def simulation(theta, sim_args):
    
    # Pull out data
    auxiliary_data = sim_args[0]
    L = sim_args[1]
    
    # Signal
    dl = dlofz(theta, auxiliary_data)
        
    # Noise
    noise = np.dot(L, np.random.normal(0, 1, len(L)))
    
    # Return signal + noise
    return dl + noise

# Generate realisation of \mu
def simulation_seeded(theta, seed, sim_args):
    
    # Pull out data
    auxiliary_data = sim_args[0]
    L = sim_args[1]
    
    # Signal
    dl = dlofz(theta, auxiliary_data)
        
    # Noise
    np.random.seed(seed)
    noise = np.dot(L, np.random.normal(0, 1, len(L)))
    
    # Return signal + noise
    return dl + noise

def compressor(data, args):

    theta_fiducial, Finv, Cinv, dmdt, dCdt, mu, Qinv, prior_mean = args
    
    return mle(theta_fiducial, Finv, Cinv, dmdt, dCdt, mu, Qinv, prior_mean, data)

def compressor_projected(data, args):
    
    theta_fiducial, Finv, Cinv, dmdt, dCdt, mu, Qinv, prior_mean, F, P1, P2 = args

    # MOPED compress the data
    d_twidle = mle(theta_fiducial, Finv, Cinv, dmdt, dCdt, mu, Qinv, prior_mean, data)
    
    # Now do the projection
    d_twidle = np.dot(F, d_twidle - theta_fiducial - np.dot(Finv, np.dot(Qinv, prior_mean - theta_fiducial)))
    d_twidle = np.dot(Finv[0:2, 0:2], np.array([d_twidle[0] - np.dot(P1, d_twidle[2:]), d_twidle[1] - np.dot(P2, d_twidle[2:])]))
    d_twidle = d_twidle + theta_fiducial[:2] + np.dot(Finv[:2,:2], np.dot(Qinv[:2,:2], prior_mean[:2] - theta_fiducial[:2]))

    return d_twidle


