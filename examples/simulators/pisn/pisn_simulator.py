import numpy as np
import scipy.integrate as integrate
import simulators.pisn.pisn_parser as pisn_parser
from astropy.cosmology import Planck15
import h5py
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u

class PISN_Model():

    def __init__(self, pisn_data_path = '/Users/hyraland/Projects/GW/Cosmology/PISNLineCosmography/'):

        # Import data
        pisn_data, pisn_cmats = pisn_parser.b14_parse(z_min=None, z_max=None, qual_cut=False,
                                           pisn_path='/Users/hyraland/Projects/GW/Cosmology/PISNLineCosmography/')
        
        self.data = pisn_data['dl'] # pisn_data['dl_z'] # np.mean(pisn_data['dl'],axis = 1)
        delta_m_cut = 10
        self.auxiliary_data = np.column_stack([pisn_data['m90']])

        # h0, omegam, w
        self.npar = 4
        self.theta_fiducial = np.array([Planck15.H0.value, Planck15.Om0, -1.0, 52.0])

        # Covariance matrix
        self.C = pisn_cmats['dl']
        self.Cinv = np.linalg.inv(self.C)
        self.L = np.linalg.cholesky(self.C)

        # Derivative of the covariance matrix
        self.n_sn = len(self.C)
        self.dCdt = np.zeros((self.npar, self.n_sn, self.n_sn))

        # N data points
        self.ndata = len(pisn_data['dl'])

        # Compute the mean
        self.dl = self.dlofz(self.theta_fiducial)

    # luminosity distance
    def dlofz(self, theta):
        
        # Cosmological parameters
        h0 = theta[0]
        omegam = theta[1]
        om = theta[2]
        m90_0 = theta[3]
        
        # Pull out the relevant things from the data
        z = self.auxiliary_data[:,0]/m90_0-1
        
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
    def simulation(self, theta, seed):
        
        # Set the seed
        np.random.seed(seed)

        # Signal
        dls = self.dlofz(theta)
        
        # Noise
        noise = np.dot(self.L, np.random.normal(0, 1, len(self.L)))
        
        # Return signal + noise
        return dls + noise

    # Generate derivative of \mu w.r.t cosmological parameters
    def dmudt(self, theta_fiducial, h):
        
        # dmdt
        dmdt = np.zeros((self.npar, self.ndata))
        
        # Fiducial data
        d_fiducial = self.dlofz(theta_fiducial)
        
        # Loop over parameters
        for i in range(0, 4):
            
            # Step theta
            theta = np.copy(self.theta_fiducial)
            theta[i] += h[i]
            
            # Shifted data with same seed
            d_dash = self.dlofz(theta)
            
            # One step derivative
            dmdt[i,:] = (d_dash - d_fiducial)/h[i]
        
        return dmdt





