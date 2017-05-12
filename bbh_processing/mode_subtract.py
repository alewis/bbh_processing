import numpy as np
import scipy.optimize as opt
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

def sinfunc(time, param):
    A, om, B, C = param
    return A * np.sin(om * time + B) + C

def errfunc(params, time, phi):
    return phi-sinfunc(time, params)


def make_guesses(tau, dT, phi):
    A = np.sqrt(np.sum(phi**2.)/(dT*np.sqrt(2.))) 
    om = (phi[-1] - phi[0])/dT 
    B = 0
    C = (phi[-1] - phi[0])/2
    pvec = [A, om, B, C]
    return pvec

def omphi_fit(time, phi, tpks, Nsamples=200, maxiter=1):
    if len(time) != len(phi):
        raise ValueError("Inconsistent sizes in time, phi")

    #Interpolate phi onto the radial peaks  
    phispl = UnivariateSpline(time, phi, k=3, s=0)
   
    #Fit a least squares line to phi between each radial peak 
    tmids = []
    omega_phis = []
    residuals = []
    t_fits = []
    params = []
    fits=[]
    for t0, tf in zip(tpks[:-1], tpks[1:]):
        tau = np.linspace(0., tf-t0, Nsamples)
        dT = tf-t0
        tmid = t0 + dT/2.
        tmids.append(tmid)
        t_fits.append(tau+t0)
        phi_in_range = phispl(tau + t0)
        pvec = make_guesses(tau, dT, phi_in_range)
        phi_to_fit = np.copy(phi_in_range)
        for i in range(0, maxiter):
            argtup = (tau, phi_to_fit)
            pnew, pcov = opt.leastsq(errfunc, pvec, args=argtup)      
            fit = sinfunc(tau, pnew) 
            phi_to_fit -= fit
        fits.append(fit)
    return np.array(t_fits), np.array(fits)
