import numpy as np
import scipy.optimize as opt
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

# def phifunc(time, phi0, omega_phi, a, b, alpha_r, alpha_th, omega_r, omega_th):
    # phi = phi0 + omega_phi*time
    # phi += a*np.sin(omega_r*time + alpha_r)
    # phi += b*np.sin(omega_th*time + alpha_th)
    # return phi

def phifunc(time, param, omr, omth, maxk):
    phi0, omphi, a, b = param[:4]
    phases = param[4:]
    phi = np.zeros(len(time))
    phi[:] = phi0
    for k, in zip(range(1, maxk+1)):
        phase_r = phases[(k-1)]
        phase_th = phases[k]
        phi += omphi*time
        phi += a*np.sin(k*omr*time + phase_r)
        phi += b*np.sin(k*omth*time + phase_th)
    return phi

def errfunc(params, time, phi, omr, omth, maxk):
    return phi - phifunc(time, params, omr, omth, maxk)


def make_guesses(tau, dT, phi_in_range, maxk):
    omega_0i = (phi_in_range[-1] - phi_in_range[0])/dT
    phi_0i = phi_in_range[0]
    line_vals = tau*omega_0i + phi_0i
    phi_to_fit = phi_in_range - line_vals
    
    guess_a = np.sqrt(np.sum(phi_to_fit**2.)/(dT*np.sqrt(2.))) 
    guess_b = guess_a
    pvec = [phi_0i, omega_0i, guess_a, guess_b]
    for k in range(0, maxk):
        pvec.append(0)
        pvec.append(1)
    # guess_alpha_r = 0. 
    # guess_alpha_th = 1. 
    return pvec, phi_to_fit


def do_fit(pvec, maxk, tau, phi_to_fit, omega_r, omega_th):
    argtup = (tau, phi_to_fit, omega_r, omega_th, maxk)
    pvec, pcov = opt.leastsq(errfunc, pvec, args=argtup)
    #phi_to_fit -= phifunc(time, pvec, omr, omth, maxk) 
    return pvec, pcov

def omphi_fit(time, phi, omega_r, omega_th, tpks, Nsamples=200, maxk=1):
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
    for t0, tf in zip(tpks[:-1], tpks[1:]):
        tau = np.linspace(0., tf-t0, Nsamples)
        dT = tf-t0
        tmid = t0 + dT/2.
        tmids.append(tmid)
        t_fits.append(tau+t0)
        phi_in_range = phispl(tau + t0)
        pvec, phi_to_fit = make_guesses(tau, dT, phi_in_range, maxk)
        omega_phis.append(pvec[1])
        pvec, pcov = do_fit(pvec, maxk, tau, phi_to_fit, omega_r, omega_th)      

        fitval = phifunc(tau, pvec, omega_r, omega_th, maxk)
        # plt.plot(tau+t0, phi_to_fit, label="data")
        # plt.plot(tau+t0, fitval, label="fit")
        resid = phi_to_fit - fitval 
        residuals.append(resid)
        params.append(pvec)
            
        answer = pvec[0] + pvec[1] * tau
        #plt.plot(tau+t0, answer, color="black", label="final_omphi")
        #plt.legend(loc="upper left")
        #plt.show()
    omph_arr = np.transpose(np.array([tmids, omega_phis]))
    params_arr = np.array(params)
    return [t_fits, omph_arr, residuals, params_arr]

