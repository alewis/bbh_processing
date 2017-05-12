import numpy as np
import cPickle as pickle
from scipy.interpolate import UnivariateSpline
import bbh_processing.hfile_tools as hfile
import bbh_processing.peak_finders as peak
import bbh_processing.equatorial_freq_methods as eqfreq
import bbh_processing.unwrap as unwrap
import matplotlib.pyplot as plt

def krph(time, data, N, Niter=1, trim=True, minidx=200):
    ph_all = omega_phi_sliding(time, data, N, Niter=Niter, trim=trim,
                               minidx=minidx, joined=True)
    r_all = omega_r_sliding(time, data, N, Niter=Niter, trim=trim,
                            minidx=minidx, joined=True)
    k_all = r_all[:, 1] / ph_all[:, 1] 
    return k_all

def ks_equatorial(omtime, omdata, N, minidx=200):

    omph = eqfreq.omega_phi_sliding(omtime, omdata, N, minidx=minidx)
    omr = eqfreq.omega_r_sliding(omtime, omdata, N, minidx=minidx)
    pt, pv = peak.conservative_peaks(omtime, omdata, minidx=minidx)
    tt, tv = peak.conservative_peaks(omtime, omdata, minidx=minidx, pktype="trough")
    ecc = eqfreq.ecc_sqrt(pt, pv, tt, tv, omph[:, 0])
    krph = omr[:, 1] / omph[:, 1]
    outdict = {
        "time": omph[:, 0], "ecc": ecc,
        "omph": omph[:, 1], "omr": omr[:, 1], 
        "krph": krph}
    return outdict

def ks_inclined(omtime, omdata, thtime, thdata, N, minidx=200):
    pt, pv = peak.conservative_peaks(omtime, omdata, minidx=minidx)
    tt, tv = peak.conservative_peaks(omtime, omdata, minidx=minidx, pktype="trough")
    omph = eqfreq.omega_phi_from_peaks(omtime, omdata, pt, pv, N)
    omr = eqfreq.omega_r_from_peaks(omtime, omdata, pt, N)
    omth = unwrap.omega_theta_from_peaks(thtime, thdata, pt, N=1)
    ecc = eqfreq.ecc_sqrt(pt, pv, tt, tv, omph[:, 0])
    krph = omr[:, 1] / omph[:, 1]
    krth = omr[:, 1] / omth[:, 1]
    kthph = omth[:, 1] / omph[:, 1]
    
    outdict = {
        "time": omph[:, 0], "ecc": ecc,
        "omph": omph[:, 1], "omr": omr[:, 1], "omth": omth[:, 1],
        "krph": krph, "krth": krth, "kthph": kthph}
    return outdict
