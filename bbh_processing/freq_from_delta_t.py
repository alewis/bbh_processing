#freq_from_delta_t.py
#Input is a set of time coordinates t_p, usually representing extremum of some 
#time series, along with a window size N counting the number of t_p's in the
#window.
#Output is a time series of instantaneous frequencies vs. time. These are 
#computed at each extremal and half-extremal time at least N extremum
import numpy as np
from scipy import interpolate
#Construct the mean anomaly: a monotonically-increasing parameter whose value
#increases by pi at every extremum. Input is the locations of the extrema in
#time. For Keplerian ellipses this is just a linear rescaling of time by the
#orbital frequency Omega_phi: M = 2pi*(t-t_o)*Omega_phi.
def mean_anomaly(tpk,t):
  Mpk = np.pi*np.arange(0,tpk.shape[0])
  Mspl = interpolate.UnivariateSpline(tpk,Mpk,s=0,k=1)
  return Mspl(t)
  


