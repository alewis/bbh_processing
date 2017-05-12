import hfile_tools as hf
import peak_finders as pks
import cPickle as pickle
from scipy import interpolate 
import scipy.optimize as opt 
import math 
import numpy as np
import hfile_tools as hf
import matplotlib.pyplot as plt
import fit_sin_to_theta as fitsin
from scipy import signal
import bbh_processing.utilities as bbh
import frequency_tools as freq

global MEAN
global AMP
global PHASE

def get_sin_guesses(data,time,a,b,width):
  #A0 =2*np.std(data[a:b])/(2**0.5)
  interval = time[b]-time[a]
  period = interval / (2.0*width)
  #B0 = 1.0/period
  #period = time[a](time[b]-time[a]) / (2.0*width)
  B0 = (2.0*np.pi)/period
  
  #the ranges are always peak-to-peak, so the phase is either pi/2 or 3pi/2 
  C0 = 3.0*np.pi/2.0
  #mean = np.mean(data)
  return B0,C0

def sinusoid(p, x):
  return p[0]*np.sin(p[1]*x+p[2])+p[3]

#shift by 3pi/2 since we go from peak to peak
def freq_sin(p,x):
  return np.sin(p[0]*x+3.0*np.pi/2.0)

def fixed_amp_mean_sinusoid(p,x):
  return np.sin(p[0]*x+p[1])

def shift_data(data,time):
  amp, mean = get_amp_mean(data,time)
  shifted_data = np.fabs(data) - mean
  shifted_data = shifted_data / amp
  return time,shifted_data

def linchirp(p,x):
  return np.sin(p[2]*x*x+p[0]*x+p[1])

def get_amp_mean(data,time):
  mean = np.mean(data)

  #to get the amplitude, we first get the peaks.
  #then at first approximation we consider the amplitude to be the average
  #distance between peak and mean.
  turns, vals = pks.turning_points(data) 
  distances = np.fabs(mean-vals)
  amp = np.mean(distances)
  return amp,mean 

def fit_to_sinusoid(this_range, time, data,width):
  a = this_range[0]
  b = this_range[1]
  #print range
  mid = a+(b-a)/2
  midtime = time[mid]
  
  shifted_data,time = shift_data(data,time) 
  errfunc = lambda p, x, y: freq_sin(p,x) - y
  
  #errfunc = lambda p, x, y: sinusoid(p,x) - y
  #y = A * sin(t*B + C) + D
  #print "A0: ", A0, "    B0: " , B0, "     C0: ", C0, "    D0: " , D0
  B0,C0 = get_sin_guesses(shifted_data,time,a,b,width) 
  p0 = [B0]
  p1, success = opt.leastsq(errfunc, p0[:], args=(time[a:b],shifted_data[a:b]))
  return p1,midtime

def initial_fit_to_linear_chirp(this_range, time, data,width):
  a = this_range[0]
  b = this_range[1]
  #print range
  errfunc = lambda p, x, y: linchirp(p,x) - y
  
  
  #y = A * sin(t*B + C) + D
  #print "A0: ", A0, "    B0: " , B0, "     C0: ", C0, "    D0: " , D0
  amp,mean = set_mean_and_amp(data,time) 
  B0,C0 = get_sin_guesses(data,time,a,b,width) 
  p0 = [amp,B0,C0,mean, 0.0]
  p2, success = opt.leastsq(errfunc, p0[:], args=(time[a:b],data[a:b]))
  return p2
def fit_to_linear_chirp(this_range, time, data,p1):
  a = this_range[0]
  b = this_range[1]
  #print range
  errfunc = lambda p, x, y: linchirp(p,x) - y
  
  
  #y = A * sin(t*B + C) + D
  #print "A0: ", A0, "    B0: " , B0, "     C0: ", C0, "    D0: " , D0
  p0 = [p1[0], p1[1], p1[2], p1[3], 0.0,0.0]
  p2, success = opt.leastsq(errfunc, p0[:], args=(time[a:b],data[a:b]))
  return p2


#Does a polynomial fit of order 'order' to t,y
def do_fit(t,y,order):
  F = np.zeros( len(t) ) 
  if 0==order : 
    F[:] = np.mean( y )
  elif 1==order : 
    from scipy.stats import linregress
    slope, intercept, r_val, p_val, std_err = linregress(t,y)
    print "Offset computed using linear regression."
    F[:] = slope * t[:] + intercept
  elif order > 1:
    z = np.polyfit(t,y,order)
    p = np.poly1d(z)
    for i in range(0,len(F)):
      F[i] = p( t[i] )
  else:
    raise ValueError("Unanticipated value of order")
  return F

#Just a wrapper around freq_from_splines 
def omega_theta_from_splines(
    hname="Horizons.h5", 
    offset_order = 1,  #
    amp_order = 3,
    amp_smoothness = 0,
    signal_order = 3,
    signal_smoothness = 0):
  time = hf.t(hname)
  theta = hf.theta_r_chi(hname)
  data = np.zeros((len(time),2))
  data[:,0] = np.copy(time)
  data[:,1] = np.copy(theta)
  return freq_from_splines(data,
      offset_order=offset_order,
      amp_order = amp_order,
      amp_smoothness = amp_smoothness,
      signal_order = signal_order,
      signal_smoothness = signal_smoothness
      )

def omega_theta_from_splines_theta(
    time,
    theta,
    offset_order = 1,  #
    amp_order = 3,
    amp_smoothness = 0,
    signal_order = 3,
    signal_smoothness = 0):
  data = np.zeros((len(time),2))
  data[:,0] = np.copy(time)
  data[:,1] = np.copy(theta)
  return freq_from_splines(data,
      offset_order=offset_order,
      amp_order = amp_order,
      amp_smoothness = amp_smoothness,
      signal_order = signal_order,
      signal_smoothness = signal_smoothness
      )

#Frequency of the theta_r_chi as 1/the time between successive peaks.
def omega_theta_from_peaks(width=2,hname="Horizons.h5",pname="Omega.pkl"):
  time = hf.t(hname)
  theta = hf.theta_r_chi(hname)
  return freq(time,theta,width=width)

def omega_theta_from_peaks_avgd(t,theta,minwidth=2,maxwidth=8):
  th_arr = np.zeros((len(t),2))
  th_arr[:,0] = t
  th_arr[:,1] = theta
  return bbh.mean_windowed_series(freq.freq_from_ranges,th_arr,
      minwidth,maxwidth,2)


def omega_theta_from_fits(width=4,hname="Horizons.h5",pname="Omega.pkl"):
  #Frequency of the theta_r_chi as a function of time.
  time = hf.t(hname)
  theta = hf.theta_r_chi(hname)
  return freq.freq_from_fits(t,theta)

def omega_theta_from_fits_theta(time, theta, width=2):
  #Frequency of the theta_r_chi as a function of time.
  return freq.freq_from_fits(time, theta)


def omega_theta_from_fits_avgd(t,theta,minwidth=2,maxwidth=8):
  th_arr = np.zeros((len(t),2))
  th_arr[:,0]=t
  th_arr[:,1]=theta
  return bbh.mean_windowed_series(freq.freq_from_fits,th_arr,
      minwidth,maxwidth,2)
