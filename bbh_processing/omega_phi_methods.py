#Here we store various methods of computing omega_phi
import peak_finders as pks
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import frequency_tools as freq
import bbh_processing.utilities as bbh
#Here we fit envelopes between the peaks and troughs of the frequency curve
#and take omega_phi to be the average of these ('method 1' in the notebook). 
def omega_phi_from_envelopes(OmegaMag):
  errcode, left,right,peaks,troughs = pks.monotonic_peaks2(OmegaMag[:,0],
      OmegaMag[:,1])

  if errcode==1:
    raise ValueError("Not enough peaks to fit envelopes.")

  #we now require the envelopes interpolated to the OmegaMag[:,0] time series
  topenv = peaks(OmegaMag[:,0])
  botenv = troughs(OmegaMag[:,0])
  omega_phi = 0.5*(topenv+botenv)
  omphi_array = np.zeros(OmegaMag.shape)
  omphi_array[:,0] = np.copy(OmegaMag[:,0])
  omphi_array[:,1] = np.copy(omega_phi)
  return omphi_array

#Take peak-to-trough windows of the function and compute their averages.
#width is the number of half peaks in each window; i.e. for width=1 
#each window is half an orbit. Data are assigned to the midpoint of the window,
#then interpolated back to the original time-series.
def omega_phi_from_time_averages(OmegaMag,width=2,stride=1,append_time=True):
  turns, vals = pks.turning_points(OmegaMag[:,1]) 
  print width
  peakrng = zip(turns[:-width:stride],turns[width::stride])
  
  om_phi_list  = [np.average(OmegaMag[a:b,1]) for a,b in peakrng]
  t = [OmegaMag[a,0]+0.5*(OmegaMag[b,0]-OmegaMag[a,0]) for a,b in peakrng]

   
  #k=3 -> cubic spline
  #s=0 -> smoothing factor = 0; i.e. demand the spline passes through all
  #points in the list (at the possible expense of its smoothness)
  omphispline = interpolate.UnivariateSpline(t,om_phi_list,k=3,s=0) 
  omphi=np.array(omphispline(OmegaMag[:,0])) 
  if append_time: 
    omphi_array = np.copy(OmegaMag) 
    omphi_array[:,1] = omphi
    return omphi_array
  else:
    return omphi
  #plt.plot(OmegaMag[:,0],omr)
  #plt.plot(OmegaMag[:,0],OmegaMag[:,1])
  #plt.show()

def omega_phi_expanding_mean(OmegaMag,append_time=True):
  import pandas as pds
  #turns, vals = pks.turning_points(OmegaMag[:,1]) 
  #turnrng = zip(turns[:-width:width],turns[width::width])
  om_phi = pds.expanding_mean(OmegaMag[:,1])
  #om_phi_list  = [np.average(OmegaMag[a:b,1]) for a,b in turnrng]
  #t = [OmegaMag[a,0]+0.5*(OmegaMag[b,0]-OmegaMag[a,0]) for a,b in turnrng]

  #k=3 -> cubic spline
  #s=0 -> smoothing factor = 0; i.e. demand the spline passes through all
  #points in the list (at the possible expense of its smoothness)
  #omphispline = interpolate.UnivariateSpline(t,om_phi_list,k=3,s=0) 
  #omphi=np.array(omphispline(OmegaMag[:,0])) 
  if append_time: 
    omphi_array = np.copy(OmegaMag) 
    omphi_array[:,1] = om_phi
    return omphi_array
  else:
    return omphi
  #plt.plot(OmegaMag[:,0],omr)
  #plt.plot(OmegaMag[:,0],OmegaMag[:,1])
  #plt.show()


def extrapolated_time_averages(OmegaMag,minwidth=2,maxwidth=10,
    compute_variances=False):
  return bbh.extrapolated_windowed_series(omega_phi_from_time_averages,
      OmegaMag,minwidth,maxwidth,stride=2)
  return means, stddevs


#Compute the time averaged omega_phi from window sizes from minwidth to
#maxwidth, then return the mean of each at each timestep. Also return
#the standard deviation (for error bars).
def mean_time_averages(OmegaMag,minwidth=2,maxwidth=10,
    compute_variances=False):
  return bbh.mean_windowed_series(omega_phi_from_time_averages,
      OmegaMag,minwidth,maxwidth,2)
  return means, stddevs


def omega_phi_from_phi(t,phi,minwidth=2,maxwidth=8):
  phi_arr = np.zeros((len(t),2))
  phi_arr[:,0] = t
  phi_arr[:,1] = phi 
  return bbh.mean_windowed_series(freq.freq_from_ranges,phi_arr,
      minwidth,maxwidth,2)

def omega_phi_from_phi_fits(t,phi,minwidth=2,maxwidth=8):
  phi_arr = np.zeros((len(t),2))
  phi_arr[:,0] = t
  phi_arr[:,1] = phi 
  return bbh.mean_windowed_series(freq.freq_from_fits,phi_arr,
      minwidth,maxwidth,2)
#Do mean_time_averages for a range of minimum and maximum window sizes.
#At each point, return the result with the minimum standard deviation.
#def minimum_stddev_average(OmegaMag,lowmin=2,highmin=6,
#    highmax=12,increment=2):
#  
#  beststdev=np.ones(len(OmegaMag[:,0]))
#  bestmean = np.zeros(len(OmegaMag[:,0]))
#  for thismin in range(lowmin,highmin,increment):
#    for thismax in range(thismin+4,highmax,increment):
#      print (thismin,thismax)
#      means,stddevs = mean_time_averages(OmegaMag,thismin,thismax,
#            compute_variances=False)
#      for i in len(means):
#        if stddevs[i]<beststddev:
#          beststdev[i]=stddevs[i]
#          bestmean[i] =means[i]
#  return bestmean,beststdev  

