import peak_finders as pks
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import fit_sin_to_theta as fitsin
#compute dot(omega_phi)/omega_phi^2
def adiabaticity(delta_t,omega):
  omega_dot= np.fabs(np.gradient(omega,delta_t))
  omega_squared = np.square(omega)
  return np.divide(omega_dot,omega_squared)

#Define omega as 1/ the time between successive peaks
def freq_from_ranges(array, width=2):
  
  time = array[:,0] 
  data = array[:,1]
  turns, vals = pks.turning_points(data) 
  turnrng = zip(turns[:-width:width],turns[width::width])
  tlist = []
  omlist = []
  ncycles = width/2.0
  for a,b in turnrng:
    t = time[a]+0.5*(time[b]-time[a])
    f = ncycles/(time[b]-time[a]) 
    omlist.append(f)
    tlist.append(t)
  #now interpolate
  omspline = interpolate.UnivariateSpline(tlist,omlist,k=3,s=0) 
  om=omspline(time) 
  om_array = np.zeros((len(time),2)) 
  om_array[:,0] = np.copy(time)
  om_array[:,1] = 2.*np.pi*om
  #plt.plot(om_array[:,0],om_array[:,1])
  #plt.show()
  return om_array

#fit sinusoids and return their frequencies
def freq_from_fits(time, theta, width=2):
  rngs, ps,midtimes = fitsin.fit_sin_from_guessed_params(
      time, theta,
      width=width)
  freqs = [p[1] for p in ps]
  #f = interp1d(om_t_times,om_t,kind='cubic')
  omspline = interpolate.UnivariateSpline(midtimes,freqs,k=3,s=0) 
  omega = np.zeros((th_arr.shape))
  omega[:,0] = np.copy(th_arr[:,0])
  omega[:,1] = 2*np.pi*omspline(th_arr[:,0])
  return omega 

def peak_to_peak_fft(array,width=6):
  time = array[:,0] 
  data = array[:,1]
  a = pks.turning_points(data) 
  print a
  turnrng = zip(turns[:-width:width],turns[width::width])
  tlist = []
  omlist = []
  ncycles = width/2.0
  
  omfreqs=[]
  for a,b in turnrng:
    #window = np.hanning(len(time[a:b]))
    #windowed = data[a:b]*window
    ffreqs = 2*np.pi*np.fft.rfftfreq(len(time[a:b]),0.5)
    omfreqs.append(ffreqs[1])
    #hs = np.abs(np.fft.rfft(windowed))
    #print ffreqs[1]
    #plt.plot(time[a:b],data[a:b],label="data")
    #plt.plot(time[a:b],windowed, label="windowed")
    #plt.legend()
    #plt.plot(time[a:b],window)
    #plt.show()
    #plt.plot(ffreqs,hs)
    #plt.xlim((0,0.06))
    #plt.ylim((0,5))
    #plt.show()
    
    t = time[a]+0.5*(time[b]-time[a])
    #f = ncycles/(time[b]-time[a]) 
    #omlist.append(f)
    tlist.append(t)
  #now interpolate
  omspline = interpolate.UnivariateSpline(tlist,omfreqs,k=3,s=0) 
  om=omspline(time) 
  om_array = np.zeros((len(time),2)) 
  om_array[:,0] = np.copy(time)
  om_array[:,1] = om
  #plt.plot(om_array[:,0],om_array[:,1])
  #plt.show()
  return om_array

  
