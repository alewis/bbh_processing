import numpy as np
import HorizonsTools 

class TimeSeries:
  def __init__(self,time_series):
    self.ts = time_series

class ScrubRadialDecay:
  def __init__(self,time_series):
    t = time_series[0,:]
    r = time_series[1,:]
    decay = HorizonsTools.MonotonicFit(r,t)[1,:]
    rscrub = r - decay
    self.monotonic = np.vstack((t,rscrub))

class Plotter(TimeSeries):
  def plt(self):
    import matplotlib.pyplot as matplot
    matplot.plot(self.ts[0,:],self.ts[1,:])

def unstack(time_series):
  return time_series[0,:], time_series[1,:]


##############################################################################
#WINDOW MAKERS
##############################################################################
class WindowMaker:
  def make(self,time_series):
    return -1
  
  #def __init__(self, time_series):
  #  self.windows = turning_points(time_series)
  #  self.n = len(windows)

class ExtremaWindows(WindowMaker):
  #indices of the time_series' extrema
  def peaks(self,time_series,peaks_in_range):
    turn_ids=[]
    for pt in range(1,time_series.shape[1]-1):
        #maxima
        if (time_series[1,pt] > time_series[1,pt-1]):
          if (time_series[1,pt] > time_series[1,pt+1]):
            turn_ids.append(pt)

        #minima
        elif (time_series[1,pt] < time_series[1,pt-1]):
          if (time_series[1,pt] < time_series[1,pt+1]):
            turn_ids.append(pt)
    p = peaks_in_range
    return zip(turn_ids[0::p],turn_ids[p::p])
  
  def make(self,time_series):
    return self.peaks(time_series,2)

#Find, for each point, a window 
class AllPointsWindows(ExtremaWindows):
  #indices of the time_series' extrema
  def make(self,time_series):
    peaks = self.peaks(time_series,1)
    windows=[]
    #for w in range(0,len(onepeaks)-2):
    for w in zip(peaks,peaks[2:]):
      #location of the leftmost points of our ranges
      left_a = w[0][0]
      left_b = w[0][1]
      
      #location of the rightmost points: two peaks to right
      right_a = w[1][0]
      right_b = w[1][1]
      #print w[0], w[1]
      left_series = time_series[1,left_a:left_b]
      right_series = time_series[1,right_a:right_b]
      #Loop over the points in the left window. Find a point in the 
      #right window with the same value.
      #for i in range(left_a, left_b):
      for i in range(0,len(left_series)):
        val = left_series[i]
        shifted_series = (right_series - val)**2
        try:
          minimum_idx = right_a+np.argmin(shifted_series)
        except ValueError:
          print "Warning: series was empty."
        #print "Window: " ,i+left_a, minimum_idx
        
        print i, val,i+left_a, minimum_idx
        windows.append((i+left_a,minimum_idx)) 
    #print windows
    return windows



##############################################################################
#DIFFERENTIATORS
##############################################################################
from abc import ABCMeta, abstractmethod
class AbstractFirstDerivative:
  __metaclass__ = ABCMeta

  def __init__(self,time_series):
    f_t = TimeSeries(time_series).ts
    self.dfdt = self.first_deriv(f_t)
  
  @abstractmethod
  def first_deriv(self,time_series,dh):
    return NotImplemented

#2nd order central in volume; 1st one-sided at boundaries.
#Returned array is same shape as input. Assumes constant time-stepping.
class NumpyGradient(AbstractFirstDerivative):
  def first_deriv(self, f_t):
    dt = f_t[0,1] - f_t[0,0]
    f = f_t[1,:]
    df = np.gradient(f,dt)
    return np.vstack((f_t[0,:],df)) 
  
  


##############################################################################
#ANALYTIC SIGNAL APPROACH
##############################################################################

class AnalyticSignal:
  def __init__(self,time_series):
    from scipy.signal import hilbert as hilbert
    ts = TimeSeries(time_series).ts
    t,f = unstack(ts)
    hilb = np.imag(hilbert(f))
    f_a = f + 1j*hilb
    self.t = t
    self.f_a = f_a


class InstantaneousPhase:
  def __init__(self,time_series):
    analytic = AnalyticSignal(time_series)
    t = analytic.t
    f_a = analytic.f_a
    phi = np.angle(f_a)
    self.phase = np.vstack((t,phi))

class InstantaneousFrequency:
  def __init__(self,time_series):
    phase = InstantaneousPhase(time_series).phase
    self.freq = NumpyGradient(phase).dfdt

class InstantaneousFrequencyNoPhase:
  def __init__(self,time_series):
    analytic = AnalyticSignal(time_series)
    t = analytic.t
    f_a = analytic.f_a
    
    u = f_a.real
    v = f_a.imag
    mag = np.absolute(f_a)

    uts = np.vstack((t,u))
    vts = np.vstack((t,v))
    dudt = NumpyGradient(uts).dfdt[1,:]
    dvdt = NumpyGradient(vts).dfdt[1,:]

    freq = (u*dvdt - dudt*v)/(2*np.pi*mag)
    self.freq = np.vstack((t,freq))
    
    
  
##############################################################################
#FREQ MAKERS
##############################################################################
class FreqMaker:
  def __init__(self,window_maker):
    self.window_maker = window_maker
  def make(self,time_series):
    return -1

class DFTFrequencyMaker(FreqMaker): 
  import matplotlib.pyplot as plt
  def GetStrongestDFTPeak(self,data,window,n):
    sp = np.fft.fft(data)
    freq = np.fft.fftfreq(n)
    sp_mag = np.sqrt(sp.real**2 + sp.imag**2)
    max_idx = np.argmax(sp_mag)
    max_amp = np.amax(sp_mag)**0.5 / 2.*np.pi
    max_freq = 2.*np.pi*freq[max_idx]
    return max_amp, max_freq
  
  def DFTOverOneWindow(self,time_series,window):
    time = time_series[0,window[0]:window[1]]
    data = time_series[1,window[0]:window[1]]
    n = window[1]-window[0]
    amp = 0
    freq = 0
    while(True):
      plt.plot(time,data)
      plt.show()
      amp, freq = self.GetStrongestDFTPeak(data,window,n)
      if(freq!=0): break
      data = data - amp
      print amp, freq
    return amp, freq

class WindowedDFTFrequencyMaker(DFTFrequencyMaker):
  def DFTOverWindows(self,time_series):
    import matplotlib.pyplot as plt
    windows = self.window_maker.make(time_series)
    #zerothmodes
    times=[]
    amps=[]
    freqs=[]
    for window in windows:
      print window
      wsize = window[1]-window[0]
      this_time = time_series[0,window[0]+wsize/2]
      #plt.plot(time_series[0,:],time_series[1,:])
      #plt.plot(time_series[0,window[0]:window[1]],
      #        time_series[1,window[0]:window[1]])
      #plt.show()
      amp,freq = self.DFTOverOneWindow(time_series,window)
      print amp
      amps.append(amp)
      freqs.append(freq)
      times.append(this_time)
    plt.plot(times,freqs)
    plt.show()
    return freqs

  def DFTOverRollingWindows(self,time_series):
    import matplotlib.pyplot as plt
    print "Rolling..."
    windows = self.window_maker.make(time_series)
    #zerothmodes
    times=[]
    amps=[]
    freqs=[]
    for window in windows:
      size = window[1]-window[0]
      for i in range(window[0]+1,window[1]):
        this_window = (i,i+size)
        this_time = time_series[0,i+size/2]
        this_amp, this_freq = self.DFTOverOneWindow(time_series,this_window) 
        print "Time: ", this_time, "Freq: " , this_freq, " Amp: ", this_amp
        times.append(this_time)
        amps.append(this_amp)
        freqs.append(this_freq)
    plt.plot(times,freqs)
    plt.show()
    return freqs







##############################################################################
#TIME SERIES AND FREQS
##############################################################################
class TimeSeriesAndFreqs:
  """Represents a time series f(t) and its instantaneous frequencies"""
  """w(t)."""
  
  #feed me a numpy array A with A[:,0] the time and A[:,1] the data.
  def __init__(self, time,data,window_maker,freq_maker, name="TimeSeries"):
    self.time_series = np.vstack((time,data))
    windows = window_maker.make(time_series)
    self.freq_series = freq_maker.make(time_series,windows)
    self.name=name

  def plot_time_series(self):
    import matplotlib.pyplot as plt
    plt.plot(self.time_series[:,0],self.time_series[:,1],label=name)
    plt.show()
  
  def plot_freq_series(self):
    import matplotlib.pyplot as plt
    plt.plot(self.freq_series[:,0],self.freq_series[:,1],label="Omega_"+name)
    plt.show()


