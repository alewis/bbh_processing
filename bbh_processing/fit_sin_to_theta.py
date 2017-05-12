import peak_finders as pks
#from scipy import interpolate 
import scipy.optimize as opt 
import math 
import numpy as np
#import matplotlib.pyplot as plt

#Here we assume the signal is of the 'general chirp' form
#    y(t) = A(t) * sin[ phi(t) ] + M(t)
#with A(t) the amplitude, phi(t) the phase, and M(t) the offset.
#The instantaneous frequency is then
#    f = 1/2pi * d(phi)/dt
#Taking Q(t) = [y(t) - M(t)]/A(t) we then have
#    Q(t) = sin[ phi(t) ] 
#so that
#    d(phi)/dt = d/dt{ arcsin [ Q(t) ] }
#    f(t)      = 1/2pi * d/dt{ arcsin [ Q(t) ] }.
#  
#M(t) is estimated by performing an order offset_order regression on the data.
#  If offset_order=0, it is taken to be simply the mean.
#
#A(t) is estimated by fitting splines of order amp_order and smoothness
#  amp_smoothness to the peaks and troughs of the data, then computing at each
#  time t
#    A(t) = 0.5*{ |peak(t) - M(t)|/|trough(t) - M(t)| }.
#
#Q(t) can then be computed directly by interpolating M(t) and A(t) to time and
#  then appealing to its definiton. 
#Finally, we compute arcsin[Q(t)], fit a spline of order signal_order and 
#  smoothness signal_smoothness, and extract that spline's derivative.
#signal_order should be at least 3 to ensure a well-behaved derivative. If a 
# well-behaved second derivative is desired (to compute the stationary phase
# condition, for example) it should be at least 6.
def freq_from_splines(
    data,
    offset_order = 1,  #
    amp_order = 3,
    amp_smoothness = 0,
    signal_order = 5,
    signal_smoothness = 1):
  t = data[:,0]
  y = data[:,1]

  #Enforce correct types and assumptions
  if not isinstance(data, np.ndarray):
    raise TypeError("data should be a numpy array")
  if data.ndim != 2:
    raise ValueError("data should be 2 dimensional; first axis time")
  
  if not isinstance(offset_order, int):
    raise TypeError("offset_order must be of integer type")
  if (offset_order < 0):
    raise ValueError("offset_order must be nonnegative.")
  
  if not isinstance(amp_order, int):
    raise TypeError("amp_order must be of integer type")
  if (amp_order < 0):
    raise ValueError("amp_order must be nonnegative.")
  
  if not isinstance(amp_smoothness, int):
    raise TypeError("amp_smoothness must be of integer type")
  if (amp_smoothness < 0):
    raise ValueError("amp_smoothness must be nonnegative.")
  
  if not isinstance(signal_order, int):
    raise TypeError("signal_order must be of integer type")
  if (signal_order< 3):
    raise ValueError("signal_order must be at least 3.")
  
  if not isinstance(signal_smoothness, int):
    raise TypeError("signal_smoothness must be of integer type")
  if (signal_smoothness < 0):
    raise ValueError("signal_smoothness must be nonnegative.")
  
  #Compute the offset M
  M = do_fit(t,y,offset_order) 
  #Get peaks and troughs
  #The indices and values of the peaks.
  turns,turn_vals = pks.turning_points(y,0) 
  
  #Distinguish the peaks from the troughs 
  peaks, peak_vals, troughs, trough_vals = pks.peaks_and_troughs(
                                                turns,turn_vals)
  #Spline interpolations of both (this seems to work better than polyfits)
  pspl = interpolate.UnivariateSpline(t[peaks],peak_vals,
      k=amp_order, s= amp_smoothness)
  trspl = interpolate.UnivariateSpline(t[troughs],trough_vals,
      k=amp_order, s= amp_smoothness)
  P = np.zeros( len(t) )
  P = pspl(t)
  Tr = np.zeros( len(t) )
  Tr = trspl(t)
  
  #The amplitude is half the distanc between the peaks and trough
  #at each point.
  A = np.zeros( len(t) )
  A[:] = 0.5* ( P - Tr ) 

  #Now we can form Q
  Q = np.zeros( len(t) )
  Q[:] = (y[:] - M[:]) / A[:]
  
  hsig = signal.hilbert(Q)
  #print hsig
  #print np.absolute(hsig)
  #plt.plot(t,hsig.real,color="green",label="real")
  #print hsig.imag
  #plt.plot(t,hsig.imag,color="grey",label="imag")
  #We need Qdot
  #interpolate and take derivative
  Qspl = interpolate.UnivariateSpline(t,Q,
      k=signal_order,s=signal_smoothness)
  Qdot_spl = Qspl.derivative()
  Qdot = Qdot_spl(t)
  #plt.plot(t,Qdot/hsig.imag)
  #plt.plot(t,Qspl(t),label="Qspl")
  #plt.plot(t,Qdot_spl(t),label="Qdotspl")
  
  #We also need Q^2
  #restrict to [-1,1] and square
  np.clip(Q,-1.0,1.0,out=Q) 
  Qsq = np.square(Q)
  denom = np.sqrt(1.0-Qsq) 

  freq = np.zeros( len(t) )
  freq[:] = (1./(2.*np.pi))*(Qdot/denom)

  #spline interpolate the arcin to extract its derivative
  asQ = np.arcsin(Q)
  asQspl = interpolate.UnivariateSpline(t, asQ,
      k=signal_order, s= signal_smoothness)
  asQdot_spl = asQspl.derivative()
  #freq = np.zeros( len(t) )
  #freq[:] = (1./(2.*np.pi))*Qdot_spl( t )
  #plt.plot(t,freq,label="freq")
  
  #plt.plot(t,y,color = "blue", label="data")
  #plt.plot(t, asQ, color="violet", label="asQ")
  #plt.plot( t, Q, color="blue", label="Q" )
  #plt.plot( t, M, color = "black", label="offset" )
  #plt.plot( t, P, color = "red", label="peaks" )
  #plt.plot( t, Tr, color = "red", label="troughs" )

  plt.legend() 
  plt.show()
  return 0

def get_amp_mean(data,time):
  #to get the amplitude, we first get the peaks.
  #then at first approximation we consider the amplitude to be the average
  #distance between peak and mean.
  ptimes, pvals, ttimes, tvals = pks.turning_points(data) 
  vals = pvals + tvals
  mean = np.mean(data)
  distances = mean - vals 
  amp = np.mean(distances)
  return amp,mean 

#Sine wave in opt friendly language.
def sin_errfunc(p,x,y):
  return full_sin(p[0],p[1],p[2],p[3],x)-y
  
#Sine wave in opt friendly language.
def sin_only_freq(f,AMP,PHASE,OFFSET,x,y):
  return full_sin(AMP,f,PHASE,OFFSET,x)-y

#Sine wave in opt friendly language.
def sin_freq_phase(p,AMP,OFFSET,x,y):
  return full_sin(AMP,p[0],p[1],OFFSET,x)-y
def full_sin(AMP,FREQ,PHASE,OFFSET,x):
  return AMP*np.sin(2*np.pi*FREQ*(x+PHASE))+OFFSET

def full_cos(AMP,FREQ,PHASE,OFFSET,x):
  return AMP*np.cos(2*np.pi*FREQ*(x+PHASE))+OFFSET

def cos_phase(AMP,FREQ,OFFSET,y,x):
  arg = (y-OFFSET)/AMP
  if arg>1.0: 
    arg=1.0
  if arg<(-1.0):
    arg=-1.0
  phs = (np.arccos(arg) / (2.0*np.pi*FREQ)) - x
  return phs 

def sin_phase(AMP,FREQ,OFFSET,y,x):
  arg = (y-OFFSET)/AMP
  if arg>1.0: 
    arg=1.0
  if arg<-1.0:
    arg=-1.0
  phs = (np.arcsin(arg) / (2.0*np.pi*FREQ)) - x
  
  return phs


#def quadratic_chirp(x, AMP)
#  return AMP*np.sin(LIN*t + QUAD*t**2) + OFFSET

#def fit_quadratic_chirp(time,data):
#  amp,mean,rngs,freqs,phases=fit_sinusoids(time,data,width=2)
#  freq = np.mean(freqs)
#  phase = np.mean(phases)


#Creates guesses for the amplitude frequency, phase, and offset of windowed
#sine curves to time vs. theta data from an eccentric inspiral. 
#Width is the number of peaks within each window; thus width/2 cycles are covered
#in each range. Width = 1 (i.e. a half cycle) seems to work well. 
#The frequency is measured in 1/M, so you'll need a factor of 2pi in the sine
#function.
def fit_sinusoids(time,data,width=2):
  
  #plot the data
  #plt.plot(time,data,color="red",label="data")

  #get the amplitude and mean of the data (approx.constant)
  amp,mean = get_amp_mean(data,time) 
  
  #take peak-to-peak ranges of data
  ptimes, pvals, ttimes, tvals = pks.turning_points(data) 
  rngs = zip(turns[:-width:width],turns[width::width])
  
  #number of cycles each rng covers
  ncycles = width/2.0
  
  #set this to 2*pi if you want an angular frequency
  fscale = 1.0
  freqs=[]
  for a,b in rngs:
    #the frequency in radians per second
    period = (time[b]-time[a])/ncycles 
    freqs.append(fscale/period) 
 
  #set the phase to ensure the fitted curve will pass through the data at the 
  #initial point
  phases=[]
  for rng,freq in zip(rngs,freqs):
    t = time[rng[0]]
    th = data[rng[0]]
    phases.append(sin_phase(amp,freq,mean,th,t))
  
  return amp, mean, rngs, freqs, phases 


#This refines the guesses from fit_sinusoids by performing nonlinear fits
def fit_sin_from_guessed_params(time,data,width=2):
  #plt.plot(time,data,color="red") 
  amp,mean,rngs,freqs,phases = fit_sinusoids(time,data,width=width)
  z = zip(rngs,freqs,phases) 

  ps = []
  midtimes = []
  
  for rng, freq, phase in z:
    a = rng[0]
    b = rng[1]
    p0 = [amp,freq,phase,mean]
    t = time[a:b] 
    th = data[a:b]
    midtimes.append(time[a] + 0.5*(time[b]-time[a])) 
    
    p1, success = opt.leastsq(sin_errfunc,p0[:],args=(t,th)) 
    ps.append(p1) 
  return rngs,ps,midtimes



