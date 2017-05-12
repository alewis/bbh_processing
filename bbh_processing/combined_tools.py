
#Finds peaks and takes the frequency as 2*pi/(peak-to-peak distance)
def frequencies_from_peak_finder(time,data,left=0,right=3):
  import bbh_processing.utilities as bbh
  import numpy as np
  
  turns,turn_vals,rngs = bbh.turning_points(data,
      first_left=left,first_right=right, stride=1)
  
  mid_times = [time[a+(b-a)/2] for a,b in rngs]
  mid_freqs = [1./(time[b]-time[a]) for a,b in rngs] 
  
  trimmed_time = bbh.trim_array(time,mid_times)
  
  freqs = bbh.interp_onto(mid_freqs,mid_times,trimmed_time,
                                        method_string='cubic')
  return trimmed_time, 2*np.pi*freqs


#def ecc_from_e_l(pname="Omega.pkl", fname="Fluxes.dat"):
#  import cPickle as pickle
#  import bbh_processing.utilities as bbh
#  f = open(pname)
#  Omega,OmegaMag = pickle.load(f)
#  f.close()
  

#Omega r: frequency of omega total oscillations
def omega_r(pname="Omega.pkl",method="ranges",width=2):
  import omega_r_methods as om
  import cPickle as pickle
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  if method=="ranges":
    omr = om.omega_r_from_ranges(OmegaMag)
  return omr[:,0],omr[:,1] 
   
#Omega phi: monotonic component of omega total
def omega_phi(pname="Omega.pkl",method="time_averages",width=2):
  import omega_phi_methods as om
  import cPickle as pickle
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  if method=="envelopes":
    omphi = om.omega_phi_from_envelopes(OmegaMag)
  elif method=="time_averages":
    omphi = om.omega_phi_from_time_averages(OmegaMag,width)
  return omphi[:,0],omphi[:,1] 



def omega_phi_with_fit(pname="Omega.pkl"):
  import cPickle as pickle
  import bbh_processing.utilities as bbh
  import scipy.optimize as opt 
  import numpy as np
  #fit to a monotonic power law
  def freq_func(p,t):
    return p[0]*(p[1]-t)**p[2]
  
  def error_func(p, t, data):
    return p[0]*(p[1]-t)**p[2] \
          -data
  
  f = open(pname)
  Omega,OmegaMag = pickle.load(f)
  f.close()
  t = OmegaMag[:,0]
  om = OmegaMag[:,1]
  
  guess_time,guess_phi = omega_phi(pname)

  p1_0 = np.max(guess_time) #shift
  p2_0 = (-3./2.) #kepler
  p0_0 = guess_phi[10]/(p1_0-guess_time[10])**p2_0
  p0=[p0_0,p1_0,p2_0] 
  

  p=[0,0,0]
  p, success = opt.leastsq(error_func, p0[:], args=(t,om))

  om_phi = np.zeros(t.shape)
  om_phi[:] = freq_func(p,t)

  return t, om_phi



#Frequency of the theta_r_chi as a function of time.
def omega_theta(hname="Horizons.h5"):
  import omega_theta_methods as omth
  omega_theta_arr = omth.omega_theta_from_fit(hname)
  return omega_theta_arr[:,0],omega_theta_arr[:,1]
  
