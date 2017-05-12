
import omega_theta_methods as omth
import omega_r_methods as omr 
import omega_phi_methods as omphi
#import bbh_processing.plot_makers as plot_makers
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np

#compute the 'angle variables' - W_i = Omega_i * t
def W_i(i="theta",hname="Horizons.h5",pname="Omega.pkl"):
  if i=="theta":
    om_i = omth.omega_theta_from_fits(width=4,hname=hname)
  elif i=="r":
    with open(pname) as f:
      Omega,OmegaMag = pickle.load(f)
    om_i = omr.omega_r_from_ranges(OmegaMag,width=4)
  elif i=="phi":
    with open(pname) as f:
      Omega,OmegaMag = pickle.load(f)
    om_i = omphi.omega_phi_from_time_averages(OmegaMag,width=4)
  else:
    raise ValueError("i must be 'theta', 'r', or 'phi'")
  
  W = np.zeros( (om_i.shape) )
  W[:,0] = np.copy( (om_i[:,0]) )
  P = om_i[:,0] * om_i[:,1]
  #impose 2pi periodicity
  W[:,1] = P #% (2.0*np.pi)
  return W

def torus_plot_around_resonance_th_r(ti,tf):
  from scipy import interpolate as interp
  W_th = W_i(i="theta")
  #plt.plot(W_th[:,0],W_th[:,1])
  W_r  = W_i(i="r")
  Wthspl = interp.UnivariateSpline(W_th[:,0],W_th[:,1],k=3,s=0)
  Wrspl = interp.UnivariateSpline(W_r[:,0],W_r[:,1],k=3,s=0)
  sz = 2.*(tf-ti)
  traj = np.zeros((sz,2)) 
  trng = np.linspace(ti,tf,sz) 
  traj[:,0]=Wthspl(trng)%(2.0*np.pi)
  traj[:,1]=Wrspl(trng)%(2.0*np.pi)
  #t = np.linspace(0,1000,10000)
  #y1 = 4.0*t % 2.0*np.pi
  #y2 = 3.0*t % 2.0*np.pi
  #segs = segmented_lines(y1,y2,5)
  #for s in segs:
  #  plt.plot(s[:,0],s[:,1])
  #plt.show()
  return segmented_lines(traj[:,0],traj[:,1],5)

def segmented_lines(W1,W2,tol):
  #find the line breaks
  pos1 = np.where(np.abs(np.diff(W1[:]))>=tol)[0]+1
  pos2 = np.where(np.abs(np.diff(W2[:]))>=tol)[0]+1
  pos = np.sort(np.append(pos1,pos2))
  i=0 
  #if there are no breaks, just return the whole array 
  segments=[]
  if pos.size==0:
    thisarr=np.zeros((len(W1),2))
    thisarr[:,0]=np.copy(W1[:])
    thisarr[:,1]=np.copy(W2[:])
    segments.append(thisarr)
  else:
    for p in pos:
      thisarr=np.zeros((len(W1[i:p]),2))
      thisarr[:,0]=np.copy(W1[i:p])
      thisarr[:,1]=np.copy(W2[i:p])
      segments.append(thisarr)
      i=p
  return segments 
def torus_plot_th_r(hname="Horizons.h5",pname="Omega.pkl"):
  W_th = W_i(i="theta")
  W_th[:,1] = W_th[:,1]%(2.0*np.pi)
  #plt.plot(W_th[:,0],W_th[:,1])
  W_r  = W_i(i="r")
  W_r[:,1] = W_r[:,1]%(2.0*np.pi)
  #find the line breaks
  return segmented_lines(W_th[:,1],W_r[:,1],5) 
  
  #rat = [5./4.,3./2.,6./5.]
  #col = ["red","blue","green"]
  #for s in segments: 
  #  import math
  #  from scipy import stats
  #  m,b,r_val,p_val,stderr = stats.linregress(s[:,0],s[:,1])
  #  if (not math.isnan(m)):
  #    for r,c in zip(rat,col):
  #      if np.fabs(r-m) < 0.1:
  #        plt.plot(s[:,0],s[:,1],color=c)
  #insert nan at the discontinuities so pylab doesn't try to join them
  #W_r_t = np.copy(W_r[:,0])
  #plt.plot(W_r[:,1],W_th[:,1])
  #plt.xlim([0,2*np.pi])
  #plt.ylim([0,2*np.pi])
  #plt.legend()
  #plt.show()
