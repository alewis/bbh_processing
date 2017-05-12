import numpy as np
from scipy.interpolate import interp1d
#Behaves similarly to the matlab 'interp' function.

def interp(ys, mul):
  #import pdb; pdb.set_trace()
  #linear extrapolation for last (mul-1) points
  ys = list(ys)
  ys.append(2*ys[-1] - ys[-2])
  #make interpolation function
  xs = np.arange(len(ys))
  print "Up to interp1d"
  fn = interp1d(xs,ys,kind='slinear')
  #call it on desired data points
  new_xs = np.arange(len(ys)-1,step=1./mul)
  print "Exiting interp"
  return fn(new_xs)
