#ndAngle between two vectors A and B.
"Utility functions not directly invoking a particular kind of file."
import numpy as np
import numpy.linalg as la
from scipy import optimize
from numpy.lib import stride_tricks
import peak_finders as peaks
#returns the angle between two (arrays of) vectors
def angle(A,B):
  return np.arccos(cosang(A,B))

#returns cos(theta) from two (arrays of) vectors
def cosang(A,B):

  if A.shape != B.shape:
    raise ValueError("Arrays had inconsistent sizes.")
  if len(A.shape)!= 2:
    raise ValueError("Input should be an array of vectors.")
  dotAB = np.einsum('ij,ij->i',A,B)
  dotAA = np.einsum('ij,ij->i',A,A)
  dotBB = np.einsum('ij,ij->i',B,B)
  return np.divide(dotAB, np.sqrt(dotAA*dotBB))

# def cosang2d(A, B):
  # """Handle the 2D case.
  # """
  # dotAB = np.einsum('i,i->i', A, B)   

#returns sin(theta) from two (arrays of) vectors
def sinang3D(A,B):
  crossAB = la.norm(np.cross(A,B),axis=1)
  dotAA = np.einsum('ij,ij->i',A,A)
  dotBB = np.einsum('ij,ij->i',B,B)
  return np.divide(crossAB,np.sqrt(dotAA*dotBB))  


def unwrapped_ang(ang):
    peaktimes, peakvals, troughtimes, troughvals = peaks.turning_points(ang)  
    unwrapped = np.zeros(ang.shape)
    grad = np.gradient(ang, 0.5)
    unwrapped = ang[0]+np.cumsum(np.fabs(grad)) 
    return unwrapped  

#Returns the 'unwrapped phase' between A and B (in the plane): similar to 
#the angle between them, but expressed as a monotonically increasing rather
#than a periodic function.
def phase(A, B):
    ang = angle(A, B)
    return np.unwrap(ang)
def cosphase(A, B):
    ang = cosang(A, B)
    return np.unwrap(ang)
#def phase(A,B):
#  Acomplex = A[:,0] + 1.0j*A[:,1]
#  Bcomplex = B[:,0] + 1.0j*B[:,1]
#  ang = np.angle(Acomplex) - np.angle(Bcomplex) 
#  #return ang
#  return np.unwrap(ang)


def one_over_r(x,a,b):
  return a*1./x + b
#compute the results of 'function' on 'data' using window sizes from minwidth to
#maxwidth by stride. Plot results as a function of window size.
def extrapolated_windowed_series(function,data,minwidth=1,maxwidth=10,stride=2,
    axis=1):
  import matplotlib.pyplot as plt
  windowrange = range(minwidth,maxwidth,stride)
  outdata = [function(data,N)[:,1] for N in windowrange]
  out_array = np.vstack(outdata)
 

  #extrap_list=[]
  #for i in range(0,len(out_array[0,:])):
    #fit,cov = optimize.curve_fit(one_over_r,windowrange,out_array[:,i])
    #print fit
    #extrap_list.append(one_over_r(10000,fit[0],fit[1]))
  means = np.average(out_array,axis=0)
  step=15000
  print "means:", means[step]
  thismean = means[step]
  thismeanarr=np.zeros(len(windowrange))
  thismeanarr[:]=thismean
  print "thismeanarr:", thismeanarr
  plt.plot(windowrange,out_array[:,step])
  plt.plot(windowrange,thismeanarr)
  plt.show()
  #stddevs = np.std(out_array,axis=0)
  return np.array(extrap_list) 

#compute the results of 'function' on 'data' using window sizes from minwidth to
#maxwidth by stride. Average the results.
def mean_windowed_series(function,data,minwidth=2,maxwidth=10,stride=1,
    axis=1):
  windowrange = range(minwidth,maxwidth,stride)
  outdata = [function(data,N)[:,1] for N in windowrange]
  out_array = np.vstack(outdata)
  
  means = np.average(out_array,axis=0)
  stddevs = np.std(out_array,axis=0)
  return means, stddevs

#Dot product of vectors in a time series. Returns a vector of dot products,
#one for each timestep. The arrays are assumed to be columnwise; i.e. A[:,0] 
#refers to the 0th time series. The arrays passed in should *not* 
#themselves include the times.
def timedot(A,B):
  if (A.shape != B.shape):
    raise ValueError("Inputs to dot product had different shapes.")
  return np.einsum('ij,ij->i',A,B) #this does the dot product


#Given an independent array y, a dependent array x (representing a 
#discretized function y=f(x)), and a value y0, return all x_y0 such that
#y0 = f(x_y0).
def invert(y, x, y0, interp_type="slinear"):
  if (x.shape != y.shape):
    raise ValueError("Inputs x and y to 'invert' had different shapes.")
  
  #first interpolate the function
  from scipy import interpolate, optimize
  y_p = np.abs(y-y0)
  f_x = interpolate.interp1d(x,y_p,kind=interp_type)
  


  def y_as_function(x):
    return f_x(x)

  print optimize.fsolve(y_as_function,2)
  
  
  #delta = np.empty(y.shape)
  #delta.fill(1.e-3)
  #print y_p
  #print np.minimum(y_p,delta)
   
  #now find crossings of the x axis
  x_i = []
  #lastx = x[0]
  #lasty = y[0]
  #for xi,yi in zip(x[1:],y_p[1:]):
  #  if (yi==0):
  #    x_i.append(xi)
  #  elif (lasty * yi < 0):
  #    x_i.append(0.5*(xi-lastx))
  #  lastx = xi
  #  lasty = yi
  #print "\n", x_i, "***\n" 
  return x_i

#Returns a unit vector parallel to each A in a time series.
def unit_vector(A):
  invmagA = 1./la.norm(A,axis=1)
  return scaled_vector(invmagA, A)
      #multiplies every element in A by the element of invmagA in its same row

def scaled_vector(c, A):
  return np.einsum('i,ij->ij', c, A)

#The coordinates of the separation vector projected into the orbital plane
def r_in_orbital_plane(r,chi,t=None):
  #Project r onto the plane normal to chi
  r_onto_chi = projected_onto(t,r,chi) 
  r_onto_plane = r-r_onto_chi[:,1:]
  if t is None:
    return r_onto_plane
  else:
    return np.column_stack([t,r_onto_plane])

#Generate a unit vector orthonormal to both vector(t) and the input
def orthonormal_axis(vec1,vec2):
  for v in (vec1,vec2):
    if v.shape[1] != 3:
      raise ValueError("input must be 3D")
    if len(vec1[:,0])!=len(vec2[:,0]):
      raise ValueError("input must have same length as series")
  orthogonal_v = np.cross(vec1,vec2)
  orthogonal_mag = la.norm(orthogonal_v,axis=1)
  vx = np.divide(orthogonal_v[:,0],orthogonal_mag)
  vy = np.divide(orthogonal_v[:,1],orthogonal_mag)
  vz = np.divide(orthogonal_v[:,2],orthogonal_mag)
  orthonormal_v = np.column_stack([vx,vy,vz])
  return orthonormal_v


#Calculate and 'a#ccumulating' phi variable; i.e. one which is not
#periodic.
def accumulating_phi(r_array):
  rho = r_array[:,:2]
  rho_0 = rho[0,:]
  #zeros = np.zeros(rho.shape[0]) 
  #rho_3 = np.column_stack([rho,zeros])
  rho_ref = stride_tricks.as_strided(rho_0,strides=(0,1*8),shape=rho.shape)
  return phase(rho,rho_ref) 


#Projection of A onto B.
def one_value_projected_onto(A,B):
  b_hat = unit_vector(B)
  scalar_projection = np.dot(A,b_hat)
  return [scalar_projection*b for b in b_hat] 

#A and B should be the spatial parts only
def projected_onto(t,A,B):
  b_hat = unit_vector(B)
  #need element-wise dot product - surprisingly, no native numpy
  #function does this.
  dotproducts = np.sum(A*B,axis=1) 
  C_space = np.zeros(B.shape)
  for i in range(0,C_space.shape[1]):
    C_space[:,i] =  dotproducts[:] * b_hat[:,i] 
  C = np.column_stack([np.copy(t),C_space])
  return C

def projected_onto(A,B):
  b_hat = unit_vector(B)
  #need element-wise dot product - surprisingly, no native numpy
  #function does this.
  dotproducts = np.sum(A*B,axis=1) 
  C = np.zeros(B.shape)
  for i in range(0,C.shape[1]):
    C[:,i] =  dotproducts[:] * b_hat[:,i] 
  return C

#Here b is just one vector instead of a series
def projected_onto_const(A,B):
  b_hat = unit_vector(B)
  #need element-wise dot product - surprisingly, no native numpy
  #function does this.
  dotproducts = np.sum(A*B,axis=1) 
  C = np.zeros(B.shape)
  for i in range(0,C_space.shape[1]):
    C[:,i] =  dotproducts[:] * b_hat[:,i] 
  return C

def interp_onto(data,oldt,newt,method_string='cubic'):
  from scipy.interpolate import interp1d
  dims = 0
  try:
    dims = len(data.shape)
  except AttributeError: #thrown if 'data' is a list
    dims = 1

  if dims==1:
    C = np.zeros(len(newt))
    f = interp1d(oldt,data,kind=method_string)
    C[:] = f(newt)  
  else:
    nrows = len(newt)
    ncols = data.shape[1]
    C = np.zeros((nrows,ncols), dtype=np.float64)
    for i in np.arange(0,ncols):
      #print oldt.shape, data[:,i].shape
      f = interp1d(oldt,data[:,i],kind=method_string)
      C[:,i] = f(newt)
  return C


#Trims long_array so its range is the same as short_array. Both are assumed
#to be 1D and monotonic - normally they will be the time in a time series. 
#Useful when interpolating.
def trim_array(long_array,short_array):
  min_datum = short_array[0]
  max_datum = short_array[-1]
  min_idx = np.where(long_array==min_datum)[0][0]
  max_idx = np.where(long_array==max_datum)[0][0]
  return long_array[min_idx:max_idx]

def trim_array_left(long_array,short_array):
  min_datum = short_array[0]
  min_idx = np.where(long_array==min_datum)[0][0]
  return long_array[min_idx:-1]

def trim_array_right(long_array,short_array):
  max_datum = short_array[-1]
  max_idx = np.where(long_array==max_datum)[0][0]
  return short_array[0:max_idx]

def trim_time_with_data(long_t,long_data,short_t,short_data):
  trimmed_time = trim_array(long_t,short_t)
  trimmed_long_data = long_data[trimmed_time[0]:trimmed_time[-1]]
  return trimmed_time,trimmed_long_data

#Get all the extrema of a curve "data". Turns stores the indices of the
#extrema, while data stores their values. rngs will store tuples of 
#the first_left-th to the first_right-th extremum, skipping by stride. 
#The actual code has been moved to peak_finders.py; this is just here so as
#not to break older scripts.
def turning_points(data, first_left=0, first_right=4, stride=1):
  from peak_finders import turning_points_with_range
  return turning_points_with_range(data,first_left=first_left, first_right=first_right,
      stride=stride)
