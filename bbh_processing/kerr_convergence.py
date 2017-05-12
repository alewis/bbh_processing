from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def omega_phi_kerr_converge(time,data,N=1, minN=0):
  splrep = UnivariateSpline(time,data,s=0,k=4)
  tpks   = splrep.derivative().roots()[1:]
  deltat = tpks[N]-tpks[minN]
  defint = splrep.integral(tpks[minN],tpks[N])
  omega_phi = defint/deltat
  return omega_phi

def omega_r_kerr_converge(time,data,N=1, minN=0):
  splrep = UnivariateSpline(time,data,s=0,k=4)
  tpks   = splrep.derivative().roots()[1:]
  deltat = tpks[N]-tpks[minN]
  omega_r = 2.0*np.pi*(N-minN)/(2*deltat)
  return omega_r

def guess(y1,ymax,y2=None):
  c = ymax
  if y2==None:
    y0guess = 0.0
  else:
    y0guess = 2.0*y1-y2
  a = y0guess-c
  if y1==ymax:
    b=0
  else:
    b = -1. * np.log((y1 - c)/a) 
  return a,b,c

#Nrange is a range of window sizes used to calculate series.
#We will be using these as x,y ranges to fit a y= a*exp(-b*x)+c
#law. c is thus the predicted value for N=inf.
def extrap_exp(nrange,series):
  try:
    a,b,c = guess(series[0],series[-1],series[1]) 
  except IndexError:
    a,b,c = guess(series[0],series[-1])
  
  def myexp(x,d):
    return a*np.exp(-b*x) + d
  print len(nrange), len(series) 
  popt,pcov = opt.curve_fit(myexp,
                       nrange,
                       series,
                       c
                      )
  return popt
