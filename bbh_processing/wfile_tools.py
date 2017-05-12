import numpy as np

#These process data from Abdul's script.
def t(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(0,))

def omega_phi(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(1,))

def omega_r(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(2,))

def k(wname="W1.6.dat"):
  return np.loadtxt(wname,usecols=(3,))
