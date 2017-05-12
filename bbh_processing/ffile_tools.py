#These process data from fluxes.dat files (the outputs from extract E and J.py)
import numpy as np

def t(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(0,))
def de_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(1,))
def djx_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(2,))
def djy_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(3,))
def djz_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(4,))
def dpx_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(5,))
def dpy_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(6,))
def dpz_dt(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(7,))
def e(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(8,))
def jx(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(9,))
def jy(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(10,))
def jz(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(11,))
def px(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(12,))
def py(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(13,))
def pz(fname="Fluxes.dat"):
  return np.loadtxt(fname,usecols=(14,))
  
def dj_dt(fname="Fluxes.dat"):
  return np.vstack((djx_dt(fname),djy_dt(fname),djz_dt(fname))).transpose()
def j(fname="Fluxes.dat"):
  return np.vstack((jx(fname),jy(fname),jz(fname))).transpose()
def dp_dt(fname="Fluxes.dat"):
  return np.vstack((dpx_dt(fname),dpy_dt(fname),dpz_dt(fname))).transpose()
def p(fname="Fluxes.dat"):
  return np.vstack((px(fname),py(fname),pz(fname))).transpose()

