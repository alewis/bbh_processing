import cPickle as pickle
import bbh_processing.hfile_tools as hf
import bbh_processing.peak_finders as peaks
import numpy as np

#We'll need the Omega.pkl and r data

def first_peak_before_t(
    hname="Horizons.h5",
    pname="Omega.pkl"):
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  t = hf.t(hname)
  r = hf.r_mag(hname)
  turns,turn_vals=  peaks.turning_points(OmegaMag[:,1])
  
  peakt = t[turns]
  peakr = r[turns]
  peakom = OmegaMag[turns,1]

  for a,b,c in zip(peakt,peakr,peakom):
    print a,b,c
if __name__=="__main__":
  first_peak_before_t()
