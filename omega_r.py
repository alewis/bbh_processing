import bbh_processing.omega_r_methods as omr
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
def omega_r(pklfile = "Omega.pkl", method="time_averages"):
  with open(pklfile) as f:
    Omega,OmegaMag = pickle.load(f)

  if method=="envelopes":
    om_r = omr.omega_r_from_envelopes(OmegaMag)
  elif method=="time_averages":
    om_r = omr.omega_r_from_time_averages(OmegaMag,width=2)
    om_r1 = omr.omega_r_from_time_averages(OmegaMag,width=1)
    om_r3 = omr.omega_r_from_time_averages(OmegaMag,width=3)
    om_r4 = omr.omega_r_from_time_averages(OmegaMag,width=4)
    t = OmegaMag[:,0]
    plt.plot(t,OmegaMag[:,1])
    plt.plot(t,om_r1,label="w=1")
    plt.plot(t,om_r,label="w=2")
    plt.plot(t,om_r3,label="w=3")
    plt.plot(t,om_r4,label="w=4")
  return om_r

if __name__ == "__main__":
  om_r = omega_r(method="time_averages")
  om_re = omega_r(method="envelopes")
  plt.plot(om_re[:,0],om_re[:,1],label="envelopes")
  plt.legend()
  plt.show()

