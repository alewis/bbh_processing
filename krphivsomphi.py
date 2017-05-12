#Generate a dat file plotting omega_r / omega_phi vs omega_phi 

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import bbh_processing.ecctools as ecctools
import bbh_processing.hfile_tools as hfile
import bbh_processing.omega_phi_methods as omphi
import bbh_processing.omega_r_methods as omr
import bbh_processing.ks_from_datfile as ks
import cPickle as pickle

def krphi_vs_omphi():
  hname = "Horizons.h5"
  pname = "Omega.pkl"

  #Need the orbital frequency to get omega phi
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  
  #Get omegas
  om_phi = omphi.omega_phi_from_time_averages(OmegaMag,2)
  om_r = omr.omega_r_from_ranges(OmegaMag,2)
  
  t, krphi = ks.ks_from_input(om_r,om_phi,equatorial=True)
  return krphi,om_phi 

if __name__=="__main__":
  krphi,om_phi = krphi_vs_omphi()
  with open("krphi_vs_omphi.dat","w") as f:
    f.write("""#
      # [1] = omega_phi 
      # [2] = omega_r/omega_phi
      """)
    np.savetxt(f,np.transpose([om_phi[:,1],krphi]))
  print "Done!"
