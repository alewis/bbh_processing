#Generate a dat file plotting omega_r / omega_phi vs omega_phi 

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import bbh_processing.ecctools as ecctools
import bbh_processing.hfile_tools as hfile
import bbh_processing.omega_phi_methods as omphi
import bbh_processing.omega_r_methods as omr
import bbh_processing.combined_tools as ct
import bbh_processing.ks_from_datfile as ks
import cPickle as pickle

def kthr_vs_omphi():
  hname = "Horizons.h5"
  pname = "Omega.pkl"

  #Need the orbital frequency to get omega phi
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  
  #Get omegas
  om_phi = omphi.omega_phi_from_time_averages(OmegaMag,2)
  om_r = omr.omega_r_from_ranges(OmegaMag,2)
  tht,om_th = ct.omega_theta(hname)
  om_th_arr = np.zeros((len(tht),2))
  om_th_arr[:,0] = tht
  om_th_arr[:,1] = om_th
  t, krphi,krth,kthphi = ks.ks_from_input(om_r,om_phi,
      om_th_arr,equatorial=False)
  return 1./krth,om_phi 

if __name__=="__main__":
  kthr,om_phi = kthr_vs_omphi()
  with open("kthr_vs_omphi.dat","w") as f:
    f.write("""#
      # [1] = omega_phi 
      # [2] = omega_th/omega_r
      """)
    np.savetxt(f,np.transpose([om_phi[:,1],kthr]))
  print "Done!"
