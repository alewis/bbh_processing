#Generate a dat file plotting eccentricity against omega_phi

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import bbh_processing.ecctools as ecctools
import bbh_processing.hfile_tools as hfile
import bbh_processing.omega_phi_methods as omphi
import cPickle as pickle

def ecc_vs_omphi():
  hname = "Horizons.h5"
  pname = "Omega.pkl"

  #We need the separation to get the eccentricity
  rmag = hfile.r_mag(hname)
  ht = hfile.t(hname)
  #And the orbital frequency to get omega phi
  with open(pname) as f:
    Omega,OmegaMag = pickle.load(f)
  #Get the eccentricity
  errode,peakt,peakecc = ecctools.ecc_peaks_from_data(ht,rmag)
  #Interpolate onto the orbital frequency times
  eccspline = interp.UnivariateSpline(peakt,peakecc,k=3,s=0)
  ecc = eccspline(OmegaMag[:,0]) 

  #Get omega phi
  om_phi = omphi.omega_phi_from_time_averages(OmegaMag,2)
  return om_phi[:,1],ecc

if __name__=="__main__":
  om_phi,ecc = ecc_vs_omphi()
  with open("ecc_vs_omphi.dat","w") as f:
    f.write("""#
      # [1] = omega_phi
      # [2] = ecc
      """)
    np.savetxt(f,np.transpose([om_phi,ecc]))
  print "Done!"
