import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import bbh_processing.ffile_tools as ffile

#Return the energy flux interpolated to OmegaMag and OmegaMag
def make_edot_vs_omega(path):
  fluxpath = path+"/Fluxes.dat"
  ft = ffile.t(fluxpath)
  edot = ffile.de_dt(fluxpath)
  
  omegapath = path+"/Omega.pkl"
  with open(omegapath) as f:
    Omega,OmegaMag = pickle.load(f)

  edotspline = interp.UnivariateSpline(ft,edot,k=3,s=0) 
  edotinterp = edotspline(OmegaMag[:,0])
  return OmegaMag,edotinterp

#Return a list of (eta, OmegaMag, edot)
def make_list_of_fluxes():
  dirs = ["q5","q7","q7.5","q8.5","q9","q9.5","q10"]
  qs = [5.,7.,7.5,8.5,9.,9.5,10.]
  etas = [q/((q+1.)**2) for q in qs]
  omega_list = []
  edot_list = []
  for p,n in zip(dirs,etas):
    this_om,this_edot = make_edot_vs_omega(p)
    omega_list.append(this_om)
    edot_list.append(this_edot)
  return etas,omega_list,edot_list

#Plot the energy flux as a function of eta for constant omega
def flux_at_constant_om(omega):
  fname = "Omega"+str(omega)+".dat"
  etas,oms,edots = make_list_of_fluxes()
   
  fluxlist = []
  for om,edot in zip(oms,edots):
    edotspline = interp.UnivariateSpline(om[:,1],edot,k=3)
    fluxlist.append(float(edotspline(omega)))
  return etas,fluxlist,fname

if __name__=="__main__":
  omegas=[0.03,0.035,0.04,0.045,0.05,0.055,0.06]
  for om in omegas:
    eta,flux,fname = flux_at_constant_om(om)
    with open(fname,"w") as f:
      f.write("""#
        # [1] = eta
        # [2] = edot
        """)
      np.savetxt(f,np.transpose([eta,flux]))
  print "Done!"    
  
