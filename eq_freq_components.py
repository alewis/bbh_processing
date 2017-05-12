""
"This generates plots of Omega_phi vs. N for an equatorial run."
""
import numpy as np
import cPickle as pickle
import bbh_processing.equatorial_freq_methods as bbhfreq
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def gen_components(path):
    with open(path+"/Omega.pkl") as f:
      Omega,OmegaMag = pickle.load(f)
    omr_extrap,omr_errs = bbhfreq.extrapolated_omega(
                                    "r",OmegaMag[:,0],OmegaMag[:,1],
                                    range(1,30,1))
    #valcol = "blue"
    #errcol = "red"
    #plt.plot(OmegaMag[:,0],omr_extrap,color=valcol)
    #plt.plot(OmegaMag[:,0],omr_extrap+omr_errs,color=errcol)
    #plt.plot(OmegaMag[:,0],omr_extrap-omr_errs,color=errcol)
    omph_extrap,omph_errs = bbhfreq.extrapolated_omega(
                                    "phi",OmegaMag[:,0],OmegaMag[:,1],
                                    range(1,30,1))
    with open(path+"/frequency_components.pkl","wb") as g:
        pickle.dump([omr_extrap,omr_errs,omph_extrap,omph_errs],g)
    #valcol = "blue"
    #errcol = "green"
    #plt.plot(OmegaMag[:,0],omph_extrap,color=valcol)
    #plt.plot(OmegaMag[:,0],omph_extrap+omph_errs,color=errcol)
    #plt.plot(OmegaMag[:,0],omph_extrap-omph_errs,color=errcol)
    #plt.show()

if __name__ == "__main__":
    p1 = "/Users/adamlewis/d2/AaronsEccSeries2/Ecc_pt2_Equatorial"
    p2 = "/Users/adamlewis/d2/AaronsEccSeries2/Ecc_pt3_Equatorial"
    p3 = "/Users/adamlewis/d2/AaronRuns/Equatorial"
    p4 = "/Users/adamlewis/d2/AaronRuns/EccEqPt3"
    lst = [p1,p2,p3,p4]
    for p in lst:
      gen_components(p)
  



