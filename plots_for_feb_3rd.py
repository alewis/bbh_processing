import bbh_processing.plot_makers as plot_makers
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser

prs=OptionParser()
prs.add_option("--name",type="string",help="name of top level directory for plot")
#prs.add_option("--eq",action="store_true",default=False,help="equatorial run")
(opts,args) = prs.parse_args()


def make_plots(hname="Horizons.h5",fname="Fluxes.dat",
    pname="Omega.pkl", wname="W1.6.dat",
    outname="plots.pdf"):
  rows = 2
  cols = 4
  i=1

  #theta(t)
  plt.subplot(rows,cols,i)
  plot_makers.theta_t(hname)
  i+=1 
  
  #chi dot L
  plt.subplot(rows,cols,i)
  plot_makers.L_dot_chi(hname,pname)
  i+=1

  #ecc_from_ap/peri
  plt.subplot(rows,cols,i)
  plot_makers.ecc_from_app_peri(hname)
  i+=1

  pathname = "/home/adlewis/m/ResonanceAnalysis/Notes/GenericBBH/Notes/plots_from_aarons_runs/feb3rd"
  outname = opts.name+".pdf" 
  fullname = pathname+outname

  fig = plt.gcf()
  fig.set_size_inches(18.5,10.5)
  fig.savefig(fullname)
  #plt.show()

make_plots()
