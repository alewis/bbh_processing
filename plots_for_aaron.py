import bbh_processing.plot_makers as plot_makers
import matplotlib.pyplot as plt
import numpy as np
from optparse import OptionParser

prs=OptionParser()
prs.add_option("-d",action="store_true",default=False,help="display only")
prs.add_option("--name",type="string",help="name of top level directory for plot")
prs.add_option("--eq",action="store_true",default=False,help="equatorial run")
prs.add_option("--restimes",action="store_true",default=False,help="display times of resonant crossings")
(opts,args) = prs.parse_args()


def make_plots(hname="Horizons.h5",fname="Fluxes.dat",
    pname="Omega.pkl", wname="W1.6.dat",
    outname="plots.pdf"):
  rows = 2
  cols = 4
  
  i=1
    
  #plt.subplot(rows,cols,i)
  #plot_makers.j(fname)
  #i+=1 
  
  plt.subplot(rows,cols,i)
  plot_makers.j_on_s(hname,fname,opts.eq)
  i+=1 
  
  #plt.subplot(rows,cols,i)
  #plot_makers.p(fname)
  #i+=1 
  
  plt.subplot(rows,cols,i)
  plot_makers.p_on_s(hname,fname,opts.eq)
  i+=1 
  
  plt.subplot(rows,cols,i)
  plot_makers.s_t(hname,opts.eq)
  i+=1 
  
  plt.subplot(rows,cols,i)
  plot_makers.energy(fname)
  i+=1 
  
  plt.subplot(rows,cols,i)
  plot_makers.edot(fname)
  i+=1 


  #The frequencies
  plt.subplot(rows,cols,i)
  plot_makers.omegas(pname,hname,opts.eq)
  i+=1 
  plt.ylim([0.01,0.06])

  #plt.subplot(rows,cols,i)
  #plot_makers.omegas_from_abdul(wname,hname,pname)
  #i+=1 
  
  plt.subplot(rows,cols,i)
  plot_makers.ks(pname,hname,opts.eq,False)
  i+=1
  plt.ylim([0.5,1])
  
  plt.subplot(rows,cols,i)
  plot_makers.ks(pname,hname,opts.eq,True)
  i+=1
  plt.ylim([0.5,1])
  
  #plt.subplot(rows,cols,i)
  #plot_makers.theta_t(hname)
  #i+=1
  
  #add resonant lines
  if (opts.restimes):
    if (opts.eq):
      rphi_t = plot_makers.get_resonant_times(pname,hname,
          opts.eq)
      for i in range(1,7):
        plt.subplot(rows,cols,i)
        for r in rphi_t:
          print r
          plt.axvline(r)
    else:
      rphi_t,rth_t,thphi_t = plot_makers.get_resonant_times(pname,hname,
          opts.eq)
      for i in range(1,7):
        plt.subplot(rows,cols,i)
        for r in rphi_t:
          plt.axvline(r)
        for r in rth_t:
          plt.axvline(r)
        for r in thphi_t:
          plt.axvline(r)

  if (opts.d):
    fig = plt.gcf()
    plt.show()
  else:
    pathname = "/home/adlewis/m/ResonanceAnalysis/Notes/GenericBBH/Notes/plots_from_aarons_runs/"
    outname = opts.name+".pdf" 
    fullname = pathname+outname

    fig = plt.gcf()
    fig.set_size_inches(18.5,10.5)
    fig.savefig(fullname)
    plt.show()

make_plots()
