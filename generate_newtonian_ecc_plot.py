#Make the ecc at apastron frequency Omega_0 vs apastron radius r_o 
#colour plot.

#Omega_a^2 = (1-e) / r_a^3
#e = Omega_a^2 r_a^3 - 1
import matplotlib.pyplot as plt
import matplotlib.colors as pltcolors
import numpy as np

  
#PARAMETERS
def make_newtonian_colourmap(
  OMEGA_MIN = 0.,      #minimum frequency 
  OMEGA_MAX = 0.03,     #maximum frequency
  NOMEGA    = 1000,     #number of omega points 
  R_MIN     = 0.,      #minimum radius 
  R_MAX     = 50.,    #maximum radius
  NR        = 1000,     #number of r points
  SAVEFILE  = True,    #whether to save a .dat file of the results
  SHOWPLOT  = True):   #whether to show the colour map immediately

  #make numpy arrays representing omega and r
  r = np.linspace(R_MIN,R_MAX,NR)
  omega = np.linspace(OMEGA_MIN,OMEGA_MAX,NOMEGA)
  r_v, om_v = np.meshgrid(r,omega)

  e = 1-np.power(r_v,3) * np.power(om_v,2)
  unbound_indices = e>1.
  e[unbound_indices]=1.
  neg_indices = e<0.
  e[neg_indices]=0.
  #x and y are bounds, so e should be the value *inside* those bounds.
  #therefore, remove the last value from the e array
  #e = r_v*om_v - 1
  etrim = e[:-1,:-1]
  print e
  #do threshold filter on e (set everything > 1 to 1)
  #iprint e.shape

  fig, ax = plt.subplots(1)
  
  #r=np.log(ri)
  #e=np.log(e)
  p = ax.pcolorfast(r,omega,etrim,cmap=plt.get_cmap('Oranges'),
      vmax=1.,vmin=0.)
  fig.colorbar(p)
  levels = np.arange(0,1.0,0.1)
  c = ax.contour(r,omega,e,levels)
  ax.clabel(c,inline=1,fontsize=10)
  
  #ax.set_xscale("symlog")
  #ax.set_yscale("symlog")
  #ax.set_xticks(r)
  #ax.set_yticks(omega)
  plt.show()

make_newtonian_colourmap()
