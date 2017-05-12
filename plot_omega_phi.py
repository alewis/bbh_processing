#import bbh_processing.plot_makers as plot_makers
from bbh_processing.combined_tools import omega_phi
import matplotlib.pyplot as plt
import numpy as np
def make_plot(om_name="Omega.pkl"):
  #plot omega_phi
  t, om = omega_phi(om_name)
  
  out=file("omega_phi.dat",'w')
  out.write("""#
    # [1] = t
    # [2] = omega_phi 
    """)
  np.savetxt(out,np.transpose([t,om]))

def helloworld():
  print "helloworld"

def show_plot(om_name="Omega.pkl"):
  print "Pkl file: ", om_name
  sp1 = plt.subplot(1,1,1)
  t, om = omega_phi(om_name)
  return plt.plot(t,om)

if __name__ == "__main__":
  make_plot()
