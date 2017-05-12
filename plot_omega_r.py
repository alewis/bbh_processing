from bbh_processing.combined_tools import omega_r
#import bbh_processing.plot_makers as plot_makers
import matplotlib.pyplot as plt
import numpy as np
def make_plot(
    om_name="Omega.pkl",
    outname="plots.pdf"):
  #plot omega_r 
  sp1 = plt.subplot(1,1,1)
  t, om = omega_r(om_name)
  
  out=file("omega_r.dat",'w')
  out.write("""#
    # [1] = t
    # [2] = Omega_r 
    """)
  np.savetxt(out,np.transpose([t,om]))

def show_plot(om_name="Omega.pkl"):
  print "Pkl file: ", om_name
  sp1 = plt.subplot(1,1,1)
  t, om = omega_r(om_name)
  return plt.plot(t,om)

if __name__ == "__main__":
  make_plot()
