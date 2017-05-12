import bbh_processing.omega_theta_methods as omth
import bbh_processing.fit_sin_to_theta as fitsin
#import bbh_processing.plot_makers as plot_makers
import matplotlib.pyplot as plt
import numpy as np
def make_plot(hname="Horizons.h5"):
  #plot omega_theta 
  #sp1 = plt.subplot(1,1,1)
  om_theta = omth.omega_theta(hname)
  print om 
  out=file("omega_theta.dat",'w')
  out.write("""#
    # [1] = t
    # [2] = omega_theta 
    """)
  np.savetxt(out,om_theta)

def show_plot(hname="Horizons.h5"):
  print "Horizons file: ", hname
  #sp1 = plt.subplot(2,1,1)
  #fitsin.fit_sinusoids(width=1,hname=hname) 
  #fitsin.fit_sin_from_guessed_params(width=4) 
  #om_theta_fit4 = omth.omega_theta_from_fits(width=4,hname=hname)
  #om_theta_rng = omth.omega_theta_from_peaks(width=2,hname=hname) 
  #om_theta_chirp = omth.omega_theta_from_chirps(width=4,hname=hname) 
  #plt.legend()
  
  om_theta_spline = omth.omega_theta_from_splines(hname=hname,
      offset_order = 3,
      amp_order = 3,
      amp_smoothness = 0,
      signal_order = 3,
      signal_smoothness = 0)
  
  #sp2=plt.subplot(2,1,2) 
  #plt.plot(om_theta_fit2[:,0],om_theta_fit2[:,1],label="fit2")
  #plt.plot(om_theta_fit3[:,0],om_theta_fit3[:,1],label="fit3")
  #plt.plot(om_theta_fit4[:,0],om_theta_fit4[:,1],label="fit4")
  #plt.plot(om_theta_chirp[:,0],om_theta_chirp[:,1],label="chirp")
  #plt.plot(om_theta_fit5[:,0],om_theta_fit5[:,1],label="fit5")
  #plt.plot(om_theta_fit6[:,0],om_theta_fit6[:,1],label="fit6")
  #plt.plot(om_theta_rng[:,0],om_theta_rng[:,1],label="rng")
  #plt.legend()
if __name__ == "__main__":
  show_plot()
  plt.show()
  #make_plot()
