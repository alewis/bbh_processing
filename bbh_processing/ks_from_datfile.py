import numpy as np
import matplotlib.pyplot as plt

def ks_from_input(
    omr,omphi,
    omtheta=None,
    equatorial=False
    ): 
    
    if not np.array_equal(omr[:,0],omphi[:,0]):
      raise ValueError("Times didn't line up")
    #r_over_phi=np.zeros(omr.shape)
    #r_over_phi[:,0]=omr[:,0]
    t = omr[:,0]
    r_over_phi=omr[:,1]/omphi[:,1]
    #plt.plot(t,r_over_phi) 
    if not equatorial:
      #interpolate omtheta onto the same time series as the others
      from scipy import interpolate
      omthspline = interpolate.UnivariateSpline(
          omtheta[:,0],omtheta[:,1],k=3,s=0)
      newomth = omthspline(t)
      r_over_theta = omr[:,1]/newomth
      theta_over_phi = newomth/omphi[:,1]
      return t, r_over_phi, r_over_theta, theta_over_phi     
    return t, r_over_phi

def ks_from_datfile(
    equatorial=False,
    omega_r="omega_r.dat",
    omega_phi="omega_phi.dat",
    omega_theta="omega_theta.dat"):
    
    omr = np.loadtxt(omega_r)
    omphi = np.loadtxt(omega_phi)
    if not np.array_equal(omr[:,0],omphi[:,0]):
      raise ValueError("Times didn't line up")
    #r_over_phi=np.zeros(omr.shape)
    #r_over_phi[:,0]=omr[:,0]
    t = omr[:,0]
    r_over_phi=omr[:,1]/omphi[:,1]
    #plt.plot(t,r_over_phi) 
    if not equatorial:
      omtheta = np.loadtxt(omega_theta)
       
      #trim the left side
      #shift = omr[0,0]
      #l = omtheta[0,0]
      #lidx = 2*(l-shift)
      
      #t = t[lidx:-1]  
      #trimr = omr[lidx:-1,1]
      #trimphi = omphi[lidx:-1,1]
      #r_over_phi = r_over_phi[lidx:-1]
  
      #trim the right side
      #shift2 = t[0]
      #ridx = 2*(omr[-1,0]-shift2)
      #trimtheta = omtheta[:ridx,1] 
      
      #interpolate omtheta onto the same time series as the others
      from scipy import interpolate
      omthspline = interpolate.UnivariateSpline(
          omtheta[:,0],omtheta[:,1],k=3,s=0)
      newomth = omthspline(t)
      r_over_theta = omr[:,1]/newomth
      theta_over_phi = newomth/omphi[:,1]
     
      f = open("ks.dat","w")
      f.write("""#
        # [1] = t
        # [2] = r/phi
        # [3] = r/theta
        # [4] = theta/phi
        """)
      np.savetxt(f,np.transpose([t,r_over_phi,r_over_theta,theta_over_phi]))
    
    else:  
      f = open("ks.dat","w")
      f.write("""#
        # [1] = t
        # [2] = r/phi
        """)
      np.savetxt(f,np.transpose([t,r_over_phi]))
       
if __name__=="__main__":
  ks_from_datfile(False)

