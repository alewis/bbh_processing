import numpy as np
import HorizonsTools  
import cPickle as pickle
import matplotlib.pyplot as plt

#radial frequency
r,rmag=HorizonsTools.SeparationVector()
rtime =HorizonsTools.Times()
om_r_t,om_r,rfit = HorizonsTools.RollingPowerLawSinusoid(rmag,rtime)

#theta frequency
theta=HorizonsTools.AngleWithSpin()
thtime =HorizonsTools.Times()
om_t_t,om_t,tfit = HorizonsTools.RollingSinusoid(theta,thtime)

#orbital frequency
f=open("Omega.pkl",'r')
Omega,OmegaMag = pickle.load(f)
f.close()

OmegaMagReduce = OmegaMag[:,1] / (2.0*np.pi)
res = [t*3/4 for t in om_t]

plt.plot(om_r_t,rfit)
plt.plot(rtime,rmag)
plt.show()
plt.plot(om_t_t,tfit)
plt.plot(thtime,theta)
plt.show()
plt.plot(OmegaMag[:,0],OmegaMagReduce,label="Orbital Frequency / 2.0Pi")
plt.plot(om_r_t,om_r,label="Radial frequency")
plt.plot(om_t_t,om_t,label="Theta frequency")
plt.plot(om_t_t,res,label="Resonant radial frequency")
plt.legend()
plt.show()

