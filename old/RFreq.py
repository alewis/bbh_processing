import numpy as np
import matplotlib.pyplot as plt
import HorizonsTools  

data=HorizonsTools.SeparationVector()
time =HorizonsTools.Times()
om_t_times,om_t = HorizonsTools.RollingSinusoid(data,time)
 
plt.ylabel("Ordinary frequency")
plt.xlabel("Time (CM frame)")
#plt.title("Frequency of spin vector oscillations")
#plt.plot(OmegaMag[:,0],OmegaMag[:,1]/(2.0*pi), label = 'Orbital frequency')
#plt.plot(time,data, label = 'data')
plt.plot(om_t_times,om_t, label = 'r frequency (Om_t)')
#plt.legend()
plt.show()
