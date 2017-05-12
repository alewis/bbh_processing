import numpy as np
import scipy
import matplotlib.pyplot as plt
import HorizonsTools  
import QuasiPeriodicTimeSeries as TS
import wigner_ville
import stft_and_plotting
r,rmag=HorizonsTools.SeparationVector()
#theta = HorizonsTools.AngleWithSpin()
t =np.asarray(HorizonsTools.Times())

#maxt = 100
#t = np.linspace(0,maxt,maxt*100)

def monotonic(t):
  #parabolic decay
  r0 = 1000. #initial separation
  tf = 100. #time at merger
  a = -r0/(tf*tf)
  return a*t*t + r0

def oscillating(t):
  f0 = 0.1 #initial frequency
  k = 0.01 #chirp param.
  A = 100. #initial amplitude
  c = 0.   #amplitude param.
  return (A-c*t)*np.sin(2*np.pi*(f0+k*t)*t)

def full(t):
  return monotonic(t) + oscillating(t)

#r = full(t)
#plt.plot(t,r)
#plt.show()
#f0 = 20
#Ts = [f0, f0/2., f0/3., f0/4.] #the shifts for the imbedding plots
#shifted_rs = [full(t+T) for T in Ts]  #r data, shifted by T

shift = 50
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')
plt.plot(rmag[:-2*shift],rmag[shift:-shift],rmag[2*shift:])
plt.show()
#four subplots; row and column sharing
#f, ((ax1,ax2), (ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
#ax1.set_title('T='+str(Ts[0]))
#ax1.plot(r,shifted_rs[0])

#ax2.set_title('T='+str(Ts[1]))
#ax2.plot(r,shifted_rs[1])

#ax3.set_title('T='+str(Ts[2]))
#ax3.plot(r,shifted_rs[2])

#ax4.set_title('T='+str(Ts[3]))
#ax4.plot(r,shifted_rs[3])

#plt.show()
