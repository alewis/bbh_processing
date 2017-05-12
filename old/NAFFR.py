import numpy as np
import scipy
import matplotlib.pyplot as plt
import HorizonsTools  
import QuasiPeriodicTimeSeries as TS
import wigner_ville
import stft_and_plotting
r,rmag=HorizonsTools.SeparationVector()
theta = HorizonsTools.AngleWithSpin()
t =np.asarray(HorizonsTools.Times())[0:10000]

x = np.asarray(rmag)[0:10000]


#sinusoid
f0 = 440 #440 hz sinusoid
fs = 8000 #sampled at 8 kHz
T = 5 #for 5 seconds
framesz = 0.05 #frame size of 50ms
hop = 0.010 #hop size of 20 ms
#t = scipy.linspace(0,T,T*fs,endpoint=False)
#x = scipy.sin(2*scipy.pi*f0*t)
#X = stft_and_plotting.stft(x,framesz)

plt.imshow(scipy.absolute(X.T),origin='lower',aspect='auto',interpolation='nearest')
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.show()



#N = t.shape[-1]
#sp=np.fft.fft(f)
#fr = np.fft.fftfreq(N,tstep)
#plt.subplot(212)
#plt.plot(fr,(2/N)*np.abs(sp))
#plt.show()
#print "Entering PWD"
#import cProfile
#cProfile.run('freq,y=PseudoWignerDistribution.PseudoWigner(t,f,2*len(f)+1)')
#freq,y=PseudoWignerDistribution.PseudoWigner(t,f, 2*len(f)+1)
print "Done"


#import cPickle as pickle
#fname = "PickledPWD.pkl"
#f = open(fname,'w')
#pickle.dump([freq,y],f)
#f.close()


#plotting
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111,projection='3d')
print t.shape
print freq.shape
print y.shape
#ax.plot_wireframe(t,freq,np.abs(y))
#plt.pcolormesh(t,freq,np.abs(y))
#plt.colorbar()
#plt.ylabel("Frequency")
#plt.xlabel("Time")
#plt.show()

