import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
g = open("PWD.pkl",'r')
t,freq,y = pickle.load(g)
print "Loaded!"
plt.imshow(y)
#plt.pcolormesh(t,freq,y)
plt.colorbar()
plt.ylabel("Frequency")
plt.xlabel("Time")
plt.show()
