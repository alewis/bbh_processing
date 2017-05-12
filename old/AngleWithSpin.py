import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import cPickle as pickle
f = open("Omega.pkl","r")
Omega, OmegaMag = pickle.load(f)

plt.plot(OmegaMag[:,0],OmegaMag[:,1])
plt.show()


