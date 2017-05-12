import matplotlib.pyplot as plt
"Test my angle script"

#Circle of radius 1 and angular velocity 1.
import numpy as np
theta = np.linspace(0,np.pi,1000)
coord = [(z[0],z[1]) for z in zip(np.cos(theta),np.sin(theta))]
x_axis = [(1,0) for t in theta]
from bbh_processing.utilites import time_series_angle
new_theta = time_series_angle(coord,x_axis)
plt.plot(theta,theta-new_theta)
plt.show()
