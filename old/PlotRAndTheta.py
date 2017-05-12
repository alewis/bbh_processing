import numpy as np
import HorizonsTools  
import matplotlib.pyplot as plt

theta=HorizonsTools.AngleWithSpin()
r,rmag = HorizonsTools.SeparationVector()
time =HorizonsTools.Times()

plt.plot(time,theta,label="theta")
plt.plot(time,rmag,label="r")
plt.legend()
plt.show()
