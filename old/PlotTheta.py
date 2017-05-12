import numpy as np
import HorizonsTools  
import matplotlib.pyplot as plt
theta=HorizonsTools.AngleWithSpin()
time =HorizonsTools.Times()

plt.plot(time,theta)
plt.show()
