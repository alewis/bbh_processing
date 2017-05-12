import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import bbh_processing.utilites as bbh
import mpl_toolkits.mplot3d.axes3d as p3

fig = plt.figure()
ax = p3.Axes3D(fig)
ax.autoscale(enable=False)
ax.set_xlim3d([-1.0,1.0])
ax.set_ylim3d([-1.0,1.0])
ax.set_zlim3d([-1.0,1.0])
#ax = fig.add_subplot(111,autoscale_on=False,projection='3d',
#    xlim=(-1,1),ylim=(-1,1),zlim=(-1,1))

import bbh_processing.hfile_tools as hf
time = hf.t("Horizons.h5")
rhat = hf.r_hat("Horizons.h5")
chi_hat = hf.chi_hat("Horizons.h5")

line, = ax.plot([],[],'o-',lw=2)
time_template = 'time = %.1fM'
time_text = ax.text(0.05,0.9,0.9,'',transform=ax.transAxes)
def init():
  line.set_data([],[])
  line.set_3d_properties([])
  time_text.set_text('')
  return line, #time_text


def animate_lines(i):
  this_x = [rhat[i][0],0,chi_hat[i][0]]
  this_y = [rhat[i][1],0,chi_hat[i][1]]
  this_z = [rhat[i][2],0,chi_hat[i][2]]
  line.set_data(this_x,this_y)
  line.set_3d_properties(this_z)
  time_text.set_text(time_template%(time[i]))
  return line, #time_text

ani = animation.FuncAnimation(fig,animate_lines,frames=np.arange(1,len(rhat)),
    interval=10,blit=False)
plt.show()

