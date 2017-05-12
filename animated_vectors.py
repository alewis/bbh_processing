import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import bbh_processing.utilites as bbh
fig = plt.figure()
ax = fig.add_subplot(111,autoscale_on=False,xlim=(-1,1),ylim=(-1,1))
ax.grid()

import bbh_processing.hfile_tools as hf
time = hf.t("Horizons.h5")
rhat = hf.r_hat("Horizons.h5")
chi_hat = hf.chi_hat("Horizons.h5")

line, = ax.plot([],[],'o-',lw=2)
origin=[0,0]
time_template = 'time = %.1fM'
time_text = ax.text(0.05,0.9,'',transform=ax.transAxes)
def init():
  line.set_data([],[])
  time_text.set_text('')
  return line, time_text


def animate_lines(i):
  this_x = [rhat[i][0],0,chi_hat[i][0]]
  this_y = [rhat[i][1],0,chi_hat[i][1]]
  #tihis_r = [rvec[i][0],rvec[i][1]]
  #this_chi = [chi[i][0],chi[i][1]]
  line.set_data(this_x,this_y)
  time_text.set_text(time_template%(time[i]))
  return line, time_text

ani = animation.FuncAnimation(fig,animate_lines,frames=np.arange(1,len(rhat)),
    init_func=init,interval=1,blit=True)

plt.show()

